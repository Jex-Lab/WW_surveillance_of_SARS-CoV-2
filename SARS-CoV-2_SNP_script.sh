#!/bin/sh

# assign path variable.
FILEPATH=`pwd -P`

## ----------------------------- ##
##        IMPORT METADATA        ##
## ----------------------------- ##

# delete carriage returns.
sed 's/\r$//' -i $FILEPATH/oh_names.txt
sed 's/\r$//' -i $FILEPATH/reference_files/codons.txt

# read names to arrays.
grep -o -P '(?<=fastq,).*(?=_oh1)' $FILEPATH/oh_names.txt > $FILEPATH/names.txt
readarray -t names <  $FILEPATH/names.txt

grep "oh1" $FILEPATH/oh_names.txt > $FILEPATH/oh1_names.txt
readarray -t oh1 < $FILEPATH/oh1_names.txt

grep "oh2" $FILEPATH/oh_names.txt > $FILEPATH/oh2_names.txt
readarray -t oh2 <  $FILEPATH/oh2_names.txt

## ----------------------------- ##
##        REMOVE OVERHANGS       ##
## ----------------------------- ##

# DEMULTIPLEXED FILES MUST BE IN DIRECTORY CALLED "demulti".

mkdir $FILEPATH/trim

for i in $(ls $FILEPATH/demulti/)
do
cutadapt -g GTGACCTATGAACTCAGGAGTC $FILEPATH/demulti/${i} | cutadapt -a GCTGCGATGTGCAAGTCTCAG - -o $FILEPATH/trim/${i}
done

## ----------------------------- ##
##         GET OH1 READS         ##
## ----------------------------- ##

mkdir $FILEPATH/oh1

for i in ${oh1[@]}
do
cp $FILEPATH/trim/"$i" oh1
done

## ----------------------------- ##
##        GET OH1 5' READS       ##
## ----------------------------- ##

mkdir $FILEPATH/oh1_5

for i in $(ls $FILEPATH/oh1/)
do
seqkit grep $FILEPATH/oh1/${i} -s -P -r -p ^TTTTCCAATGTTACTTGGTTCCAT --out-file $FILEPATH/oh1_5/${i}
done

## ----------------------------- ##
##       GET OH1 3' READS        ##
## ----------------------------- ##

mkdir $FILEPATH/oh1_3

for i in $(ls $FILEPATH/oh1/)
do
seqkit grep $FILEPATH/oh1/${i} -s -P -r -p ^CAAATCGCTCCAGGGCAAAC --out-file $FILEPATH/oh1_3/${i}
done

## ----------------------------- ##
##   TRIM & FILTER OH1 5' READS  ##
## ----------------------------- ##

mkdir $FILEPATH/oh1_5_trim

for i in $(ls $FILEPATH/oh1_5/)
do
trimmomatic SE -threads 24 -phred33 $FILEPATH/oh1_5/${i} $FILEPATH/oh1_5_trim/${i} HEADCROP:24 CROP:144 SLIDINGWINDOW:5:28 MINLEN:144
done

## ----------------------------- ##
##   TRIM & FILTER OH1 3' READS  ##
## ----------------------------- ##

mkdir $FILEPATH/oh1_3_trim

for i in $(ls $FILEPATH/oh1_3/)
do
trimmomatic SE -threads 24 -phred33 $FILEPATH/oh1_3/${i} $FILEPATH/oh1_3_trim/${i} HEADCROP:20 CROP:134 SLIDINGWINDOW:5:28 MINLEN:134
done

## ----------------------------- ##
##         GET OH2 READS         ##
## ----------------------------- ##

mkdir $FILEPATH/oh2

for i in ${oh2[@]}
do
cp $FILEPATH/trim/"$i" oh2
done

## ----------------------------- ##
##        GET OH2 5' READS       ##
## ----------------------------- ##

mkdir $FILEPATH/oh2_5

for i in $(ls $FILEPATH/oh2/)
do
seqkit grep $FILEPATH/oh2/${i} -s -P -r -p ^TGCAATTATTCGCACTAGAATAAAC --out-file $FILEPATH/oh2_5/${i}
done

## ----------------------------- ##
##        GET OH2 3' READS       ##
## ----------------------------- ##

mkdir $FILEPATH/oh2_3

for i in $(ls $FILEPATH/oh2/)
do
seqkit grep $FILEPATH/oh2/${i} -s -P -r -p ^AGTACTACTACTCTGTATGGTTGGT.AC --out-file $FILEPATH/oh2_3/${i}
done

## ----------------------------- ##
##     REMOVE PDs FROM OH2 3'    ##
## ----------------------------- ##

mkdir $FILEPATH/oh2_3_pd1

for i in $(ls $FILEPATH/oh2_3/)
do
seqkit grep $FILEPATH/oh2_3/${i} -s -p TAACACACTGACTAGAGACTAGTGG -v --out-file $FILEPATH/oh2_3_pd1/${i}
done

mkdir $FILEPATH/oh2_3_pd2

for i in $(ls $FILEPATH/oh2_3_pd1/)
do
seqkit grep $FILEPATH/oh2_3_pd1/${i} -s -p AGTAAAGCAGAGATCATTTAATTTAGTAGG -v --out-file $FILEPATH/oh2_3_pd2/${i}
done

## ----------------------------- ##
##   TRIM & FILTER OH2 5' READS  ##
## ----------------------------- ##

mkdir $FILEPATH/oh2_5_trim

# the rev primer is partially retained to enable identification of DEL156-157

for i in $(ls $FILEPATH/oh2_5/)
do
trimmomatic SE -threads 24 -phred33 $FILEPATH/oh2_5/${i} $FILEPATH/oh2_5_trim/${i} HEADCROP:16 CROP:143 SLIDINGWINDOW:5:28 MINLEN:143
done

## ----------------------------- ##
##   TRIM & FILTER OH2 3' READS  ##
## ----------------------------- ##

mkdir $FILEPATH/oh2_3_trim

for i in $(ls $FILEPATH/oh2_3_pd2/)
do
trimmomatic SE -threads 24 -phred33 $FILEPATH/oh2_3_pd2/${i} $FILEPATH/oh2_3_trim/${i} HEADCROP:28 CROP:134 SLIDINGWINDOW:5:28 MINLEN:134
done

## ------------------------------##
##       REV COMP OH2 READS      ##
## ------------------------------##

mkdir $FILEPATH/oh2_5_rev

for i in $(ls $FILEPATH/oh2_5_trim/)
do
seqkit seq $FILEPATH/oh2_5_trim/${i} -r -p --out-file $FILEPATH/oh2_5_rev/${i}
done

mkdir $FILEPATH/oh2_3_rev

for i in $(ls $FILEPATH/oh2_3_trim/)
do
seqkit seq $FILEPATH/oh2_3_trim/${i} -r -p --out-file $FILEPATH/oh2_3_rev/${i}
done

## ----------------------------- ##
##    COMBINE OH1 & OH2 READS    ##
## ----------------------------- ##

mkdir $FILEPATH/combi

sampNos1=`echo ${#names[@]}`
let "sampNos2=$sampNos1-1"

for i in $(seq 0 1 $sampNos2)
do
cat $FILEPATH/oh1_5_trim/${oh1[i]} $FILEPATH/oh1_3_trim/${oh1[i]} $FILEPATH/oh2_5_rev/${oh2[i]} $FILEPATH/oh2_3_rev/${oh2[i]} > $FILEPATH/combi/${names[i]}.fastq
done

## ----------------------------- ##
##        MAP READS TO REF       ##
## ----------------------------- ##

mkdir $FILEPATH/bams

for i in ${names[@]}
do
bwa-mem2 mem -t 24 -R '@RG\tID:'${i}'\tSM:'${i}'' $FILEPATH/reference_files/s_gene.fa $FILEPATH/combi/${i}.fastq | samtools view -bo $FILEPATH/bams/${i}.bam -
done

## ----------------------------- ##
##   SORT & INDEX MAPPED READS   ##
## ----------------------------- ##

mkdir $FILEPATH/sorted

for i in ${names[@]}
do
samtools sort $FILEPATH/bams/${i}.bam -o $FILEPATH/sorted/${i}.sorted.bam
done

for i in ${names[@]}
do
samtools index $FILEPATH/sorted/${i}.sorted.bam
done

## ----------------------------- ##
##     MASK NON-TARGET REGIONS   ##
## ----------------------------- ##

mkdir $FILEPATH/masks

for i in ${names[@]}
do
bedtools maskfasta -fi $FILEPATH/reference_files/s_gene.fa -bed $FILEPATH/reference_files/s_gene_mask.bed -fo $FILEPATH/masks/${i}.masked.fa
done

## ----------------------------- ##
##     ADD INDEL QUALITY INFO    ##
## ----------------------------- ##

mkdir $FILEPATH/indelq

for i in ${names[@]}
do
lofreq indelqual $FILEPATH/sorted/${i}.sorted.bam --dindel -f $FILEPATH/masks/${i}.masked.fa -o $FILEPATH/indelq/${i}.indel.bam
done

## ----------------------------- ##
##            RE-INDEX           ##
## ----------------------------- ##

for i in ${names[@]}
do
samtools index $FILEPATH/indelq/${i}.indel.bam
done

## ----------------------------- ##
##        CALL VARIANTS          ##
## ----------------------------- ##

mkdir $FILEPATH/vcfs

N=20

for i in ${names[@]}
do
((j=j%N))
((j++==0)) && wait
lofreq call --call-indels -f $FILEPATH/masks/${i}.masked.fa -o $FILEPATH/vcfs/${i}.vcf $FILEPATH/indelq/${i}.indel.bam &
done

## ----------------------------- ##
##      FILTER VARIANT CALLS     ##   
## ----------------------------- ##

mkdir $FILEPATH/filtered

for i in ${names[@]}
do
lofreq filter -i $FILEPATH/vcfs/${i}.vcf -o $FILEPATH/filtered/${i}.filtered.vcf --af-min 0.02 --cov-min 400
done

## ----------------------------- ##
##    SUMMARISE MUTATION DATA    ##
## ----------------------------- ##

mkdir $FILEPATH/summary

for i in ${names[@]}
do
gatk VariantsToTable -F POS -F REF -F ALT -F AF -F DP -V $FILEPATH/filtered/${i}.filtered.vcf -O $FILEPATH/summary/${i}.txt
done

## ----------------------------- ##
##     REFORMAT MUTATION DATA    ##
## ----------------------------- ##

# update nucleotide numbering.
for i in ${names[@]}
do
awk '$1+=21562' $FILEPATH/summary/${i}.txt > $FILEPATH/summary/${i}_snps.txt
done

# remove header.
for i in ${names[@]}
do
sed '1d' $FILEPATH/summary/${i}_snps.txt > $FILEPATH/summary/${i}_snps1.txt
done

# remove leading space.
for i in ${names[@]}
do
sed -r 's/\s//1' $FILEPATH/summary/${i}_snps1.txt > $FILEPATH/summary/${i}_snps2.txt
done

# insert ">" character.
for i in ${names[@]}
do
sed -r 's/\s/\>/1' $FILEPATH/summary/${i}_snps2.txt > $FILEPATH/summary/${i}_snps3.txt
done

# add file names.
for i in ${names[@]}
do
sed '1 i\'${i}'' $FILEPATH/summary/${i}_snps3.txt > $FILEPATH/summary/${i}_snps4.txt
done

# replace line breaks with spaces.
for i in ${names[@]}
do
tr "\n" " " < $FILEPATH/summary/${i}_snps4.txt > $FILEPATH/summary/${i}_snps5.txt
done

# delete trailing space.
for i in ${names[@]}
do
sed -r 's/\s$//' $FILEPATH/summary/${i}_snps5.txt > $FILEPATH/summary/${i}_snps6.txt
done

# concatenate files.
awk 'FNR==1{print ""}1' $FILEPATH/summary/*snps6.txt > $FILEPATH/summary/all_samples.txt

# remove odd numbered lines.
sed '1d; n; d' -i $FILEPATH/summary/all_samples.txt

## ----------------------------- ##
##         REPLACE NAMES         ##
## ----------------------------- ##

# create associative array.
declare -A codons

# populate array (add key-value pair data).
while IFS=',' read a b; do codons["$a"]="$b"; done < $FILEPATH/reference_files/codons.txt

# replace nucleotide names with codon names.
while IFS=' ' read line; do
for i in ${!codons[@]}; do
line=${line//$i/${codons[$i]}}
done
echo $line
done < $FILEPATH/summary/all_samples.txt > $FILEPATH/summary/all_mutations.txt

# replace every third space with forward slash.
sed '-es/ /\//'{600..3..3} -i $FILEPATH/summary/all_mutations.txt

# move final summary file into job directory.
mv $FILEPATH/summary/all_mutations.txt $FILEPATH/

wait

## ----------------------------- ##
##     TIDY UP JOB DIRECTORY     ##
## ----------------------------- ##

rm -r $FILEPATH/summary/*snps*.txt
rm -r $FILEPATH/demulti
rm -r $FILEPATH/trim
rm -r $FILEPATH/oh1
rm -r $FILEPATH/oh1_5
rm -r $FILEPATH/oh1_3
rm -r $FILEPATH/oh1_5_trim
rm -r $FILEPATH/oh1_3_trim
rm -r $FILEPATH/oh2
rm -r $FILEPATH/oh2_5
rm -r $FILEPATH/oh2_3
rm -r $FILEPATH/oh2_3_pd1
rm -r $FILEPATH/oh2_3_pd2
rm -r $FILEPATH/oh2_5_trim
rm -r $FILEPATH/oh2_3_trim
rm -r $FILEPATH/oh2_5_rev
rm -r $FILEPATH/oh2_3_rev
rm -r $FILEPATH/bams
rm -r $FILEPATH/masks
rm -r $FILEPATH/indelq
rm -r $FILEPATH/filtered
rm -r $FILEPATH/summary

## ----------------------------- ##
## ---------- THE END ---------- ##
## ----------------------------- ##

