## ----------------------------- ##
## ------- SNP ANALYSIS -------- ##
## ----------------------------- ##

# assign path variable.
FILEPATH=`pwd -P`

# delete carriage returns in metadata files.
sed 's/\r$//' -i $FILEPATH/meta/sample_names.txt
sed 's/\r$//' -i $FILEPATH/meta/oh1_names.txt
sed 's/\r$//' -i $FILEPATH/meta/oh2_names.txt
sed 's/\r$//' -i $FILEPATH/../sgene_VOC_codons/codons.txt

# create sample and file name reference files.
readarray -t names < $FILEPATH/meta/sample_names.txt
readarray -t oh1 < $FILEPATH/meta/oh1_names.txt
readarray -t oh2 < $FILEPATH/meta/oh2_names.txt

# define 'for loop' value (total number of samples - 1).
sampNos1=`grep -c '^[^ ]' $FILEPATH/meta/sample_names.txt`
let "sampNos2=$sampNos1-1"

## ----------------------------- ##
##          DEMULTIPLEX          ##
## ----------------------------- ##

rm -r $FILEPATH/demulti
mkdir $FILEPATH/demulti

cutadapt \
-j 12 \
-e 0 \
--no-indels \
-g file:$FILEPATH/meta/fwd_indexes.fasta \
-G file:$FILEPATH/meta/rev_indexes.fasta \
-o $FILEPATH/demulti/{name1}_{name2}.fastq \
-p $FILEPATH/demulti/{name1}_{name2}_rev.fastq \
$FILEPATH/raw_seqs/*R1*.fastq.gz \
$FILEPATH/raw_seqs/*R2*.fastq.gz

# delete redundant files.
rm -r $FILEPATH/demulti/*_rev.fastq
rm -r $FILEPATH/demulti/*unknown*

wait

## ----------------------------- ##
##        REMOVE OVERHANGS       ##
## ----------------------------- ##

rm -r $FILEPATH/trim
mkdir $FILEPATH/trim

for i in $(ls $FILEPATH/demulti/)
do
cutadapt -g GTGACCTATGAACTCAGGAGTC $FILEPATH/demulti/${i} | cutadapt -a GCTGCGATGTGCAAGTCTCAG - -o $FILEPATH/trim/${i}
done

wait

## ----------------------------- ##
##         GET OH1 READS         ##
## ----------------------------- ##

rm -r $FILEPATH/oh1
mkdir $FILEPATH/oh1

for i in ${oh1[@]}
do
cp $FILEPATH/trim/"$i" oh1
done

wait

## ----------------------------- ##
##        GET OH1 5' READS       ##
## ----------------------------- ##

rm -r $FILEPATH/oh1_5
mkdir $FILEPATH/oh1_5

for i in $(ls $FILEPATH/oh1/)
do
seqkit grep $FILEPATH/oh1/${i} -s -P -r -p ^TTTTCCAATGTTACTTGGTTCCAT --out-file $FILEPATH/oh1_5/${i}
done

wait

## ----------------------------- ##
##       GET OH1 3' READS        ##
## ----------------------------- ##

rm -r $FILEPATH/oh1_3
mkdir $FILEPATH/oh1_3

for i in $(ls $FILEPATH/oh1/)
do
seqkit grep $FILEPATH/oh1/${i} -s -P -r -p ^CAAATCGCTCCAGGGCAAAC --out-file $FILEPATH/oh1_3/${i}
done

wait

## ----------------------------- ##
##   TRIM & FILTER OH1 5' READS  ##
## ----------------------------- ##

rm -r $FILEPATH/oh1_5_trim
mkdir $FILEPATH/oh1_5_trim

for i in $(ls $FILEPATH/oh1_5/)
do
trimmomatic SE -threads 12 -phred33 $FILEPATH/oh1_5/${i} $FILEPATH/oh1_5_trim/${i} HEADCROP:24 CROP:144 SLIDINGWINDOW:5:28 MINLEN:144
done

wait

## ----------------------------- ##
##   TRIM & FILTER OH1 3' READS  ##
## ----------------------------- ##

rm -r $FILEPATH/oh1_3_trim
mkdir $FILEPATH/oh1_3_trim

for i in $(ls $FILEPATH/oh1_3/)
do
trimmomatic SE -threads 12 -phred33 $FILEPATH/oh1_3/${i} $FILEPATH/oh1_3_trim/${i} HEADCROP:20 CROP:134 SLIDINGWINDOW:5:28 MINLEN:134
done

wait

## ----------------------------- ##
##         GET OH2 READS         ##
## ----------------------------- ##

rm -r $FILEPATH/oh2
mkdir $FILEPATH/oh2

for i in ${oh2[@]}
do
cp $FILEPATH/trim/"$i" oh2
done

wait

## ----------------------------- ##
##        GET OH2 5' READS       ##
## ----------------------------- ##

rm -r $FILEPATH/oh2_5
mkdir $FILEPATH/oh2_5

# extra bases are included in the primer sequence to detect/exclude primer-dimers

for i in $(ls $FILEPATH/oh2/)
do
seqkit grep $FILEPATH/oh2/${i} -s -P -r -p ^TGCAATTATTCGCACTAGAATAAAC --out-file $FILEPATH/oh2_5/${i}
done

wait

## ----------------------------- ##
##        GET OH2 3' READS       ##
## ----------------------------- ##

rm -r $FILEPATH/oh2_3
mkdir $FILEPATH/oh2_3

for i in $(ls $FILEPATH/oh2/)
do
seqkit grep $FILEPATH/oh2/${i} -s -P -r -p ^AGTACTACTACTCTGTATGGTTGGT.AC --out-file $FILEPATH/oh2_3/${i}
done

wait

## ----------------------------- ##
##     REMOVE PDs FROM OH2 3'    ##
## ----------------------------- ##

rm -r $FILEPATH/oh2_3_pd1
mkdir $FILEPATH/oh2_3_pd1

# primer dimers are rarely generated but it's better to be safe than sorry.

for i in $(ls $FILEPATH/oh2_3/)
do
seqkit grep $FILEPATH/oh2_3/${i} -s -p TAACACACTGACTAGAGACTAGTGG -v --out-file $FILEPATH/oh2_3_pd1/${i}
done

wait

rm -r $FILEPATH/oh2_3_pd2
mkdir $FILEPATH/oh2_3_pd2

for i in $(ls $FILEPATH/oh2_3_pd1/)
do
seqkit grep $FILEPATH/oh2_3_pd1/${i} -s -p AGTAAAGCAGAGATCATTTAATTTAGTAGG -v --out-file $FILEPATH/oh2_3_pd2/${i}
done

wait

## ----------------------------- ##
##   TRIM & FILTER OH2 5' READS  ##
## ----------------------------- ##

rm -r $FILEPATH/oh2_5_trim
mkdir $FILEPATH/oh2_5_trim

# the rev primer is partially retained to enable identification of DEL156-157

for i in $(ls $FILEPATH/oh2_5/)
do
trimmomatic SE -threads 12 -phred33 $FILEPATH/oh2_5/${i} $FILEPATH/oh2_5_trim/${i} HEADCROP:16 CROP:143 SLIDINGWINDOW:5:28 MINLEN:143
done

wait

## ----------------------------- ##
##   TRIM & FILTER OH2 3' READS  ##
## ----------------------------- ##

rm -r $FILEPATH/oh2_3_trim
mkdir $FILEPATH/oh2_3_trim

for i in $(ls $FILEPATH/oh2_3_pd2/)
do
trimmomatic SE -threads 12 -phred33 $FILEPATH/oh2_3_pd2/${i} $FILEPATH/oh2_3_trim/${i} HEADCROP:28 CROP:134 SLIDINGWINDOW:5:28 MINLEN:134
done

wait

## ------------------------------##
##       REV COMP OH2 READS      ##
## ------------------------------##

rm -r $FILEPATH/oh2_5_rev
mkdir $FILEPATH/oh2_5_rev

for i in $(ls $FILEPATH/oh2_5_trim/)
do
seqkit seq $FILEPATH/oh2_5_trim/${i} -r -p --out-file $FILEPATH/oh2_5_rev/${i}
done

wait

rm -r $FILEPATH/oh2_3_rev
mkdir $FILEPATH/oh2_3_rev

for i in $(ls $FILEPATH/oh2_3_trim/)
do
seqkit seq $FILEPATH/oh2_3_trim/${i} -r -p --out-file $FILEPATH/oh2_3_rev/${i}
done

wait

## ----------------------------- ##
##    COMBINE OH1 & OH2 READS    ##
## ----------------------------- ##

rm -r $FILEPATH/combi
mkdir $FILEPATH/combi

for i in $(seq 0 1 $sampNos2)
do
cat $FILEPATH/oh1_5_trim/${oh1[i]} $FILEPATH/oh1_3_trim/${oh1[i]} $FILEPATH/oh2_5_rev/${oh2[i]} $FILEPATH/oh2_3_rev/${oh2[i]} > $FILEPATH/combi/${names[i]}.fastq
done

wait

## ----------------------------- ##
##        MAP READS TO REF       ##
## ----------------------------- ##

rm -r $FILEPATH/bams
mkdir $FILEPATH/bams

for i in ${names[@]}
do
bwa-mem2 mem -t 12 -R '@RG\tID:'${i}'\tSM:'${i}'' $FILEPATH/../sgene_ref/s_gene.fa $FILEPATH/combi/${i}.fastq | samtools view -bo $FILEPATH/bams/${i}.bam -
done

wait

## ----------------------------- ##
##   SORT & INDEX MAPPED READS   ##
## ----------------------------- ##

rm -r $FILEPATH/sorted
mkdir $FILEPATH/sorted

for i in ${names[@]}
do
samtools sort $FILEPATH/bams/${i}.bam -o $FILEPATH/sorted/${i}.sorted.bam
done

wait

for i in ${names[@]}
do
samtools index $FILEPATH/sorted/${i}.sorted.bam
done

wait

## ----------------------------- ##
##     MASK NON-TARGET REGIONS   ##
## ----------------------------- ##

rm -r $FILEPATH/masks
mkdir $FILEPATH/masks

for i in ${names[@]}
do
bedtools maskfasta -fi $FILEPATH/../sgene_ref/s_gene.fa -bed $FILEPATH/../sgene_ref/s_gene_mask.bed -fo $FILEPATH/masks/${i}.masked.fa
done

wait

## ----------------------------- ##
##     ADD INDEL QUALITY INFO    ##
## ----------------------------- ##

rm -r $FILEPATH/indelq
mkdir $FILEPATH/indelq

for i in ${names[@]}
do
lofreq indelqual $FILEPATH/sorted/${i}.sorted.bam --dindel -f $FILEPATH/masks/${i}.masked.fa -o $FILEPATH/indelq/${i}.indel.bam
done

wait

## ----------------------------- ##
##            RE-INDEX           ##
## ----------------------------- ##

for i in ${names[@]}
do
samtools index $FILEPATH/indelq/${i}.indel.bam
done

wait

## ----------------------------- ##
##        CALL VARIANTS          ##
## ----------------------------- ##

rm -r $FILEPATH/vcfs
mkdir $FILEPATH/vcfs

# number of samples to batch process as an alternative to call-parallel which doesn't always behave.
N=10

for i in ${names[@]}
do
((j=j%N))
((j++==0)) && wait
lofreq call --call-indels -f $FILEPATH/masks/${i}.masked.fa -o $FILEPATH/vcfs/${i}.vcf $FILEPATH/indelq/${i}.indel.bam &
done

wait

## ----------------------------- ##
##      FILTER VARIANT CALLS     ##   
## ----------------------------- ##

rm -r $FILEPATH/filtered
mkdir $FILEPATH/filtered

for i in ${names[@]}
do
lofreq filter -i $FILEPATH/vcfs/${i}.vcf -o $FILEPATH/filtered/${i}.filtered.vcf --af-min 0.02 --cov-min 400
done

wait

## ----------------------------- ##
##    SUMMARISE MUTATION DATA    ##
## ----------------------------- ##

rm -r $FILEPATH/summary
mkdir $FILEPATH/summary

for i in ${names[@]}
do
gatk VariantsToTable -F POS -F REF -F ALT -F AF -F DP -V $FILEPATH/filtered/${i}.filtered.vcf -O $FILEPATH/summary/${i}.txt
done

wait

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
while IFS=',' read a b; do codons["$a"]="$b"; done < $FILEPATH/../sgene_VOC_codons/codons.txt

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
