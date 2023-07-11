library(readxl)
library(plyr)
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(writexl)

# import sample metadata.
SAMPLEMETA <- read_xlsx("suppl_table3.xlsx", col_types = c("numeric","text","text","date","date",
                                                               "text","text","text","text","text","text",
                                                               "numeric","numeric","numeric","numeric","numeric","numeric","numeric"))
# sort.
SAMPLEMETA <- SAMPLEMETA[order(SAMPLEMETA$Sample_No),]

# subset.
samplemeta <- SAMPLEMETA[,c(1:5)]

# import SNP data  - data for samples that failed sequencing will be absent.
allmut <- read.table(file = "all_mutations.txt", header = FALSE, sep = " ")

# set column name.
colnames(allmut)[1] <- "Sample_No"

# edit sample names.
allmut$Sample_No <- gsub("Repeat-", "", allmut$Sample_No)
allmut$Sample_No <- gsub("Repeat_", "", allmut$Sample_No)
allmut$Sample_No <- gsub("_2$", "", allmut$Sample_No)

# add names of samples that failed sequencing to SNP dataframe.
`%ni%` <- Negate(`%in%`)
seqfail <- data.frame(samplemeta$Sample_No[samplemeta$Sample_No %ni% allmut$Sample_No])
colnames(seqfail) <- c("Sample_No")
allmut <- rbind.fill(allmut, seqfail)

# set column selection variables.
x <- c("V3","V5","V7","V9","V11","V13","V15","V17","V19","V21","V23","V25","V27","V29","V31","V33","V35","V37")
y <- c("V2","V4","V6","V8","V10","V12","V14","V16","V18","V20","V22","V24","V26","V28","V30","V32","V34","V36")

# remove counts.
allmut <- allmut %>% mutate(across(all_of(x), ~substr(. , 1, 8)))

# multiply % by 100.
allmut <- allmut %>% mutate(across(all_of(x), ~(as.numeric(.) *100)))

# round %.
allmut <- allmut %>% mutate(across(all_of(x), ~round(. , 1)))

# remove incidental snps (synonymous, STOP, uninterpreted).
allmut <- allmut %>% mutate(across(all_of(y), ~gsub("\\^", NA, .)))
allmut <- allmut %>% mutate(across(all_of(y), ~gsub("STOP", NA, .)))
allmut <- allmut %>% mutate(across(all_of(y), ~gsub(">", NA, .)))

# combine snp names and %.
snps <- cbind(allmut[1], mapply(paste0, allmut[, seq(2, 37, 2)], sep = "|", allmut[, seq(3, 37, 2)]))

# get characteristic snps.
char <- snps %>% mutate(across(all_of(y), ~gsub("\\*", NA, .)))
char <- char %>% mutate(across(all_of(y), ~gsub("NA", NA, .)))

# get characteristic snp names.
charn <- char %>% mutate(across(all_of(y), ~gsub("\\|.*", "", .)))

# get non-characteristic snps.
nonc <- snps %>% mutate(across(all_of(y), ~gsub("[A-Z]\\|", NA, .)))
nonc <- nonc %>% mutate(across(all_of(y), ~gsub("[0-9]\\|", NA, .)))
nonc <- nonc %>% mutate(across(all_of(y), ~gsub("\\|NA", NA, .)))

# get non-characteristic snp names.
noncn <- nonc %>% mutate(across(all_of(y), ~gsub("\\*\\|.*", "", .)))

# get non-characteristic snp %.
noncp <- nonc %>% mutate(across(all_of(y), ~gsub(".*\\|", "", .)))

# combine characteristic snp names.
charn <- cbind(charn[,1], (unite(charn[,2:19], `SNPs detected`, sep = ", ", na.rm = TRUE)))
colnames(charn)[1] <- "Sample_No"

# combine data.
combined <- merge(samplemeta, charn, by="Sample_No", all.y = FALSE)

# combine non-characteristic % names.
noncn <- cbind(noncn[,1], (unite(noncn[,2:19], `Non-characteristic SNPs`, sep = ", ", na.rm = TRUE)))
colnames(noncn)[1] <- "Sample_No"

# combine data.
combined <- merge(combined, noncn, by="Sample_No", all.y = FALSE)

# combine non-characteristic snp %.
noncp <- cbind(noncp[,1], (unite(noncp[,2:19], `Non-characteristic SNPs (%)`, sep = ", ", na.rm = TRUE)))
colnames(noncp)[1] <- "Sample_No"

# combine data.
combined <- merge(combined, noncp, by="Sample_No", all.y = FALSE)

# assess results.
combined$Result <- ""
combined$Result[grep("T95I|G142D|E156G|Del157-158|L452R|T478K", combined$`SNPs detected`)] <- "Delta VOC detected"
combined$Result[grep("A67V|Del143-145|N440K|G446S|S477N|E484A", combined$`SNPs detected`)] <- "Omicron and Delta VOC detected"
combined$Result[intersect(grep("K417N", combined$`SNPs detected`),grep("N501Y", combined$`SNPs detected`))] <- "Omicron and Delta VOC detected"
# length(grep("Omicron and Delta VOC detected", combined$Result))
combined$Result[which(str_count(combined$`SNPs detected`, ",") < 1)] <- "Sequencing failure"

# add footnotes.
combined$Footnote <- ""

combined$Footnote[which(SAMPLEMETA$`Mean Cq` > 37)] <- "2" # Quantitative virus levels below optimal performance threshold required for sequencing

aaa <- gsub("\\D"," ", combined$`SNPs detected`)
bbb <- grep("67", aaa) # NTD 5' SNPs
ccc <- grep("95|142|143|145|156|157|158", aaa) # NTD 3' SNPs
ddd <- grep("440|446|452", aaa) # RBM 5' SNPs
eee <- grep("477|478|484|493|496|498", aaa) # RBM 3' SNPs
fff <- union(union(ccc,ddd),eee) # NTD 3' | RBM 5' | RBM 3'
ggg <- setdiff(bbb,fff) # NTD 5' only
hhh <- union(union(bbb,ddd),eee) # NTD 5' | RBM 5' | RBM 3'
iii <- setdiff(ccc,hhh) # NTD 3' only
jjj <- union(union(bbb,ccc),eee) # NTD 5' | NTD 3' | RBM 5'
kkk <- setdiff(ddd,jjj) # RBM 5' only
lll <- union(union(bbb,ccc),ddd) # NTD 5' | NTD 3' | RBM 5'
mmm <- setdiff(eee,lll) # RBM 3' only
nnn <- union(union(union(ggg,iii),kkk),mmm) # all single amplicon detections
ooo <- which(combined$Result == "Sequencing failure")
ppp <- setdiff(nnn,ooo) # single amplicon detections sans seq failures
combined$Footnote[ppp] <- "1" # A valid VOC interpretation has been made from partial sequence information

# output data to Excel file.
write_xlsx(combined, path = "suppl_table2.xlsx")

# ------------------------- #

# Table 2 values

# replace Cq>45
SAMPLEMETA$`N Cq (Rep1)` <- gsub("Cq>45", NA, SAMPLEMETA$`N Cq (Rep1)`)
SAMPLEMETA$`N Cq (Rep2)` <- gsub("Cq>45", NA, SAMPLEMETA$`N Cq (Rep2)`)
SAMPLEMETA$`Orf1 Cq (Rep1)` <- gsub("Cq>45", NA, SAMPLEMETA$`Orf1 Cq (Rep1)`)
SAMPLEMETA$`Orf1 Cq (Rep2)` <- gsub("Cq>45", NA, SAMPLEMETA$`Orf1 Cq (Rep2)`)

dualPos <- sum(SAMPLEMETA$`Target positive` == "Dual")
singlePosN <- sum(na.omit(SAMPLEMETA$`Single loci` == "N"))
singlePosOrf1 <- sum(na.omit(SAMPLEMETA$`Single loci` == "Orf"))
totalSamples <- dualPos + singlePosN + singlePosOrf1

dualPosPercent <- round(dualPos/totalSamples*100, 1)
singlePosNPercent <- round(singlePosN/totalSamples*100, 1)
singlePosOrf1Percent <- round(singlePosOrf1/totalSamples*100, 1)
totalPercent <- dualPosPercent + singlePosNPercent + singlePosOrf1Percent

dualPosNCqMean <- round(mean(na.omit(SAMPLEMETA$`N Cq (AVG)`[which(SAMPLEMETA$`Target positive` == "Dual")])), 1)
dualPosOrf1CqMean <- round(mean(na.omit(SAMPLEMETA$`Orf1 Cq (AVG)`[which(SAMPLEMETA$`Target positive` == "Dual")])), 1)
singlePosNCqMean <- round(mean(na.omit(SAMPLEMETA$`N Cq (AVG)`[which(SAMPLEMETA$`Target positive` != "Dual")])), 1)
singlePosOrf1CqMean <- round(mean(na.omit(SAMPLEMETA$`Orf1 Cq (AVG)`[which(SAMPLEMETA$`Target positive` != "Dual")])), 1)

dualPosNCqRange <- round(range(na.omit(SAMPLEMETA$`N Cq (AVG)`[which(SAMPLEMETA$`Target positive` == "Dual")])), 1)
dualPosOrf1CqRange <- round(range(na.omit(SAMPLEMETA$`Orf1 Cq (AVG)`[which(SAMPLEMETA$`Target positive` == "Dual")])), 1)
singlePosNCqRange <- round(range(na.omit(SAMPLEMETA$`N Cq (AVG)`[which(SAMPLEMETA$`Target positive` != "Dual")])), 1)
singlePosOrf1CqRange <- round(range(na.omit(SAMPLEMETA$`Orf1 Cq (AVG)`[which(SAMPLEMETA$`Target positive` != "Dual")])), 1)

dualPosCopRxn <- round(mean(na.omit(SAMPLEMETA$`Interpolated (copies/rxn)`[which(SAMPLEMETA$`Target positive` == "Dual")])), 0)
singlePosNCopRxn <- round(mean(na.omit(SAMPLEMETA$`Interpolated (copies/rxn)`[which(SAMPLEMETA$`Single loci` == "N")])), 0)
singlePosOrf1CopRxn <- round(mean(na.omit(SAMPLEMETA$`Interpolated (copies/rxn)`[which(SAMPLEMETA$`Single loci` == "Orf")])), 0)

dualPosCopSampler <- round(mean(na.omit(SAMPLEMETA$`Interpolated (copies/sampler)`[which(SAMPLEMETA$`Target positive` == "Dual")])), 0)
singlePosNCopSampler <- round(mean(na.omit(SAMPLEMETA$`Interpolated (copies/sampler)`[which(SAMPLEMETA$`Single loci` == "N")])), 0)
singlePosOrf1CopSampler <- round(mean(na.omit(SAMPLEMETA$`Interpolated (copies/sampler)`[which(SAMPLEMETA$`Single loci` == "Orf")])), 0)

# ------------------------- #

# Table 3 values

# combine SNP and RT-qPCR data
SAMPLEMETA <- SAMPLEMETA[order(SAMPLEMETA$Sample_No),]
combined <- combined[order(combined$Sample_No),]
combined2 <- cbind(SAMPLEMETA, combined[,6:10])

# # output data to Excel file.
# write_xlsx(combined2, path = "supplementary_SNP_data2.xlsx")

rrr <- which(combined2$`Target positive` == "Dual")
sss <- which(combined2$`Target positive` == "Single")

ttt <- which(combined2$Result == "Omicron and Delta VOC detected")
uuu <- which(combined2$Result == "Delta VOC detected")
vvv <- which(combined2$Result == "Sequencing failure")

# Omicron & Delta detected, Dual/Single RT-PCR positive, NTD+RBM sequenced
omiDeltaDualPos <- length(setdiff((intersect(ttt,rrr)),ppp))
omiDeltaSingPos <- length(setdiff((intersect(ttt,sss)),ppp))
sub1 <- omiDeltaDualPos + omiDeltaSingPos
omiDeltaDualPosCq <- round(mean(combined2$`Mean Cq`[setdiff((intersect(ttt,rrr)),ppp)]),1)
omiDeltaDualPosRange <- range(combined2$`Mean Cq`[setdiff((intersect(ttt,rrr)),ppp)])

# Delta detected, Dual/Single RT-PCR positive, NTD+RBM sequenced
deltaDualPos <- length(setdiff((intersect(uuu,rrr)),ppp))
deltaSingPos <- length(setdiff((intersect(uuu,sss)),ppp))
sub2 <- deltaDualPos + deltaSingPos
DeltaDualPosCq <- round(mean(combined2$`Mean Cq`[setdiff((intersect(uuu,rrr)),ppp)]),1)
DeltaDualPosRange <- range(combined2$`Mean Cq`[setdiff((intersect(uuu,rrr)),ppp)])
DeltaSingPosCq <- round(mean(combined2$`Mean Cq`[setdiff((intersect(uuu,sss)),ppp)]),1)
DeltaSingPosRange <- range(combined2$`Mean Cq`[setdiff((intersect(uuu,sss)),ppp)])

# Unable to sequence, Dual/Single RT-PCR positive
unableDualPos <- length(intersect(vvv,rrr))
unableSingPos <- length(intersect(vvv,sss))
sub3 <- unableDualPos + unableSingPos
unableDualPosCq <- round(mean(combined2$`Mean Cq`[setdiff((intersect(vvv,rrr)),ppp)]),1)
unableDualPosRange <- range(combined2$`Mean Cq`[setdiff((intersect(vvv,rrr)),ppp)])
unableSingPosCq <- round(mean(combined2$`Mean Cq`[setdiff((intersect(vvv,sss)),ppp)]),1)
unableSingPosRange <- range(combined2$`Mean Cq`[setdiff((intersect(vvv,sss)),ppp)])

# Dual/Single RT-PCR positive, Either NTD or RBM sequenced
onlyNTDorRBMdualPos <- length(intersect(ppp,rrr))
onlyNTDorRBMsingPos <- length(intersect(ppp,sss))
sub4 <- onlyNTDorRBMdualPos + onlyNTDorRBMsingPos
onlyNTDorRBMDualPosCq <- round(mean(combined2$`Mean Cq`[intersect(ppp,rrr)]),1)
onlyNTDorRBMDualPosRange <- range(combined2$`Mean Cq`[intersect(ppp,rrr)])
onlyNTDorRBMSingPosCq <- round(mean(combined2$`Mean Cq`[intersect(ppp,sss)]),1)
onlyNTDorRBMSingPosRange <- range(combined2$`Mean Cq`[intersect(ppp,sss)])

sub1 + sub2 + sub3 + sub4

# ------------------------- #

# Table 4 values

dualPos
singlePosN + singlePosOrf1  

dualPosSeqSuccess <- length(intersect(union(ttt,uuu),rrr))
singPosSeqSuccess <- length(intersect(union(ttt,uuu),sss))

dualPosSeqSuccess + singPosSeqSuccess

# ------------------------- #

# Table 5 values

xxx <- which(combined2$`Mean Cq` <= 36)
yyy <- which(combined2$`Mean Cq` > 36)

dualPosCqLT36 <- length(intersect(xxx,rrr))
singPosCqLT36 <- length(intersect(xxx,sss))

dualPosCqGT36 <- length(intersect(yyy,rrr))
singPosCqGT36 <- length(intersect(yyy,sss))

dualPosCqLT36SeqSuccess <- length(intersect(intersect(union(ttt,uuu),rrr),xxx))
  
dualPosCqGT36SeqSuccess <- length(intersect(intersect(union(ttt,uuu),rrr),yyy))
singPosCqGT36SeqSuccess <- length(intersect(intersect(union(ttt,uuu),sss),yyy))
  
# ------------------------- #

# value check for manuscript

# NTD Omi == A67V, Del143-145
# RBM Omi == N440K, G446S, S477N, E484A, Q493R, G496S, Q498R | K417N & N501Y 
omiNTD <- grep("67|143|145", aaa) # Omi in NTD detected
omiRBM1 <- grep("440|446|477|484|493|496|498",aaa) # Omi in RBM detected
omiRBM2 <- intersect(grep("K417N", combined$`SNPs detected`),grep("N501Y", combined$`SNPs detected`))
omiRBM <- union(omiRBM1, omiRBM2)
combined2$Sample_No[setdiff(omiNTD,omiRBM)] # Omi detected in NTD but not RBM
combined2$Sample_No[setdiff(omiRBM,omiNTD)] # Omi detected in RBM but not in NTD

deltaPosTotal <- length(uuu)
deltaOmiTotal <- length(ttt)

# number of samples between 15 Nov 2021 to 30 Nov 2021.
length(which(combined2$Retrieved >= "2021-11-15" & combined2$Retrieved < "2021-11-30"))
