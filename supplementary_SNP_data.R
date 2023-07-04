library(readxl)
library(plyr)
library(dplyr)
library(tidyr)
library(lubridate)
library(writexl)

# import sample metadata.
samplemeta <- read_xlsx("sample_metadata.xlsx", col_types = c("numeric","text","text","date","date"))

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

# output data to Excel file.
write_xlsx(combined, path = "supplementary_SNP_data.xlsx")
