library(readxl)
library(plyr)
library(dplyr)
library(tidyr)
library(lubridate)
library(writexl)

# get SNP data  - data for samples that failed sequencing will be absent.
run001 <- read.table(file = "all_mutations_001.txt", header = F, fill = T, na.strings = "", col.names = paste0("V", seq_len(max(count.fields("all_mutations_001.txt")))))
run002 <- read.table(file = "all_mutations_002.txt", header = F, fill = T, na.strings = "", col.names = paste0("V", seq_len(max(count.fields("all_mutations_002.txt")))))
run003 <- read.table(file = "all_mutations_003.txt", header = F, fill = T, na.strings = "", col.names = paste0("V", seq_len(max(count.fields("all_mutations_003.txt")))))
run004 <- read.table(file = "all_mutations_004.txt", header = F, fill = T, na.strings = "", col.names = paste0("V", seq_len(max(count.fields("all_mutations_004.txt")))))
run005 <- read.table(file = "all_mutations_005.txt", header = F, fill = T, na.strings = "", col.names = paste0("V", seq_len(max(count.fields("all_mutations_005.txt")))))
run006 <- read.table(file = "all_mutations_006.txt", header = F, fill = T, na.strings = "", col.names = paste0("V", seq_len(max(count.fields("all_mutations_006.txt")))))
run007 <- read.table(file = "all_mutations_007.txt", header = F, fill = T, na.strings = "", col.names = paste0("V", seq_len(max(count.fields("all_mutations_007.txt")))))
run008 <- read.table(file = "all_mutations_008.txt", header = F, fill = T, na.strings = "", col.names = paste0("V", seq_len(max(count.fields("all_mutations_008.txt")))))
run009 <- read.table(file = "all_mutations_009.txt", header = F, fill = T, na.strings = "", col.names = paste0("V", seq_len(max(count.fields("all_mutations_009.txt")))))
run010 <- read.table(file = "all_mutations_010.txt", header = F, fill = T, na.strings = "", col.names = paste0("V", seq_len(max(count.fields("all_mutations_010.txt")))))
run011 <- read.table(file = "all_mutations_011.txt", header = F, fill = T, na.strings = "", col.names = paste0("V", seq_len(max(count.fields("all_mutations_011.txt")))))
run012 <- read.table(file = "all_mutations_012.txt", header = F, fill = T, na.strings = "", col.names = paste0("V", seq_len(max(count.fields("all_mutations_012.txt")))))
run013 <- read.table(file = "all_mutations_013.txt", header = F, fill = T, na.strings = "", col.names = paste0("V", seq_len(max(count.fields("all_mutations_013.txt")))))
run014 <- read.table(file = "all_mutations_014.txt", header = F, fill = T, na.strings = "", col.names = paste0("V", seq_len(max(count.fields("all_mutations_014.txt")))))
run015 <- read.table(file = "all_mutations_015.txt", header = F, fill = T, na.strings = "", col.names = paste0("V", seq_len(max(count.fields("all_mutations_015.txt")))))
run016 <- read.table(file = "all_mutations_016.txt", header = F, fill = T, na.strings = "", col.names = paste0("V", seq_len(max(count.fields("all_mutations_016.txt")))))
run017 <- read.table(file = "all_mutations_017.txt", header = F, fill = T, na.strings = "", col.names = paste0("V", seq_len(max(count.fields("all_mutations_017.txt")))))

# combine SNP data.
runALL <- rbind.fill(run001,run002,run003,run004,run005,run006,run007,run008,run009,run010,run011,run012,run013,run014,run015,run016,run017)

# set column name.
colnames(runALL)[1] <- "Sample_No"

# edit repeat sample naming.
runALL$Sample_No <- gsub("Repeat-", "", runALL$Sample_No)
runALL$Sample_No <- gsub("Repeat_", "", runALL$Sample_No)
runALL$Sample_No <- gsub("_2$", "", runALL$Sample_No)

# get reportable sample names.
samples <- read_xlsx("sample_names_reportable.xlsx", col_types = c("numeric","text","text","date","date"))

# retain reportable samples.
report <- runALL[runALL$Sample_No %in% samples$Sample_No,]

# sort by Sample_No.
report <- report[order(report$Sample_No),]

# reset row numbering.
row.names(report) <- NULL

# # identify replicates.
# report$Sample_No[which(duplicated(report$Sample_No))]
# which(duplicated(report$Sample_No))

# remove replicates.
report <- report[-c(155,177,180,197,216,218,224,227,231,234,235,237,239,242,243,245,247,249,263,265,278,294), ]

# set column selection variables.
x <- c("V3","V5","V7","V9","V11","V13","V15","V17","V19","V21","V23","V25","V27","V29","V31","V33","V35","V37")
y <- c("V2","V4","V6","V8","V10","V12","V14","V16","V18","V20","V22","V24","V26","V28","V30","V32","V34","V36")

# remove counts.
report <- report %>% mutate(across(all_of(x), ~substr(. , 1, 8)))

# multiply % by 100.
report <- report %>% mutate(across(all_of(x), ~(as.numeric(.) *100)))

# round %.
report <- report %>% mutate(across(all_of(x), ~round(. , 1)))

# remove incidental snps (synonymous, STOP, uninterpreted).
report <- report %>% mutate(across(all_of(y), ~gsub("\\^", NA, .)))
report <- report %>% mutate(across(all_of(y), ~gsub("STOP", NA, .)))
report <- report %>% mutate(across(all_of(y), ~gsub(">", NA, .)))

# combine snp names and %.
snps <- cbind(report[1], mapply(paste0, report[, seq(2, 37, 2)], sep = "|", report[, seq(3, 37, 2)]))

# get characteristic snps.
char <- snps %>% mutate(across(all_of(y), ~gsub("\\*", NA, .)))
char <- char %>% mutate(across(all_of(y), ~gsub("NA", NA, .)))

# get characteristic snp names.
charn <- char %>% mutate(across(all_of(y), ~gsub("\\|.*", "", .)))

# get non-characteristic snps.
nonc <- snps %>% mutate(across(all_of(y), ~gsub("[A-Z]\\|", NA, .)))
nonc <- nonc %>% mutate(across(all_of(y), ~gsub("[0-9]\\|", NA, .)))

# get non-characteristic snp names.
noncn <- nonc %>% mutate(across(all_of(y), ~gsub("\\*\\|.*", "", .)))

# get non-characteristic snp %.
noncp <- nonc %>% mutate(across(all_of(y), ~gsub(".*\\|", "", .)))

# combine characteristic snp names.
charn <- cbind(charn[,1], (unite(charn[,2:19], `SNPs detected`, sep = ", ", na.rm = TRUE)))
colnames(charn)[1] <- "Sample_No"

# combine data.
table1 <- merge(samples, charn, by="Sample_No", all.y = FALSE)

# combine non-characteristic % names.
noncn <- cbind(noncn[,1], (unite(noncn[,2:19], `Non-characteristic SNPs`, sep = ", ", na.rm = TRUE)))
colnames(noncn)[1] <- "Sample_No"

# combine data.
table1 <- merge(table1, noncn, by="Sample_No", all.y = FALSE)

# combine non-characteristic snp %.
noncp <- cbind(noncp[,1], (unite(noncp[,2:19], `Non-characteristic SNPs (%)`, sep = ", ", na.rm = TRUE)))
colnames(noncp)[1] <- "Sample_No"

# combine data.
table1 <- merge(table1, noncp, by="Sample_No", all.y = FALSE)

# combine data - add metadata for samples that failed sequencing.
`%ni%` <- Negate(`%in%`)
table2 <- rbind.fill(table1, (samples[samples$Sample_No %ni% table1$Sample_No,]))

# sort by Sample_No.
table2 <- table2[order(table2$Sample_No),]

# reset row numbering.
row.names(table2) <- NULL

# output SNP data to Excel file.
write_xlsx(table2, path = "supplementary_SNP_data.xlsx")
