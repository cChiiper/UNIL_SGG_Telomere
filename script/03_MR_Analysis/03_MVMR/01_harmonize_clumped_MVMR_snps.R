################################################################################
### Harmonize clumped SNPs                                                   ###
### Author: Samuel Moix                                                      ###
### Date: 26.10.2022                                                         ###
################################################################################

################################################
### Libraries ##################################
library(dplyr)

################################################
### Set traits of interest #####################
EXPOSURES <- c("APOB","BMI","BFM","CHOLEST","LDL","MDD","TRI","EDUAGE","AGELASTB")
OUTCOME <- c("TELOMERE")

################################################
### Set directories ############################
summary_stats_folder <- "../"
snps_folder <- "../"
output_folder <- "../"

################################################
### Merge clumped SNPs files ################### 

### Get file location (in sub-folder)
all_files <- c()
for (EXPOSURE in EXPOSURES) {
  file_name <- paste0(EXPOSURE, "_to_", OUTCOME, "_clumped_snps.txt") 
  file_name <- list.files(path = snps_folder, pattern = file_name, recursive = TRUE)
  all_files <- c(all_files, file.path(snps_folder, file_name[1]))
  print(paste("Found",length(file_name),"file(s) for clumped SNPs of", EXPOSURE))
  if(length(file_name) == 0){
    stop(paste("File for",EXPOSURE,"missing"))
  }
}

### Load snps data
snp_list <- lapply(all_files,read.table) # Read snp data
snps <- dplyr::bind_rows(snp_list) # Merge snp data
snps <- unique(snps) # Take unique SNPs
names(snps) <- "SNP"
print(paste("Total number of SNPs:", nrow(snps)))
rm(snp_list)
 
################################################
### Prepare MVMR data ########################## 
 
### Load outcome data
file_name <- paste0(OUTCOME,"_gwas_summary_uk10kck.ma")
outcome_df <- read.table(file.path(summary_stats_folder, file_name), header = T)
print("Outcome df is in memory..")

counter_exposure <- length(EXPOSURES)
for (EXPOSURE in EXPOSURES) {
  file_name <- paste0(EXPOSURE,"_gwas_summary_uk10kck.ma")
  expo_df <- read.table(file.path(summary_stats_folder, file_name), header = T)
  print(paste(EXPOSURE, "as exposure df is in memory."))
  
  expo_df <- expo_df[expo_df$SNP %in% snps$SNP,]
  print("Exposure df is subsetted to clumped SNPs.")
  
  if(counter_exposure == length(EXPOSURES)){
    df <- merge(outcome_df, expo_df, by = 'SNP', suffixes = c('_OUT', '_EXPO'))
    rm(outcome_df)
    first_exposure <- FALSE
  }else{
    df <- merge(df, expo_df, by = 'SNP', suffixes = c(paste0('_',prev_EXP), '_EXPO'))
    df <- df %>%
      rename(A1_EXPO = A1) %>%
      rename(A2_EXPO = A2)
  }
  print(paste0("Number of SNPs before allele harmonizing: ", nrow(df)))
  
  
  ## harmonize alleles
  # Change columns data type 
  df$A1_EXPO <- as.character(df$A1_EXPO)
  df$A2_EXPO <- as.character(df$A2_EXPO)
  df$A1_OUT <- as.character(df$A1_OUT)
  df$A2_OUT <- as.character(df$A2_OUT)
  df$Freq_OUT <- as.numeric(df$Freq_OUT)
  df$b_OUT <- as.numeric(df$b_OUT)
  
  # Flip flipped alleles
  mask_flip <- ((df$A1_EXPO == df$A2_OUT) & (df$A2_EXPO == df$A1_OUT))
  df[mask_flip, "A1_OUT"] <- df[mask_flip, "A1_EXPO"]
  df[mask_flip, "A2_OUT"] <- df[mask_flip, "A2_EXPO"]
  df[mask_flip, "Freq_OUT"] <- 1-df[mask_flip, "Freq_OUT"] 
  df[mask_flip, "b_OUT"] <- -df[mask_flip, "b_OUT"]
  
  ## check if remaining alleles are correct
  mask_correct <- ((as.character(df$A1_EXPO) == as.character(df$A1_OUT)) & (as.character(df$A2_EXPO) == as.character(df$A2_OUT)))
  df <- df[mask_correct,]
  
  print(paste0("Number of SNPs of allele harmonizing: ", nrow(df)))
  df <- select(df, -c("A1_EXPO","A2_EXPO"))
  if(counter_exposure > 1){
    colnames(df) = gsub("_EXPO", "", colnames(df))
    prev_EXP <- EXPOSURE
    counter_exposure <- counter_exposure-1
  }else{
    colnames(df) = gsub("_EXPO", paste0("_",EXPOSURE), colnames(df))
  }
}

### Rename allele columns
df <- df %>% 
  rename(A1 = A1_OUT) %>%
  rename(A2 = A2_OUT) 

################################################
### Export results ############################# 

write.table(df, file.path(output_folder, "mv_harmonized_data_TL.txt"), row.names = F, quote = F)






