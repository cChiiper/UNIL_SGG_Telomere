#################################################################################################################
### Add new traits manually to the initial dataframe                                                          ###
### Author: Samuel Moix                                                                                       ###
### Date: 30.06.2022                                                                                          ###
#################################################################################################################

################################################
### Libraries ##################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(matrixStats)

################################################
### Working directories ########################

data_folder = "../"
export_folder = "../"


################################################
### STEP 1: Add traits to add ##################

### SET YOUR TRAITS TO ADD HERE <---------------
FieldID	<- c(22192, 6177, 6153,
             1588, 1578, 1608, 5364,
             1568, 1598,
             46, 47,
             6138, 3872)
PHENO <- c("TL_add_check", "cholest_medication_M", "cholest_medication_F",
           "alcool_01_bc", "alcool_02_cw", "alcool_03_fw", "alcool_04_o",
           "alcool_05_rw", "alcool_06_s",
           "hand_grip_l", "hand_grip_r",
           "qualifications","age_birth_primi")
par_categorical_list <- c(TRUE, TRUE, TRUE,
                          FALSE, FALSE, FALSE, FALSE,
                          FALSE, FALSE,
                          FALSE, FALSE,
                          TRUE, TRUE)

# Dataframe of field ID and short name for the corresponding phenotype
pheno_list <- as.data.frame(cbind(FieldID, PHENO))

#################################################
### STEP 2: Order raw phenotype files ###########

# As phenotypic data is spread throughout multiple files downloaded from the UKBB portal, we start by listing all phenotype files available and order them by date
raw_pheno_file <- file.info(list.files("../", pattern = "^ukb.*csv", full.names = T, recursive = F))
raw_pheno_file <- data.frame(File = rownames(raw_pheno_file), Date = raw_pheno_file[,4])
raw_pheno_file <- raw_pheno_file[order(raw_pheno_file$Date, decreasing = T), ]
print(paste0("There are ", nrow(raw_pheno_file), " raw phenotype files"))

#################################################
### STEP 3: Identify most recent location #######
# As mentioned in STEP 2, there are multiple phenotype files, with some phenotypes being in several files
# For the purpose of the analysis, we use phenotype information originating from the most recent phenotype file

pheno_list$File <- NA
pheno_list$File <- as.character(pheno_list$File) 
pheno_list$Col <- NA
pheno_list$Col <- as.character(pheno_list$Col) 

# Loop over all phenotypes to detect the most recent location
for (p in 1:nrow(pheno_list)) {
  # Define the phenotype
  pheno <- pheno_list[p, "PHENO"] 
  ID <- pheno_list[p, "FieldID"] 
  print(paste0("Extracting most recent file for ", pheno))
  counter <- 1
  
  # Determine the most recent location --> loop through raw files
  while (is.na(pheno_list[p, "File"])) {
    
    header <- fread(as.character(raw_pheno_file[counter, "File"]), nrow = 0)
    col <- grep(paste0("^", ID, "-"), names(header))
    if (length(col) > 0) {
      pheno_list[p, "File"] <- as.character(raw_pheno_file[counter, "File"])
      print(paste0("Most recent location: ", sub(".*/", "", pheno_list[p, "File"])))
      pheno_list[p, "Col"] <- paste(col, collapse = "_")}
    counter <- counter +1
  }
} 
rm(p, pheno, ID, counter, header, col)

#################################################
### STEP 4: Extract phenotypes             ######

# This file contains the full dataframe
phenotypes <- as.data.frame(fread(file.path(data_folder, "complete_DF.txt") , header = T))

# Loop over all phenotypes to extract the data and calculate the average over different measurement instances
for(p in 1:nrow(pheno_list)) {
  # Whether categorical or not
  par_categorical <- par_categorical_list[p]
  
  # Define the phenotype
  pheno <- as.character(pheno_list[p, "PHENO"]) 
  ID <- pheno_list[p, "FieldID"] 
  file <- pheno_list[p, "File"] 
  col <- pheno_list[p, "Col"] 
  print(paste0("Extracting ", pheno, " from ", sub(".*/", "", file)))
  
  # Extract columns corresponding to the defined phenotype
  temp_pheno <- as.data.frame(fread(file, header = T, select = c(1, as.numeric(str_split(col, pattern = "_")[[1]]))))
  
  # Average (if more than one column present)
  if(par_categorical){
    if (ncol(temp_pheno) > 2) {
      temp_pheno <- data.frame(eid = temp_pheno[, 1], pheno = temp_pheno[, 2])}
  } else{
    if (ncol(temp_pheno) > 2) {
      # Set traits that can have negative values
      traits_with_negatives <- c("TDI", "TL") 
      if(!pheno %in% traits_with_negatives){
        for(i in 2:length(temp_pheno)){
          temp_pheno[,i][which(temp_pheno[,i] < 0)] <- NA
        }
      }
      temp_pheno <- data.frame(eid = temp_pheno[, 1], pheno = rowMeans(temp_pheno[, 2:ncol(temp_pheno)], na.rm = T))}
  }
  
  # Merge
  colnames(temp_pheno) <- c("eid", pheno)
  phenotypes <- merge(phenotypes, temp_pheno, by="eid", all.x = TRUE)
  
}
rm(pheno, ID, file, col, temp_pheno)
print(paste0("Dimensions of phenotype table: ", ncol(phenotypes), " pheno x ", nrow(phenotypes), " eids"))

# Drop eids as sensitive data
phenotypes <- select(phenotypes, -c("eid"))

# Save new dataframe
fwrite(phenotypes, file.path(export_folder,  "complete_DF_ladd_noeid.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

