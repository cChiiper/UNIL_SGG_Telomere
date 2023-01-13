################################################################################
### Check rsids in GWAS catalog                                              ###
### Author: Samuel Moix                                                      ###
### Date: 14.09.2022                                                         ###
################################################################################


################################################
### Libraries ##################################
library(dplyr)
library(data.table)

################################################
### Working directories ########################
data_folder = "../"
path_to_SNPs <- "../"

################################################
### Load data ##################################

### Gwas catalog data
GWAS_cat <- as.data.frame(fread(file = file.path(data_folder,  
                                        "gwas_catalog_v1.0-associations_e107_r2022-10-08.tsv"), 
                       sep = '\t', header = TRUE, quote = ""))

### Telomere IVs
TLIVs <- as.data.frame(fread(file = file.path(data_folder,  
                                                 "TELOMERE_IVs.txt"), 
                                sep = '\t', header = TRUE))

### Loci to remove
HLA_rsids <- fread(file.path(path_to_SNPs,  "HLA_rsids.txt"), header = F); colnames(HLA_rsids) <- "SNP"
HBB_rsids <- fread(file.path(path_to_SNPs,  "HBB_rsids.txt"), header = F); colnames(HBB_rsids) <- "SNP"

################################################
### Filter data ################################

GWAS_cat <- select(GWAS_cat, c("SNPS","P-VALUE","DISEASE/TRAIT", "MAPPED_GENE"))
GWAS_cat$`P-VALUE` <- as.numeric(GWAS_cat$`P-VALUE`)
### Keep TL IVs SNPs found in the GWAS catalog
GWAS_cat <- filter(GWAS_cat, SNPS %in% TLIVs$SNP)

print(paste("There are",length(TLIVs$SNP[which(!TLIVs$SNP %in% GWAS_cat$SNPS)]),
            "SNPs that were not found in the GWAS catalog"))

####

TL_snps <- GWAS_cat %>% 
  filter(grepl("Telo|telo",`DISEASE/TRAIT`)) %>%
  distinct(SNPS) %>%
  pull(SNPS)

TL_snps_df <- GWAS_cat %>%
  filter(SNPS %in% TL_snps)

col_names <- c("SNP","Trait","NB_trait","Gene")
filtered_df <- data.frame(matrix(ncol = length(col_names), nrow = length(TL_snps)))
colnames(filtered_df) <- col_names
i <- 1
for(snp in TL_snps){
  temp_df <- TL_snps_df[which(TL_snps_df$SNPS == snp),]
  temp <- temp_df[,"DISEASE/TRAIT"]
  trait <- paste(unique(temp), collapse = '; ')
  filtered_df[i,] <- c(snp, trait, length(unique(temp)), unique(temp_df[,"MAPPED_GENE"]))
  i <- i + 1
}
rm(i)
GWAS_cat$MAPPED_GENE

smallest_P <- GWAS_cat %>%
  group_by(SNPS) %>%
  top_n(1, `P-VALUE`) 


################################################################################
### Sensitivity analysis on restrained SNPs                                  ###
### Author: Samuel Moix                                                      ###
### Date: 14.09.2022                                                         ###
################################################################################

################################################
### Libraries ##################################
library(dplyr)
library(data.table)
library(phenoscanner)

################################################
### Working directories ########################
data_folder = ".."
path_to_SNPs <- "../"

################################################
### Load data ##################################

### Telomere IVs
TLIVs <- as.data.frame(fread(file = file.path(data_folder,  
                                              "TELOMERE_IVs.txt"), 
                             sep = '\t', header = TRUE))

### Loci to remove
HLA_rsids <- fread(file.path(path_to_SNPs,  "HLA_rsids.txt"), header = F); colnames(HLA_rsids) <- "SNP"
HBB_rsids <- fread(file.path(path_to_SNPs,  "HBB_rsids.txt"), header = F); colnames(HBB_rsids) <- "SNP"
TLIVs <- TLIVs$SNP[which(!TLIVs$SNP %in% c(HLA_rsids$SNP, HBB_rsids$SNP))]

### Phenoscanner query
# Divide in multiple queries as phenoscanner only allows 100 queries at the time
queries <- list(TLIVs[1:100],TLIVs[101:200],TLIVs[201:300],
                                                  TLIVs[301:length(TLIVs)])
# Create result dataframe
df_res_IV_trait = data.frame(matrix(nrow = 0, ncol = 2))
colnames(df_res_IV_trait) <- c("rsid","trait")

for (query in queries) {
  # Query
  res <- phenoscanner::phenoscanner(snpquery=query)
  snp_df <- res$results[,c("rsid","trait")]
  # Keeping distinct rows
  snp_df <- dplyr::distinct(snp_df)
  # Append query results
  df_res_IV_trait <- rbind(df_res_IV_trait, snp_df)
}

fwrite(df_res_IV_trait, file.path(path_to_SNPs,  "TL_IVs_phenoscanner.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
df_res_IV_trait <- as.data.frame(fread(file.path(path_to_SNPs,  "TL_IVs_phenoscanner.txt"), header = T))

non_pleio_snp <- names(table(df_res_IV_trait$rsid))[which(table(df_res_IV_trait$rsid) == 1)]
TL_snp_query <- df_res_IV_trait$rsid[which(df_res_IV_trait$trait %in% c("Telomere length",
                                                                        "Leukocyte telomere length",
                                                                        "Leukocyte telomere length in females"))]


### Only non-pleiotropic found SNP == rs8105767



################################################################################
### P-value method to remove "pleiotropic" SNPs                              ###
### Author: Samuel Moix                                                      ###
### Date: 12.12.2022                                                         ###
################################################################################   

################################################
### Working directories ########################
data_folder = "../"
path_to_SNPs <- "../"

################################################
### Load data ##################################

### Telomere IVs
TLIVs <- as.data.frame(fread(file = file.path(data_folder,  
                                              "TELOMERE_IVs.txt"), 
                             sep = '\t', header = TRUE))

### Load p-values (!Changed to tab separated as fread doesn't like sep = ' ')
df_pval <- as.data.frame(fread(file = file.path(data_folder, "TL_IV_snp_pvalues_tab.txt"), 
                                      header = TRUE, sep = "\t"))
  
# Remove TL and put rownames as rsids
row.names(df_pval) <- df_pval$SNP
df_pval <- df_pval %>% select(-c("TELOMERE","SNP"))

# Get minimal pvalue
df_min_pval <- apply(df_pval, 1, FUN = min, na.rm = TRUE)
df_min_pval <- as.data.frame(df_min_pval)
colnames(df_min_pval) <- "pval"
df_filtered <- df_min_pval %>% 
  filter(pval > 0.0005)

df_pval[rownames(df_filtered),]
GWAS_cat_sub <- filter(GWAS_cat, SNPS %in% rownames(df_filtered))
