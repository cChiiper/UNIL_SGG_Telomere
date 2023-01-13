################################################################################
### Harmonize clumped SNPs                                                   ###
### Author: Samuel Moix                                                      ###
### Date: 26.10.2022                                                         ###
################################################################################

################################################
### Libraries ##################################
library(TwoSampleMR)
library(dplyr)
library(gridExtra)
library(data.table)

################################################
### Set traits of interest #####################
#c("APOB","BMI","BFM","CHOLEST","LDL","MDD","TRI","EDUAGE","AGELASTB")
EXPOSURES <- c("APOB","TRI","BMI","AGELASTB")
OUTCOME <- c("TELOMERE")

################################################
### Working directories ########################
data_folder <- "../"
export_folder <- "../"
file_path_MR <- "../TRAIT_on_TELOMERE_mr_results.txt"
path_to_SNPs <- ".."

################################################ mv_harmonized_data_ALB_AEE_BMI.txt
### Load data ################################## mv_harmonized_data_TL_lipid_ses.txt
dat <- as.data.frame(fread(file.path(data_folder, "mv_harmonized_data_TL_lipid_ses.txt"), header = T))
print(paste0("Number of SNPs: ", nrow(dat)))

### Remove SNPs from HLA locus (chr6 25MB-37MB)
HLA_rsids <- fread(file.path(path_to_SNPs,
                             "HLA_rsids.txt"), header = F)
dat <- filter(dat, !SNP %in% HLA_rsids$V1)
print(paste0("Number of SNPs without HLA locus: ", nrow(dat)))

### Remove specific SNPS
#snps_remove <- c("rs6857","rs142042446","rs11642015","rs56094641")
#dat <- filter(dat, !SNP %in% snps_remove)

### Remove SNPs from not-wanted exposures
list_exp_snp <- c()
for (exp_i in EXPOSURES) {
  # Search file
  file_name <- paste0(exp_i, "_to_", OUTCOME, "_MR_data.tsv") 
  # Adapt file name if in sub-directory
  file_name <- list.files(path = path_to_SNPs, pattern = file_name, recursive = TRUE)
  # Read file
  print(file_name[1])
  exp_clumped_snps <- fread(file.path(paste0(path_to_SNPs,"/", file_name[1])), header = T)
  # Add SNPs to list
  list_exp_snp <- c(list_exp_snp,exp_clumped_snps$SNP)
}
# Filter SNPs 
dat <- filter(dat, SNP %in% list_exp_snp)
rm(exp_clumped_snps)
rm(list_exp_snp)

### Allele frequency filter 
for (exp_i in EXPOSURES) {
  dat <- dat[abs(dat[, paste0("Freq_",exp_i)] - dat$Freq_OUT) < 0.05,]
}

print(paste0("Number of SNPs after frequency check: ", nrow(dat)))

### Steiger filter: remove SNPs with larger outcome than exposure effects
for (exp_i in EXPOSURES) {
  dat$zval_steiger <- (abs(dat[, paste0("b_",exp_i)])-abs(dat$b_OUT))/sqrt(dat[, paste0("se_",exp_i)]**2 + dat$se_OUT**2)
  dat <- dat[dat$zval_steiger > -1.96, ]
}

print(paste0("Number of SNPs after steiger filter: ", nrow(dat)))


################################################
### Format into TwoSampleMR ####################

### Read outcome data
dat$Phenotype <- OUTCOME
outcome_dat  <- format_data(
  dat,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "b_OUT",
  se_col = "se_OUT",
  eaf_col = "Freq_OUT",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "p_OUT",
  samplesize_col = "N_OUT",
)

### Read exposure data
exposures_df <- list()
index <- 1
for (TRAIT in EXPOSURES) {
  dat$Phenotype <- TRAIT
  exposures_df[[index]] <- format_data(
    dat,
    type = "exposure",
    snp_col = "SNP",
    beta_col = paste0("b_",TRAIT),
    se_col = paste0("se_",TRAIT),
    eaf_col = paste0("Freq_",TRAIT),
    effect_allele_col = "A1",
    other_allele_col = "A2",
    pval_col = paste0("p_",TRAIT),
    samplesize_col = paste0("N_",TRAIT),
  )
  index <- index + 1
}
rm(index)

exposure_dat <- dplyr::bind_rows(exposures_df)
rm(exposures_df)

### Re-harmonize data

mvdat <- mv_harmonise_data(exposure_dat, outcome_dat, harmonise_strictness = 2)



################################################
### MVMR #######################################

mv_basic(mvdat, pval_threshold = 5e-08)
mv_ivw_results <-mv_ivw(mvdat, pval_threshold = 5e-08)

mv_plots <- gridExtra::grid.arrange(grobs = c(mv_ivw_results$plots[1],
                                              mv_ivw_results$plots[2],
                                              mv_ivw_results$plots[3]),
                                    ncol=3)

mvmr_res <- mv_multiple(mvdat)$result
print(mvmr_res)

################################################
### Compare with ivw result ####################

df_res_MR <- read.table(file_path_MR, sep = '\t', header = TRUE)
df_res_MR <- df_res_MR %>% 
  filter(exposure %in% EXPOSURES) %>%
  filter(method == "Inverse variance weighted")

for (trait_compare in EXPOSURES) {
  b.x <- df_res_MR$b[which(df_res_MR$exposure == trait_compare)]
  se.x <- df_res_MR$se[which(df_res_MR$exposure == trait_compare)]
  CI_sup.x <- b.x + 1.96 * se.x
  CI_inf.x <- b.x - 1.96 * se.x
  
  b.y <- mvmr_res$b[which(mvmr_res$exposure == trait_compare)]
  se.y <- mvmr_res$se[which(mvmr_res$exposure == trait_compare)]
  pval.y <- mvmr_res$pval[which(mvmr_res$exposure == trait_compare)]
  CI_sup.y <- b.y + 1.96 * se.y
  CI_inf.y <- b.y - 1.96 * se.y
  
  # Compute pdiff
  pdiff <- 2*pnorm(-abs((b.x-b.y)/sqrt(se.x**2+se.y**2)), mean = 0, sd = 1)
  #Format p-value for printing
  pdiffformat <- unlist(base::strsplit(as.character(pdiff),"e"))
  pvalformat <- unlist(base::strsplit(as.character(pval.y),"e"))
  #Print
  print(trait_compare)
  print(paste("From",round(b.x,4),"se",round(se.x,4),#paste0(" CI (",round(CI_inf.x,4)," - ",round(CI_sup.x,4),")"),
              "to", round(b.y,4),"se",round(se.y,4), #paste0(" CI (",round(CI_inf.y,4)," - ",round(CI_sup.y,4),")"),
              "with pval:",ifelse(length(pvalformat) == 1, 
                                  as.character(round(as.numeric(pvalformat),4)),
                                  paste0(round(as.numeric(pvalformat[1]),1)," x 10",pvalformat[2])),
              "and with pdiff:", ifelse(length(pdiffformat) == 1, 
                                    as.character(round(as.numeric(pdiffformat),4)),
                                    paste0(round(as.numeric(pdiffformat[1]),1)," x 10",pdiffformat[2]))))
}



MVMR_merge <- mvmr_res[,c("exposure","b","se","pval")]
MVMR_merge$method <- "MVMR"
df_res_MR$method <- "IVW"
plot_dataframe <- rbind(df_res_MR[,c("exposure","b","se","pval","method")],MVMR_merge)

require(ggforestplot)
ggforestplot::forestplot(
  df = plot_dataframe,
  name = exposure,
  estimate = b,
  pvalue = pval,
  psignif = 0.05,
  xlab = expression(""*alpha*" (CI95%)"),
  title = "",
  colour = method,  se = se
)
