################################################################################
### Observation against MR results for other traits                          ###
### Author: Samuel Moix                                                      ###
### Date: 02.08.2022                                                         ###
################################################################################


################################################
### Libraries ##################################
library(dplyr)
library(ggplot2)
require(ggrepel)

################################################
### Set trait of interest ######################
st_TRAIT <- "HDL"

################################################
### Working directories ########################

data_folder = "../"
data_folder_obs = "../"
export_folder = "../"
MR_results_file <- paste0("../","/TRAIT_on_",st_TRAIT,"_mr_results.txt")

################################################
### Parameters  ################################
indep_trait <- 135

################################################
### Load Data and prepare data #################

### Load FieldID ref
df_code <- read.csv(file = file.path(data_folder,  "Phenotypes_Selection_Table/TL_metadata.csv"), sep = ',', header = TRUE)
df_code <- rename(df_code, exposure = "MR_name")

################################################
### Load MR results ############################

### Load MR results file
df_res_MR <- read.table(file = file.path(data_folder,  MR_results_file), sep = '\t', header = TRUE)


### Select method
df_res_MR <- filter(df_res_MR, method %in% c("Inverse variance weighted"))

### Add CI to MR data frame
df_res_MR <- df_res_MR %>% 
  mutate(Beta_CI_inf_MR = b - 1.96*se) %>%
  mutate(Beta_CI_sup_MR = b + 1.96*se)

### Select columns of interest
df_res_MR <- df_res_MR[,c("exposure","b","Beta_CI_inf_MR","Beta_CI_sup_MR","pval")]

### Add "pheno" names
df_res_MR <- merge(df_res_MR, df_code[,c("exposure","pheno")],by="exposure", all.x = TRUE)

################################################
### Calculating regression coef ################

#... RUN code from 03_regression_analysis_TL_corrected.R
df <- readRDS(file.path(data_folder_obs, "df_scale_cor.rds"))

### Regression coefficients on df
df_trait <- df_code$pheno[which(df_code$exposure == st_TRAIT)]
df <- select(df, c(df_res_MR$pheno, all_of(df_trait))) ### Select traits of interest

x <- c("pheno", "mean", "sd", "rcoef", "CI1", "CI2", "p_value", "se")

df_res <- data.frame(matrix(ncol = length(x), nrow = 0))
colnames(df_res) <- x

### Calculate p-value function
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

for(i in 1:ncol(df)) {       # for-loop over columns
  reg_TRAIT_pheno <- lm(df[,df_trait] ~ df[, i])
  df_res[nrow(df_res) + 1,] <- c(colnames(df)[i], 
                                 if(class(df[,i]) == 'integer' | class(df[,i]) == 'numeric'){
                                   mean(df[, i], na.rm = T)
                                 } else NA,
                                 if(class(df[,i]) == 'integer' | class(df[,i]) == 'numeric'){
                                   sd(df[, i], na.rm = T)
                                 } else NA,
                                 reg_TRAIT_pheno$coefficients[[2]],
                                 confint(reg_TRAIT_pheno)[2,1],
                                 confint(reg_TRAIT_pheno)[2,2],
                                 lmp(reg_TRAIT_pheno),
                                 coef(summary(reg_TRAIT_pheno))[2, 2])
  print(colnames(df)[i])
  print(reg_TRAIT_pheno)
}
rm(reg_TRAIT_pheno)


# Put data in df_res as numeric and factor
df_res$pheno <- as.factor(df_res$pheno)
df_res[,2:ncol(df_res)] <- lapply(df_res[,2:ncol(df_res)], as.numeric)

# Arrange in descending order of rcoef
df_res <- arrange(df_res, desc(rcoef))


################################################
### Merge and plot data ########################

### Merge observational and MR data
df_both <- merge(df_res_MR, df_res[,c("pheno","rcoef","CI1","CI2","p_value")],by="pheno", all.x = TRUE)

### Merge forward data and association
df_res_MR$type_an <- "a_r"
### Rename observation columns to match MR names for plotting
df_res_OBS <- df_both %>% select("rcoef","CI1","CI2","p_value","pheno") %>%
  rename(b = rcoef, Beta_CI_inf_MR = CI1, Beta_CI_sup_MR = CI2, pval = p_value)
df_res_OBS$type_an <- "b"

test <- rbind(select(df_res_MR, -c("exposure")), df_res_OBS)

test %>%
  mutate(pheno = forcats::fct_reorder(pheno, desc(b))) %>%
  ggplot(aes(pheno, b, colour = type_an)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(position = position_dodge(width=0.75), 
             shape = ifelse(test$pval < 0.05/indep_trait, 16,1)) +
  geom_errorbar(aes(ymax = Beta_CI_sup_MR, ymin = Beta_CI_inf_MR, width = .1, color = type_an), 
                position = position_dodge(width=0.75),
                alpha = 0.5) +
  xlab("") +
  coord_flip() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))


obs_vs_mr_plot <- ggplot(df_both, aes(x=rcoef, y=b)) +
  geom_errorbar(aes(ymax = Beta_CI_sup_MR, ymin = Beta_CI_inf_MR, width = 0), color = "red", alpha = 0.1) +
  geom_errorbar(aes(xmax = CI2, xmin = CI1, width = 0), color = "blue", alpha = 0.1) +
  geom_point(shape=1) + # Show dots
  geom_abline(slope = 1, color = "grey", alpha = 0.5) +
  ggrepel::geom_text_repel(
    label=df_both$pheno, 
    size = 2
  ) +
  xlab("rcoef") +
  ylab("MR beta (ivw)")


################################################################################
################### Look for confounder by product #############################
################################################################################

### Libraries
library(tidyr)
library(forcats)


### Data folder
data_folder <- "../"
metadata_folder <- "../"
metadata <- read.csv(file.path(metadata_folder, "TL_metadata.csv"))

### Load data
# Choose studied traits
st_TRAIT <- "TRI"
# Traits with effect on TELOMERE
df_res_MR_tl <- read.table(file = file.path(data_folder,  "TRAIT_on_TELOMERE_mr_results.txt"), sep = '\t', header = TRUE)
df_res_MR_tl <- filter(df_res_MR_tl, exposure %in% c("MDD","ALCOWF","HYPERTENSION","BFM","BMI",
                                                     "LYMPH","EOSINO","URATE","WBC",
                                                     "CREAC","MCH","HEIGHT","IGF1","APOB","RBC",
                                                     "CHOLEST","LDL","LIPIDDISORDER",
                                                     "EDUAGE","LSMOKE","AGEFIRSTB","AGELASTB",
                                                     "TRI"))
# Traits with effect on studied trait
df_res_MR_st <- read.table(file = file.path(data_folder,  
                                            paste0("../TRAIT_on_",st_TRAIT,"_mr_results.txt")), 
                           sep = '\t', header = TRUE)

df_res_MR <- rbind(df_res_MR_tl,df_res_MR_st); rm(df_res_MR_st); rm(df_res_MR_tl)
# Select "ivw"
df_res_MR <- filter(df_res_MR, method == "Inverse variance weighted")
# Remove traits with missing effects
traits_mis <- names(which(table(df_res_MR$exposure) == 1))
df_res_MR <- df_res_MR[which(!df_res_MR$exposure %in% traits_mis),]
# Change labels 
names_coding <- select(metadata, c("pheno","MR_name"))
colnames(names_coding) <- c("pheno", "exposure")
df_res_MR <- merge(df_res_MR,names_coding,by="exposure", all.x = TRUE)


### Summarize data to product of traits
# and get CI of the product


### Compute estimate product and standard error according to
# https://stats.stackexchange.com/questions/498224/standard-error-of-estimated-sum-or-product-mean 
df_res_prod <- df_res_MR %>%
  mutate(direct = ifelse(outcome == st_TRAIT, "TR", "TL")) %>%
  select(c("pheno", "b", "se", "direct")) %>%
  tidyr::pivot_wider(names_from = direct, values_from = c("b","se")) %>%
  mutate(prod_b = b_TR*b_TL) %>%
  mutate(se = sqrt(((b_TR**2)*(se_TL**2))+((b_TL**2)*(se_TR**2))+ ((se_TR**2)*(se_TL**2)))) %>%
  arrange(desc(abs(prod_b)))


# ### Plot
# df_res_prod %>%
#   mutate(pheno = forcats::fct_reorder(pheno, (prod_b))) %>%
#   ggplot(aes(x = pheno,y = prod_b)) +
#   geom_bar(stat="identity", fill = "#95BF74") +
#   geom_errorbar(
#     aes(x=pheno, 
#         ymin = prod_b - 1.96*se, # 
#         ymax = prod_b + 1.96*se), 
#     color = "red"
#   ) +
#   #geom_text(aes(label = pheno), nudge_y = 0.004, size = 2) +
#   ylab("\u03b1 (ivw) product (CI95%)") +
#   xlab(paste("Product from",st_TRAIT)) +
#   coord_flip() + 
#   theme(
#     #axis.text.y = element_blank(),
#     axis.ticks = element_blank())


### Add p-value (double check that formula)
df_res_prod <- df_res_prod %>%
  mutate(z = abs(prod_b/se)) %>%
  mutate(pval = 2*pnorm(z, mean = 0, sd = 1, lower.tail = FALSE)) # or -abs(z)

### Put descriptive title
df_res_prod <- merge(df_res_prod, metadata[,c("pheno","description")], by="pheno")
df_res_prod <- df_res_prod %>%
  select(-pheno) %>%
  rename(pheno = description)

### Reorder by beta
df_order <- df_res_prod %>%
  arrange(prod_b)

df_res_prod <- df_res_prod %>%
  group_by(pheno) %>%  
  arrange(factor(pheno, levels = df_order$pheno)) %>%
  ungroup()

### Plot
ggforestplot::forestplot(
  df = df_res_prod,
  name = pheno,
  estimate = prod_b,
  pvalue = pval,
  logodds = FALSE,
  psignif = 0.05,
  xlab = expression(paste(""*alpha*""["product"]," (95% CI)")),
  title = "",
  se = se
) +
  ggplot2::coord_cartesian(xlim = c(-0.06,0.02)) 

#TRI c(-0.06,0.02)
#HDL c(-0.02,0.1)
#WEIGHT c(-0.15,0.01)
  


