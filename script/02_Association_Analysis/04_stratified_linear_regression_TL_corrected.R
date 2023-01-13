#################################################################################################################
### Stratified regression coefficients with multiple phenotypes                                                          ###
### Author: Samuel Moix                                                                                       ###
### Date: 04.04.2022                                                                                          ###
#################################################################################################################


################################################
### Libraries ##################################

library(data.table) # Read data 
library(ggplot2) # Plotting basics
library(dplyr)
library(forcats)
require(ggrepel)
require(tibble)

################################################
### Parameters  ################################

### Scale all data
par_scale <- TRUE
###  Use only continuous traits
par_only_cont <- FALSE
###  Use only factorial traits
par_only_fact <- FALSE
### Set number of independent traits
indep_trait <- 135 
### Whether to remove outliers from data
par_outlier_filter <- TRUE
### Correct traits (cholesterol for medication)
par_corr_traits <- TRUE
par_keep_not_cor <- FALSE # Whether to keep traits not corrected 
### Whether to remove patients with blood cancer
par_rm_blood_cancer <- TRUE

####### LOAD !!!!!!!!!

data_folder = "../"

df <-  as.data.frame(fread(file.path(data_folder,  "filtered_DF.txt"), header = T))

### Remove single-sex traits
df <- select(df, -c(menarche, menopause, BC, OC, endometriosis, menstruation, nb_birth,
                      birth_weight_first_child, age_first_birth, age_last_birth, PC, facial_hair,
                    balding, reproductive_lp))

### Remove blood cancer patients 
if(par_rm_blood_cancer){
  # Remove patients with blood cancer
  df <- df[which(df$blood_cancer == 0 | is.na(df$blood_cancer)),]
  df <- select(df, -c("blood_cancer"))
}

### Correct traits
if(par_corr_traits){
  # Correct cholesterol for cholesterol medication
  # Changes in mmol/L from https://pubmed.ncbi.nlm.nih.gov/22354153/
  # and from Simvastatin https://bmcprimcare.biomedcentral.com/articles/10.1186/1471-2296-4-18 
  
  df <- df %>% add_column(cholesterol_cor = ifelse(df$cholest_medication == 0, df$cholesterol,
                                                   df$cholesterol + 1.64)
                          , .after = "cholesterol" )
  
  df <- df %>% add_column(HDL_cor = ifelse(df$cholest_medication == 0, df$HDL,
                                           df$HDL - 0.1)
                          , .after = "HDL" )
  
  df <- df %>% add_column(LDL_cor = ifelse(df$cholest_medication == 0, df$LDL,
                                           df$LDL + 1.7)
                          , .after = "LDL" )
  
  df <- df %>% add_column(TG_cor = ifelse(df$cholest_medication == 0, df$TG,
                                          df$TG + 0.4) # All doses
                          , .after = "TG" )
  
  if(!par_keep_not_cor){
    df$cholesterol <- df$cholesterol_cor
    df$HDL <- df$HDL_cor
    df$LDL <- df$LDL_cor
    df$TG <- df$TG_cor
    df <- select(df, -c("cholesterol_cor","HDL_cor","LDL_cor","TG_cor"))
  }
}

################################################
### Data type segregation ######################

df_M <- df[which(df$sex == 1), -4] # -4 to remove sex
df_F <- df[which(df$sex == 2), -4]


df_M$TLc <- residuals(lm(df_M$Zadj_TS ~ df_M$age + df_M$age2 + df_M$array))
df_F$TLc <- residuals(lm(df_F$Zadj_TS ~ df_F$age + df_F$age2 + df_F$array))


if (par_only_cont){
  ### Create df only containing continuous traits ###!
  contvec <- c(1:which(colnames(df_M) == "age_cancer"), which(colnames(df_M) == "TLc"))
  df_M <- df_M[,contvec]
  df_F <- df_F[,contvec]
  df_M <- select(df_M, -c(blood_cancer))
  df_F <- select(df_F, -c(blood_cancer))
}

if (par_only_fact){
  ### Create df only containing factorial traits
  contvec <- c(which(colnames(df_M) == "blood_cancer"),
               which(colnames(df_M) == "household_income"):which(colnames(df_M) == "TLc"))
  df_M <- df_M[,contvec]
  df_F <- df_F[,contvec]
}


# Outliers defined as �5SDs
contvec <- c(1:which(colnames(df_M) == "age_cancer"), which(colnames(df_M) == "TLc"))
if(par_outlier_filter){
  for (i in contvec) {
    upper <- mean(df_M[,i], na.rm = T) + sd(df_M[,i], na.rm = T)*5
    lower <- mean(df_M[,i], na.rm = T) - sd(df_M[,i], na.rm = T)*5
    df_M[,i][which(df_M[,i] > upper | df_M[,i] < lower)] <- NA
    rm(upper)
    rm(lower)
  }
}

if(par_outlier_filter){
  for (i in contvec) {
    upper <- mean(df_F[,i], na.rm = T) + sd(df_F[,i], na.rm = T)*5
    lower <- mean(df_F[,i], na.rm = T) - sd(df_F[,i], na.rm = T)*5
    df_F[,i][which(df_F[,i] > upper | df_F[,i] < lower)] <- NA
    rm(upper)
    rm(lower)
  }
}
rm(contvec)
################################################
### Total scaling sex-specific #################

if (par_scale){
  #### To put all traits as numeric and scale 
  df_M[,1:length(df_M)] <- lapply(df_M[,1:length(df_M)], as.numeric) 
  df_M <- as.data.frame(scale(df_M)) 
  df_F[,1:length(df_F)] <- lapply(df_F[,1:length(df_F)], as.numeric) 
  df_F <- as.data.frame(scale(df_F)) 
}

################################################
### Calculating regression coef ################

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


### Regression coefficients on df for males
df <- df_M
x <- c("pheno", "mean", "sd", "rcoef", "CI1", "CI2", "p_value", "se")
df_res <- data.frame(matrix(ncol = length(x), nrow = 0))
colnames(df_res) <- x
for(i in 1:ncol(df)) {       # for-loop over columns
  reg_TLc_pheno <- lm(df$TLc ~ df[, i])
  df_res[nrow(df_res) + 1,] <- c(colnames(df)[i], 
                                 if(class(df[,i]) == 'integer' | class(df[,i]) == 'numeric'){
                                   mean(df[, i], na.rm = T)
                                 } else NA,
                                 if(class(df[,i]) == 'integer' | class(df[,i]) == 'numeric'){
                                   sd(df[, i], na.rm = T)
                                 } else NA,
                                 reg_TLc_pheno$coefficients[[2]],
                                 confint(reg_TLc_pheno)[2,1],
                                 confint(reg_TLc_pheno)[2,2],
                                 lmp(reg_TLc_pheno),
                                 coef(summary(reg_TLc_pheno))[2, 2])
  print(colnames(df)[i])
  print(reg_TLc_pheno)
}
# Put data in df_res as numeric and factor
df_res$pheno <- as.factor(df_res$pheno)
df_res[,2:ncol(df_res)] <- lapply(df_res[,2:ncol(df_res)], as.numeric)
# Remove TL data (too correlated)
df_res <- filter(df_res, !pheno %in% c("TLc", "Zadj_TS"))
# Arrange in descending order of rcoef
df_res_M <- arrange(df_res, desc(rcoef))


### Regression coefficients on df for females
df <- df_F
x <- c("pheno", "mean", "sd", "rcoef", "CI1", "CI2", "p_value", "se")
df_res <- data.frame(matrix(ncol = length(x), nrow = 0))
colnames(df_res) <- x
for(i in 1:ncol(df)) {       # for-loop over columns
  reg_TLc_pheno <- lm(df$TLc ~ df[, i])
  df_res[nrow(df_res) + 1,] <- c(colnames(df)[i], 
                                 if(class(df[,i]) == 'integer' | class(df[,i]) == 'numeric'){
                                   mean(df[, i], na.rm = T)
                                 } else NA,
                                 if(class(df[,i]) == 'integer' | class(df[,i]) == 'numeric'){
                                   sd(df[, i], na.rm = T)
                                 } else NA,
                                 reg_TLc_pheno$coefficients[[2]],
                                 confint(reg_TLc_pheno)[2,1],
                                 confint(reg_TLc_pheno)[2,2],
                                 lmp(reg_TLc_pheno),
                                 coef(summary(reg_TLc_pheno))[2, 2])
  print(colnames(df)[i])
  print(reg_TLc_pheno)
}
# Put data in df_res as numeric and factor
df_res$pheno <- as.factor(df_res$pheno)
df_res[,2:ncol(df_res)] <- lapply(df_res[,2:ncol(df_res)], as.numeric)
# Remove TL data (too correlated)
df_res <- filter(df_res, !pheno %in% c("TLc", "Zadj_TS"))
# Arrange in descending order of rcoef
df_res_F <- arrange(df_res, desc(rcoef))

### Merge dataframes
df_res_M$sex <- "M"
df_res_F$sex <- "F"
df_res <- bind_rows(df_res_M, df_res_F)
df_res$sex <- as.factor(df_res$sex)


################################################
### Plotting ###################################

df_res_all <- left_join(df_res_M, df_res_F, by = "pheno")
df_res_all <- mutate(df_res_all, 
               pdiff = 2*pnorm(-abs((rcoef.x-rcoef.y)/sqrt(se.x**2+se.y**2)), mean = 0, sd = 1))

# Remove which are not significantly different 
df_res_plot <- df_res_all

### Put descriptive title
df_code <- read.csv(file = file.path("../",  
                                     "../TL_metadata.csv"), 
                    sep = ',', header = TRUE)
df_res_plot <- merge(df_res_plot, df_code[,c("pheno","description")], by="pheno")

# Remove non-wanted labels
df_res_plot$description[which(df_res_plot$pdiff > 0.05)] <- ""



# Plot
comp_plot <- ggplot(df_res_plot, aes(x=rcoef.x, y=rcoef.y)) +
  geom_errorbar(aes(ymax = CI2.y, ymin = CI1.y, width = 0), color = "#E9AfA3",
                alpha = ifelse(df_res_plot$pdiff < 0.05,0.8,0)) +
  geom_errorbar(aes(xmax = CI2.x, xmin = CI1.x, width = 0), color = "#99B2DD", 
                alpha = ifelse(df_res_plot$pdiff < 0.05,0.8,0)) +
  geom_point(pch=21, 
             fill =  ifelse(df_res_plot$pdiff < 0.05/indep_trait,"#8338ec","#CBBDDD"), 
             size = 1.5, 
             colour = "black",
             alpha = ifelse(df_res_plot$pdiff < 0.05,0.9,0.2)) + # Show dots
  geom_abline(slope = 1, color = "#685044", alpha = 0.2, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  ggrepel::geom_text_repel(
    label=df_res_plot$description, 
    size = 2.8,
    #fontface = "bold"
    color = ifelse(df_res_plot$pdiff < 0.05/indep_trait, "red", "black"),
    max.overlaps = 55,
    box.padding = 1.2,
    seed = 1996
    )+
  xlab(expression(""*beta*" males (95% CI)")) + 
  ylab(expression(""*beta*" females (95% CI)"))
  
mylim <- 0.075
comp_plot + 
  theme_light() +
  scale_x_continuous(breaks = seq(-0.06, 0.06, by = 0.02), limits = c(-mylim, mylim)) +
  scale_y_continuous(breaks = seq(-0.06, 0.06, by = 0.02), limits = c(-mylim, mylim)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.title=element_text(size=12)) 
  
### Get information
df_res_sub_info <- df_res_all[which(df_res_all$pdiff < 0.05),]

for(pheno_analysed in df_res_sub_info$pheno){
  index <- which(df_res_sub_info$pheno == pheno_analysed)
  pvalformat_m <- unlist(base::strsplit(as.character(df_res_sub_info$p_value.x[index]),"e"))
  pvalformat_f <- unlist(base::strsplit(as.character(df_res_sub_info$p_value.y[index]),"e"))
  print(paste0(pheno_analysed, " (\u03B2male = ", round(df_res_sub_info$rcoef.x[index],4), #"; 95%CI ",
               #round(df_res_sub_info$CI1.x[index],4)," - ",round(df_res_sub_info$CI2.x[index],4),
               "; pmale = ", round(as.numeric(pvalformat_m[1]),1)," x 10",pvalformat_m[2],
        " | ",   "\u03B2female = ", round(df_res_sub_info$rcoef.y[index],4), #"; 95%CI ",
               #round(df_res_sub_info$CI1.x[index],4)," - ",round(df_res_sub_info$CI2.x[index],4),
               "; pfemale = ", round(as.numeric(pvalformat_f[1]),1)," x 10",pvalformat_f[2],")"))
}
