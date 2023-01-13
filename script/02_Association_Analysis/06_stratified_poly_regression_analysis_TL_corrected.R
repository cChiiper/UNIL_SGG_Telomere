#################################################################################################################
### Stratified regression coefficients with multiple phenotypes                                                          ###
### Author: Samuel Moix                                                                                       ###
### Date: 31.05.2022                                                                                          ###
#################################################################################################################


################################################
### Libraries ##################################

library(data.table) # Read data 
library(ggplot2) # Plotting basics
library(dplyr)
library(forcats)

################################################
### Parameters  ################################

### Scale all data
par_scale <- FALSE
### Set number of independent traits
indep_trait <- 135


####### LOAD !!!!!!!!!

data_folder = "../"

df <-  as.data.frame(fread(file.path(data_folder,  "filtered_DF.txt"), header = T))

### Remove single-sex traits
df <- select(df, -c(menarche, menopause, BC, OC, endometriosis, menstruation, nb_birth,
                      birth_weight_first_child, age_first_birth, age_last_birth, PC, facial_hair,
                    balding, reproductive_lp))


################################################
### Data type segregation ######################

df_M <- df[which(df$sex == 1), -4] # -4 to remove sex
df_F <- df[which(df$sex == 2), -4]


df_M$TLc <- residuals(lm(df_M$Zadj_TS ~ df_M$age + df_M$age2 + df_M$array))
df_F$TLc <- residuals(lm(df_F$Zadj_TS ~ df_F$age + df_F$age2 + df_F$array))



### Create df only containing continuous traits ###!
contvec <- c(1:which(colnames(df_M) == "age_cancer"), which(colnames(df_M) == "TLc"))
df_M <- df_M[,contvec]
df_F <- df_F[,contvec]
df_M <- select(df_M, -c(blood_cancer))
df_F <- select(df_F, -c(blood_cancer))


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

### Regression coefficients on df for males
df <- df_M
x <- c("pheno", "rcoef", "CI1", "CI2", "p_value", "se")
df_res <- data.frame(matrix(ncol = length(x), nrow = 0))
colnames(df_res) <- x
for(i in 1:ncol(df)) {       # for-loop over columns
  reg_TLc_pheno <- lm(df$TLc[which(!is.na(df[, i]))] ~ df[, i][which(!is.na(df[, i]))] + 
                        I(df[, i][which(!is.na(df[, i]))]**2))
  df_res[nrow(df_res) + 1,] <- c(colnames(df)[i], 
                                 reg_TLc_pheno$coefficients[[3]],
                                 confint(reg_TLc_pheno)[3,1], #3,1 if rcoef2
                                 confint(reg_TLc_pheno)[3,2],
                                 coef(summary(reg_TLc_pheno))[, 4][[3]],
                                 coef(summary(reg_TLc_pheno))[3, 2])
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
x <- c("pheno", "rcoef", "CI1", "CI2", "p_value", "se")
df_res <- data.frame(matrix(ncol = length(x), nrow = 0))
colnames(df_res) <- x
for(i in 1:ncol(df)) {       # for-loop over columns
  reg_TLc_pheno <- lm(df$TLc[which(!is.na(df[, i]))] ~ df[, i][which(!is.na(df[, i]))] + 
                        I(df[, i][which(!is.na(df[, i]))]**2))
  df_res[nrow(df_res) + 1,] <- c(colnames(df)[i], 
                                 reg_TLc_pheno$coefficients[[3]],
                                 confint(reg_TLc_pheno)[3,1], #3,1 if rcoef2
                                 confint(reg_TLc_pheno)[3,2],
                                 coef(summary(reg_TLc_pheno))[, 4][[3]],
                                 coef(summary(reg_TLc_pheno))[3, 2])
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

###
df_res_M$sex <- "M"
df_res_F$sex <- "F"
df_res <- bind_rows(df_res_M, df_res_F)
df_res$sex <- as.factor(df_res$sex)


### Plot 2 !!!

df_res_all <- left_join(df_res_M, df_res_F, by = "pheno")
df_res_all <- mutate(df_res_all, 
               pdiff = 2*pnorm(-abs((rcoef.x-rcoef.y)/sqrt(se.x**2+se.y**2)), mean = 0, sd = 1))

# Remove which are not significantly different 
df_res_all_sub <- df_res_all[which(df_res_all$pdiff < 0.05/indep_trait),]

ggplot(df_res_all_sub, aes(x=rcoef.x, y=rcoef.y)) +
  geom_errorbar(aes(ymax = CI2.y, ymin = CI1.y, width = 0), color = "red", alpha = 0.1) +
  geom_errorbar(aes(xmax = CI2.x, xmin = CI1.x, width = 0), color = "blue", alpha = 0.1) +
  geom_point(shape=1) + # Show dots
  geom_abline(slope = 1, color = "grey", alpha = 0.5) +
  geom_text(
    label=df_res_all_sub$pheno, 
    nudge_y = 0.0002, #0.0002 when scaled 0.0001 if not
    check_overlap = T,
    size = 2
  ) +
  xlab("rcoef males") +
  ylab("rcoef females") +
  ylim(0,0.035)
                
