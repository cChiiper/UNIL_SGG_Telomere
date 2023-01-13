################################################################################
### Analysis to test whether pregnancy has an effect on TL                   ###
### Author: Samuel Moix                                                      ###
### Date: 11.12.2022                                                         ###
################################################################################

################################################
### Libraries ##################################
library(data.table) # Read data 
library(ggplot2) # Plotting basics
library(dplyr)

################################################
### Working directories ########################
data_folder = "../"
export_folder = "../"

################################################
### Parameters  ################################

###  Number of independent traits
indep_trait <- 135
### Whether to remove patients with blood cancer
par_rm_blood_cancer <- TRUE
### Remove outliers
par_outlier_filter <- TRUE

################################################
### Load Data ##################################

df <- as.data.frame(fread(file.path(data_folder,  "filtered_DF.txt"), header = T)) 
metadata <- read.csv(file = "../TL_metadata.csv", sep = ',', header = TRUE)


### Remove blood cancer patients 
if(par_rm_blood_cancer){
  # Remove patients with blood cancer
  df <- df[which(df$blood_cancer == 0 | is.na(df$blood_cancer)),]
  df <- select(df, -c("blood_cancer"))
}

### Linear regression model TL corrected 
df$TLc <- residuals(lm(df$Zadj_TS ~ df$age + df$sex + df$age:df$sex + df$age2 + df$age2:df$sex + df$array + df$array:df$sex))

df <- df %>%
  select(c("Zadj_TS","TLc","sex","age","nb_birth", "age_last_birth", "age_first_birth"))

# Outliers defined as �5SDs
if(par_outlier_filter){
  print("Outliers filtered")
  for (i in c("Zadj_TS","TLc","nb_birth", "age_last_birth", "age_first_birth")) {
    upper <- mean(df[,i], na.rm = T) + sd(df[,i], na.rm = T)*5
    lower <- mean(df[,i], na.rm = T) - sd(df[,i], na.rm = T)*5
    df[,i][which(df[,i] > upper | df[,i] < lower)] <- NA
    rm(upper)
    rm(lower)
  }
}

### Remove TL outliers
df <- df[which(!is.na(df$Zadj_TS)),]

################################################
### Scale data and add column ##################

df$Zadj_TS <- as.vector(scale(df$Zadj_TS))
df$TLc <- as.vector(scale(df$TLc))
df$age_last_birth_scaled <- as.vector(scale(df$age_last_birth))


### Add ever had a pregnancy
df$ever_pregnant <- ifelse(df$nb_birth > 0, 1, 0)

################################################
### Calculating regression coef ################

### P-value function
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

### Regress against ever had kids or not
df$ever_pregnant_scaled <- as.vector(scale(df$ever_pregnant))
model <- lm(df$TLc ~ df$ever_pregnant_scaled)
print(lmp(model))
print(model)
print(confint(model)[2,1])
print(confint(model)[2,2])

# Welch two-sample t-test
wtst_test <- t.test(TLc ~ ever_pregnant, data = df)
print(wtst_test)
wtst_mean_0 <- wtst_test$estimate[[1]]
wtst_mean_1 <- wtst_test$estimate[[2]]
(wtst_mean_0-wtst_mean_1)/wtst_mean_1*100


################################################
### TL ~ age_last_birth by age #################

df_young_female <- df[which(df$sex == 2 & df$age <= 50),]
df_older_female <- df[which(df$sex == 2 & df$age > 50),]

### Check all together
lm(df$TLc ~ df$age_last_birth_scaled)

### Compare models
young_reg <- lm(df_young_female$TLc ~ df_young_female$age_last_birth_scaled)
young_b <- young_reg$coefficients[[2]]
young_se <- coef(summary(young_reg))[2, 2]

older_reg <- lm(df_older_female$TLc ~ df_older_female$age_last_birth_scaled)
older_b <- older_reg$coefficients[[2]]
older_se <- coef(summary(older_reg))[2, 2]

### Show results
print("TLc ~ age last birth by age group")
pdiff <- 2*pnorm(-abs((young_b-older_b)/sqrt(young_se**2+older_se**2)),mean = 0, sd = 1)
print(paste("Young b:", round(young_b,3), "Young se:", round(young_se,4),
            "Older b:", round(older_b,3), "Older se:", round(older_se,4),
            "P_diff:", pdiff), sep=" ")

rm(df_young_female) 
rm(df_older_female)

################################################
### TL ~ age by sex (corrected by age last birth

df_M <- df[which(df$sex == 1),]
df_F <- df[which(df$sex == 2),]

### Adjust female length by age_last_birth

df_F <- df_F[which(!is.na(df_F$age_last_birth_scaled)),]
df_F$Zadj_TS_corrected <- residuals(lm(df_F$Zadj_TS ~ df_F$age_last_birth_scaled))

### Scale TL
df_M$Zadj_TS <- as.vector(df_M$Zadj_TS)
df_F$Zadj_TS <- as.vector(df_F$Zadj_TS)
df_F$Zadj_TS_corrected <- as.vector(df_F$Zadj_TS_corrected)

### Compare models
male_reg <- lm(df_M$Zadj_TS ~ df_M$age)
male_b <- male_reg$coefficients[[2]]
male_se <- coef(summary(male_reg))[2, 2]

female_reg <- lm(df_F$Zadj_TS_corrected ~ df_F$age)
#female_reg <- lm(df_F$Zadj_TS ~ df_F$age)
female_b <- female_reg$coefficients[[2]]
female_se <- coef(summary(female_reg))[2, 2]

### Show results
print("Female adjusted for age_last_birth")
pdiff <- 2*pnorm(-abs((male_b-female_b)/sqrt(male_se**2+female_se**2)),mean = 0, sd = 1)
print(paste("Male b:", round(male_b,3), "male se:", round(male_se,4),
            "Female b:", round(female_b,3), "female se:", round(female_se,4),
            "P_diff:", pdiff), sep=" ")


################################################
### TL ~ age by sex (without pregnant women) ###

df_M <- df[which(df$sex == 1),]
df_F <- df[which(df$sex == 2 & df$ever_pregnant == 0),]

### Scale TL
df_M$Zadj_TS <- as.vector(df_M$Zadj_TS)
df_F$Zadj_TS <- as.vector(df_F$Zadj_TS)

### Compare models
male_reg <- lm(df_M$Zadj_TS ~ df_M$age)
male_b <- male_reg$coefficients[[2]]
male_se <- coef(summary(male_reg))[2, 2]

female_reg <- lm(df_F$Zadj_TS ~ df_F$age)
female_b <- female_reg$coefficients[[2]]
female_se <- coef(summary(female_reg))[2, 2]

### Show results
print("No pregnant women")
pdiff <- 2*pnorm(-abs((male_b-female_b)/sqrt(male_se**2+female_se**2)),mean = 0, sd = 1)
print(paste("Male b:", round(male_b,3), "male se:", round(male_se,4),
            "Female b:", round(female_b,3), "female se:", round(female_se,4),
            "P_diff:", pdiff), sep=" ")


################################################
### Add time since last birth ##################

df$time_since_last_birth <- df$age - df$age_last_birth
df$time_since_last_birth_scaled <- as.vector(scale(df$time_since_last_birth))

model <- lm(df$TLc ~ df$time_since_last_birth_scaled)
print(lmp(model))
print(model)


################################################
### TL ~ age had kids vs never had hids ########

df_F_kid <- df[which(df$ever_pregnant == 1),]
#df_F_kid <- df_F_kid[sample(1:141000,30000),]
df_F_nokid <- df[which(df$ever_pregnant == 0),]

### Scale TL
df_F_nokid$Zadj_TS <- as.vector((df_F_nokid$Zadj_TS))
df_F_kid$Zadj_TS <- as.vector((df_F_kid$Zadj_TS))

### Compare models
nokid_reg <- lm(df_F_nokid$Zadj_TS ~ df_F_nokid$age)
nokid_b <- nokid_reg$coefficients[[2]]
nokid_se <- coef(summary(nokid_reg))[2, 2]

kid_reg <- lm(df_F_kid$Zadj_TS ~ df_F_kid$age)
kid_b <- kid_reg$coefficients[[2]]
kid_se <- coef(summary(kid_reg))[2, 2]

### Show results
print("No kids vs kids")
pdiff <- 2*pnorm(-abs((nokid_b-kid_b)/sqrt(nokid_se**2+kid_se**2)),mean = 0, sd = 1)
print(paste("No kid b:", round(nokid_b,4), "No kid se:", round(nokid_se,4),
            "Kid b:", round(kid_b,4), "kid se:", round(kid_se,4),
            "P_diff:", pdiff), sep=" ")
