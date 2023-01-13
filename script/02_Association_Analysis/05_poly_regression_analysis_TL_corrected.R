################################################################################
### Polynomial regression coefficients with multiple phenotypes              ###
### Author: Samuel Moix                                                      ###
### Date: 26.05.2022                                                         ###
################################################################################

################################################
### Libraries ##################################
library(data.table) # Read data 
library(ggplot2) # Plotting basics
library(dplyr)
library(forcats)
require(gridExtra) # Putting multiple plots together
library(tibble)

################################################
### Working directories ########################
data_folder = "../"
export_folder = "../"

################################################
### Parameters  ################################

### Scale all data
par_scale <- TRUE
###  Calculate number of independent traits automatically if FALSE set (indep_trait <- XXX)
par_calc_indep_traits <- FALSE
indep_trait <- 135 # (See script 03_regression_analysis...)
### TRUE to keep non-signficative traits on the regression plot
par_keep_not_sign <- FALSE
### Weather to check if exponential regressions fit better than linear/polynomial 
par_exp_reg_analysis <- TRUE
# True if linear, False if compare with polynomial rsquared
par_linear <- FALSE
### Weather to plot the fitted lines
par_plot_fit <- FALSE
### Weather to remove outliers from data
par_outlier_filter <- TRUE
### Whether to remove patients with blood cancer
par_rm_blood_cancer <- TRUE
### Correcting parameters
### Correct traits (cholesterol for medication and female traits for SES)
par_corr_traits <- FALSE 
par_keep_not_cor <- FALSE # Whether to keep traits not corrected 

################################################
### Load Data ##################################

df <- as.data.frame(fread(file.path(data_folder,  "filtered_DF.txt"), header = T)) 

### Linear regression model TL corrected
df$TLc <- residuals(lm(df$Zadj_TS ~ df$age + df$sex + df$age:df$sex + df$age2 + 
                         df$age2:df$sex + df$array + df$array:df$sex))


### Remove blood cancer patients 
if(par_rm_blood_cancer){
  # Remove patients with blood cancer
  df <- df[which(df$blood_cancer == 0 | is.na(df$blood_cancer)),]
  df <- select(df, -c("blood_cancer"))
}


################################################
### Correct Data ###############################


### Correct traits
if(par_corr_traits){
  # Correct cholesterol for cholesterol medication
  # Changes in mmol/L from https://pubmed.ncbi.nlm.nih.gov/22354153/
  # and from Simvastatin https://bmcprimcare.biomedcentral.com/articles/10.1186/1471-2296-4-18 
  
  df <- df %>% tibble::add_column(cholesterol_cor = ifelse(df$cholest_medication == 0, df$cholesterol,
                                                   df$cholesterol + 1.64)
                          , .after = "cholesterol" )
  
  df <- df %>% tibble::add_column(HDL_cor = ifelse(df$cholest_medication == 0, df$HDL,
                                           df$HDL - 0.1)
                          , .after = "HDL" )
  
  df <- df %>% tibble::add_column(LDL_cor = ifelse(df$cholest_medication == 0, df$LDL,
                                           df$LDL + 1.7)
                          , .after = "LDL" )
  
  if(!par_keep_not_cor){
    df$cholesterol <- df$cholesterol_cor
    df$HDL <- df$HDL_cor
    df$LDL <- df$LDL_cor
    df <- select(df, -c("cholesterol_cor","HDL_cor","LDL_cor"))
  }
  
  # Correct female trait for SES
  trait_to_correct_SES <- c("age_first_birth", "age_last_birth", "birthweight",
                            "menopause", "reproductive_lp", "menstruation",
                            "birth_weight_first_child","nb_birth")
  
  if(par_keep_not_cor){
    df <- add_column(df, "age_first_birth_cor_SES" = NA, .after = "age_first_birth") 
    df <- add_column(df, "age_last_birth_cor_SES" = NA, .after = "age_last_birth") 
    df <- add_column(df, "birthweight_cor_SES" = NA, .after = "birthweight") 
    df <- add_column(df, "menopause_cor_SES" = NA, .after = "menopause") 
    df <- add_column(df, "reproductive_lp_cor_SES" = NA, .after = "reproductive_lp") 
    df <- add_column(df, "menstruation_cor_SES" = NA, .after = "menstruation") 
    df <- add_column(df, "birth_weight_first_child_cor_SES" = NA, .after = "birth_weight_first_child") 
    df <- add_column(df, "nb_birth_cor_SES" = NA, .after = "nb_birth") 
  }
  
  temp_name <- ""
  if(par_keep_not_cor){
    temp_name <- "_cor_SES"
  }
  
  for (trait in trait_to_correct_SES) { # Correct for all traits in list
    df[paste0(trait, temp_name)] <- residuals(lm(df[,which(colnames(df) == trait)] ~ df$TDI + df$household_income +
                                                   df$age_end_education, na.action=na.exclude))
  }
  rm(trait)
  rm(temp_name)
}

################################################
### Total scaling and data segregation by type #

### Create df only containing continuous traits ###
contvec <- c(1:which(colnames(df) == "age_cancer"), which(colnames(df) == "TLc"))
df_continuous <- df[,contvec]
df_continuous <- select(df_continuous, -contains(c("blood_cancer", "sex")))

### Remove all outliers from data (defined as �5SDs)
if(par_outlier_filter){
  for (i in 1:length(df_continuous)) {
    upper <- mean(df_continuous[,i], na.rm = T) + sd(df_continuous[,i], na.rm = T)*5
    lower <- mean(df_continuous[,i], na.rm = T) - sd(df_continuous[,i], na.rm = T)*5
    df_continuous[,i][which(df_continuous[,i] > upper | df_continuous[,i] < lower)] <- NA
    rm(upper)
    rm(lower)
    }
}

if (par_scale){
  ### Scale continuous data and replace in data frame ###
  df_continuous <- as.data.frame(scale(df_continuous)) ###
}

#### Change input df
df <- df_continuous
rm(df_continuous)



################################################
### Calculating number of independent traits ###

if (par_calc_indep_traits){
  ### Create a correlation matrix, do not consider NA
  pheno_cor <- cor(df, use="pairwise.complete.obs")
  ### Set missing values to 0 to allow svd() ###
  pheno_cor[which(is.na(pheno_cor))] <- 0
  ### Calculate the eigenvalues ###
  # svd() computes the singular-value decomposition of a rectangular matrix, with $d being a vector containing the singular values of the decomposition
  pheno_EV <- svd(pheno_cor)$d
  ### Calculate Neff ###
  # Neff is defined as the number of eigenvalues required to explain 99.5% of the variation from the CNV data 
  sum_EV <- 0
  count <- 0
  while(sum_EV/sum(pheno_EV) < 0.995) {
    count <- count + 1
    sum_EV <- sum_EV + pheno_EV[count]
  }
  ### Set number of independent traits
  indep_trait <- count
}

################################################
### Calculating regression coef ################

### Function to calculate p-value
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

x <- c("pheno", "rcoef", "CI1", "CI2", "p_value", "rsquared")

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
                                 summary(reg_TLc_pheno)$adj.r.squared)
  print(colnames(df)[i])
  print(reg_TLc_pheno)
}

# Put data in df_res as numeric and factor
df_res$pheno <- as.factor(df_res$pheno)
df_res[,2:ncol(df_res)] <- lapply(df_res[,2:ncol(df_res)], as.numeric)

# Remove TL data (too correlated)
df_res <- filter(df_res, !pheno %in% c("TLc", "Zadj_TS"))

# Arrange in descending order of rcoef
df_res <- arrange(df_res, desc(rcoef))

# Subset significant traits
df_res_sub <- df_res[which(df_res$p_value < 0.05/indep_trait),]
if (par_keep_not_sign){
  df_res_sub <- df_res
}

### Plotting 
reg_plot_all_scaled <- df_res_sub %>%
  mutate(pheno = fct_reorder(pheno, desc(rcoef))) %>%
  ggplot(aes(x=pheno, y=rcoef)) +
  geom_errorbar(aes(ymax = CI2, ymin = CI1, width = .3), color = "grey", alpha = 0.8) +
  geom_point(colour = ifelse(df_res_sub$p_value < 0.05/indep_trait, "black", "grey")) +
  coord_flip() +
  #ylim(round(min(df_res_sub$CI1),2), round(max(df_res_sub$CI2),2)) +
  xlab("") +
  ylab(expression("Beta"[2]*" (95% CI)")) +
  #labs(title = "TL corrected ~ trait (sd = 1, u = 0)")+
  theme(axis.text.y = element_text(size=8))  # CHANGE LABEL SIZE

reg_plot_all_scaled + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

### Export plot
export_file_name <- "reg_plot_pol2"

if(par_scale){
  export_file_name <- paste0(export_file_name, "_scaled")
}

if(par_outlier_filter){
  export_file_name <- paste0(export_file_name, "_no_outliers")
}

if(par_corr_traits){
  export_file_name <- paste0(export_file_name, "_cor_traits")
}

if(par_keep_not_cor){
  export_file_name <- paste0(export_file_name, "_and_NC")
}

if(par_rm_blood_cancer){
  export_file_name <- paste0(export_file_name, "_no_bc")
}

export_file_name <- paste0(export_file_name, ".pdf")
ggsave(file.path(export_folder, export_file_name),
       reg_plot_all_scaled, width=3, height=3, units="in", scale=3)

### Plot fitted model line over scatter plot 
# Add descriptive name
metadata <- read.csv(file = "../TL_metadata.csv", sep = ',', header = TRUE)
df_res_sub <- merge(df_res_sub, metadata[,c("pheno","description")], by = "pheno", all.x = TRUE)
df_res_sub$description[which(df_res_sub$pheno == "physical_activity")] <- "Physical activity"
df_res_sub$description[which(df_res_sub$pheno == "WBC_count")] <- "White blood cell count"
df_res_sub$description[which(df_res_sub$pheno == "MCH")] <- "MCH"

# Plot
if(par_plot_fit){
  dfs <- df[sample(nrow(df), 10000), ] # Subset for plotting <-------
  list_traits_pol <- as.vector(df_res_sub$pheno)
  plot_list <- list()
  for (i in 1:length(list_traits_pol)) {
    plot_list[[i]]<-ggplot(dfs, aes_string(x = list_traits_pol[i], y = "TLc")) + 
      geom_point(alpha = 0.2)+ 
      ylab("") +
      xlab(df_res_sub$description[which(df_res_sub$pheno == list_traits_pol[i])]) +
      stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
      theme(axis.title.x = element_text(size = 8))
  }
  grid_plot_pol_fit <- gridExtra::grid.arrange(grobs = plot_list, ncol=5,top="TLc ~ x + x**2")
  
  ### Save plot
  export_file_name <- "grid_plot_fit_pol"
  
  if(par_outlier_filter){
    export_file_name <- paste0(export_file_name, "_no_outliers")
  }
  
  if(par_rm_blood_cancer){
    export_file_name <- paste0(export_file_name, "_no_bc")
  }
  
  export_file_name <- paste0(export_file_name, ".pdf")
  ggsave(file.path(export_folder, export_file_name),
         grid_plot_pol_fit, width=3, height=3, units="in", scale=3)
}


################################################
### Exponential regression comparison ##########

if (par_exp_reg_analysis){
  ### Subset only significant in bionmial 
  df <- select(df, c(as.character(df_res_sub$pheno), "TLc"))

    if(par_linear){
    
    ### Calculating linear regression coef
    ### Results dataframe
    x <- c("pheno", "rcoef", "CI1", "CI2", "p_value", "rsquared")
    
    df_res <- data.frame(matrix(ncol = length(x), nrow = 0))
    colnames(df_res) <- x
    
    
    for(i in 1:ncol(df)) {       # for-loop over columns
      reg_TLc_pheno <- lm(df$TLc ~ df[, i])
      df_res[nrow(df_res) + 1,] <- c(colnames(df)[i], 
                                     reg_TLc_pheno$coefficients[[2]],
                                     confint(reg_TLc_pheno)[2,1],
                                     confint(reg_TLc_pheno)[2,2],
                                     lmp(reg_TLc_pheno),
                                     summary(reg_TLc_pheno)$adj.r.squared)
      print(colnames(df)[i])
      print(reg_TLc_pheno)
    }
    
    # Put data in df_res as numeric and factor
    df_res$pheno <- as.factor(df_res$pheno)
    df_res[,2:ncol(df_res)] <- lapply(df_res[,2:ncol(df_res)], as.numeric)
    # Remove TL data (too correlated)
    df_res <- filter(df_res, !pheno %in% c("TLc", "Zadj_TS"))
    # Arrange in descending order of rcoef
    df_res <- arrange(df_res, desc(rcoef))
    # No need to subset as already subset beforhand 
    df_res_sub <- df_res
    
  }
  
  ### Add constant (min value) to avoid having y <= 0 for log
  df$TLc <- df$TLc + ceiling(abs(min(df$TLc, na.rm = TRUE)))
  
  df_res_exp <- data.frame(matrix(ncol = length(x), nrow = 0))
  colnames(df_res_exp) <- x
  
  for(i in 1:ncol(df)) {       # for-loop over columns
    reg_TLc_pheno <- lm(log(df$TLc) ~ df[, i])
    df_res_exp[nrow(df_res_exp) + 1,] <- c(colnames(df)[i], 
                                           reg_TLc_pheno$coefficients[[2]],
                                           confint(reg_TLc_pheno)[2,1], 
                                           confint(reg_TLc_pheno)[2,2],
                                           lmp(reg_TLc_pheno),
                                           summary(reg_TLc_pheno)$adj.r.squared)
    print(colnames(df)[i])
    print(reg_TLc_pheno)
  }
  
  # Put data in df_res as numeric and factor
  df_res_exp$pheno <- as.factor(df_res_exp$pheno)
  df_res_exp[,2:ncol(df_res_exp)] <- lapply(df_res_exp[,2:ncol(df_res_exp)], as.numeric)
  # Remove TL data (too correlated)
  df_res_exp <- filter(df_res_exp, !pheno %in% c("TLc"))
  # Arrange in descending order of rcoef
  df_res_exp <- arrange(df_res_exp, desc(rcoef))
  # Subset relevant traits
  df_res_exp_sub <- df_res_exp[which(df_res_exp$p_value < 0.05/indep_trait),]
  
  ### See which have a better r-squared
  exp_better_rsquare <- df_res_exp_sub %>%
    select(c("pheno", "rsquared")) %>%
    rename(rsquared_exp = "rsquared") %>%
    inner_join(select(df_res_sub, c("pheno", "rsquared")), by = "pheno") %>%
    filter(rsquared_exp > rsquared)
  
  print(paste("There are", length(exp_better_rsquare$pheno), 
              "phenotypes with a higher rsquared with an expoential regression model"))
  
  ### Plotting 
  if(par_linear){ # Remove traits where R.adj not bigger as linear
    df_res_exp_sub_radj <- filter(df_res_exp_sub, pheno %in% exp_better_rsquare$pheno)
    
    reg_plot_exp <- df_res_exp_sub_radj %>%
      mutate(pheno = fct_reorder(pheno, desc(rcoef))) %>%
      ggplot(aes(x=pheno, y=rcoef)) +
      geom_errorbar(aes(ymax = CI2, ymin = CI1, width = .3), color = "grey", alpha = 0.8) +
      geom_point(colour = ifelse(df_res_exp_sub_radj$p_value < 0.05/indep_trait, "black", "grey")) +
      coord_flip() +
      #ylim(round(min(df_res_sub$CI1),2), round(max(df_res_sub$CI2),2)) +
      xlab("") +
      ylab("Regression coefficient") +
      #labs(title = "TL corrected ~ trait (sd = 1, u = 0)")+
      theme(axis.text.y = element_text(size=8))  # CHANGE LABEL SIZE
    
    reg_plot_exp
    
    rm(df_res_exp_sub_radj)
    
    ### Export plot
    export_file_name <- "reg_plot_exp"
    
    if(par_scale){
      export_file_name <- paste0(export_file_name, "_scaled")
    }
    
    if(par_outlier_filter){
      export_file_name <- paste0(export_file_name, "_no_outliers")
    }
    
    if(par_corr_traits){
      export_file_name <- paste0(export_file_name, "_cor_traits")
    }
    
    if(par_keep_not_cor){
      export_file_name <- paste0(export_file_name, "_and_NC")
    }
    
    if(par_rm_blood_cancer){
      export_file_name <- paste0(export_file_name, "_no_bc")
    }
    
    export_file_name <- paste0(export_file_name, ".pdf")
    ggsave(file.path(export_folder, export_file_name),
           reg_plot_exp, width=3, height=3, units="in", scale=3)
  }
}


### Plot fitted model line over scatter plot 
if(par_exp_reg_analysis){
  if(length(exp_better_rsquare$pheno != 0)){
    if(par_plot_fit){
      dfs <- df[1:10000,] # Subset for plotting <-------
      list_traits_exp <- as.vector(exp_better_rsquare$pheno)
      plot_list <- list()
      for (i in 1:length(list_traits_exp)) {
        log.model <-lm(log(TLc) ~ dfs[,which(colnames(dfs) == list_traits_exp[i])], 
                       data = dfs,  na.action=na.exclude)
        plot_list[[i]]<-ggplot(dfs, aes_string(x = list_traits_exp[i], y = "TLc")) + 
          geom_point(alpha = 0.2) +
          geom_line(data = dfs, aes_string(x = list_traits_exp[i], exp(fitted(log.model)))
                    , size = 1, linetype = 1, col ="blue")+
          ylab("") 
      }
      
      grid_plot_exp_fit <- gridExtra::grid.arrange(grobs = plot_list, ncol=1,top="log(TLc + c) ~ x")
      
      ### Save plot
      export_file_name <- "grid_plot_fit_exp"
      
      if(par_outlier_filter){
        export_file_name <- paste0(export_file_name, "_no_outliers")
      }
      
      if(par_rm_blood_cancer){
        export_file_name <- paste0(export_file_name, "_no_bc")
      }
      
      export_file_name <- paste0(export_file_name, ".pdf")
      ggsave(file.path(export_folder, export_file_name),
             grid_plot_exp_fit, width=3, height=3, units="in", scale=3)
    }
  }
}
