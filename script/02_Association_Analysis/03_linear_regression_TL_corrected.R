################################################################################
### Regression coefficients with multiple phenotypes                         ###
### Author: Samuel Moix                                                      ###
### Date: 04.04.2022                                                         ###
################################################################################

################################################
### Libraries ##################################
library(data.table) # Read data 
library(ggplot2) # Plotting basics
library(ggExtra) # Extra plotting XP
library(ggpubr) # Plotting regression stuff
library(dplyr)
library(tibble)
library(forcats)
require(pheatmap)
library(ggforestplot) # Forest plot

################################################
### Working directories ########################
data_folder = "../"
export_folder = "../"

################################################
### Parameters  ################################

### Scale all data
par_scale <- TRUE
###  Use only continuous traits
par_only_cont <- FALSE
###  Use only factorial traits
par_only_fact <- FALSE
###  Calculate number of independent traits automatically if FALSE set (indep_trait <- XXX)
par_calc_indep_traits <- FALSE
indep_trait <- 135
### TRUE to keep non-significant traits on the regression plot
par_keep_not_sign <- FALSE
### Generate correlation plots
par_cor_plot <- FALSE
###  Generate correlation plot with PCs 
par_check_cor_PC <- FALSE
### Whether to remove outliers from data
par_outlier_filter <- TRUE
### Whether to remove patients with blood cancer
par_rm_blood_cancer <- TRUE
### Correcting parameters
### Correct traits (cholesterol for medication and female traits for SES)
par_corr_traits <- FALSE
par_keep_not_cor <- FALSE # Whether to keep traits not corrected 
### Compare two results
par_compare <- FALSE # Has to be on true to generate rds
par_save_rds <- FALSE # Set first TRUE (save first results) then FALSE for comparison
par_analysis_names <- c("Not corrected", "Corrected") # First for first analysis
### Finalizing parameters
par_final_corr <- FALSE
par_final_plot <- FALSE


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

# Subset to traits for correction analysis
if(par_final_corr){
  df <- select(df,c("TLc","age_first_birth", "age_last_birth",
                        "menopause", "reproductive_lp", "menstruation",
                        "birth_weight_first_child","nb_birth","cholesterol", "HDL", 
                        "LDL","TG", "TDI", "cholest_medication", "household_income",
                        "age_end_education"))
}

################################################
### Correct Data ###############################

### Correct traits
if(par_corr_traits){
  # Correct cholesterol for cholesterol medication
  # Changes in mmol/L from https://pubmed.ncbi.nlm.nih.gov/22354153/
  # and from Simvastatin https://bmcprimcare.biomedcentral.com/articles/10.1186/1471-2296-4-18 
 
  df <- df %>% add_column(cholesterol_cor = ifelse(df$cholest_medication == 0, df$cholesterol,
                                                   df$cholesterol + 1.64) # Patients without CV
                    , .after = "cholesterol" )
  
  df <- df %>% add_column(HDL_cor = ifelse(df$cholest_medication == 0, df$HDL,
                                                   df$HDL - 0.1) # All doses
                          , .after = "HDL" )
  
  df <- df %>% add_column(LDL_cor = ifelse(df$cholest_medication == 0, df$LDL,
                                                   df$LDL + 1.4) # All doses
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
  
  # Correct female trait for SES
  trait_to_correct_SES <- c("age_first_birth", "age_last_birth",
                            "menopause", "reproductive_lp", "menstruation",
                            "birth_weight_first_child","nb_birth")
  
  if(par_keep_not_cor){
    df <- add_column(df, "age_first_birth_cor_SES" = NA, .after = "age_first_birth") 
    df <- add_column(df, "age_last_birth_cor_SES" = NA, .after = "age_last_birth") 
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
  
  if(FALSE){
    # Correct traits for BMI
    trait_to_correct_BMI <- c("hip_cir") #Not WHR, WHRadjBMI
    
    if(par_keep_not_cor){
      df <- add_column(df, "hip_cir_cor_BMI" = NA, .after = "hip_cir") 
    }
    
    temp_name <- ""
    if(par_keep_not_cor){
      temp_name <- "_cor_BMI"
    }
    
    for (trait in trait_to_correct_BMI) { # Correct for all traits in list
      # Subset females for which both traits don't have any NA's
      df_F <- subset(df, sex == 2 & !is.na(df[,which(colnames(df) == trait)]) & !is.na(df$BMI))
      # Regression
      res_F <- residuals(lm(df_F[,which(colnames(df) == trait)] ~ df_F$BMI, na.action=na.exclude))
      # Add residuals to the main data frame
      df[which(df$sex == 2 & !is.na(df[,which(colnames(df) == trait)]) & !is.na(df$BMI)), paste0(trait, temp_name)] <- res_F
      rm(df_F);rm(res_F)
      # Repeat with males
      df_M <- subset(df, sex == 1 & !is.na(df[,which(colnames(df) == trait)]) & !is.na(df$BMI))
      res_M <- residuals(lm(df_M[,which(colnames(df) == trait)] ~ df_M$BMI, na.action=na.exclude))
      df[which(df$sex == 1 & !is.na(df[,which(colnames(df) == trait)]) & !is.na(df$BMI)), paste0(trait, temp_name)] <- res_M
      rm(df_M);rm(res_M)
      # Put NA's on trait's that contained NA's (needed because of the loop syntaxe)
      df[is.na(df[,which(colnames(df) == trait)]) | is.na(df$BMI),which(colnames(df) == paste0(trait, temp_name))] <- NA
    }
    rm(trait)
    rm(temp_name)
  }
  
}


################################################
### Segregation by type ########################

### Create df only containing continuous traits ###
# Make vector of continous traits
contvec <- c(1:which(colnames(df) == "age_cancer"), which(colnames(df) == "TLc"))
# Remove sex and BC
contvec <- contvec[! contvec %in% which(colnames(df) %in% c("blood_cancer", "sex"))] 

if (par_only_cont){
  df_continuous <- df[,contvec]

  if (par_scale){
    ### Scale continuous data and replace in data frame ###
    df_continuous <- as.data.frame(scale(df_continuous)) ###
  }
  
  #### Change input df
  df <- df_continuous
  rm(df_continuous)
}


if (par_only_fact){
  ### Create df only containing factorial traits
  factvec <- c(which(colnames(df) == "blood_cancer"), which(colnames(df) == "sex"),
               which(colnames(df) == "alcool_frq"):which(colnames(df) == "TLc"))
  df <- df[,factvec]
}

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
### Remove continuous traits outliers from data #

# Outliers defined as �5SDs
if(par_outlier_filter & !par_only_fact){
  print("Outliers filtered")
  for (i in contvec) {
    upper <- mean(df[,i], na.rm = T) + sd(df[,i], na.rm = T)*5
    lower <- mean(df[,i], na.rm = T) - sd(df[,i], na.rm = T)*5
    df[,i][which(df[,i] > upper | df[,i] < lower)] <- NA
    rm(upper)
    rm(lower)
  }
}
rm(contvec)

################################################
### Total scaling ##############################

if (par_scale){
  #### To put all traits as numeric
  df[,1:length(df)] <- lapply(df[,1:length(df)], as.numeric) ###
  ### Scale whole data
  df <- as.data.frame(scale(df)) ###
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


### Regression coefficients on df
x <- c("pheno", "mean", "sd", "rcoef", "CI1", "CI2", "p_value", "se")

df_res <- data.frame(matrix(ncol = length(x), nrow = ncol(df)))
colnames(df_res) <- x

for(i in 1:ncol(df)) {       # for-loop over columns
  reg_TLc_pheno <- lm(df$TLc ~ df[, i])
  df_res[i,] <- c(colnames(df)[i], 
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

# Add metadata
df_res <- merge(df_res, metadata[,c("pheno","description","Category")], 
              by = "pheno", all.x = TRUE)


# Arrange in descending order of rcoef
df_res <- arrange(df_res, (rcoef))

# Subset relevant traits
#df_res_sub <- filter(df_res, rcoef < -0.08 | rcoef > 0.01)
df_res_sub <- df_res[which(df_res$p_value < 0.05/indep_trait),]

if (par_keep_not_sign){
  df_res_sub <- df_res
}

# Remove some categories
if(par_final_plot){
  df_res_sub$Category[which(!df_res_sub$Category %in% c("Blood composition","Cardiovascular",
                             "Reproductive","Pulmonary","Cancer",
                             "Lipoproteins","Morphological","Renal"))] <- "Not assigned"
}

### Plotting 

reg_plot_all_scaled <- ggforestplot::forestplot(
  df = df_res_sub,
  name = description,
  estimate = rcoef,
  pvalue = p_value,
  psignif = 0.05/indep_trait,
  xlab = expression(""*beta*" (95% CI)"),
  colour = Category,
  se = se,
) + 
  theme(axis.text.y = element_text(size = 9))
reg_plot_all_scaled

if(par_final_plot){
  reg_plot_all_scaled + #
    scale_color_manual(values=c("#de4f33", "#81ba20", "#bb292c","#FFE900","#fe9b00",
                                "#000000", "#1e90ff", "#499999", "#BF40BF")) +
    guides(color = guide_legend(override.aes = list(size = 1.5)))
}


################################################
### Save plot ##################################

export_file_name <- "reg_plot"

if(par_only_cont){
  export_file_name <- paste0(export_file_name, "_cont")
}

if(par_only_fact){
  export_file_name <- paste0(export_file_name, "_fact")
}

if(par_scale){
  export_file_name <- paste0(export_file_name, "_scaled")
}

if(par_corr_traits){
  export_file_name <- paste0(export_file_name, "_cor_traits")
}

if(par_keep_not_cor){
  export_file_name <- paste0(export_file_name, "_and_NC")
}

if(par_outlier_filter){
  export_file_name <- paste0(export_file_name, "_no_outliers")
}

if(par_rm_blood_cancer){
  export_file_name <- paste0(export_file_name, "_no_bc")
}

export_file_name <- paste0(export_file_name, ".pdf")

ggsave(file.path(export_folder, export_file_name),
       reg_plot_all_scaled, width=5, height=3, units="in", scale=3)


################################################
### Comparison depending on correction #########
################################################

### Compare one results vs another

if(par_compare){
  if(par_save_rds){
    saveRDS(df_res, file.path(data_folder, "df_res_1.rds"))
  }
  
  if(!par_save_rds){
    ### Compare first analysis with current analysis
    df_res_1 <- readRDS(file.path(data_folder, "df_res_1.rds"))
    # Keep same columns
    df_res <- select(df_res, colnames(df_res_1))
    df_res_comp <- left_join(df_res_1, df_res, by = "pheno")
    df_res_comp <- mutate(df_res_comp, 
                          pdiff = 2*pnorm(-abs((rcoef.x-rcoef.y)/sqrt(se.x**2+se.y**2)), 
                                          mean = 0, sd = 1))
    df_res_comp_sub <- df_res_comp[which(df_res_comp$pdiff < 0.05),] # 0.05/indep_trait
    df_res_comp_sub <- select(df_res_comp_sub, c("pheno","pdiff","rcoef.x","rcoef.y"))
    print(as.vector(df_res_comp_sub$pheno))
    
    ### Plotting difference 
    df_res_1$analysis <- par_analysis_names[1]
    df_res$analysis <- par_analysis_names[2]
    
    df_res_comp_plot <- rbind(df_res_1, df_res)
    df_res_comp_plot <- df_res_comp_plot[which(df_res_comp_plot$pheno %in% df_res_comp_sub$pheno),]
    
    df_res_comp_plot %>%
      mutate(pheno = fct_reorder(pheno, desc(rcoef))) %>%
      ggplot(aes(pheno, rcoef, colour = analysis)) +
      geom_point(position = position_dodge(width=0.75), 
                 shape = ifelse(df_res_comp_plot$p_value < 0.05/indep_trait, 16, 1),
                 size = 2) +
      geom_errorbar(aes(ymax = CI2, ymin = CI1, width = .1, color = analysis), 
                    position = position_dodge(width=0.75),
                    alpha = 0.5) +
      theme(legend.position="top",
            legend.title=element_blank())+
      coord_flip()
    
    ### Make the final clean plot
    if(par_final_corr){
      # Merge before and after
      df_res_comp_plot <- rbind(df_res_1, df_res)
      # Set list of traits to consider
      traits_to_look_at <- c("age_first_birth", "age_last_birth",
                             "menopause", "reproductive_lp", "menstruation",
                             "birth_weight_first_child","nb_birth","cholesterol", 
                             "HDL","LDL","TG")
      trait_to_correct_SES <- c("age_first_birth", "age_last_birth",
                                "menopause", "reproductive_lp", "menstruation",
                                "birth_weight_first_child","nb_birth")
      df_res_comp_plot <- df_res_comp_plot[which(df_res_comp_plot$pheno %in% traits_to_look_at),]
      # Set type of correction
      df_res_comp_plot$type_corr <- "Cholesterol medication"
      df_res_comp_plot$type_corr[which(df_res_comp_plot$pheno %in% trait_to_correct_SES)] <- "Socioeconomic status"
      # Add shape for significance
      df_res_comp_plot$shape_factor <- ifelse(df_res_comp_plot$pheno %in% df_res_comp_sub$pheno,
                                              ifelse(df_res_comp_plot$p_value < 0.05/indep_trait, "p < 0.05/135; p diff. < 0.05", "p > 0.05/135; p diff. < 0.05"),
                                              ifelse(df_res_comp_plot$p_value < 0.05/indep_trait, "p < 0.05/135; p diff. > 0.05", "p > 0.05/135; p diff. > 0.05"))
      df_res_comp_plot$shape_factor <- as.factor(df_res_comp_plot$shape_factor)
      
      ###Plot
      plot <- ggforestplot::forestplot(
        df = df_res_comp_plot,
        name = description,
        estimate = rcoef,
        pvalue = p_value,
        psignif = 0.05/indep_trait,
        xlab = expression(""*beta*" (95% CI)"),
        title = "",
        colour = analysis,
        shape = shape_factor,
        se = se
      )
      
      plot + facet_grid(vars(type_corr), space = "free", scales = "free") +
        theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
              strip.text.y = element_text(size = 8.5),
              axis.text.y=element_text(size=12),
              axis.title.x = element_text(size = 14),
              strip.background=element_rect(fill="grey"),
              legend.position="bottom",
              legend.box = "horizontal",
              legend.direction = "vertical",
              legend.title=element_blank(),
              legend.key.height=unit(0.5, "cm"),
              legend.text=element_text(size=13)) +
        scale_shape_manual(name="shape_factor",values = c(16,15,0,0)) +
        guides(color = guide_legend(override.aes = list(size = 2))) +
        guides(shape = guide_legend(override.aes = list(size = 2))) +
        scale_color_manual(values=c("#a30000", "#F07B15"))
    }
  }
}

################################################
### Correlation ################################
################################################

if (par_cor_plot){
  ################################################
  ### Correlation heatmap ########################
  
  if (par_only_cont){
    # Clustered heatmap
    cormat <- cor(df, use="pairwise.complete.obs")
    diag(cormat) <- NA
    
    library(pheatmap)
    mtscaled = as.matrix(cormat)
    H = pheatmap(mtscaled)
    hm <- pheatmap(mtscaled[H$tree_row$order,H$tree_col$order],cluster_rows = T,cluster_cols = F, 
                   main = "Pearson correlation heatmap", 
                   fontsize = 11, fontsize_row = 4, fontsize_col = 4)
    
    ggsave(file.path(export_folder, "cor_plot_cont.pdf"),
           hm, width=3, height=3, units="in", scale=3)
  }
  
  
  ################################################
  ### Summary df #################################
  
  if (par_only_cont){
    ### Regression coefficients on df
    x <- c("pheno", "mean", "median", "sd", "min", "max", "nb")
    df_summary <- data.frame(matrix(ncol = length(x), nrow = 0))
    colnames(df_summary) <- x
    
    for(i in 1:ncol(df)) {       # for-loop over columns
      reg_TLc_pheno <- lm(df$TLc ~ df[, i])
      df_summary[nrow(df_summary) + 1,] <- c(colnames(df)[i], 
                                             mean(df[, i], na.rm = T), 
                                             median(df[, i], na.rm = T),
                                             sd(df[, i], na.rm = T),
                                             min(df[, i], na.rm = T),
                                             max(df[, i], na.rm = T),
                                             sum(!is.na(df[, i]))
      )
    }
    
    # Put data as numeric and factor
    df_summary$pheno <- as.factor(df_summary$pheno)
    df_summary[,2:ncol(df_summary)] <- lapply(df_summary[,2:ncol(df_summary)], as.numeric)
  }
  
  
  
  ################################################
  ### Correlation with factors ###################
  
  
  if (par_only_fact){
    # Set 1 for NA's in diseases (because way we built the )
    # Only works if diseases are all together (should be updated)
    colid1 <- which(colnames(df) == "BC") # First disease
    colid2 <- which(colnames(df) == "menstruation") # Last disease
    df[,colid1:colid2][is.na(df[,colid1:colid2])] <- 1
    # Set 0 for blood cancer
    df$blood_cancer[is.na(df$blood_cancer)] <- 0
    
    # Remove sex specific factors
    df <- select(df, -c(BC, OC, endometriosis, menstruation, PC, balding))
    
    ### Rename columns
    rename_dict <- metadata[,c("pheno","description")]
    rename_dict <- rbind(rename_dict, c("TLc","TLc"))
    
    rename_df <- rename_dict %>%
      filter(pheno %in% colnames(df)) %>%
      group_by(pheno) %>%  
      arrange(factor(pheno, levels = colnames(df))) %>%
      ungroup()
    
    colnames(df) <- rename_df$description
    
    ### Compute cor. matrix
    cormat <- cor(df, use="pairwise.complete.obs")
    diag(cormat) <- NA
    
    mtscaled <- as.matrix(cormat)
    H <- pheatmap(mtscaled)
    hm <- pheatmap(mtscaled[H$tree_row$order,H$tree_col$order],cluster_rows = T,cluster_cols = F, 
                   main = "Pearson correlation heatmap", 
                   fontsize = 11, fontsize_row = 6, fontsize_col = 6)
    
    ggsave(file.path(export_folder, "cor_plot_fact.pdf"),
           hm, width=3, height=3, units="in", scale=3)
  }
  
  
  ################################################
  ### Control correlation with PCs ###############
  
  
  if (par_check_cor_PC){ # Attention depending on the names .x or without
    selvec <- c("Zadj_TS.x", "age.x", "sex", "blood_cancer_TF", "array")
    for(i in 1:40){
      selvec <- c(selvec, paste("PC", i, sep = ""))
    }
    
    df_control <- as.data.frame(fread(file.path(data_folder,  "complete_DF_noeid.txt"), header = T, select = selvec)) # Load data
    # Rename if needed
    df_control <- df_control %>% 
      rename(Zadj_TS = Zadj_TS.x) %>%
      rename(age = age.x)
  
    df_control <- df_control %>% add_column(age2 = (df_control$age)**2 , .after = "age") 
    
    df_control$array <- as.factor(df_control$array)
    
    df_control$TLc <- residuals(lm(df_control$Zadj_TS ~ df_control$age + 
                                     df_control$sex + df_control$age:df_control$sex + 
                                     df_control$age2 + df_control$age2:df_control$sex + df_control$array + 
                                     df_control$array:df_control$sex))
    
    ### Uncomment and remove from selection to add array in correlation plot ###
    df_control$array <- as.character(df_control$array)
    df_control$array[which(df_control$array == "UKBL")] <- 0
    df_control$array[which(df_control$array == "UKBB")] <- 1
    df_control$array <- as.integer(df_control$array)
    
    
    df_map <- select(df_control, -c(blood_cancer_TF, age2, Zadj_TS))
    
    # Clustered heatmap
    cormat <- cor(df_map, use="pairwise.complete.obs")
    diag(cormat) <- NA
    
    
    mtscaled = as.matrix(cormat)
    H = pheatmap(mtscaled)
    hm <- pheatmap(mtscaled[H$tree_row$order,H$tree_col$order],cluster_rows = T,cluster_cols = F, 
                   main = "Pearson correlation heatmap", 
                   fontsize = 13, fontsize_row = 9, fontsize_col = 9)
  }
  
}


