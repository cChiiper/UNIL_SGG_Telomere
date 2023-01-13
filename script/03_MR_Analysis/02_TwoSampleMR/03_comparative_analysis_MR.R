################################################################################
### Compare MR results with Codd and al. results                             ###
### Author: Samuel Moix                                                      ###
### Date: 11.07.2022                                                         ###
################################################################################


################################################
### Libraries ##################################
library(data.table)
library(ggplot2)
library(dplyr)
require(stringr)
require(forcats)
require(ggforestplot)

################################################
### Working directories ########################
data_folder = "../"
export_folder = "../"

################################################
### Load data ##################################
df_codd <- read.csv(file = file.path(data_folder,  "MR_compare/CompareWithCodd.csv"), sep = ';', header = TRUE)
df_res_MR<- read.csv(file = file.path(data_folder,  "MR_plots/TELOMERE_on_TRAIT_mr_results.txt"), sep = '\t', header = TRUE)
lm_res <- read.csv(file = file.path(data_folder,  "Phenotypes_Selection_Table/TL_metadata.csv"), sep = ',', header = TRUE)

##########
indep_trait <- 135

################################################
### Adjust data format #########################

### Adjust association beta 
df_codd[c('Beta', 'Beta_CI')] <- stringr::str_split_fixed(df_codd$Beta..95..CI., ' ', 2)
df_codd[c('Beta_CI_inf', 'Beta_CI_sup')] <- stringr::str_split_fixed(df_codd$Beta_CI, ' ', 2)
df_codd$Beta_CI_inf <- gsub("\\(", "", df_codd$Beta_CI_inf)
df_codd$Beta_CI_inf <- gsub(",", "", df_codd$Beta_CI_inf)
df_codd$Beta_CI_sup <- gsub("\\)", "", df_codd$Beta_CI_sup)

### Adjust MR beta
df_codd[c('Beta_MR', 'Beta_CI_MR')] <- stringr::str_split_fixed(df_codd$Beta..95..CI..MR, ' ', 2)
df_codd[c('Beta_CI_inf_MR', 'Beta_CI_sup_MR')] <- stringr::str_split_fixed(df_codd$Beta_CI_MR, ' ', 2)
df_codd$Beta_CI_inf_MR <- gsub("\\(", "", df_codd$Beta_CI_inf_MR)
df_codd$Beta_CI_inf_MR <- gsub(",", "", df_codd$Beta_CI_inf_MR)
df_codd$Beta_CI_sup_MR <- gsub("\\)", "", df_codd$Beta_CI_sup_MR)

### Adjust P values
df_codd$P.value <- gsub("<", "", df_codd$P.value)
df_codd$P.value.MR <- gsub("<", "", df_codd$P.value.MR)
df_codd$P.value <- gsub(",", ".", df_codd$P.value)
df_codd$P.value.MR <- gsub(",", ".", df_codd$P.value.MR)

### Select and rename columns of interest
df_codd <- df_codd[c('Trait','Field.Code', 'Beta', 'Beta_CI_inf','Beta_CI_sup', 'P.value',
          'Beta_MR','Beta_CI_inf_MR','Beta_CI_sup_MR','P.value.MR')]

df_codd <- df_codd %>% 
  mutate(se = (as.numeric(Beta_CI_sup_MR) - as.numeric(Beta_CI_inf_MR))/3.92)

### Match names
df_codd <- rename(df_codd, FieldID = Field.Code)
df_codd <- rename(df_codd, b = Beta_MR)
df_codd <- rename(df_codd, pval = P.value.MR)

### Add field ID's
lm_res[which(lm_res$pheno == "WHR"),"FieldID"] <- "48 � 49"
lm_res[which(lm_res$pheno == "GRIP"),"FieldID"] <- "47"
lm_res <- rename(lm_res, outcome = "MR_name")

df_res_MR <- merge(df_res_MR, lm_res[,c("pheno","outcome","FieldID")], by="outcome", all.x = T)
df_res_MR[which(df_res_MR$outcome == "WHR"), "FieldID"] <- "48 � 49"
df_res_MR[which(df_res_MR$outcome == "GRIP"), "FieldID"] <- "47"
df_res_MR <- filter(df_res_MR, !is.na(FieldID))
df_res_MR <- filter(df_res_MR, method == "Inverse variance weighted")


### Change numeric variables to numeric
df_codd[,-c(1,2)] <- data.frame(lapply(df_codd[,-c(1,2)], function(x) as.numeric(as.character(x)))) 

### Combine my results and theirs !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Codd_tc <- df_codd[c("FieldID","b","pval","se")]
Our_tc <- df_res_MR[c("FieldID","b","pval", "se")]

Codd_tc$author <- "Codd"
Our_tc$author <- "Auwerx_Moix"
test <- rbind(Codd_tc, Our_tc)

### Adding phenotype names
temp <- lm_res[,1:2]
temp$FieldID <- as.character(temp$FieldID)
#test <- test[which(test$FieldID %in% keep),]
test <- merge(test,temp,by="FieldID", all.x = FALSE)

### Put descriptive title
test <- merge(test, lm_res[,c("pheno","description")], by="pheno")
test <- test %>%
  select(-pheno) %>%
  rename(pheno = description)

# Sort
test <- arrange(test, (b))


################################################
### Plotting ###################################

### Keep only phenotypes sign. in at least one study
sig_beta <- test$pheno[which(test$pval < 0.05)]
test2 <- filter(test, pheno %in% sig_beta)
test2 <- test2 %>%
  filter(pheno %in% names(which(table(test2$pheno) != 1))) %>%
  mutate(pshape = ifelse(pval < 0.05/indep_trait,"multi","simple")) 

### Plot
ggforestplot::forestplot(
  df = test2,
  name = pheno,
  estimate = b,
  pvalue = pval,
  psignif = 0.05,
  xlab = expression("MR ivw "*alpha*" (CI95%)"),
  title = "",
  colour = author,
  shape = pshape,
  se = se
) + 
  ggplot2::scale_shape_manual(
    values = c(21L, 22L),
    labels = c("p < 0.05/135", "p < 0.05"),
    name = ""
  ) + 
  ggplot2::scale_color_manual(
    values = c("#a30000","#F07B15"),
    labels = c("Moix & Auwerx",  expression(paste("Codd ", italic("et al.")))),
    name = ""
  ) + 
  theme(axis.text.y = element_text(size = 9),
        legend.position="bottom") +            # x-axis text size
      guides(color = guide_legend(override.aes = list(size = 1.5))) +
      guides(shape = guide_legend(override.aes = list(size = 1.5)))





################################################################################
################### Compare disease outcome ####################################
################################################################################

################################################
### Load data ##################################
df_codd <- read.csv(file = file.path(data_folder,  "../disease_compare.csv"), sep = ';', header = TRUE)
df_res_MR<- read.csv(file = file.path(data_folder,  "../TELOMERE_on_TRAIT_mr_results.txt"), sep = '\t', header = TRUE)
lm_res <- read.csv(file = file.path(data_folder,  "../TL_metadata.csv"), sep = ',', header = TRUE)


################################################
### Adjust data format #########################

### Adjust association beta 
df_codd[c('Beta', 'Beta_CI')] <- stringr::str_split_fixed(df_codd$Hazard.Ratio..95..CI..Obs, ' ', 2)
df_codd[c('Beta_CI_inf', 'Beta_CI_sup')] <- stringr::str_split_fixed(df_codd$Beta_CI, ', ', 2)
df_codd$Beta_CI_inf <- gsub("\\(", "", df_codd$Beta_CI_inf)
df_codd$Beta_CI_sup <- gsub("\\)", "", df_codd$Beta_CI_sup)

### Adjust MR beta
df_codd[c('Beta_MR', 'Beta_CI_MR')] <- stringr::str_split_fixed(df_codd$Odds.Ratio..95..CI..MR, ' ', 2)
df_codd[c('Beta_CI_inf_MR', 'Beta_CI_sup_MR')] <- stringr::str_split_fixed(df_codd$Beta_CI_MR, ',', 2)
df_codd$Beta_CI_inf_MR <- gsub("\\(", "", df_codd$Beta_CI_inf_MR)
df_codd$Beta_CI_sup_MR <- gsub("\\)", "", df_codd$Beta_CI_sup_MR)

### Adjust pvalues
df_codd$P.value.Obs <- gsub(",", ".", df_codd$P.value.Obs)
df_codd$P.value.MR <- gsub(",", ".", df_codd$P.value.MR)


### Select and rename columns of interest
df_codd <- df_codd[c('pheno', 'Beta', 'Beta_CI_inf','Beta_CI_sup', 'P.value.Obs',
                     'Beta_MR','Beta_CI_inf_MR','Beta_CI_sup_MR','P.value.MR')]

### Change numeric variables to numeric
df_codd[,-1] <- data.frame(lapply(df_codd[,-1], function(x) as.numeric(as.character(x)))) 

### Put as beta (not odds ratio)
df_codd$Beta_MR <- log(df_codd$Beta_MR)
df_codd$Beta_CI_inf_MR <- log(df_codd$Beta_CI_inf_MR)
df_codd$Beta_CI_sup_MR <- log(df_codd$Beta_CI_sup_MR)


### Add MR standard error
df_codd <- df_codd %>% 
  mutate(se = ((as.numeric(Beta_CI_sup_MR) - as.numeric(Beta_CI_inf_MR))/3.92))

### Match names
df_codd <- rename(df_codd, b = Beta_MR)
df_codd <- rename(df_codd, pval = P.value.MR)



### Adjust our results
df_res_MR <- filter(df_res_MR, method == "Inverse variance weighted")
df_res_MR <- rename(df_res_MR, MR_name = outcome)
df_res_MR <- merge(df_res_MR, lm_res[,c("pheno","MR_name")], by="MR_name")
disease_names <- lm_res$pheno[which(!is.na(lm_res$Disease_cases))]
df_res_MR <- df_res_MR[which(df_res_MR$pheno %in% disease_names),] 

### Get data of interest
Codd_tc <- df_codd[c("pheno","b","pval","se")]
Our_tc <- df_res_MR[c("pheno","b","pval", "se")]

Codd_tc <- Codd_tc[which(!is.na(Codd_tc$pheno)),]
Our_tc <- Our_tc[which(Our_tc$pheno %in% Codd_tc$pheno),]

Codd_tc$author <- "Codd"
Our_tc$author <- "Auwerx_Moix"

test <- rbind(Codd_tc, Our_tc)



### Add shape for plotting
test <- test %>%
  mutate(pshape = ifelse(pval < 0.05/indep_trait,"multi","simple")) 

test <- arrange(test, (b))

### Plot
ggforestplot::forestplot(
  df = test,
  name = pheno,
  estimate = b,
  pvalue = pval,
  psignif = 0.05,
  xlab = "MR ivw OR or HR(CI95%)",
  title = "",
  colour = author,
  shape = pshape,
  se = se,
  logodds = TRUE
) + 
  ggplot2::scale_shape_manual(
    values = c(21L, 22L),
    labels = c(expression(paste("p < 3.7 �", 10^{-4})), "p < 0.05"),
    name = ""
  ) + 
  ggplot2::scale_color_manual(
    values = c("#a30000","#F07B15"),
    labels = c("Moix & Auwerx", "Codd et al."),
    name = ""
  ) + 
  theme(axis.text.y = element_text(size = 8),
        legend.position="bottom")

################################################
### As scatter plot (because we can't really get OR and HR)

test2 <- merge(Our_tc, Codd_tc, by = "pheno", suffixes = c("_our","_codd"))
test2 <- test2 %>%
  select(-c("author_our","author_codd")) %>%
  mutate(CI_sup_our = b_our + 1.96*se_our) %>%
  mutate(CI_inf_our = b_our - 1.96*se_our) %>%
  mutate(CI_sup_codd = b_codd + 1.96*se_codd) %>%
  mutate(CI_inf_codd = b_codd - 1.96*se_codd) %>%
  filter(pval_our < 0.05 | pval_codd < 0.05) %>%
  mutate(shape_p = ifelse((pval_our < 0.05/indep_trait & pval_codd < 0.05/indep_trait),
         21,22)) %>%
  mutate(fill_p = ifelse(pval_our < 0.05/indep_trait,
         ifelse(pval_codd < 0.05/indep_trait,"black","#a30000"),
         ifelse(pval_codd < 0.05/indep_trait,"#F07B15","grey"))) %>%
  mutate(b_codd = exp(b_codd)) %>%
  mutate(CI_inf_codd = exp(CI_inf_codd)) %>%
  mutate(CI_sup_codd = exp(CI_sup_codd))
  

comp_plot <- ggplot(test2, aes(x=b_our, y=b_codd,fill = fill_p)) +
  geom_errorbar(aes(ymax = CI_sup_codd, ymin = CI_inf_codd, width = 0), color = "#F07B15", alpha = 0.7) +
  geom_errorbar(aes(xmax = CI_sup_our, xmin = CI_inf_our, width = 0), color = "#a30000", alpha = 0.7) +
  geom_point(#shape=21,
             shape=test2$shape_p, 
             #fill =test2$fill_p, 
             size = 1.5, 
             colour = "black") + # Show dots
  geom_abline(slope = 1, color = "#685044", alpha = 0.2, linetype = "dashed") +
  geom_vline(xintercept = 0, color = "black", alpha = 0.4, linetype = "dashed") +
  geom_hline(yintercept = 1, color = "black", alpha = 0.4, linetype = "dashed") +
  ggrepel::geom_text_repel(
    label=test2$pheno, 
    size = 3,
    #fontface = "bold"
  )+
  xlab(expression("Moix & Auwerx : MR ivw "*alpha*" (CI95%)")) + 
  ylab(expression(paste("Codd ", italic("et al.")," : MR ivw OR (95% CI)"))) +
  ggplot2::coord_cartesian(ylim = c(0.6,2.15),xlim = c(-0.2,0.2))


comp_plot + 
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.title=element_text(size=11),
        legend.position="bottom",
        legend.title=element_blank()) +
  scale_fill_manual(values=c("#a30000","#F07B15", "black", "grey"),
                    labels=c(
                      "Moix & Auwerx p < 0.05/135",
                      expression(paste("Codd ", italic("et al.")," p < 0.05/135")),
                      "Both p < 0.05/135",
                      "Both/Either p < 0.05"
                    ))+
  guides(fill = guide_legend(override.aes = list(shape=21, size = 2.5)))

### Put descriptive title
test2 <- merge(test2, lm_res[,c("pheno","description")], by="pheno")
test2 <- test2 %>%
  select(-pheno) %>%
  rename(pheno = description)

