################################################################################
### Observation against MR results                                           ###
### Author: Samuel Moix                                                      ###
### Date: 02.08.2022                                                         ###
################################################################################


################################################
### Libraries ##################################
library(dplyr)
library(ggplot2)
require(ggrepel)
library(ggforestplot) # Forest plots

################################################
### Working directories ########################

data_folder = "../"
export_folder = "../"
# File names
MR_results_files <- c("/TELOMERE_on_TRAIT_mr_results.txt", "MR_plots/TRAIT_on_TELOMERE_mr_results.txt")


################################################
### Parameters  ################################
par_only_sign <- FALSE ### For now keep FALSE for test !!!!!!
indep_trait <- 135

################################################
### Load Data and prepare data #################

### Load results
#df_res <- ADD file written to avoid relaunching all everytime !!!!!!!!!!!!
df_res <- readRDS(file.path("../df_res_1.rds"))

### Load metadata
df_code <- read.csv(file = file.path(data_folder,  "Phenotypes_Selection_Table/TL_metadata.csv"), sep = ',', header = TRUE)
disease_traits <- df_code$pheno[which(!is.na(df_code$Disease_cases))]

par_foward <- TRUE
for (file_MRR in MR_results_files) {
  ### Load MR results files
  df_res_MR <- read.table(file = file.path(data_folder,  file_MRR), sep = '\t', header = TRUE)
  # Choose preferred MR method
  # Inverse variance weighted     Simple mode     Weighted mode     MR Egger
  # Weighted median
  df_res_MR <- filter(df_res_MR, method == "Inverse variance weighted")
  
  ### Add CI to MR data frame
  df_res_MR <- df_res_MR %>% 
    mutate(Beta_CI_inf_MR = b - 1.96*se) %>%
    mutate(Beta_CI_sup_MR = b + 1.96*se)
  
  ### Select column of interest
  if(par_foward){
    df_code <- rename(df_code, outcome = "MR_name")
    df_res_MR <- merge(df_res_MR, df_code[,c("pheno","outcome")],by="outcome", all.x = TRUE)
    df_res_MR <- df_res_MR[,c("pheno","outcome","b","Beta_CI_inf_MR","Beta_CI_sup_MR","pval","se")]
    
  }else{
    df_code <- rename(df_code, exposure = "outcome")
    df_res_MR <- merge(df_res_MR, df_code[,c("pheno","exposure")],by="exposure", all.x = TRUE)
    df_res_MR <- df_res_MR[,c("pheno","exposure","b","Beta_CI_inf_MR","Beta_CI_sup_MR","pval","se")]
  }
  
  ### Add linear regression coefficients
  df_both <- merge(df_res_MR, df_res[,c("pheno","rcoef","CI1","CI2","p_value")],by="pheno", all.x = TRUE)
  if(par_only_sign){
    df_both <- filter(df_both, pval < 0.05/indep_trait)
  }
  
  ################################################
  ### Plotting ###################################
  
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
    ylab(paste("MR beta (ivw)", ifelse(par_foward, "forward", "reverse")))
  
  if(par_foward){
    obs_vs_mr_plot_forward <- obs_vs_mr_plot
    df_res_MRF <- df_res_MR
  }else{
    obs_vs_mr_plot_reverse <- obs_vs_mr_plot
  }
  rm(obs_vs_mr_plot)
  par_foward <- FALSE
}

# obs_vs_mr_plot_forward
# obs_vs_mr_plot_reverse

################################################################################
### MR results plotting                                                      ###
################################################################################

### Merge forward data to reverse and association
df_res_MR$type_an <- "reverse"
df_res_MRF$type_an <- "forward"
### Rename observation columns to match MR names for plotting (comment last for all)
df_res_OBS <- df_res %>% select("rcoef","CI1","CI2","p_value","pheno","se") %>%
  rename(b = rcoef, Beta_CI_inf_MR = CI1, Beta_CI_sup_MR = CI2, pval = p_value) %>%
  filter(pheno %in% df_res_MR$pheno | pheno %in% df_res_MRF$pheno)
df_res_OBS$type_an <- "b"

### Merge forward and reverse results
plot_df_full <- rbind(select(df_res_MR, 
                     c("b","Beta_CI_inf_MR","Beta_CI_sup_MR","pval","pheno","type_an","se")),
              select(df_res_MRF, 
                     c("b","Beta_CI_inf_MR","Beta_CI_sup_MR","pval","pheno","type_an","se")), 
              df_res_OBS)

### Recode direction
plot_df_full$type_an <- as.factor(plot_df_full$type_an)
plot_df_full$type_an <- relevel(plot_df_full$type_an , "forward")
plot_df_full$type_an <- relevel(plot_df_full$type_an , "reverse")

### Add significance for plotting
plot_df_full <- plot_df_full %>%
  mutate(pshape = ifelse(pval < 0.05/indep_trait,"multi","simple")) 


### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### 

### Reset plotting dataframe
test <- plot_df_full

### Parameters 
direction <- "forward"
direction <- "reverse"
par_trait <- TRUE
par_trait <- FALSE


### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### 

### Setting x-axis
if(direction == "forward"){
  if(par_trait){
    par_xlim_plot <- c(-0.2, 0.2)
  } else {
    par_xlim_plot <- c(-0.3, 0.4)
  }
} else{
  if(par_trait){
    par_xlim_plot <- c(-0.21, 0.28)
  } else {
    par_xlim_plot <- c(-0.28, 0.1)
  }
}

### Keep only traits of interest
# Select new traits
if(direction == "forward"){
  select_traits <- df_code$pheno[which(df_code$Codd == 2)]
} else{
  select_traits <- df_code$pheno
}


# Direction and significance
select_traits <- select_traits[which(select_traits %in% 
                                       unique(c(test$pheno[which(test$pval < 0.05 & 
                                                                   test$type_an == direction)])))]

### Subset dataframe
test <- test %>% 
  filter(pheno %in% select_traits)

# Filter diseases 
if(par_trait){
  test <- filter(test, !pheno %in% disease_traits)
}else{
  test <- filter(test, pheno %in% disease_traits) 
}


### Put descriptive title
test <- merge(test, df_code[,c("pheno","description")], by="pheno")
test <- test %>%
  select(-pheno) %>%
  rename(pheno = description)

### Reorder by beta
test_order <- test %>%
  filter(type_an == direction) %>%
  arrange(b)

test <- test %>%
  group_by(pheno) %>%  
  arrange(factor(pheno, levels = test_order$pheno)) %>%
  ungroup()

### Plot
ggforestplot::forestplot(
  df = test,
  name = pheno,
  estimate = b,
  pvalue = pval,
  logodds = FALSE,
  psignif = 0.05,
  xlab = expression(""*alpha*" or "*beta*" (CI95%)"),
  title = "",
  colour = type_an,
  shape = pshape,
  se = se
) + 
  theme(axis.text.y = element_text(size = 10), # x-axis text size
        legend.position="bottom",
        legend.title=element_blank(),
        legend.box = "horizontal",
        legend.direction = "vertical", #) +
        legend.justification = c(1, 0)) +            
  ggplot2::coord_cartesian(xlim = par_xlim_plot) +
  ggplot2::scale_color_manual(
    values = c("#004777","#a30000","#8EA604"),
    labels = c("MR (ivw) reverse", "MR (ivw) forward", expression("Observational ("*beta*")")),
    name = ""
  ) +
  ggplot2::scale_shape_manual(
    values = c(21L, 22L),
    labels = c("p < 0.05/135", "p < 0.05"), # expression(paste("p < 3.7 �", 10^{-4}))
    name = ""
  ) +
  guides(color = guide_legend(override.aes = list(size = 1.5), ncol=2)) +
  guides(shape = guide_legend(override.aes = list(size = 1.5)))



### Reverse+Forward plot
df_res_MR_FR <- merge(df_res_MR, df_res_MRF, by = "pheno")
df_res_MR_FR <- df_res_MR_FR %>% mutate(b = b.x + b.y, 
                                        se.x = (Beta_CI_sup_MR.x-b.x)/1.96, 
                                        se.y = (Beta_CI_sup_MR.y-b.y)/1.96) %>%
  mutate(se = sqrt(se.x**2+se.y**2)) %>% # Calculate SE from CI
  mutate(Beta_CI_inf_MR = b - 1.96*se,
         Beta_CI_sup_MR = b + 1.96*se) %>%
  rename(se_mr = se) %>%
  select("pheno","exposure","b","Beta_CI_inf_MR","Beta_CI_sup_MR", "se_mr")

### Add linear regression coefficients
df_both <- merge(df_res_MR_FR, df_res[,c("pheno","rcoef","CI1","CI2","p_value", "se")],by="pheno", all.x = TRUE)
# Add pval significantly diff
df_both <- df_both %>%  
  mutate(pval = 2*pnorm(-abs((rcoef - b)/sqrt(se**2+se_mr**2)), mean = 0, sd = 1))

if(par_only_sign){
  df_both <- filter(df_both, pval < 0.05/indep_trait)
}

################################################
### Plotting ###################################

# Put description
df_both <- merge(df_both, df_code[,c("pheno","description")], by="pheno")
df_both <- df_both %>%
  select(-pheno) %>%
  rename(pheno = description)

obs_vs_mr_plot_FR <- ggplot(df_both, aes(x=rcoef, y=b)) +
  geom_errorbar(aes(ymax = Beta_CI_sup_MR, ymin = Beta_CI_inf_MR, width = 0), color = "#D18EF0", alpha = 0.7) +
  geom_errorbar(aes(xmax = CI2, xmin = CI1, width = 0), color = "#8EA604", alpha = 0.7) +
  geom_point(shape=16) + # Show dots
  geom_abline(slope = 1, color = "grey", alpha = 0.5) +
  ggrepel::geom_text_repel(
    label=df_both$pheno, 
    size = 2
  ) +
  xlab(expression(""*beta*" (95% CI)")) +
  ylab(expression(paste(""*alpha*""["forward + reverse"]," (95% CI)")))



################################################
### Adjust Plotting ############################

# obs_vs_mr_plot_FR + 
#   geom_smooth(method='lm')
# obs_vs_mr_plot_forward  + 
#   geom_smooth(method='lm')
# obs_vs_mr_plot_reverse + 
#   geom_smooth(method='lm')

require(ggforce)
mylim <- 0.1
obs_vs_mr_plot_forward + 
  ggforce::facet_zoom(ylim = c(-mylim,0.05), xlim = c(-mylim, 0.05))
mylim <- 0.1
obs_vs_mr_plot_reverse +
  ggforce::facet_zoom(ylim = c(-mylim,mylim), xlim = c(-mylim, mylim))
mylim <- 0.2
obs_vs_mr_plot_FR +
  ggforce::facet_zoom(ylim = c(-mylim,mylim), xlim = c(-mylim, mylim))


obs_vs_mr_plot_FR +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.9) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
