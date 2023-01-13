#################################################################################################################
### Plotting of Telomere length with age, sex and cancer                                                      ###
### Author: Samuel Moix                                                                                       ###
### Date: 24.03.2022                                                                                          ###
#################################################################################################################


################################################
### Libraries ##################################
library(data.table) # Read data 
library(ggplot2) # Plotting basics
library(ggExtra) # Extra plotting XP
library(ggpubr) # Plotting regression stuff

################################################
### Working directories ########################
data_folder = "../"


################################################
### Data loading and filtering #################
df <- as.data.frame(fread(file.path(data_folder,  "filtered_DF.txt"), header = T)) 
df$Zadj_TS <- scale(df$Zadj_TS)
df$sex <- as.factor(df$sex) # Set sex as factor


################################################
### Add corrected phenotypes ###################

df$TLc_age = residuals(lm(df$Zadj_TS ~ df$age))
df$TLc_age_sex <- residuals(lm(df$Zadj_TS ~ df$age * df$sex))


################################################
### Plotting ###################################

### Boxplot TL ~ age regression line ########### ########### ########### ########### ########### ###########

lm_plot_age <- ggplot(df, aes(x = age, y = Zadj_TS)) +
  geom_point(alpha = 0) + # Has to be set to fit regression line alpha = 0  to not show points
  geom_boxplot(data = df, aes(x = age, y = Zadj_TS, group = factor(age)), 
               alpha = 0.7, lwd = 0.7, outlier.alpha = 0.05, fill = "#E3E4B4") + 
  labs(title = "Telomere length vs age",  # Editing labs
       x = "Age [years]", 
       y = "Z-adjusted T/S log [-]",
       fill = "Sex") +
  scale_x_continuous(breaks = seq(40, 70, by = 2))+ # Editing x scale
  geom_smooth(method = "lm", aes(x = age, y = Zadj_TS), formula = y ~ x, color = "#CB0000") + # Draw lm
  stat_regline_equation( # Add lm equation
    aes(x = age, y = Zadj_TS, label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
    formula = y ~ x, show.legend = FALSE, color = "#CB0000") +
  theme(legend.position="left",  # Change layout
        plot.title = element_text(size=20),
        axis.text = element_text(size = 12),
        axis.title=element_text(size=14),
        legend.key.height= unit(1, 'cm'), 
        legend.key.width= unit(1, 'cm'), 
        legend.title = element_text(face = "bold") 
  ) +
  ylim(-5,5)

ggMarginal(lm_plot_age, margins = "both", type = "density", size = 10, fill = "#E3E4B4")


# plot(df$age, df$TLc_age, ylab="Residuals", xlab="Age", main="TL ~ age") 
# abline(0, 0)  
  
### Boxplot TL ~ age by sex regression line ########### ########### ########### ########### ########### ####

lm_plot_age <- ggplot(df, aes(x = age, y = Zadj_TS, group = factor(sex), color = factor(sex))) +
  geom_point(alpha = 0) + # Has to be set to fit regression line alpha = 0  to not show points
  geom_boxplot(data = df, aes(x = age, y = Zadj_TS, group = interaction(age, factor(sex)), fill = factor(sex)), 
               alpha = 0.7, lwd = 0.7, outlier.alpha = 0.05) +
  labs(title = "Telomere length vs age by sex",  # Editing labs
       x = "Age [years]", 
       y = "Z-adjusted T/S log [-]",
       fill = "Sex") +
  scale_x_continuous(breaks = seq(40, 70, by = 2))+ # Editing x scale
  geom_smooth(method = "lm", aes(x = age, y = Zadj_TS, group=sex, color=sex), formula = y ~ x) + # Draw lm
  stat_regline_equation( # Add lm equation
    aes(x = age, y = Zadj_TS, color = sex, label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
    formula = y ~ x, show.legend = FALSE) +
  scale_fill_manual(values=c("#94CDE1", "#EA6262"), labels = c("M", "F")) + #Set colors
  scale_color_manual(values=c("#2471A3","#C0392B")) +
  guides(colour = "none") + # Takes legend of geom_point away
  theme(legend.position="left",  # Change layout
        plot.title = element_text(size=20),
        axis.text = element_text(size = 12),
        axis.title=element_text(size=14),
        legend.key.height= unit(1, 'cm'), 
        legend.key.width= unit(1, 'cm'), 
        legend.title = element_text(face = "bold") 
  ) +
  ylim(-5,5)

# Add marginal plots (distribution)
ggMarginal(lm_plot_age, margins = "both", type = "density", groupColour = T, groupFill = T, size = 10)


### Density curves plot by age group

df$age_group <- NA
df$age_group[which(df$age >= 40 & df$age <= 45)] <- "40-45"
df$age_group[which(df$age >= 46 & df$age <= 50)] <- "46-50"
df$age_group[which(df$age >= 51 & df$age <= 55)] <- "51-55"
df$age_group[which(df$age >= 56 & df$age <= 60)] <- "56-60"
df$age_group[which(df$age >= 61 & df$age <= 65)] <- "61-65"
df$age_group[which(df$age >= 66 & df$age <= 70)] <- "66-70"


library(plyr)

mu <- ddply(df, "age_group", summarise, grp.mean=mean(Zadj_TS))
ggplot(data=df, aes(x=Zadj_TS, group=age_group, fill=age_group)) +
  geom_density(adjust=1.5, alpha=.2) +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=age_group),
             linetype="dashed") +
  xlim(-4,4)


### Boxplot TL ~ age by presence of blood malignancy ########### ########### ########### ########### ######

df$blood_cancer_TF[which(df$blood_cancer_TF == TRUE)] <- 1
df$blood_cancer_TF[which(df$blood_cancer_TF == FALSE)] <- 2
df$blood_cancer_TF <- factor(df$blood_cancer_TF)

lm_plot_age <- ggplot(df, aes(x = age, y = Zadj_TS, group = factor(blood_cancer_TF), color = factor(blood_cancer_TF))) +
  geom_point(alpha = 0) + # Has to be set to fit regression line alpha = 0  to not show points
  geom_boxplot(data = df, aes(x = age, y = Zadj_TS, group = interaction(age, factor(blood_cancer_TF)), fill = factor(blood_cancer_TF)), 
               alpha = 0.7, lwd = 0.7, outlier.alpha = 0.05) +
  labs(title = "Telomere length vs age by presence of blood malignancy",  # Editing labs
       x = "Age [years]", 
       y = "Z-adjusted T/S log [-]",
       fill = "Blood malignancy") +
  scale_x_continuous(breaks = seq(40, 70, by = 2))+ # Editing x scale
  geom_smooth(method = "lm", aes(x = age, y = Zadj_TS, group=blood_cancer_TF, color=blood_cancer_TF), formula = y ~ x) + # Draw lm
  stat_regline_equation( # Add lm equation
    aes(x = age, y = Zadj_TS, color = blood_cancer_TF, label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
    formula = y ~ x, show.legend = FALSE) +
  scale_fill_manual(values=c("#A3CF6A", "#F0B642"), labels = c("W/o BM", "BM")) + #Set colors
  scale_color_manual(values=c("#628339","#BF8A1F")) +
  guides(colour = "none") + # Takes legend of geom_point away
  theme(legend.position="left",  # Change layout
        plot.title = element_text(size=20),
        axis.text = element_text(size = 12),
        axis.title=element_text(size=14),
        legend.key.height= unit(1, 'cm'), 
        legend.key.width= unit(1, 'cm'), 
        legend.title = element_text(face = "bold") 
  ) +
  ylim(-5,5)

lm_plot_age
# Add marginal plots (distribution)
ggMarginal(lm_plot_age, margins = "both", type = "density", groupColour = T, groupFill = T, size = 10)



### Boxplot TL ~ blood cancer ##################
library(grid)



p1 <- ggplot(df, aes(x=blood_cancer_TF, y=Zadj_TS)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) + 
  labs(title = "", x = "Blood cancer", y = "Z-adjusted T/S log [-]") +
  ylim(-1,1)

p2 <- ggplot(df, aes(x=blood_cancer_TF, y=TLc_age)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) + 
  labs(title = "", x = "Blood cancer", y = "T/S corrected for age [-]")+
  ylim(-1,1)

require(gridExtra)
gridExtra::grid.arrange(p1, p2, ncol = 2,
             top=textGrob("TL distribution in function of presence of BM",
                                           gp=gpar(fontsize=15,font=3)))


### Welch two sample t-test for group with and without cancer
t.test(Zadj_TS ~ blood_cancer, data = df) 
t.test(TLc_age ~ blood_cancer, data = df)




################################################
### Statistics #################################

### TL ~ age all regression ####################
TL_age_model = lm(df$Zadj_TS ~ df$age)

# residu <- residuals(TL_age_model)
# par(mfrow=c(1,3))
# plot(df$age, residu) ; abline(lm(residu~df$age))
# plot(df$age, abs(residu)) ; abline(lm(abs(residu)~df$age))
# qqnorm(residu) ; qqline(residu)
# cor.test(abs(residu), df$age, method="spearman") # p-value = 0.03148 but already log transformed

par(mfrow=c(2,2)) 
plot(TL_age_model)


### TL ~ age Male Female regression ############

df$age = scale(df$age)
lm_male = lm(df$Zadj_TS[which(df$sex == 1)] ~ df$age[which(df$sex == 1)])
lm_female = lm(df$Zadj_TS[which(df$sex == 2)] ~ df$age[which(df$sex == 2)])

print(lm_male)
print(lm_female)

### Pearson's cor test between age and TL #####
cor.test(df$age, df$Zadj_TS)

### Average TL between male and female


ggplot(df, aes(x=as.factor(sex), y=Zadj_TS)) + 
  geom_boxplot(fill="slateblue", alpha=0.2,outlier.shape = NA) + 
  coord_cartesian(ylim = c(-3,3)) +
  xlab("Sex") +
  ylab("TL") +
  scale_x_discrete(labels = c('M','F'))


t.test(Zadj_TS ~ sex, data = df) 

################################################
### Analysis of lifestyle confounders ##########
df <- a
par_vege <- TRUE

if(par_vege){
  df <- df %>%
    filter(!is.na(alcool_frq)) %>%
    filter(!is.na(alcool_weekly)) %>%
    filter(!is.na(smoking_status)) %>%
    filter(!is.na(fruit)) %>%
    filter(!is.na(beef)) %>%
    filter(!is.na(vegetable))
}

if(!par_vege){
  df <- df %>%
    filter(!is.na(testosterone)) %>%
    filter(!is.na(SHBG))
}

df$Zadj_TS <- scale(df$Zadj_TS)
df$age <- scale(df$age)


if(par_vege){
  ### Add corrected trait
  df$TLc <- residuals(lm(df$Zadj_TS ~ df$alcool_frq + df$alcool_weekly +
                                df$smoking_status + df$fruit + df$beef +
                                df$vegetable))
}

if(!par_vege){
  df$TLc <- residuals(lm(df$Zadj_TS ~ df$testosterone + df$SHBG))
}

df$TLc <- scale(df$TLc)


df_M_ls <- filter(df, sex == 1) %>% select(-sex)
df_F_ls <- filter(df, sex == 2) %>% select(-sex)


### Regression
# Non-adjusted regression
male_reg <- lm(df_M_ls$Zadj_TS ~ df_M_ls$age)
male_b <- male_reg$coefficients[[2]]
male_se <- coef(summary(male_reg))[2, 2]

female_reg <- lm(df_F_ls$Zadj_TS ~ df_F_ls$age)
female_b <- female_reg$coefficients[[2]]
female_se <- coef(summary(female_reg))[2, 2]

pdiff_1 <- 2*pnorm(-abs((male_b-female_b)/sqrt(male_se**2+female_se**2)),mean = 0, sd = 1)
print(paste("Male b:", round(male_b,3), "male se:", round(male_se,4),
            "Female b:", round(female_b,3), "female se:", round(female_se,4),
            "P_diff:", pdiff_1), sep=" ")


# Adjusted regression
male_reg_c <- lm(df_M_ls$TLc ~ df_M_ls$age)
male_b_c <- male_reg_c$coefficients[[2]]
male_se_c <- coef(summary(male_reg_c))[2, 2]

female_reg_c <- lm(df_F_ls$TLc ~ df_F_ls$age)
female_b_c <- female_reg_c$coefficients[[2]]
female_se_c <- coef(summary(female_reg_c))[2, 2]

pdiff_2 <- 2*pnorm(-abs((male_b_c-female_b_c)/sqrt(male_se_c**2+female_se_c**2)),mean = 0, sd = 1)
print(paste("Male b:", round(male_b_c,3), "male se:", round(male_se_c,4),
            "Female b:", round(female_b_c,3), "female se:", round(female_se_c,4),
            "P_diff:", pdiff_2), sep=" ")

table_1 <- data.frame(male_b  = c(round(male_b,3),round(male_b_c,3)),
                      male_se = c(round(male_se,4),round(male_se_c,4)),
                      female_b = c(round(female_b,3),round(female_b_c,3)),
                      female_se = c(round(female_se,4),round(female_se_c,4)),
                      pdiff = c(pdiff_1,pdiff_2))
rownames(table_1) <- c("TL","TLc")


### Difference between adjusted and non-adjusted
pdiff_m <- 2*pnorm(-abs((male_b_c-male_b)/sqrt(male_se_c**2+male_se**2)),mean = 0, sd = 1)
pdiff_f <- 2*pnorm(-abs((female_b_c-female_b)/sqrt(female_se_c**2+female_se**2)),mean = 0, sd = 1)

table_1 ; print(paste("pdiff males:", pdiff_m)) ; print(paste("pdiff females:", pdiff_f)) ; print(paste("N:", nrow((df))))

