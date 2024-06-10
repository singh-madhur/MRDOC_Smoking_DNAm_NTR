## Examine the associations between Smoking Phenotypes and Smoking PRSs from GSCAN2

## clear workspace
rm(list = ls()); gc()

## Directories #####################

baseDir <- "/base_dir/" 
prsDir  <- paste0(baseDir,"data/PRS_Smk/GSCAN2/")
phenDir <- paste0(baseDir,"data/phenos/")
outDir  <- paste0(baseDir,"out/pheno_descr/")


# Libraries
library(tidyverse)
library(polycor)      # for polyserial correlation of binary outcome
library(ggrepel)
library(gee)          # for R2 for continuous traits
library(fastDummies)  # for creating dummy variables for R2 calculations


# Phenotypic Data

## Fisnumber	= Person ID. Use this ID to merge to genotype data (polygenic scores)
load(file = paste0(phenDir,"twin_covariates_long_data.Rdata"))
dim(long_twin_data)
glimpse(head(long_twin_data))
table(long_twin_data$zygo, useNA = "ifany")


## Smoking pheno variables
smk_phenos <- c("everSmk","currentSmk","formerSmk","smkCess","hFTNDmax_o",
                "cotin","packyears","nyearstop","ncigday", "hFTNDmax")




# Standardize continuous variables ---------------------------------------------

psych::describe(long_twin_data$DNA_BloodSampling_Age)

long_twin_data <- long_twin_data |> 
  mutate(age_sc = scale(DNA_BloodSampling_Age),
         cotin_sc = scale(cotin),
         packyears_sc = scale(packyears),
         nyearstop_sc = scale(nyearstop),
         ncigday_sc = scale(ncigday),
         hFTNDmax_sc = scale(hFTNDmax)
  )




# Smoking Initiation PRS --------------------------------------------------------

## Source: GSCAN2. Saunders et al, 2022

## Primary trait: Binary - Ever Smoker (= current or former) / Never Smoker
PRS_SmokingInitiation <- haven::read_sav(paste0(prsDir,"NTR-DSR-2962_Smoking_SmokingInitiation_PMID36477530_MRG16_PedMergedWithScores.sav"))

dim(PRS_SmokingInitiation)
glimpse(head(PRS_SmokingInitiation))

## Check the IDs across datasets
table(is.element(long_twin_data$fisnumber, PRS_SmokingInitiation$FISnumber))

## `FISnumber` is the same variable as long_twin_data$fisnumber
PRS_SmokingInitiation$fisnumber <- PRS_SmokingInitiation$FISnumber
PRS_SmokingInitiation$FISnumber <- NULL

## Add PRS to the main dataset
phenos_prs_smk_init <- long_twin_data |> 
  left_join(PRS_SmokingInitiation, by = "fisnumber")

dim(phenos_prs_smk_init)

## Filter EUR ancestry ---------------------------------------------------------

phenos_prs_smk_init |> 
  count(EUR_1KG_Outlier)

twin_long_prs_EUR <- phenos_prs_smk_init |> 
  filter(EUR_1KG_Outlier == 0)

twin_long_prs_EUR |> 
  count(zygo)

twin_long_prs_EUR |> 
  count(Smoking)

twin_long_prs_EUR |> 
  count(id_sex) |> 
  mutate(prop = n / sum(n))

psych::describe(twin_long_prs_EUR$DNA_BloodSampling_Age)



## Correlations -----------------------------------------------------------------

## result dataframe
prs <- twin_long_prs_EUR |> 
  select(contains("SCORE_Smoking_SmokingInitiation_MRG16_LDp1")) |> 
  names()

prs_cor_smk_int <- as.data.frame(cbind(prs, 
                                       "causal_fr" = NA_real_,
                                       "r" = NA_real_,
                                       "lb" = NA_real_,
                                       "ub" = NA_real_))

## polyserial correlation of a continuous variable (prs) and a binary variable (everSmk)

twin_long_prs_EUR <- data.frame(twin_long_prs_EUR) ## polycor cannot handle tibble dataframes

for (i in 1:length(prs)) {
  prs_cor_smk_int$prs[i] <- prs[i]
  prs_cor_smk_int$causal_fr[i] <- substr(prs[i], start = 1, stop = 6)
  test <- polyserial(twin_long_prs_EUR[,prs[i]],
                     twin_long_prs_EUR[,"everSmk"], 
                     std.err = T)
  prs_cor_smk_int$r[i] <- test$rho
  se <- sqrt(test$var)
  prs_cor_smk_int$ub[i] <- test$rho + qnorm(p = 0.975)*se
  prs_cor_smk_int$lb[i] <- test$rho - qnorm(p = 0.975)*se
}


prs_cor_smk_int$r <- as.numeric(prs_cor_smk_int$r)
prs_cor_smk_int$lb <- as.numeric(prs_cor_smk_int$lb)
prs_cor_smk_int$ub <- as.numeric(prs_cor_smk_int$ub)

prs_cor_smk_int

pdf(paste0(outDir,"prs_gscan2_smkInit_corr.pdf"), width = 8, height = 6)
ggplot(prs_cor_smk_int, aes(x=causal_fr, y=r, label=round(r,4))) +
  geom_col(fill = "steelblue", color = "steelblue4", alpha = 0.75) +
  geom_errorbar(aes(ymin = lb, ymax = ub, 
                    width = 0.25)) +
  geom_label_repel(min.segment.length = 0,
                   position = position_nudge_repel(y=0.005)) +
  labs(title = "Correlation between GSCAN2 PRS and Smoking Initiation in the NTR",
       x="Causal Fraction in LDpred", y="Polyserial Correlation (95% Confidence Interval)",
       caption = "Crude correlations, not accounting for twin clustering") +
  coord_cartesian(ylim = c(0.0, 0.25)) +
  theme_light(12) +
  theme(axis.text.x = element_text(face="bold"))
dev.off()


## Correlation with Other Smoking Phenotypes ###############

## Using p<0.1 causal fraction (with the strongest correlation with Smk Initiation)
twin_long_prs_EUR |> 
  select(P_0_1_SCORE_Smoking_SmokingInitiation_MRG16_LDp1, all_of(smk_phenos)) |> 
  glimpse()


cor_prs_smk_int <- as.data.frame(cbind(smk_phenos,
                                       "r" = NA_real_,
                                       "lb" = NA_real_,
                                       "ub" = NA_real_))
# First 5 smk phenos are categorical (binary or ordinal)
for (i in 1:5) {
  test <- polyserial(twin_long_prs_EUR$P_0_1_SCORE_Smoking_SmokingInitiation_MRG16_LDp1,
             twin_long_prs_EUR[,smk_phenos[i]], std.err = T) 
  cor_prs_smk_int$r[i] <- test$rho
  se <- sqrt(test$var)
  cor_prs_smk_int$ub[i] <- test$rho + qnorm(p = 0.975)*se
  cor_prs_smk_int$lb[i] <- test$rho - qnorm(p = 0.975)*se
}


## Rest smk phenos are continuous
for (i in 6:length(smk_phenos)) {
  test <- cor.test(twin_long_prs_EUR$P_0_1_SCORE_Smoking_SmokingInitiation_MRG16_LDp1,
                   twin_long_prs_EUR[,smk_phenos[i]])
  cor_prs_smk_int$r[i] <- test$estimate
  cor_prs_smk_int$lb[i] <- test$conf.int[1]
  cor_prs_smk_int$ub[i] <- test$conf.int[2]
}

str(cor_prs_smk_int)

cor_prs_smk_int$r <- as.numeric(cor_prs_smk_int$r)
cor_prs_smk_int$lb <- as.numeric(cor_prs_smk_int$lb)
cor_prs_smk_int$ub <- as.numeric(cor_prs_smk_int$ub)

cor_prs_smk_int$direction <- factor(ifelse(cor_prs_smk_int$r>0, "pos", "neg"))

cor_prs_smk_int <- cor_prs_smk_int |> 
  mutate(smk_phenos = fct_reorder(as_factor(smk_phenos), -r))

cor_prs_smk_int

pdf(paste0(outDir,"prs_gscan2_smkInit_corr_across_smk_phenos.pdf"), width = 9, height = 6)
ggplot(cor_prs_smk_int, aes(x=smk_phenos, y=r,
                            label=paste(smk_phenos,round(r,3),sep=": "))) +
  geom_col(aes(fill = direction, color = direction), alpha = 0.75) +
  geom_errorbar(aes(ymin = lb, ymax = ub, 
                    width = 0.25)) +
  geom_label_repel(min.segment.length = 0, 
                   position = position_nudge_repel(y=-0.005)) +
  geom_hline(yintercept = 0, color = "grey50") +
  scale_fill_manual(values = c("indianred","steelblue")) +
  scale_color_manual(values = c("indianred4","steelblue4")) +
  labs(title = "Correlation between PRS of Smk Initiation and Smk Phenotypes",
       x="Smoking Trait", y="Polyserial Correlation (95% Confidence Interval)",
       caption = "Crude correlations, not accounting for twin clustering") +
  coord_cartesian(ylim = c(-0.2, 0.3)) +
  theme_light(12) +
  theme(axis.text.x = element_blank(),
        legend.position = "none")
dev.off()


# Save Pheno+PRS Data for EUR ancestry ----------------------------------------------------

dim(twin_long_prs_EUR)


# Regression of Current v Never Smoking ---------------------------------------------------

twin_long_prs_EUR |> 
  count(PLATFORM)

twin_long_prs_EUR |> 
  count(id_sex)

# standardize PCs

twin_long_prs_EUR <- twin_long_prs_EUR |> 
  mutate(female = ifelse(id_sex == "vrouw", 1, 0), 
         PC1_sc = scale(PC1_1KG),
         PC2_sc = scale(PC2_1KG),
         PC3_sc = scale(PC3_1KG),
         PC4_sc = scale(PC4_1KG),
         PC5_sc = scale(PC5_1KG),
         PC6_sc = scale(PC6_1KG),
         PC7_sc = scale(PC7_1KG),
         PC8_sc = scale(PC8_1KG),
         PC9_sc = scale(PC9_1KG),
         PC10_sc = scale(PC10_1KG),
         Platform = factor(PLATFORM)
  )


# create dummy variables for `Platform`
twin_long_prs_EUR <- dummy_cols(twin_long_prs_EUR,
                                select_columns = c("Platform"), 
                                ignore_na = TRUE, 
                                remove_most_frequent_dummy = TRUE)

table(twin_long_prs_EUR$Platform_2, useNA = "ifany")
table(twin_long_prs_EUR$Platform_3, useNA = "ifany")


## Save variables as objects
currentSmk = as.factor(twin_long_prs_EUR$currentSmk)
age = as.numeric(twin_long_prs_EUR$age_sc)
female = as.numeric(twin_long_prs_EUR$female)

PC1 = as.numeric(twin_long_prs_EUR$PC1_sc)
PC2 = as.numeric(twin_long_prs_EUR$PC2_sc)
PC3 = as.numeric(twin_long_prs_EUR$PC3_sc)
PC4 = as.numeric(twin_long_prs_EUR$PC4_sc)
PC5 = as.numeric(twin_long_prs_EUR$PC5_sc)
PC6 = as.numeric(twin_long_prs_EUR$PC6_sc)
PC7 = as.numeric(twin_long_prs_EUR$PC7_sc)
PC8 = as.numeric(twin_long_prs_EUR$PC8_sc)
PC9 = as.numeric(twin_long_prs_EUR$PC9_sc)
PC10 = as.numeric(twin_long_prs_EUR$PC10_sc)

Platform = as.factor(twin_long_prs_EUR$Platform)
Platform_2 = as.numeric(twin_long_prs_EUR$Platform_2)
Platform_3 = as.numeric(twin_long_prs_EUR$Platform_3)

prs_everSmk_p_0_01 = scale(as.numeric(twin_long_prs_EUR$P_0_01_SCORE_Smoking_SmokingInitiation_MRG16_LDp1))
prs_everSmk_p_0_05 = scale(as.numeric(twin_long_prs_EUR$P_0_05_SCORE_Smoking_SmokingInitiation_MRG16_LDp1))
prs_everSmk_p_0_1 = scale(as.numeric(twin_long_prs_EUR$P_0_1_SCORE_Smoking_SmokingInitiation_MRG16_LDp1))
prs_everSmk_p_0_2 = scale(as.numeric(twin_long_prs_EUR$P_0_2_SCORE_Smoking_SmokingInitiation_MRG16_LDp1))
prs_everSmk_p_0_3 = scale(as.numeric(twin_long_prs_EUR$P_0_3_SCORE_Smoking_SmokingInitiation_MRG16_LDp1))
prs_everSmk_p_0_5 = scale(as.numeric(twin_long_prs_EUR$P_0_5_SCORE_Smoking_SmokingInitiation_MRG16_LDp1))
prs_everSmk_p_inf = scale(as.numeric(twin_long_prs_EUR$P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1))


## Baseline Regression (only the covariates)
reg_currentSmk=as.data.frame(
  glm(currentSmk ~ female + age + I(age^2) + Platform + 
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
      family=binomial(link='logit'), na.action=na.omit)$coefficients)

vp2_currentSmk=var(   
  reg_currentSmk["female",]*female+
    reg_currentSmk["age",]*age+
    reg_currentSmk["I(age^2)",]*(age^2)+
    reg_currentSmk["Platform2",]*Platform_2+
    reg_currentSmk["Platform3",]*Platform_3+
    reg_currentSmk["PC1",]*PC1+
    reg_currentSmk["PC2",]*PC2+
    reg_currentSmk["PC3",]*PC3+
    reg_currentSmk["PC4",]*PC4+
    reg_currentSmk["PC5",]*PC5+
    reg_currentSmk["PC6",]*PC6+
    reg_currentSmk["PC7",]*PC7+
    reg_currentSmk["PC8",]*PC8+
    reg_currentSmk["PC9",]*PC9+
    reg_currentSmk["PC10",]*PC10,
  na.rm = TRUE
)     #variance predictors       
vare1 = pi^2 / 3     # residual (homoskedastic) variance
R2logist_currentSmk = vp2_currentSmk / (vp2_currentSmk + vare1)   
R2logist_currentSmk



### p<0.01 -----------------------------------------------------------------------

reg_currentSmk_p_0_01 = as.data.frame(
  glm(currentSmk ~ prs_everSmk_p_0_01 + female + age + I(age^2) + Platform + 
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
      family=binomial(link='logit'), na.action=na.omit)$coefficients)

vp2_currentSmk_p_0_01 = var(   
  reg_currentSmk_p_0_01["prs_everSmk_p_0_01",]*prs_everSmk_p_0_01+
    reg_currentSmk_p_0_01["female",]*female+
    reg_currentSmk_p_0_01["age",]*age+
    reg_currentSmk_p_0_01["I(age^2)",]*(age^2)+
    reg_currentSmk_p_0_01["Platform2",]*Platform_2+
    reg_currentSmk_p_0_01["Platform3",]*Platform_3+
    reg_currentSmk_p_0_01["PC1",]*PC1+
    reg_currentSmk_p_0_01["PC2",]*PC2+
    reg_currentSmk_p_0_01["PC3",]*PC3+
    reg_currentSmk_p_0_01["PC4",]*PC4+
    reg_currentSmk_p_0_01["PC5",]*PC5+
    reg_currentSmk_p_0_01["PC6",]*PC6+
    reg_currentSmk_p_0_01["PC7",]*PC7+
    reg_currentSmk_p_0_01["PC8",]*PC8+
    reg_currentSmk_p_0_01["PC9",]*PC9+
    reg_currentSmk_p_0_01["PC10",]*PC10,
  na.rm = TRUE
)     #variance predictors       
vare1 = pi^2 / 3     # residual (homoskedastic) variance
R2logist_currentSmk_p_0_01 = vp2_currentSmk_p_0_01 / (vp2_currentSmk_p_0_01 + vare1)   
R2logist_currentSmk_p_0_01

# Diff. in the R2
r2_currentSmk_p_0_01 <- matrix(100*(R2logist_currentSmk_p_0_01 - R2logist_currentSmk)) # % variance explained by PRS
rownames(r2_currentSmk_p_0_01) <- "p<0.01"
r2_currentSmk_p_0_01



### p<0.05 -----------------------------------------------------------------------

reg_currentSmk_p_0_05 = as.data.frame(
  glm(currentSmk ~ prs_everSmk_p_0_05 + female + age + I(age^2) + Platform + 
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
      family=binomial(link='logit'), na.action=na.omit)$coefficients)

vp2_currentSmk_p_0_05 = var(   
  reg_currentSmk_p_0_05["prs_everSmk_p_0_05",]*prs_everSmk_p_0_05+
    reg_currentSmk_p_0_05["female",]*female+
    reg_currentSmk_p_0_05["age",]*age+
    reg_currentSmk_p_0_05["I(age^2)",]*(age^2)+
    reg_currentSmk_p_0_05["Platform2",]*Platform_2+
    reg_currentSmk_p_0_05["Platform3",]*Platform_3+
    reg_currentSmk_p_0_05["PC1",]*PC1+
    reg_currentSmk_p_0_05["PC2",]*PC2+
    reg_currentSmk_p_0_05["PC3",]*PC3+
    reg_currentSmk_p_0_05["PC4",]*PC4+
    reg_currentSmk_p_0_05["PC5",]*PC5+
    reg_currentSmk_p_0_05["PC6",]*PC6+
    reg_currentSmk_p_0_05["PC7",]*PC7+
    reg_currentSmk_p_0_05["PC8",]*PC8+
    reg_currentSmk_p_0_05["PC9",]*PC9+
    reg_currentSmk_p_0_05["PC10",]*PC10,
  na.rm = TRUE
)     #variance predictors       
R2logist_currentSmk_p_0_05 = vp2_currentSmk_p_0_05 / (vp2_currentSmk_p_0_05 + vare1)   
R2logist_currentSmk_p_0_05

# Diff. in the R2
r2_currentSmk_p_0_05 <- matrix(100*(R2logist_currentSmk_p_0_05 - R2logist_currentSmk)) # % variance explained by PRS
rownames(r2_currentSmk_p_0_05) <- "p<0.05"
r2_currentSmk_p_0_05


### p<0.1 -----------------------------------------------------------------------

reg_currentSmk_p_0_1 = as.data.frame(
  glm(currentSmk ~ prs_everSmk_p_0_1 + female + age + I(age^2) + Platform + 
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
      family=binomial(link='logit'), na.action=na.omit)$coefficients)

vp2_currentSmk_p_0_1 = var(   
  reg_currentSmk_p_0_1["prs_everSmk_p_0_1",]*prs_everSmk_p_0_1+
    reg_currentSmk_p_0_1["female",]*female+
    reg_currentSmk_p_0_1["age",]*age+
    reg_currentSmk_p_0_1["I(age^2)",]*(age^2)+
    reg_currentSmk_p_0_1["Platform2",]*Platform_2+
    reg_currentSmk_p_0_1["Platform3",]*Platform_3+
    reg_currentSmk_p_0_1["PC1",]*PC1+
    reg_currentSmk_p_0_1["PC2",]*PC2+
    reg_currentSmk_p_0_1["PC3",]*PC3+
    reg_currentSmk_p_0_1["PC4",]*PC4+
    reg_currentSmk_p_0_1["PC5",]*PC5+
    reg_currentSmk_p_0_1["PC6",]*PC6+
    reg_currentSmk_p_0_1["PC7",]*PC7+
    reg_currentSmk_p_0_1["PC8",]*PC8+
    reg_currentSmk_p_0_1["PC9",]*PC9+
    reg_currentSmk_p_0_1["PC10",]*PC10,
  na.rm = TRUE
)     #variance predictors       
R2logist_currentSmk_p_0_1 = vp2_currentSmk_p_0_1 / (vp2_currentSmk_p_0_1 + vare1)   
R2logist_currentSmk_p_0_1

# Diff. in the R2

r2_currentSmk_p_0_1 <- matrix(100*(R2logist_currentSmk_p_0_1 - R2logist_currentSmk)) # % variance explained by PRS
rownames(r2_currentSmk_p_0_1) <- "p<0.1"
r2_currentSmk_p_0_1


### p<0.2 -----------------------------------------------------------------------

reg_currentSmk_p_0_2 = as.data.frame(
  glm(currentSmk ~ prs_everSmk_p_0_2 + female + age + I(age^2) + Platform + 
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
      family=binomial(link='logit'), na.action=na.omit)$coefficients)

vp2_currentSmk_p_0_2 = var(   
  reg_currentSmk_p_0_2["prs_everSmk_p_0_2",]*prs_everSmk_p_0_2+
    reg_currentSmk_p_0_2["female",]*female+
    reg_currentSmk_p_0_2["age",]*age+
    reg_currentSmk_p_0_2["I(age^2)",]*(age^2)+
    reg_currentSmk_p_0_2["Platform2",]*Platform_2+
    reg_currentSmk_p_0_2["Platform3",]*Platform_3+
    reg_currentSmk_p_0_2["PC1",]*PC1+
    reg_currentSmk_p_0_2["PC2",]*PC2+
    reg_currentSmk_p_0_2["PC3",]*PC3+
    reg_currentSmk_p_0_2["PC4",]*PC4+
    reg_currentSmk_p_0_2["PC5",]*PC5+
    reg_currentSmk_p_0_2["PC6",]*PC6+
    reg_currentSmk_p_0_2["PC7",]*PC7+
    reg_currentSmk_p_0_2["PC8",]*PC8+
    reg_currentSmk_p_0_2["PC9",]*PC9+
    reg_currentSmk_p_0_2["PC10",]*PC10,
  na.rm = TRUE
)     #variance predictors       
R2logist_currentSmk_p_0_2 = vp2_currentSmk_p_0_2 / (vp2_currentSmk_p_0_2 + vare1)   
R2logist_currentSmk_p_0_2

# Diff. in the R2
r2_currentSmk_p_0_2 <- matrix(100*(R2logist_currentSmk_p_0_2 - R2logist_currentSmk)) # % variance explained by PRS
rownames(r2_currentSmk_p_0_2) <- "p<0.2"
r2_currentSmk_p_0_2



### p<0.3 ------------------------------------------------------------------------

reg_currentSmk_p_0_3 = as.data.frame(
  glm(currentSmk ~ prs_everSmk_p_0_3 + female + age + I(age^2) + Platform + 
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
      family=binomial(link='logit'), na.action=na.omit)$coefficients)

vp2_currentSmk_p_0_3 = var(   
  reg_currentSmk_p_0_3["prs_everSmk_p_0_3",]*prs_everSmk_p_0_3+
    reg_currentSmk_p_0_3["female",]*female+
    reg_currentSmk_p_0_3["age",]*age+
    reg_currentSmk_p_0_3["I(age^2)",]*(age^2)+
    reg_currentSmk_p_0_3["Platform2",]*Platform_2+
    reg_currentSmk_p_0_3["Platform3",]*Platform_3+
    reg_currentSmk_p_0_3["PC1",]*PC1+
    reg_currentSmk_p_0_3["PC2",]*PC2+
    reg_currentSmk_p_0_3["PC3",]*PC3+
    reg_currentSmk_p_0_3["PC4",]*PC4+
    reg_currentSmk_p_0_3["PC5",]*PC5+
    reg_currentSmk_p_0_3["PC6",]*PC6+
    reg_currentSmk_p_0_3["PC7",]*PC7+
    reg_currentSmk_p_0_3["PC8",]*PC8+
    reg_currentSmk_p_0_3["PC9",]*PC9+
    reg_currentSmk_p_0_3["PC10",]*PC10,
  na.rm = TRUE
)     #variance predictors       
R2logist_currentSmk_p_0_3 = vp2_currentSmk_p_0_3 / (vp2_currentSmk_p_0_3 + vare1)   
R2logist_currentSmk_p_0_3

# Diff. in the R2
r2_currentSmk_p_0_3 <- matrix(100*(R2logist_currentSmk_p_0_3 - R2logist_currentSmk)) # % variance explained by PRS
rownames(r2_currentSmk_p_0_3) <- "p<0.3"
r2_currentSmk_p_0_3



### p<0.5 -----------------------------------------------------------------------

reg_currentSmk_p_0_5 = as.data.frame(
  glm(currentSmk ~ prs_everSmk_p_0_5 + female + age + I(age^2) + Platform + 
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
      family=binomial(link='logit'), na.action=na.omit)$coefficients)

vp2_currentSmk_p_0_5 = var(   
  reg_currentSmk_p_0_5["prs_everSmk_p_0_5",]*prs_everSmk_p_0_5+
    reg_currentSmk_p_0_5["female",]*female+
    reg_currentSmk_p_0_5["age",]*age+
    reg_currentSmk_p_0_5["I(age^2)",]*(age^2)+
    reg_currentSmk_p_0_5["Platform2",]*Platform_2+
    reg_currentSmk_p_0_5["Platform3",]*Platform_3+
    reg_currentSmk_p_0_5["PC1",]*PC1+
    reg_currentSmk_p_0_5["PC2",]*PC2+
    reg_currentSmk_p_0_5["PC3",]*PC3+
    reg_currentSmk_p_0_5["PC4",]*PC4+
    reg_currentSmk_p_0_5["PC5",]*PC5+
    reg_currentSmk_p_0_5["PC6",]*PC6+
    reg_currentSmk_p_0_5["PC7",]*PC7+
    reg_currentSmk_p_0_5["PC8",]*PC8+
    reg_currentSmk_p_0_5["PC9",]*PC9+
    reg_currentSmk_p_0_5["PC10",]*PC10,
  na.rm = TRUE
)     #variance predictors       
R2logist_currentSmk_p_0_5 = vp2_currentSmk_p_0_5 / (vp2_currentSmk_p_0_5 + vare1)   
R2logist_currentSmk_p_0_5


# Diff. in the R2
r2_currentSmk_p_0_5 <- matrix(100*(R2logist_currentSmk_p_0_5 - R2logist_currentSmk)) # % variance explained by PRS
rownames(r2_currentSmk_p_0_5) <- "p<0.5"
r2_currentSmk_p_0_5



### p-inf -----------------------------------------------------------------------

reg_currentSmk_p_inf = as.data.frame(
  glm(currentSmk ~ prs_everSmk_p_inf + female + age + I(age^2) + Platform + 
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
      family=binomial(link='logit'), na.action=na.omit)$coefficients)

vp2_currentSmk_p_inf = var(   
  reg_currentSmk_p_inf["prs_everSmk_p_inf",]*prs_everSmk_p_inf+
    reg_currentSmk_p_inf["female",]*female+
    reg_currentSmk_p_inf["age",]*age+
    reg_currentSmk_p_inf["I(age^2)",]*(age^2)+
    reg_currentSmk_p_inf["Platform2",]*Platform_2+
    reg_currentSmk_p_inf["Platform3",]*Platform_3+
    reg_currentSmk_p_inf["PC1",]*PC1+
    reg_currentSmk_p_inf["PC2",]*PC2+
    reg_currentSmk_p_inf["PC3",]*PC3+
    reg_currentSmk_p_inf["PC4",]*PC4+
    reg_currentSmk_p_inf["PC5",]*PC5+
    reg_currentSmk_p_inf["PC6",]*PC6+
    reg_currentSmk_p_inf["PC7",]*PC7+
    reg_currentSmk_p_inf["PC8",]*PC8+
    reg_currentSmk_p_inf["PC9",]*PC9+
    reg_currentSmk_p_inf["PC10",]*PC10,
  na.rm = TRUE
)     #variance predictors       
R2logist_currentSmk_p_inf = vp2_currentSmk_p_inf / (vp2_currentSmk_p_inf + vare1)   
R2logist_currentSmk_p_inf

# Diff. in the R2

r2_currentSmk_p_inf <- matrix(100*(R2logist_currentSmk_p_inf - R2logist_currentSmk)) # % variance explained by PRS
rownames(r2_currentSmk_p_inf) <- "p-Inf"
r2_currentSmk_p_inf




# Regression of Former vs. Never Smokers -------------------------------------------------------------------

## Save variables as objects
formerSmk = as.factor(twin_long_prs_EUR$formerSmk)


## Baseline Regression (only the covariates)
reg_formerSmk=as.data.frame(
  glm(formerSmk ~ female + age + I(age^2) + Platform + 
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
      family=binomial(link='logit'), na.action=na.omit)$coefficients)

vp2_formerSmk=var(   
  reg_formerSmk["female",]*female+
    reg_formerSmk["age",]*age+
    reg_formerSmk["I(age^2)",]*(age^2)+
    reg_formerSmk["Platform2",]*Platform_2+
    reg_formerSmk["Platform3",]*Platform_3+
    reg_formerSmk["PC1",]*PC1+
    reg_formerSmk["PC2",]*PC2+
    reg_formerSmk["PC3",]*PC3+
    reg_formerSmk["PC4",]*PC4+
    reg_formerSmk["PC5",]*PC5+
    reg_formerSmk["PC6",]*PC6+
    reg_formerSmk["PC7",]*PC7+
    reg_formerSmk["PC8",]*PC8+
    reg_formerSmk["PC9",]*PC9+
    reg_formerSmk["PC10",]*PC10,
  na.rm = TRUE
)     #variance predictors       
vare1 = pi^2 / 3     # residual (homoskedastic) variance
R2logist_formerSmk = vp2_formerSmk / (vp2_formerSmk + vare1)   
R2logist_formerSmk
# 0.1654737


### p<0.01 -----------------------------------------------------------------------

reg_formerSmk_p_0_01 = as.data.frame(
  glm(formerSmk ~ prs_everSmk_p_0_01 + female + age + I(age^2) + Platform + 
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
      family=binomial(link='logit'), na.action=na.omit)$coefficients)

vp2_formerSmk_p_0_01 = var(   
  reg_formerSmk_p_0_01["prs_everSmk_p_0_01",]*prs_everSmk_p_0_01+
    reg_formerSmk_p_0_01["female",]*female+
    reg_formerSmk_p_0_01["age",]*age+
    reg_formerSmk_p_0_01["I(age^2)",]*(age^2)+
    reg_formerSmk_p_0_01["Platform2",]*Platform_2+
    reg_formerSmk_p_0_01["Platform3",]*Platform_3+
    reg_formerSmk_p_0_01["PC1",]*PC1+
    reg_formerSmk_p_0_01["PC2",]*PC2+
    reg_formerSmk_p_0_01["PC3",]*PC3+
    reg_formerSmk_p_0_01["PC4",]*PC4+
    reg_formerSmk_p_0_01["PC5",]*PC5+
    reg_formerSmk_p_0_01["PC6",]*PC6+
    reg_formerSmk_p_0_01["PC7",]*PC7+
    reg_formerSmk_p_0_01["PC8",]*PC8+
    reg_formerSmk_p_0_01["PC9",]*PC9+
    reg_formerSmk_p_0_01["PC10",]*PC10,
  na.rm = TRUE
)     #variance predictors       
vare1 = pi^2 / 3     # residual (homoskedastic) variance
R2logist_formerSmk_p_0_01 = vp2_formerSmk_p_0_01 / (vp2_formerSmk_p_0_01 + vare1)   
R2logist_formerSmk_p_0_01


# Diff. in the R2

r2_formerSmk_p_0_01 <- matrix(100*(R2logist_formerSmk_p_0_01 - R2logist_formerSmk)) # % variance explained by PRS
rownames(r2_formerSmk_p_0_01) <- "p<0.01"
r2_formerSmk_p_0_01



### p<0.05 -----------------------------------------------------------------------

reg_formerSmk_p_0_05 = as.data.frame(
  glm(formerSmk ~ prs_everSmk_p_0_05 + female + age + I(age^2) + Platform + 
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
      family=binomial(link='logit'), na.action=na.omit)$coefficients)

vp2_formerSmk_p_0_05 = var(   
  reg_formerSmk_p_0_05["prs_everSmk_p_0_05",]*prs_everSmk_p_0_05+
    reg_formerSmk_p_0_05["female",]*female+
    reg_formerSmk_p_0_05["age",]*age+
    reg_formerSmk_p_0_05["I(age^2)",]*(age^2)+
    reg_formerSmk_p_0_05["Platform2",]*Platform_2+
    reg_formerSmk_p_0_05["Platform3",]*Platform_3+
    reg_formerSmk_p_0_05["PC1",]*PC1+
    reg_formerSmk_p_0_05["PC2",]*PC2+
    reg_formerSmk_p_0_05["PC3",]*PC3+
    reg_formerSmk_p_0_05["PC4",]*PC4+
    reg_formerSmk_p_0_05["PC5",]*PC5+
    reg_formerSmk_p_0_05["PC6",]*PC6+
    reg_formerSmk_p_0_05["PC7",]*PC7+
    reg_formerSmk_p_0_05["PC8",]*PC8+
    reg_formerSmk_p_0_05["PC9",]*PC9+
    reg_formerSmk_p_0_05["PC10",]*PC10,
  na.rm = TRUE
)     #variance predictors       
R2logist_formerSmk_p_0_05 = vp2_formerSmk_p_0_05 / (vp2_formerSmk_p_0_05 + vare1)   
R2logist_formerSmk_p_0_05

# Diff. in the R2

r2_formerSmk_p_0_05 <- matrix(100*(R2logist_formerSmk_p_0_05 - R2logist_formerSmk)) # % variance explained by PRS
rownames(r2_formerSmk_p_0_05) <- "p<0.05"
r2_formerSmk_p_0_05


### p<0.1 -----------------------------------------------------------------------

reg_formerSmk_p_0_1 = as.data.frame(
  glm(formerSmk ~ prs_everSmk_p_0_1 + female + age + I(age^2) + Platform + 
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
      family=binomial(link='logit'), na.action=na.omit)$coefficients)

vp2_formerSmk_p_0_1 = var(   
  reg_formerSmk_p_0_1["prs_everSmk_p_0_1",]*prs_everSmk_p_0_1+
    reg_formerSmk_p_0_1["female",]*female+
    reg_formerSmk_p_0_1["age",]*age+
    reg_formerSmk_p_0_1["I(age^2)",]*(age^2)+
    reg_formerSmk_p_0_1["Platform2",]*Platform_2+
    reg_formerSmk_p_0_1["Platform3",]*Platform_3+
    reg_formerSmk_p_0_1["PC1",]*PC1+
    reg_formerSmk_p_0_1["PC2",]*PC2+
    reg_formerSmk_p_0_1["PC3",]*PC3+
    reg_formerSmk_p_0_1["PC4",]*PC4+
    reg_formerSmk_p_0_1["PC5",]*PC5+
    reg_formerSmk_p_0_1["PC6",]*PC6+
    reg_formerSmk_p_0_1["PC7",]*PC7+
    reg_formerSmk_p_0_1["PC8",]*PC8+
    reg_formerSmk_p_0_1["PC9",]*PC9+
    reg_formerSmk_p_0_1["PC10",]*PC10,
  na.rm = TRUE
)     #variance predictors       
R2logist_formerSmk_p_0_1 = vp2_formerSmk_p_0_1 / (vp2_formerSmk_p_0_1 + vare1)   
R2logist_formerSmk_p_0_1

# Diff. in the R2

r2_formerSmk_p_0_1 <- matrix(100*(R2logist_formerSmk_p_0_1 - R2logist_formerSmk)) # % variance explained by PRS
rownames(r2_formerSmk_p_0_1) <- "p<0.1"
r2_formerSmk_p_0_1


### p<0.2 -----------------------------------------------------------------------

reg_formerSmk_p_0_2 = as.data.frame(
  glm(formerSmk ~ prs_everSmk_p_0_2 + female + age + I(age^2) + Platform + 
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
      family=binomial(link='logit'), na.action=na.omit)$coefficients)

vp2_formerSmk_p_0_2 = var(   
  reg_formerSmk_p_0_2["prs_everSmk_p_0_2",]*prs_everSmk_p_0_2+
    reg_formerSmk_p_0_2["female",]*female+
    reg_formerSmk_p_0_2["age",]*age+
    reg_formerSmk_p_0_2["I(age^2)",]*(age^2)+
    reg_formerSmk_p_0_2["Platform2",]*Platform_2+
    reg_formerSmk_p_0_2["Platform3",]*Platform_3+
    reg_formerSmk_p_0_2["PC1",]*PC1+
    reg_formerSmk_p_0_2["PC2",]*PC2+
    reg_formerSmk_p_0_2["PC3",]*PC3+
    reg_formerSmk_p_0_2["PC4",]*PC4+
    reg_formerSmk_p_0_2["PC5",]*PC5+
    reg_formerSmk_p_0_2["PC6",]*PC6+
    reg_formerSmk_p_0_2["PC7",]*PC7+
    reg_formerSmk_p_0_2["PC8",]*PC8+
    reg_formerSmk_p_0_2["PC9",]*PC9+
    reg_formerSmk_p_0_2["PC10",]*PC10,
  na.rm = TRUE
)     #variance predictors       
R2logist_formerSmk_p_0_2 = vp2_formerSmk_p_0_2 / (vp2_formerSmk_p_0_2 + vare1)   
R2logist_formerSmk_p_0_2

# Diff. in the R2

r2_formerSmk_p_0_2 <- matrix(100*(R2logist_formerSmk_p_0_2 - R2logist_formerSmk)) # % variance explained by PRS
rownames(r2_formerSmk_p_0_2) <- "p<0.2"
r2_formerSmk_p_0_2


### p<0.3 ------------------------------------------------------------------------

reg_formerSmk_p_0_3 = as.data.frame(
  glm(formerSmk ~ prs_everSmk_p_0_3 + female + age + I(age^2) + Platform + 
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
      family=binomial(link='logit'), na.action=na.omit)$coefficients)

vp2_formerSmk_p_0_3 = var(   
  reg_formerSmk_p_0_3["prs_everSmk_p_0_3",]*prs_everSmk_p_0_3+
    reg_formerSmk_p_0_3["female",]*female+
    reg_formerSmk_p_0_3["age",]*age+
    reg_formerSmk_p_0_3["I(age^2)",]*(age^2)+
    reg_formerSmk_p_0_3["Platform2",]*Platform_2+
    reg_formerSmk_p_0_3["Platform3",]*Platform_3+
    reg_formerSmk_p_0_3["PC1",]*PC1+
    reg_formerSmk_p_0_3["PC2",]*PC2+
    reg_formerSmk_p_0_3["PC3",]*PC3+
    reg_formerSmk_p_0_3["PC4",]*PC4+
    reg_formerSmk_p_0_3["PC5",]*PC5+
    reg_formerSmk_p_0_3["PC6",]*PC6+
    reg_formerSmk_p_0_3["PC7",]*PC7+
    reg_formerSmk_p_0_3["PC8",]*PC8+
    reg_formerSmk_p_0_3["PC9",]*PC9+
    reg_formerSmk_p_0_3["PC10",]*PC10,
  na.rm = TRUE
)     #variance predictors       
R2logist_formerSmk_p_0_3 = vp2_formerSmk_p_0_3 / (vp2_formerSmk_p_0_3 + vare1)   
R2logist_formerSmk_p_0_3

# Diff. in the R2
r2_formerSmk_p_0_3 <- matrix(100*(R2logist_formerSmk_p_0_3 - R2logist_formerSmk)) # % variance explained by PRS
rownames(r2_formerSmk_p_0_3) <- "p<0.3"
r2_formerSmk_p_0_3



### p<0.5 -----------------------------------------------------------------------

reg_formerSmk_p_0_5 = as.data.frame(
  glm(formerSmk ~ prs_everSmk_p_0_5 + female + age + I(age^2) + Platform + 
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
      family=binomial(link='logit'), na.action=na.omit)$coefficients)

vp2_formerSmk_p_0_5 = var(   
  reg_formerSmk_p_0_5["prs_everSmk_p_0_5",]*prs_everSmk_p_0_5+
    reg_formerSmk_p_0_5["female",]*female+
    reg_formerSmk_p_0_5["age",]*age+
    reg_formerSmk_p_0_5["I(age^2)",]*(age^2)+
    reg_formerSmk_p_0_5["Platform2",]*Platform_2+
    reg_formerSmk_p_0_5["Platform3",]*Platform_3+
    reg_formerSmk_p_0_5["PC1",]*PC1+
    reg_formerSmk_p_0_5["PC2",]*PC2+
    reg_formerSmk_p_0_5["PC3",]*PC3+
    reg_formerSmk_p_0_5["PC4",]*PC4+
    reg_formerSmk_p_0_5["PC5",]*PC5+
    reg_formerSmk_p_0_5["PC6",]*PC6+
    reg_formerSmk_p_0_5["PC7",]*PC7+
    reg_formerSmk_p_0_5["PC8",]*PC8+
    reg_formerSmk_p_0_5["PC9",]*PC9+
    reg_formerSmk_p_0_5["PC10",]*PC10,
  na.rm = TRUE
)     #variance predictors       
R2logist_formerSmk_p_0_5 = vp2_formerSmk_p_0_5 / (vp2_formerSmk_p_0_5 + vare1)   
R2logist_formerSmk_p_0_5

# Diff. in the R2
r2_formerSmk_p_0_5 <- matrix(100*(R2logist_formerSmk_p_0_5 - R2logist_formerSmk)) # % variance explained by PRS
rownames(r2_formerSmk_p_0_5) <- "p<0.5"
r2_formerSmk_p_0_5



### p-inf -----------------------------------------------------------------------


reg_formerSmk_p_inf = as.data.frame(
  glm(formerSmk ~ prs_everSmk_p_inf + female + age + I(age^2) + Platform + 
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
      family=binomial(link='logit'), na.action=na.omit)$coefficients)

vp2_formerSmk_p_inf = var(   
  reg_formerSmk_p_inf["prs_everSmk_p_inf",]*prs_everSmk_p_inf+
    reg_formerSmk_p_inf["female",]*female+
    reg_formerSmk_p_inf["age",]*age+
    reg_formerSmk_p_inf["I(age^2)",]*(age^2)+
    reg_formerSmk_p_inf["Platform2",]*Platform_2+
    reg_formerSmk_p_inf["Platform3",]*Platform_3+
    reg_formerSmk_p_inf["PC1",]*PC1+
    reg_formerSmk_p_inf["PC2",]*PC2+
    reg_formerSmk_p_inf["PC3",]*PC3+
    reg_formerSmk_p_inf["PC4",]*PC4+
    reg_formerSmk_p_inf["PC5",]*PC5+
    reg_formerSmk_p_inf["PC6",]*PC6+
    reg_formerSmk_p_inf["PC7",]*PC7+
    reg_formerSmk_p_inf["PC8",]*PC8+
    reg_formerSmk_p_inf["PC9",]*PC9+
    reg_formerSmk_p_inf["PC10",]*PC10,
  na.rm = TRUE
)     #variance predictors       
R2logist_formerSmk_p_inf = vp2_formerSmk_p_inf / (vp2_formerSmk_p_inf + vare1)   
R2logist_formerSmk_p_inf


# Diff. in the R2

r2_formerSmk_p_inf <- matrix(100*(R2logist_formerSmk_p_inf - R2logist_formerSmk)) # % variance explained by PRS
rownames(r2_formerSmk_p_inf) <- "p-Inf"
r2_formerSmk_p_inf



## END. --------------------------
