## Save the Smk and Covariate Data for MR-DOC with residualized DNAm and PRS
## Madhur Singh

# clear the workspace
rm(list = ls()); gc()

## Directories #####################

baseDir <- "/base_dir/" 
prsDir  <- paste0(baseDir,"data/PRS_Smk/GSCAN2/")
phenDir <- paste0(baseDir,"data/phenos/")
outDir  <- paste0(baseDir,"out/pheno_descr/")


# Libraries
library(tidyverse)
library(forcats)
library(polycor)      # for polyserial correlation of binary outcome


## Load the long-format covariate twin data
load(file = paste0(phenDir,"twin_covariates_long_data.Rdata"))

dim(long_twin_data)
names(long_twin_data)
head(long_twin_data$IdatName)

glimpse(head(long_twin_data))

table(long_twin_data$zygo)


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

## Filter EUR ancestry ------------------------------------

phenos_prs_smk_init |> 
  count(EUR_1KG_Outlier)

phenos_prs_smk_init |> 
  filter(EUR_1KG_Outlier == 1) |> 
  group_by(FamilyNumber) |> 
  mutate(twin_nr = row_number()) |> 
  ungroup() |> 
  count(twin_nr)

twin_long_prs <- phenos_prs_smk_init |> 
  filter(EUR_1KG_Outlier == 0)

dim(twin_long_prs)
names(twin_long_prs)
str(twin_long_prs)

twin_long_prs <- as.data.frame(twin_long_prs)

## Max R2 is with prs_p01 (5.1% for current v. never; 2.0% for former v. never smk)



## Sample Decription ------------------------------------------------------

twin_long_prs |> 
  group_by(FamilyNumber) |> 
  mutate(twin_nr = row_number()) |> 
  ungroup() |> 
  count(twin_nr)

twin_long_prs |> 
  group_by(FamilyNumber) |> 
  mutate(twin_nr = row_number()) |> 
  ungroup() |> 
  group_by(zygo) |> 
  count(twin_nr) |> 
  mutate(incomp = max(n) - min(n))

twin_long_prs |> 
  count(id_sex) |> 
  mutate(prop = n/sum(n))

psych::describe(twin_long_prs$DNA_BloodSampling_Age)





## Smoking Phenotype ------------------------------------------------------

table(twin_long_prs$Smoking)


## Residualize the PRS for PCs and genotyping platform ------------------

twin_long_prs$PLATFORM <- as.factor(twin_long_prs$PLATFORM)
class(twin_long_prs$PLATFORM)
table(twin_long_prs$PLATFORM, useNA = "ifany")

rownames(twin_long_prs) <- twin_long_prs$IdatName

modPRS <- lm(P_0_1_SCORE_Smoking_SmokingInitiation_MRG16_LDp1 ~ PLATFORM + 
               PC1_1KG + PC2_1KG + PC3_1KG + PC4_1KG + PC5_1KG + 
               PC6_1KG + PC7_1KG + PC8_1KG + PC9_1KG + PC10_1KG,
             data = twin_long_prs)
summary(modPRS)

resPRS <- as.data.frame(modPRS$residuals)
colnames(resPRS) <- "resPRS"
resPRS$IdatName <- rownames(resPRS)
summary(resPRS$resPRS)

## c.f.
mean(twin_long_prs$P_0_1_SCORE_Smoking_SmokingInitiation_MRG16_LDp1, na.rm = T)
sd(twin_long_prs$P_0_1_SCORE_Smoking_SmokingInitiation_MRG16_LDp1, na.rm = T)

resPRS$resPRS <- as.numeric(scale(resPRS$resPRS))
summary(resPRS$resPRS)
head(resPRS)
nrow(resPRS)

twin_long_prs <- merge(twin_long_prs, resPRS, by = "IdatName", all.x = T)
summary(twin_long_prs$resPRS)


## Clean the covariates --------------------------------------------

# better named variables, centered and scaled
twin_reg_data <- twin_long_prs |> 
  mutate(currentSmk = as.numeric(currentSmk) - 1,                # Needs to be a numeric dummy variable
         formerSmk = as.numeric(formerSmk) - 1, 
         smkCess = as.numeric(smkCess) - 1, 
         everSmk = as.numeric(everSmk) - 1, 
         age = as.numeric(scale(DNA_BloodSampling_Age)),
         female = ifelse(id_sex == "vrouw", 1, 0),
         PC1_sc = as.numeric(scale(PC1_1KG)),
         PC2_sc = as.numeric(scale(PC2_1KG)),
         PC3_sc = as.numeric(scale(PC3_1KG)),
         PC4_sc = as.numeric(scale(PC4_1KG)),
         PC5_sc = as.numeric(scale(PC5_1KG)),
         PC6_sc = as.numeric(scale(PC6_1KG)),
         PC7_sc = as.numeric(scale(PC7_1KG)),
         PC8_sc = as.numeric(scale(PC8_1KG)),
         PC9_sc = as.numeric(scale(PC9_1KG)),
         PC10_sc = as.numeric(scale(PC10_1KG)),
         Platform = factor(PLATFORM)
  )

# create dummy variables for `Platform`
twin_reg_data <- fastDummies::dummy_cols(twin_reg_data,
                                         select_columns = c("Platform"), 
                                         ignore_na = TRUE, 
                                         remove_most_frequent_dummy = TRUE)

table(twin_reg_data$Platform_2, useNA = "ifany")
table(twin_reg_data$Platform_3, useNA = "ifany")


## Select columns 
twin_reg_data_select <- twin_reg_data |> 
  select(
    FamilyNumber, twin_nr, IdatName, fisnumber, zygo, DNA_BloodSampling_Age, age, female,
    Smoking, everSmk, currentSmk, formerSmk, smkCess, ncigday, nyearstop, resPRS,
    Platform_2, Platform_3, Sample_Plate, Array_rownum, Neut_Perc, Mono_Perc, Eos_Perc,
    PC1_sc, PC2_sc, PC3_sc, PC4_sc, PC5_sc, PC6_sc, PC7_sc, PC8_sc, PC9_sc, PC10_sc
  ) 

head(twin_reg_data_select)

dim(twin_reg_data_select)


## Assign twin index within family - for pivoting to wide

twin_reg_data_select <- twin_reg_data_select |> 
  group_by(FamilyNumber) |> 
  mutate(twin_nr = row_number(),
         twin_nr = paste0("T",twin_nr)) |> 
  arrange(FamilyNumber, twin_nr) |> 
  ungroup()


twin_reg_data_select  |> 
  count(twin_nr)

## Pivot to wide format

twin_dat <- twin_reg_data_select |> 
  pivot_wider(id_cols = c(FamilyNumber, zygo), names_from = twin_nr, 
              values_from = c(fisnumber, IdatName, age, female, resPRS,
                              everSmk, currentSmk, formerSmk, smkCess,
                              ncigday, nyearstop,  Platform_2, Platform_3,
                              PC1_sc, PC2_sc, PC3_sc, PC4_sc, PC5_sc, 
                              PC6_sc, PC7_sc, PC8_sc, PC9_sc, PC10_sc,
                              Sample_Plate, Array_rownum, 
                              Neut_Perc, Mono_Perc, Eos_Perc), 
              names_glue = "{twin_nr}_{.value}")

dim(twin_dat)

head(names(twin_dat), 10)

table(is.na(twin_dat$T1_fisnumber))
table(is.na(twin_dat$T2_fisnumber))


## Check missingness ---------------------------------------------

covars  <- c(
  # Demographic covariates
  "age", "female", 
  # PCs
  "PC1_sc", "PC2_sc", "PC3_sc", "PC4_sc", "PC5_sc", 
  "PC6_sc", "PC7_sc", "PC8_sc", "PC9_sc", "PC10_sc",
  # Platform
  "Platform_2", "Platform_3",
  # Methylation Covariates
  "Sample_Plate", "Array_rownum", 
  "Neut_Perc", "Mono_Perc", "Eos_Perc"
) 
nCovars <- length(covars)
nCovars


defvars <- paste(rep("T2",nCovars), covars, sep = "_")
for (ii in 1:nCovars) {
  twin_dat[which(is.na(twin_dat$T2_fisnumber)), defvars[ii]] <- -99999
}


# check
twin_dat |> 
  filter(is.na(T2_fisnumber)) |> 
  count(T2_female)

# Now confirm there is no missing data on covariates (definition variables)
twin_dat |> 
  mutate(missing_def = ifelse(is.na(T1_female) | is.na(T2_female) | 
                                is.na(T1_PC1_sc) | is.na(T2_PC1_sc), 
                              TRUE, FALSE)) |> 
  count(missing_def)


## Sample Size by zygosity -----------------------------

table(twin_dat$zygo, useNA = "ifany")

twin_dat |> 
  mutate(zygo = fct_recode(as_factor(zygo), "MZ" = "1", "DZ" = "2")) |> 
  group_by(zygo) |> 
  count(is.na(T2_fisnumber))


## PRS F-stats -------------------------------------------------

## Effective sample size, given N pairs of twins
## Minica et al, https://doi.org/10.1038/mp.2014.121 
# NEmz = (2*Nmz) / (1+rMZ)
# NEdz = (2*Ndz) / (1+rDZ)


## R2 of prs_p01  = 0.05067724 for current v. never; 0.02020994 for former v. never smk

#### Current Smk ------------------------------------------------------------

## MZ twins
countMZ <- twin_long_prs |> 
  filter(zygo == 1,
         !is.na(currentSmk)) |> 
  count(FamilyNumber) |> 
  count(n)
Nmz     <- countMZ[countMZ$n==2,"nn"]
single1 <- countMZ[countMZ$n==1,"nn"]

rMZ <- polycor::polychor(twin_dat[twin_dat$zygo==1,]$T1_currentSmk,
                         twin_dat[twin_dat$zygo==1,]$T2_currentSmk)

NEmz <- (2*Nmz) / (1 + rMZ)

## DZ twins
countDZ <- twin_long_prs |> 
  filter(zygo == 2,
         !is.na(currentSmk)) |> 
  count(FamilyNumber) |> 
  count(n)
Ndz     <- countDZ[countDZ$n==2,"nn"]
single2 <- countDZ[countDZ$n==1,"nn"]

rDZ <- polycor::polychor(twin_dat[twin_dat$zygo==2,]$T1_currentSmk,
                         twin_dat[twin_dat$zygo==2,]$T2_currentSmk)

NEdz <- (2*Ndz) / (1 + rDZ)


## Total
NE = NEmz + NEdz + single1 + single2

## F
R2 <- 0.05067724
K <- 2 # (PRS + intercept)
Fstat <- ( R2 / (1-R2) ) * ( (NE-K) / (K-1) )
Fstat


#### Former Smk ------------------------------------------------------------

## MZ twins
countMZ <- twin_long_prs |> 
  filter(zygo == 1,
         !is.na(formerSmk)) |> 
  count(FamilyNumber) |> 
  count(n)
Nmz     <- countMZ[countMZ$n==2,"nn"]
single1 <- countMZ[countMZ$n==1,"nn"]

rMZ <- polycor::polychor(twin_dat[twin_dat$zygo==1,]$T1_formerSmk,
                         twin_dat[twin_dat$zygo==1,]$T2_formerSmk)

NEmz <- (2*Nmz) / (1 + rMZ)

## DZ twins
countDZ <- twin_long_prs |> 
  filter(zygo == 2,
         !is.na(formerSmk)) |> 
  count(FamilyNumber) |> 
  count(n)
Ndz     <- countDZ[countDZ$n==2,"nn"]
single2 <- countDZ[countDZ$n==1,"nn"]

rDZ <- polycor::polychor(twin_dat[twin_dat$zygo==2,]$T1_formerSmk,
                         twin_dat[twin_dat$zygo==2,]$T2_formerSmk)

NEdz <- (2*Ndz) / (1 + rDZ)


## Total
NE = NEmz + NEdz + single1 + single2

## F
R2 <- 0.02020994
K <- 2 # (PRS + intercept)
Fstat <- ( R2 / (1-R2) ) * ( (NE-K) / (K-1) )
Fstat



## Need the data in long format to merge with CpG ----------------------

twin_dat_long <-  twin_dat |>  
  pivot_longer(cols = -c(FamilyNumber, zygo),  names_to = c("Twin_nr",".value"),
               names_pattern = "T(1|2)_(.*)")
dim(twin_dat_long)

head(twin_dat_long)
table(twin_dat_long$Twin_nr, useNA = "ifany")
####

## Save ---------------------

filename <- paste0(phenDir,"twin_dat_long_mrdoc_resid.csv")
write.csv(twin_dat_long, file = filename, quote = F, row.names = F)



