## Estimate the variance in DNAm explained by the PRS/mQTL
## Madhur Singh

rm(list = ls()); gc()

# Libraries
library(data.table)
library(tidyverse)
library(gee)          # for R2 for continuous traits

# Directories
baseDir   <- "/base_dir/"
phenoDir  <- paste0(baseDir,"data/phenos/")
prsDir    <- paste0(baseDir,"out/prs_mqtls/scores/")
godmcDir  <- paste0(baseDir,"data/GoDMC/batches/")
mqtlDir   <- paste0(baseDir,"out/prs_mqtls/")
outDir    <- paste0(baseDir,"data/PRS_mQTL/Rsq/")
methDir   <- "/base_dir/Methylation/Blood_450karray_freeze2021/"


## Covariate twin data -----------------------------------
## EUR-ancestry pheno + PRS + covars
twin_covar_data <- read.csv(paste0(phenoDir,"twin_dat_long_mrdoc_resid.csv"), header = T)

## drop the missing T2
phenoDat <- twin_covar_data |> 
  filter(!is.na(fisnumber))


## DNAm Beta ####

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

cpgDat <- loadRData(paste0(methDir,"Blood.450k.BetavaluesFunNorm.062021.RData"))
methIDs <- colnames(cpgDat)
cpgDat$IdatName <- rownames(cpgDat)

phenoDat <- cpgDat |> 
  right_join(phenoDat)

dim(phenoDat)

## Need rMZ and rDZ
aceDir <- paste0(baseDir,"data/Catalogue_ACE_estimates/")
aceEst <- fread(paste0(aceDir,"41467_2016_BFncomms11115_MOESM54_ESM.txt"),
                header = T, data.table = F)

aceOut <- aceEst |> 
  filter(cgid %in% methIDs)


## PRS ####

## R2 in GEE ###########################

prsFilesRaw <- list.files(prsDir)
prsFiles    <- prsFilesRaw[grep("profile",prsFilesRaw)]

nCPG <- ncol(cpgDat)

for (jj in 1:nCPG) {          ## Start CpG Loop
  
  cpgT <- colnames(cpgDat)[jj]
  phenoDat$cpgT <- phenoDat[,cpgT]
  
  phenoFilt <- phenoDat |> 
    filter(!is.na(cpgT), !is.na(Neut_Perc), !is.na(Mono_Perc), !is.na(Eos_Perc)) |> 
    arrange(FamilyNumber)
  
  ## Effective N for Fstat -------------------------------------------------------
  ## Effective sample size, given N pairs of twins
  ## Minica et al, https://doi.org/10.1038/mp.2014.121 
  # NEmz = (2*Nmz) / (1+rMZ)
  # NEdz = (2*Ndz) / (1+rDZ)
  
  ## MZ twins
  countMZ <- phenoFilt |> 
    filter(zygo == 1) |> 
    count(FamilyNumber) |> 
    count(n)
  Nmz <- countMZ[countMZ$n==2,"nn"]
  single1 <- countMZ[countMZ$n==1,"nn"]
  
  
  ## DZ twins
  countDZ <- phenoFilt |> 
    filter(zygo == 2) |> 
    count(FamilyNumber) |> 
    count(n)
  Ndz <- countDZ[countDZ$n==2,"nn"]
  single2 <- countDZ[countDZ$n==1,"nn"]
  
  aceT <- aceOut |> 
    filter(cgid == cpgT)
  
  NEmz = (2*Nmz) / (1+aceT$rMZ)
  NEdz = (2*Ndz) / (1+aceT$rDZ)
  
  NE = NEmz + NEdz + single1 + single2
  
  
  
  ## Standardize the covariates #####
  phenoFilt <- phenoFilt |> 
    mutate(age_sc = as.numeric(scale(age)),
           Neut_Perc = as.numeric(scale(Neut_Perc)),
           Mono_Perc = as.numeric(scale(Mono_Perc)),
           Eos_Perc = as.numeric(scale(Eos_Perc)),
           Array_rownum = as.numeric(scale(Array_rownum)), 
           PC1_sc = as.numeric(scale(PC1_sc)),
           PC2_sc = as.numeric(scale(PC2_sc)),
           PC3_sc = as.numeric(scale(PC3_sc)),
           PC4_sc = as.numeric(scale(PC4_sc)),
           PC5_sc = as.numeric(scale(PC5_sc)),
           PC6_sc = as.numeric(scale(PC6_sc)),
           PC7_sc = as.numeric(scale(PC7_sc)),
           PC8_sc = as.numeric(scale(PC8_sc)),
           PC9_sc = as.numeric(scale(PC9_sc)),
           PC10_sc = as.numeric(scale(PC10_sc))
    )
  
  
  prsFilesCG <- prsFiles[grep(cpgT,prsFiles)]
  nPRS <- length(prsFilesCG)
  
  geeOut <- data.frame(matrix(NA_real_, nrow = nPRS, ncol = 8))
  colnames(geeOut) <- c("CpG","PRS","Estimate", "RobustSE", "pvalue", "R2", "R2percent", "Fstat")
  geeOut$CpG <- cpgT
  
  ## Start PRS Loop
  for (i in 1:nPRS) {          
    
    prsFileT <- prsFilesCG[i]
    prs <- fread(paste0(prsDir,prsFileT), header = T, stringsAsFactors = F, data.table = F)
    
    prsName <- gsub(pattern = paste0(cpgT,"_godmc_"), replacement = "", x = prsFileT)
    prsName <- gsub(pattern = ".profile", replacement = "", x = prsName)
    geeOut[i,"PRS"] <- prsName
    
    ## Merge with covariates
    ## Keep only the IDs with PRS
    phenoFilt <-  prs |> 
      mutate(fisnumber = as.double(IID)) |> 
      select(fisnumber, SCORE) |> 
      left_join(phenoFilt, by = "fisnumber") |> 
      arrange(FamilyNumber)
    
    
    ## Standardize continuous variables ---------------------------------------------
    
    phenoFilt <- phenoFilt |> 
      mutate(cpg_sc = as.numeric(scale(cpgT)), 
             PRS = as.numeric(scale(SCORE))
      ) |> 
      #### AGAIN ENSURE THE DATA ARE ARRANGED BY FAMILY NUMBER #########
    arrange(FamilyNumber)
    
    
    ## GEE ------------------------------------------------------------------------
    dat = phenoFilt 
    familynumber = as.numeric(dat$FamilyNumber)
    trait = as.numeric(dat$cpg_sc)
    PRS = as.numeric(dat$PRS)
    age = as.numeric(dat$age_sc)
    female = as.numeric(dat$female)              # dummy variable
    Neut_Perc = as.numeric(dat$Neut_Perc)
    Mono_Perc = as.numeric(dat$Mono_Perc)
    Eos_Perc = as.numeric(dat$Eos_Perc)
    Array_rownum = as.numeric(dat$Array_rownum)
    Sample_Plate = as.factor(dat$Sample_Plate)   # factor
    Platform_2 = as.numeric(dat$Platform_2)      # dummy variable
    Platform_3 = as.numeric(dat$Platform_3)      # dummy variable
    PC1 = as.numeric(dat$PC1_sc)
    PC2 = as.numeric(dat$PC2_sc)
    PC3 = as.numeric(dat$PC3_sc)
    PC4 = as.numeric(dat$PC4_sc)
    PC5 = as.numeric(dat$PC5_sc)
    PC6 = as.numeric(dat$PC6_sc)
    PC7 = as.numeric(dat$PC7_sc)
    PC8 = as.numeric(dat$PC8_sc)
    PC9 = as.numeric(dat$PC9_sc)
    PC10 = as.numeric(dat$PC10_sc)
    
    r1=gee(trait ~ PRS + age + female + Neut_Perc + Mono_Perc + Eos_Perc + Array_rownum + Sample_Plate + 
             Platform_2 + Platform_3 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
           id=familynumber, family=gaussian, corstr="exchangeable", maxiter=100, na.action=na.omit)
    coeff1 = summary(r1)$coeff
    estimate = coeff1["PRS","Estimate"]
    RobustSE= coeff1["PRS","Robust S.E."]  # use the robust standard error and robust z-score 
    pvalue = 2*pnorm(-abs(coeff1["PRS","Robust z"]))
    R2 = estimate^2  # proportion of variance in trait explained by PRS
    R2percent = 100*R2  # percentage of variance in trait explained by PRS
    K <- 2 # (PRS + intercept)
    Fstat <- ( R2 / (1-R2) ) * ( (NE-K) / (K-1) )
    out = c(estimate, RobustSE, pvalue, R2, R2percent, Fstat)
    
    geeOut[i,3:8] <- out
    
    phenoFilt$SCORE <- NULL
    phenoFilt$PRS <- NULL
    rm(r1,coeff1,estimate,RobustSE,pvalue,R2,R2percent,Fstat,out)
    
  }          
  ## End PRS Loop
  
  if(jj == 1) {
    prsOut <- geeOut
    rm(geeOut)
  } else {
    prsOut <- rbind(prsOut, geeOut)
    rm(geeOut)
  }
  
  if(jj == nCPG) {
    fwrite(prsOut, file = paste0(outDir,"MRDoC_CpG_cis_mQTL_Rsq.dat"), 
           col.names = T, row.names = F, quote = F)
  }
  
  print(paste("Computed R2 for CpG", jj, "of", nCPG))
  
  
}          ## End CpG Loop


## END R2 LOOP -----------

## Merge info on the number of SNPs
snp_per_prs <- fread( paste0(mqtlDir,"n_snps_p_0.05_by_clumpLDr2.dta"), header = T)

snp_per_prs <- snp_per_prs |> 
  mutate(PRS = case_when(
    clumpLDr2 == 0.5 ~ "cis_prs_LDr2_5_p.0.05",
    clumpLDr2 == 0.1 ~ "cis_prs_LDr2_1_p.0.05",
  ))

prsOut <- prsOut |> 
  left_join(snp_per_prs)

prsOut <- prsOut |> 
  mutate(nSNP = ifelse(is.na(nSNP), 1, nSNP))   # nSNP = 1 for the mQTL

## Add a "priority" variable
## If multiple scores have the same Rsq, the most parsimonious score is prioritized
prsOut <- prsOut |> 
  mutate(priority = case_when(
    clumpLDr2 == 0.5 ~ 1,
    clumpLDr2 == 0.1 ~ 2,
    is.na(clumpLDr2) ~ 3
  ))

## Preference order: mQTL > LDr<0.1 > LDr<0.5

dim(prsOut)
table(prsOut$PRS)
summary(prsOut$Estimate)
summary(prsOut$R2percent)
summary(prsOut$Fstat)

prsOut <- prsOut |> 
  group_by(CpG) |> 
  arrange(desc(R2percent), desc(priority)) |> 
  mutate( performance = row_number(),
          top = ifelse(performance == 1, TRUE, FALSE)
         ) |> 
  ungroup() |> 
  arrange(CpG)


prsOut |> 
  filter(top == TRUE) |> 
  count(PRS)

fwrite(prsOut, file = paste0(outDir,"MRDoC_CpG_cis_PRS_Rsq_sorted.dat"), 
       col.names = T, row.names = F, quote = F)

