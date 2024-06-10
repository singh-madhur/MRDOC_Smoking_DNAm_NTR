## Madhur Singh

# clear the workspace
rm(list = ls()); gc()

# Libraries
library(data.table)
library(tidyverse)

## Directories
baseDir      <- "/base_dir/"
phenoDir     <- paste0(baseDir,"data/phenos/")
mqtlDir      <- paste0(baseDir,"data/PRS_mQTL/top_score/") 
mqtlResidDir <- paste0(baseDir,"data/mQTL_resid/")
godmcDir     <- paste0(baseDir,"data/GoDMC/batches/")

## Covariate twin data -----------------------------------
## EUR-ancestry pheno + PRS + covars
twin_covar_data <- fread(paste0(phenoDir,"twin_dat_long_mrdoc_resid.csv"), data.table = F)

## drop the missing T2
phenoDat <- twin_covar_data |> 
  filter(!is.na(fisnumber)) |> 
  mutate(fisnumber = as.numeric(fisnumber))

## mQTL scores -------------------------------------------------
mqtlFilesRaw <- list.files(mqtlDir)
mqtlFiles <- mqtlFilesRaw[grep(".profile",mqtlFilesRaw)]

## CpG IDs
CpGs <- fread(paste0(phenoDir,"godmc_cis_cpg_ids.dta"), header = F)
colnames(CpGs) <- "ID"
nCpG <- length(CpGs$ID)

## Start Loop ----------------------------------------------------

for (ii in 1:nCpG) {
  
  cpgT  <- CpGs[ii]
  mqtlT <- mqtlFiles[grep(cpgT,mqtlFiles)]
  
  mQTL <- fread(paste0(mqtlDir,mqtlT), header = T)
  mQTL <- mQTL |> 
    mutate(fisnumber = as.numeric(IID),
           mQTL = SCORE) |> 
    select(fisnumber, mQTL)

  ## merge
  mQTL_info <- mQTL |> 
    inner_join(phenoDat, by = "fisnumber")
  
  ## Residualize the mQTL for PCs and genotyping platform #############
  
  mQTL_info <- data.frame(mQTL_info)
  rownames(mQTL_info) <- mQTL_info$fisnumber
  ## covariates = platform and the first 10 PCs
  modPRS <- lm(mQTL ~ Platform_2 + Platform_3 + PC1_sc + PC2_sc + PC3_sc + 
                 PC4_sc + PC5_sc + PC6_sc + PC7_sc + PC8_sc + PC9_sc + PC10_sc,
               data = mQTL_info)
  resMQTL <- as.data.frame(modPRS$residuals)
  colnames(resMQTL) <- "resMQTL"
  resMQTL$fisnumber <- rownames(resMQTL)
  resMQTL$resMQTL <- as.numeric(scale(resMQTL$resMQTL))
  
  ## Save the residualized, scaled mQTL
  filename <- paste0(mqtlResidDir,cpgT,"_godmc_resid_prs.dat")
  fwrite(resMQTL, file = filename, sep = "\t", quote = F, col.names = T, row.names = F)
  
  if (ii %% 50 == 0) {print(paste("Saved residualized score for cpg",ii,"of",nCpG))}
  
}  ## End CpG Loop

## END 

