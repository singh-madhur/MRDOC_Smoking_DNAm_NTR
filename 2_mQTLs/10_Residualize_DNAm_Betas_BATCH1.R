## Residualize the DNAm Beta score

# clear the workspace
rm(list = ls()); gc()

# Libraries
library(tidyverse)

# directories
baseDir      <- "/base_dir/"
phenoDir     <- paste0(baseDir,"data/phenos/")
outDir       <- paste0(baseDir,"data/DNAm_resid/")
methDir      <- "/base_dir/Methylation/Blood_450karray_freeze2021/"
outFile      <- paste0(outDir,"residualized_beta_blood450k_BATCH1.RData")


## Covariate twin data --------------------------
## EUR-ancestry pheno + PRS + covars
twin_covar_data <- read.csv(paste0(phenoDir,"twin_dat_long_mrdoc_resid.csv"), header = T)

## drop the missing T2
long_twin_data <- twin_covar_data |> 
  filter(!is.na(fisnumber))

long_twin_data <- long_twin_data |> 
  mutate(Sample_Plate = as_factor(Sample_Plate))

long_twin_data <- data.frame(long_twin_data)

## Methylation Data -----------------------

load(paste0(methDir,"Blood.450k.MvaluesFunNorm.062021_BATCH1.RData"))

CpGs <- colnames(tBeta_BATCH1)
head(CpGs)

nCpG <- length(CpGs)
nCpG

tBeta_BATCH1 <- data.frame(tBeta_BATCH1)
tBeta_BATCH1$IdatName <- rownames(tBeta_BATCH1)

dat <- long_twin_data |> 
  left_join(tBeta_BATCH1, by = "IdatName")

dat <- data.frame(dat)
rownames(dat) <- dat$IdatName

rm(tBeta_BATCH1); gc()

## Regression ---------------------

Nind       <- nrow(dat)
resBeta    <- as.data.frame(c(rownames(dat)))
colnames(resBeta) <- "IdatName"
nrow(resBeta) == Nind
Nind
nCpG

for (i in 1:nCpG) {
  
  cpgT <- CpGs[i]
  dat$cpgT <- dat[, colnames(dat) %in% cpgT ]
  datT <- dat |> 
    select(cpgT, female, age, Sample_Plate, Array_rownum, 
           Neut_Perc, Mono_Perc, Eos_Perc, IdatName) |> 
    filter(!is.na(cpgT), 
           ## There's some scores recorded as -Inf
           cpgT > -Inf, 
           cpgT < Inf)
  datT <- data.frame(datT)
  rownames(datT) <- datT$IdatName
  
  linMod <- lm(cpgT ~ female + age + Sample_Plate + Array_rownum + 
                 Neut_Perc + Mono_Perc + Eos_Perc, data = datT)
  
  res <- as.data.frame(cbind(residual = linMod$residuals))
  ## add variable with Idatname
  res$IdatName <- rownames(res)
  ## rename the residuals as the CpG ID
  colnames(res)[colnames(res) %in% "residual"] <- cpgT
  ## merge with the output object
  resBeta <- merge(resBeta, res, by = "IdatName", all = TRUE)
  
  dat$cpgT <- NULL; rm(datT)
  
  if(i %% 100 == 0) { print(paste("Ran CpG",i,"of",nCpG)) }
  
}


nrow(resBeta) == nrow(dat)
table(resBeta$IdatName %in% rownames(dat))
table(rownames(dat) %in% resBeta$IdatName)
rownames(resBeta) <- resBeta$IdatName
resBeta$IdatName <- NULL

ncol(resBeta) == nCpG
table(colnames(resBeta) == CpGs)
dim(resBeta)

save(resBeta, file= outFile)

rm(list = ls(all = TRUE)); gc()
