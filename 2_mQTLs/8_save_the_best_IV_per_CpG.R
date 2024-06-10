## Select the best IV and move to a common directory
## Madhur Singh

library(data.table)
setDTthreads(12)
library(tidyverse)

baseDir   <- "/base_dir/" 
rsqDir    <- paste0(baseDir,"data/PRS_mQTL/Rsq/")
outDir    <- paste0(baseDir,"data/PRS_mQTL/top_score/")
prsDir    <- paste0(baseDir,"out/prs_mqtls/scores/")

## Load the file names and R2
prsRsq <- fread(paste0(rsqDir,"MRDoC_CpG_IV_with_maxR2.dat"), header = T)
## also copy the file over to outDir
system(paste0("cp ",rsqDir,"MRDoC_CpG_IV_with_maxR2.dat ",outDir))

## PRS Files
prsFilesRaw <- list.files(prsDir)
prsFiles    <- prsFilesRaw[grep("profile",prsFilesRaw)]

nCpG <- nrow(prsRsq)

## Copy the PRS to the final directory
for (i in 1:nCpG) {
  prsFileT <- prsRsq$fileName[i]
  system(paste0("cp ",prsDir,prsFileT," ",outDir))
  if(i %% 100 == 0) print(paste("Copied file",i,"of",nCpG))
}



