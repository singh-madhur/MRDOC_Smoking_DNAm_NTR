## Harmonize the variables in the NTR Geno file and the GoDMC Sum Stats
## Madhur Singh

rm(list = ls()); gc()

## Libraries
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

## Working Directory on NTR-Compute2
baseDir <- "/base_dir/" 
dataDir <- paste0(baseDir,"data/GoDMC/") 
mqtlDir <- paste0(dataDir, "cpg_mqtls/")
qcGenoDor <- paste0(baseDir, "out/geno/")
phenDir <- paste0(baseDir,"data/phenos/")
tempDir <- paste0(baseDir,"scripts/temp/")
logDir  <- paste0(baseDir,"logs/")
genoDir <- "/base_dir/PRSScoring/"

## Path to PLINK
plinkDir <- "/base_dir/Util/plink1.9/"


## Get the IDs to keep for PRS ##################################################

#### Pheno File ####################################

phenoFile <- fread(paste0(phenDir,"twin_dat_long_mrdoc_resid.csv"))
phenoFile <- data.frame(phenoFile)
dim(phenoFile)
str(phenoFile)

mrdocIDs <- phenoFile |> 
  filter(!is.na(fisnumber), !is.na(resPRS)) |> 
  select( fisnumber) |> 
  mutate(IID = fisnumber)


#### Geno FAM file ##########################################

genoFam <- fread(paste0(genoDir,"MRG16_PRSQCed_C1-22andX_V2.fam"), 
                 sep = " ", 
                 col.names = c("FID","IID","PID","MID","sex","pheno"))
genoFam <- data.frame(genoFam)

keepIDs <- mrdocIDs |> 
  left_join(genoFam, by = "IID") |> 
  filter(!is.na(FID)) |> 
  select(FID, IID) 


write.table(keepIDs, file = paste0(phenDir,"keepIDs_for_mrdoc.dta"), 
            quote = F, row.names = F, col.names = F)


droppedIDs <- phenoFile |>
  mutate(IID = fisnumber) |> 
  left_join(genoFam, by = "IID") |> 
  filter(!is.na(fisnumber), is.na(FID)) |> 
  select(FamilyNumber, fisnumber) 


write.table(droppedIDs, file = paste0(phenDir,"droppedIDs_with_missing_geno.dta"), 
            quote = F, row.names = F, col.names = T)


## Generate Final QCed target data ###############################################

dateCheck <- gsub( pattern="-", replacement="", x=Sys.Date() )
jobName   <- paste0("generate_qc_geno_for_mqtl_prs_",dateCheck)
shellName <- paste0(tempDir,jobName,".sh",sep="")
logFile   <- paste0(logDir,jobName)

bedFile  <- paste0(genoDir,"MRG16_PRSQCed_C1-22andX_V2")
keepFile <- paste0(phenDir,"keepIDs_for_mrdoc.dta")
dropSnp  <- paste0(qcGenoDor,"ntr_mqtls.mismatch")
recodSNP <- paste0(qcGenoDor,"ntr_mqtls.a1")
outFile  <- paste0(qcGenoDor,"MRG16_PRSQCed_GoDMC_matched")

plinkCom <- paste0(
  plinkDir,"plink \\
--bfile ", bedFile, " \\
--keep ", keepFile, " \\
--exclude ", dropSnp, " \\
--a1-allele ", recodSNP," \\
--make-bed \\
--write-snplist \\
--out ", outFile
)

write(plinkCom, file = shellName)
system(paste0("chmod +x ",shellName))
system(paste0("nohup ",shellName," > ",logFile,".log 2> ",logFile,".err &"))


