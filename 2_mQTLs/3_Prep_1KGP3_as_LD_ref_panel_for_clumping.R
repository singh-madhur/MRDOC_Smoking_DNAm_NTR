## Filter 1KGP3 to keep EUR and SNPs overlapping with the GWAS Sum Stats and NTR

## GoDMC GWAS used 1KGP3 as the reference panel
## NTR used HRC as the reference panel

rm(list = ls()); gc()

## Libraries
library(data.table)
library(dplyr)
library(tidyr)

## Working Directory on NTR-Compute2
baseDir <- "/base_dir/" 
dataDir <- paste0(baseDir,"MRDoC_Smk_DNAm/data/GoDMC/")
phenoDir <- paste0(baseDir,"MRDoC_Smk_DNAm/data/phenos/")
qcGenoDor <- paste0(baseDir, "MRDoC_Smk_DNAm/out/geno/")
genoDir <- "/base_dir/PRSScoring/"
tempDir <- paste0(baseDir,"MRDoC_Smk_DNAm/scripts/temp/")
logDir  <- paste0(baseDir,"MRDoC_Smk_DNAm/logs/")

plinkDir <- paste0(baseDir,"Util/plink1.9/")

# First create a filtered 1kgp files
# Keep only EUR ids

## 1. Filter 1KGP to EUR only ###################################################


dateCheck <- gsub( pattern="-", replacement="", x=Sys.Date() )
jobName   <- paste0("filter_1kgp3_for_mqtl_prs_",dateCheck)
shellName <- paste0(tempDir,jobName,".sh",sep="")
logFile   <- paste0(logDir,jobName)

bedFile  <- paste0(genoDir,"1KG_P3V5A_B37_FixJJV2_CHR1-22X")
keepFile <- paste0(genoDir,"EUR_samples.dat")
outFile  <- paste0(qcGenoDor,"1KGP_EUR_SNPs")

plinkCom <- paste0(
  plinkDir,"plink \\
--bfile ", bedFile, " \\
--keep ", keepFile, " \\
--snps-only just-acgt \\
--maf 0.001 \\
--make-bed \\
--write-snplist \\
--out ", outFile
)

write(plinkCom, file = shellName)
system(paste0("chmod +x ",shellName))
system(paste0("nohup ",shellName," > ",logFile,".log 2> ",logFile,".err &"))


## 2. Check SNP overlap between 1KGP and GoDMC/NTR #############################

kgp_bim <- fread(paste0(outFile,".bim"),
                 sep = "\t",
                 col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"))

dim(kgp_bim)   
head(kgp_bim)

ntr_bim <- fread(paste0(qcGenoDor,"MRG16_PRSQCed_GoDMC_matched.bim"),
                 sep = "\t",
                 col.names = c("CHR", "bSNP", "bCM", "BP", "bA1", "bA2"))

dim(ntr_bim)  
head(ntr_bim)

ntr_kgp <- ntr_bim |> 
  left_join(kgp_bim, by = c("CHR","BP"))

head(ntr_kgp)

ntr_kgp <- ntr_kgp |> 
  mutate(matchSNP = ifelse(bSNP == SNP, TRUE, FALSE))

table(ntr_kgp$matchSNP, useNA = "ifany")


## 3. Match SNPs by CHR_BP_A1_A2 ##################################################

kgp_bim <- kgp_bim |> 
  # New variable = chr + bp + A1 + A2 
  tidyr::unite("snpID", c(CHR,BP,A1,A2), remove = F)

ntr_bim <- ntr_bim |> 
  # New variable = chr + bp + A1 + A2
  tidyr::unite("snpID", c(CHR,BP,bA1,bA2), 
               remove = F)

ntr_bim <- ntr_bim |> 
  # SNP id with reverse alleles = chr + bp + A2 + A1 
  tidyr::unite("snpIDrev", c(CHR,BP,bA2,bA1), 
               remove = F)

# Function for finding the complementary allele
complement <- function(x) {
  switch (
    x,
    "A" = "T",
    "C" = "G",
    "T" = "A",
    "G" = "C",
    return(NA)
  )
}

ntr_bim$C.A1 <- sapply(ntr_bim$bA1, complement)
ntr_bim$C.A2 <- sapply(ntr_bim$bA2, complement)

ntr_bim <- ntr_bim |> 
  # SNPid with complementary strand = chr + bp + C.A1 + C.A2
  tidyr::unite("snpIDcomp", c(CHR,BP,C.A1,C.A2), 
               remove = F)

ntr_bim <- ntr_bim |> 
  # SNPid with reversed complementary strand = chr + bp + C.A2 + C.A1
  tidyr::unite("snpIDcomp_rev", c(CHR,BP,C.A2,C.A1), 
               remove = F)


#### Merge by fwd SNP ##############################################################

nrow(ntr_bim) 

ntr_kgp_fwd <- ntr_bim |> 
  select(snpID, bSNP) |> 
  inner_join(kgp_bim)

glimpse(ntr_kgp_fwd)



#### Merge by rev SNP  ############################################################

ntr_kgp_rev <- ntr_bim |> 
  select(snpIDrev, bSNP) |> 
  mutate(snpID = snpIDrev) |> 
  select(snpID, bSNP) |> 
  inner_join(kgp_bim)

glimpse(ntr_kgp_rev)


#### Merge by compl SNP ########################################################

ntr_kgp_comp <- ntr_bim |> 
  select(snpIDcomp, bSNP) |> 
  mutate(snpID = snpIDcomp) |> 
  select(snpID, bSNP) |> 
  inner_join(kgp_bim)

nrow(ntr_kgp_comp)
## 0

ntr_kgp_comp_rev <- ntr_bim |> 
  select(snpIDcomp_rev, bSNP) |> 
  mutate(snpID = snpIDcomp_rev) |> 
  select(snpID, bSNP) |> 
  inner_join(kgp_bim)

nrow(ntr_kgp_comp_rev)
## 0


## Reverse A1 in KGP, where needed
## Then filter KGP for overlapping SNPs


## 4. Update the Reverse Alleles in the BIM File #################################

tmp <- ntr_kgp_rev$A1
ntr_kgp_rev$A1 <- ntr_kgp_rev$A2    # Replace A1 with A2
ntr_kgp_rev$A2 <- tmp               # Replace A2 with A1

## Update the SNP id
ntr_kgp_rev <- ntr_kgp_rev |> 
  # New variable = chr + bp + A1 + A2 
  tidyr::unite("snpID", c(CHR,BP,A1,A2), remove = F)

head(ntr_kgp_rev)


####  Save the SNP IDs with reversed A1 #######################################
fwdSNP <- ntr_kgp_fwd |> 
  select(snpID, SNP, A1)

revSNP <- ntr_kgp_rev |> 
  select(snpID, SNP, A1)

ntr_kgp_snp <- rbind(fwdSNP, revSNP)

## List of SNPs to be extracted during LD clumping
snpID <- ntr_kgp_snp |> 
  select(SNP)

fwrite( snpID, file = paste0(qcGenoDor,"kgp_godmc_match.snplist"), sep = "\t",
        quote = F, row.names = F, col.names = F )


## This file can be used to assign A1 when reading in the KGP BED files 
fwrite(ntr_kgp_snp[,c("SNP", "A1")], file = paste0(qcGenoDor,"kgp_snps.a1"),
       quote = F, row.names = F, col.names = F, sep="\t" )



## 5. Add KGP rsIDs to the GWAS Sum Stats #######################################

## Will use this rsID for clumping

## Load the GoDMC GWAS Summ Stats file
godmcRes <- fread(paste0(dataDir,"GoDMC_sumStats_noBIOS_NTR_matched_cleaned.csv"), 
                  stringsAsFactors = F, data.table = F)

head(godmcRes)

godmcSNP <- godmcRes |> 
  summarise(snpID = unique(snpID)) 
nrow(godmcSNP)


## merge by snpID
godmcRes_kgp <- kgp_bim |> 
  mutate(kgpSNP = SNP) |> 
  select(snpID, kgpSNP) |> 
  left_join(godmcRes, by = "snpID") |> 
  relocate(cpg) 



#### Save the QCed Sum Stats ############################################################

outFile <- paste0(dataDir, "GoDMC_sumStats_noBIOS_NTR_KGP_matched.dat" )
fwrite(godmcRes_kgp, file = outFile, sep = "\t", quote = F, row.names = F, col.names = T)




## 6. Make the Filtered KGP with Overlapping SNPs ##########################################

dateCheck <- gsub( pattern="-", replacement="", x=Sys.Date() )
jobName   <- paste0("filter_1kgp3_with_snp_overlap_with_godmc_",dateCheck)
shellName <- paste0(tempDir,jobName,".sh",sep="")
logFile   <- paste0(logDir,jobName)

bedFile  <- paste0(qcGenoDor,"1KGP_EUR_SNPs")
snpFile  <- paste0(qcGenoDor,"kgp_godmc_match.snplist")
snpA1    <- paste0(qcGenoDor,"kgp_snps.a1")
outFile  <- paste0(qcGenoDor,"1KGP_GoDMC_SNPs")

plinkCom <- paste0(
  plinkDir,"plink \\
--bfile ", bedFile, " \\
--extract ", snpFile, " \\
--a1-allele ", snpA1," \\
--maf 0.001 \\
--make-bed \\
--out ", outFile
)

write(plinkCom, file = shellName)
system(paste0("chmod +x ",shellName))
system(paste0("nohup ",shellName," > ",logFile,".log 2> ",logFile,".err &"))

####
