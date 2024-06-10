## Madhur Singh

library(data.table)
setDTthreads(12)
library(tidyverse)

baseDir   <- "/base_dir/" 
outDir    <- paste0(baseDir,"data/PRS_mQTL/Rsq/")
tempDir   <- paste0(baseDir,"scripts/temp/")
logDir    <- paste0(baseDir,"logs/")
mqtlDir   <- paste0(baseDir,"out/prs_mqtls/")


## Load the Rsq result files ----------------------------------------
### Keeping all CpG sites - even if no SNP had p<0.05 ---------------------------
### Will filter later by F-statistic of the allelic score

## If the minimum pval in the Sum Stats was >0.05 (the pval threshold used for C&T)
## No pval thresholding was applied
## This gives 5 types of PRS across CpGs, as PLINK changes the PRS name according to pval threshold
# 1:      cis_prs_LDr2_1_p  
# 2: cis_prs_LDr2_1_p.0.05 
# 3:      cis_prs_LDr2_5_p
# 4: cis_prs_LDr2_5_p.0.05 
# 5:    top_cis_mQTL_score 
## Keeping the first 16 characters of the PRS name gives us 3 types of PRS (as intended)


prsOut <- fread(paste0(outDir,"MRDoC_CpG_cis_PRS_Rsq.dat"))
prsOut$prs <- substr(prsOut$PRS,start = 1, stop = 16)

snp_per_prs <- fread( paste0(mqtlDir,"n_snps_p_0.05_by_clumpLDr2.dta"), header = T)
snp_per_prs <- snp_per_prs |> 
  mutate(prs = case_when(
    clumpLDr2 == 0.5 ~ "cis_prs_LDr2_5_p",
    clumpLDr2 == 0.1 ~ "cis_prs_LDr2_1_p",
  ))

prsOut <- prsOut |> 
  left_join(snp_per_prs)

prsOut <- prsOut |> 
  mutate(nSNP = ifelse(is.na(nSNP), 1, nSNP))   # nSNP = 1 for the mQTL
## nSNP is the number of clumped SNPs with p<0.05 --> those included in the PRS
## If nSNP is 0, then no pval threshold was applied --> all SNPs were included

## Examine
dim(prsOut) 
glimpse(head(prsOut))

## Add a "priority" variable
## If multiple scores have the same Rsq, the simplest score is prioritized
prsOut <- prsOut |> 
  mutate(priority = case_when(
    clumpLDr2 == 0.5 ~ 1,
    clumpLDr2 == 0.1 ~ 2,
    is.na(clumpLDr2) ~ 3
  ))

## Preference order: mQTL > LDr<0.1 > LDr<0.5

prsOut <- prsOut |> 
  group_by(CpG) |> 
  arrange(desc(R2percent), desc(priority)) |> 
  mutate( performance = row_number(),
          top = ifelse(performance == 1, TRUE, FALSE)
  ) |> 
  ungroup() |> 
  arrange(CpG)


## Save the top PRS list -------------------------------------------

prsTop <- prsOut |> 
  filter(top == TRUE) |> 
  select(CpG, PRS, Estimate, RobustSE, pvalue, R2percent, Fstat, nSNP, negEst)

## Example PRS File name = "cg00000714_godmc_cis_prs_LDr2_1_p.0.05.profile"
## Create a File Name Variable
prsTop <- prsTop |> 
  mutate(fileName = paste0(CpG,"_godmc_",PRS,".profile"))

fwrite(prsTop, file = paste0(outDir,"MRDoC_CpG_IV_with_maxR2.dat"), sep = "\t",
       col.names = T, row.names = F, quote = F)

prsTop <- fread(paste0(outDir,"MRDoC_CpG_IV_with_maxR2.dat"), data.table = F)

