## Prepare a list of CpGs with GWAS sums stats from Go-DMC
## Subset these for those with known assoc. with Smoking in Joehanes et al EWAS
## Madhur Singh,

## Libraries
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)


## Working Directory on NTR=Compute2
baseDir <- "/base_dir/" 
dataDir <- paste0(baseDir,"data/GoDMC/")
phenoDir <- paste0(baseDir,"data/phenos/")
mqtlDir <- paste0(dataDir, "cpg_mqtls/")
cisDir  <- paste0(dataDir, "cpg_mqtls/cis_mqtls/")


## GoDMC Sum Stats ###############################################################

## Load the GoDMC GWAS Summ Stats file - matched with NTR and KGP
godmcRes <- fread(paste0(dataDir,"GoDMC_sumStats_noBIOS_NTR_KGP_matched.dat"), 
                  stringsAsFactors = F, data.table = F)

cpgIDs <- godmcRes |> 
  summarise(ID = unique(cpg))

nrow(cpgIDs)
## 13,275 CpGs with valid GWAS sum stats

## Net Smk-associated autosomal CpGs used for MR-DoC = 16,940
## Of these, we have mQTLs for 13,275 CpGs (78.4%)

## Save these ids
fwrite(cpgIDs, paste0(phenoDir,"godmc_cpg_ids.dta"), 
       quote = F, sep = "\t", row.names = F, col.names = F)

cpg_cis_IDs <- godmcRes |> 
  filter(cistrans == TRUE) |> 
  summarise(ID = unique(cpg))

## save
fwrite(cpg_cis_IDs, paste0(phenoDir,"godmc_cis_cpg_ids.dta"), 
       quote = F, sep = "\t", row.names = F, col.names = F)


## Save the sumstats by CpG #######################################

glimpse(cpgIDs)
# Rows: 13,275

for (i in 1:nrow(cpgIDs)) {
  
  cpg_i <- cpgIDs$ID[i]
  mqtls <- godmcRes |> 
    filter(cpg == cpg_i) |> 
    select(kgpSNP, SNP, chr, pos, beta_a1, se, pval, samplesize,
           allele1, allele2, freq_a1, freq_se, cistrans, num_studies) |> 
    arrange(pval)
  
  nGWsnp <- nrow(mqtls)
  if(nGWsnp > 0) { 
    fwrite(mqtls, file = paste0(mqtlDir,cpg_i,"_GoDMC_sumStats_noBIOS_QC.dat"), 
           sep = "\t", quote = F, row.names = F, col.names = T) 
  }
  
  if(i %%100 == 0) print(paste("Saved Sum Stats for cpg",i,"of",nrow(cpgIDs)))
  
}


#### Cis mQTLs only ##########################

for (i in 1:nrow(cpgIDs)) {
  
  cpg_i <- cpgIDs$ID[i]
  mqtls <- godmcRes |> 
    filter(cpg == cpg_i,
           cistrans == TRUE) |> 
    select(kgpSNP, SNP, chr, pos, beta_a1, se, pval, samplesize,
           allele1, allele2, freq_a1, freq_se, cistrans, num_studies) |> 
    arrange(pval)
  
  
  nGWsnp <- nrow(mqtls)
  if(nGWsnp > 0) { 
    fwrite(mqtls, file = paste0(cisDir,cpg_i,"_GoDMC_cis_sumStats_noBIOS_QC.dat"), 
           sep = "\t", quote = F, row.names = F, col.names = T) 
  }
  
  if(i %%100 == 0) print(paste("Saved Sum Stats for cpg",i,"of",nrow(cpgIDs)))
  
}



## END ############

