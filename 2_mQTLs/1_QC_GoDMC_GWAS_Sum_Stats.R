## QC the GoDMC Sum Stats
## Madhur Singh

rm(list = ls()); gc()

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

## Working Directory on NTR-Compute2
baseDir <- "/base_dir/" 
dataDir <- paste0(baseDir,"data/GoDMC/")
phenDir <- paste0(baseDir,"data/phenos/")
genoDir <- "/base_dir/PRSScoring/"
qcGenoDor <- paste0(baseDir, "out/geno/")


## Sum Stats
gwasSumStat <- fread(paste0(dataDir, "GoDMC_sumStats_noBIOS_smkAssoc_noSNPnoXY_CpGs.csv"))
dim(gwasSumStat) 
names(gwasSumStat)

head(gwasSumStat)


## 0. Clean SNP IDs ############################################################

## The Sum Stats do not include rsIDs. 
## We will merge the GWAS sumstats with the BIM file by chr and BP position
## Create a new SNP id = chr_bp_a1_a2, which we can match in the geno bim file

gwasSumStat <- gwasSumStat |> 
  # Get the Chr and BP
  tidyr::separate(col = snp, sep = ":", 
                  into = c("chr","pos","snpType"), 
                  remove = FALSE, 
                  convert = TRUE) |> 
  # Remove "chr" prefix from the chr variable
  mutate(chr = gsub(pattern = "chr", replacement = "", x = chr)) |> 
  # New variable = chr + bp + A1 + A2
  tidyr::unite("snpID", c(chr,pos,allele1,allele2), 
               remove = F)

str(gwasSumStat)



## 1. Effect Allele ############################################################

## allele1
## beta coef. = beta_a1

## 2. Genome Build #############################################################

## Build 37
## 1000 Genomes as the reference panel 
## EUR ancestry


## 3. MAF ######################################################################

## Keeping only sites with MAF > 0.01 

summary(gwasSumStat$freq_a1)
## Min > 0.01, Max < 0.99



## 4. Get the SNP IDs ##############################################################

## This file has the sum stats for all CpGs --> Each SNP can appear multiple times
## Get the unique SNPs 
godmcSNP <- gwasSumStat |> 
  summarise(SNP = unique(snpID)) |> 
  tidyr::separate(col = SNP, sep = "_", 
                  into = c("chr","pos","allele1","allele2"), 
                  remove = FALSE, 
                  convert = TRUE)

godmcSNP <- godmcSNP |> 
  rename(snpID = SNP)

dim(godmcSNP) 

table(godmcSNP$chr)


## 5. Ambiguous SNPs #############################################################

godmcSNP <- godmcSNP |> 
  mutate( ambigSNP =
            ifelse( (allele1 == "A" & allele2 == "T") | 
                      (allele1 == "T" & allele2 == "A") | 
                      (allele1 == "G" & allele2 == "C") | 
                      (allele1 == "C" & allele2 == "G") ,
                    TRUE, FALSE)
  )

table(godmcSNP$ambigSNP)

godmcSNP_filt <- godmcSNP |> 
  filter(ambigSNP == FALSE)


## 6. Geno BIM file ################################################################

## Get similar snpID in the BIM file
## Merge with sum stats by snpID
## Save the corresponding rsID from the BIM file

## THe geno data has already been QCed

genobim <- fread(paste0(genoDir,"MRG16_PRSQCed_C1-22andX_V2.bim"),
                 sep = "\t",
                 col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"))
genobim <- data.frame(genobim)

genobim <- genobim |> 
  # New variable = chr + bp + A1 + A2
  tidyr::unite("snpID", c(CHR,BP,A1,A2), 
               remove = F)

genobim <- genobim |> 
  # SNPid with reverse direction = chr + bp + A2 + A1
  tidyr::unite("snpIDrev", c(CHR,BP,A2,A1), 
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

genobim$C.A1 <- sapply(genobim$A1, complement)
genobim$C.A2 <- sapply(genobim$A2, complement)

genobim <- genobim |> 
  # SNPid with complementary strand = chr + bp + C.A1 + C.A2
  tidyr::unite("snpIDcomp", c(CHR,BP,C.A1,C.A2), 
               remove = F)

genobim <- genobim |> 
  # SNPid with reversed complementary strand = chr + bp + C.A2 + C.A1
  tidyr::unite("snpIDcomp_rev", c(CHR,BP,C.A2,C.A1), 
               remove = F)

head(genobim)

## 7a. Merge and ID the SNPs with matching alleles in BIM File ############################

## Merge by snpID

nrow(godmcSNP_filt) 
godmcSNP_NTRgeno <- genobim |> 
  select(snpID, SNP) |> 
  inner_join(godmcSNP_filt)

glimpse(godmcSNP_NTRgeno)

nrow(godmcSNP_NTRgeno) / nrow(godmcSNP_filt)

## 7b. Merge and ID the SNPs with reverse alleles in BIM File ############################

## Merge by snpIDrev
godmcSNP_NTRgeno_rev <- genobim |> 
  select(snpIDrev, SNP) |> 
  mutate(snpID = snpIDrev) |> 
  select(snpID, SNP) |> 
  inner_join(godmcSNP_filt)

glimpse(godmcSNP_NTRgeno_rev)

nrow(godmcSNP_NTRgeno_rev) / nrow(godmcSNP_filt)


## 7c. Merge and ID the SNPs with complementary strand in BIM File ############################

## Merge by snpIDcomp

godmcSNP_NTRgeno_comp <- genobim |> 
  select(snpIDcomp, SNP) |> 
  mutate(snpID = snpIDcomp) |> 
  select(snpID, SNP) |> 
  inner_join(godmcSNP_filt)

glimpse(godmcSNP_NTRgeno_comp)


## 7d. Merge and ID the SNPs with complementary strand and opposite allele #################

## Merge by snpIDcomp_rev
godmcSNP_NTRgeno_comp_rev <- genobim |> 
  select(snpIDcomp_rev, SNP) |> 
  mutate(snpID = snpIDcomp_rev) |> 
  select(snpID, SNP) |> 
  inner_join(godmcSNP_filt)

glimpse(godmcSNP_NTRgeno_comp_rev)

## Since there are only 2 SNPs (from the same genomic locus) with flipped strands - just dropping these.

## Total
(nrow(godmcSNP_NTRgeno) + nrow(godmcSNP_NTRgeno_rev)) / nrow(godmcSNP_filt)


## 8. Update the Reverse Alleles in the BIM File #################################

# Update the recode SNPs
recode_snps <- genobim$snpIDrev %in% godmcSNP_NTRgeno_rev$snpID
table(recode_snps) 

tmp <- genobim[recode_snps,]$A1
genobim[recode_snps,]$A1 <- genobim[recode_snps,]$A2    # Replace A1 with A2
genobim[recode_snps,]$A2 <- tmp                         # Replace A2 with A1


####  Save the SNP IDs with reversed A1 ##################################

## This file can be used to assign A1 when reading in the BED files in PLINK
fwrite(genobim[,c("SNP", "A1")], file = paste0(qcGenoDor,"ntr_mqtls.a1"),
       quote = F, row.names = F, col.names = F, sep="\t" )



## 9. Merge rsIDs back with the GWAS Sum Stats #######################################

dim(godmcSNP_NTRgeno)       
dim(godmcSNP_NTRgeno_rev)   
## There should be no overlap between the two
godmcSNP_NTRgeno |> 
  select(snpID) |> 
  inner_join(godmcSNP_NTRgeno_rev)
# <0 rows>

godmcSNP_ID <- godmcSNP_NTRgeno |> 
  rbind(godmcSNP_NTRgeno_rev)


## Keep the GWAS sum stats with QCed, overlapping SNPs

godmcSumStat <- godmcSNP_ID |> 
  select(snpID, SNP) |> 
  left_join(gwasSumStat, by = "snpID") |> 
  relocate(cpg) 


#### Save the QCed Sum Stats ############################################################

outFile <- paste0(dataDir, "GoDMC_sumStats_noBIOS_NTR_matched_cleaned.csv" )
fwrite(godmcSumStat, file = outFile, quote = F, row.names = F, col.names = T)


## 10. Save the non-overlapping SNPs in NTR ########################################

## These will be excluded when reading in the BED file in PLINK

genobim <- genobim |> 
  mutate(mismatch = ifelse( ! (snpID %in% godmcSNP_ID$snpID | snpIDrev %in% godmcSNP_ID$snpID) , 
                            TRUE,
                            FALSE))

mismatch_snp <- genobim |> 
  filter(mismatch == TRUE) |> 
  select(SNP)

## save
fwrite( mismatch_snp, file = paste0(qcGenoDor,"ntr_mqtls.mismatch"), sep = "\t",
        quote = F, row.names = F, col.names = F )


# END. ===================================
