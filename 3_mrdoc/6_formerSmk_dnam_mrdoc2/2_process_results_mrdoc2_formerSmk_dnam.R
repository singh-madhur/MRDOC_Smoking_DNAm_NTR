## Process the Results for MRDoC-2 Former Smk <--> DNAm
## Madhur Singh

## Keep only the CpGs associated w/ Former Smk

# clear the workspace
rm(list = ls()); gc()

# Libraries
library(data.table)
setDTthreads(16)

#Directories
baseDir <- "/base_dir/" 
outDir  <- paste0(baseDir,"out/mrdoc2_formerSmk_dnam/")
annoDir <- paste0(baseDir,"data/annotation/")
plotDir <- paste0(outDir,"plots/")
sumDir  <- paste0(outDir,"summary/")
phenoDir <- paste0(baseDir,"data/phenos/")

## Load results #########

batches <- 1:15

for (i in batches) {
  
  outFile   <- paste0(outDir,"mrdoc2_formerSmk_DNAm_BATCH",i,".dat")
  mrdocOut  <- fread(outFile, header = T)
  print(paste0(i,": ",nrow(mrdocOut)))
  
  if(i == 1) { 
    mrdocRes <- mrdocOut 
    rm(mrdocOut); gc()
  } else {
    mrdocRes <- rbind(mrdocRes, mrdocOut)
    rm(mrdocOut); gc()
  }
}

library(ggplot2)
library(forcats)

rownames(mrdocRes) <- mrdocRes$cpg
table(mrdocRes$statusCode, useNA = "ifany")

## Missing SEs
table(rowSums(is.na(mrdocRes)) == 0)

## Check IV R2 and Fstat #############

## IV R2 File
ivR2 <- fread(paste0(baseDir,"data/PRS_mQTL/Rsq/MRDoC_CpG_IV_with_maxR2.dat"), data.table = F)

mrdocRes <- mrdocRes |> 
  left_join(ivR2, by = c("cpg" = "CpG"))

mrdocResSub <- mrdocRes |> 
  filter(
    ## strong IV
    Fstat > 10,
    ## No missing SE
    rowSums(is.na(mrdocRes)) == 0,
  )
dim(mrdocResSub)


## Current-Smk-associated CpGs Inflation Factor using {bacon} ##############################

## this is to compare the lambdas in current and former smk
library(bacon)

#### Smk --> DNAm ###################

## Setting the seed for rng will ensure reproducibility
set.seed(2023) 
bc_g1 <- bacon(teststatistics = NULL, 
               effectsizes = mrdocResSub$g1_hat, standarderrors = mrdocResSub$g1_SE) 

#### DNAm --> Smk ###################

set.seed(2023) 
bc_g2 <- bacon(teststatistics = NULL, 
               effectsizes = mrdocResSub$g2_hat, standarderrors = mrdocResSub$g2_SE, na.exclude = T) 


## Former Smoking related CpG sites ############################

## EWAS results from Joehanes et al
## CpGs associated w/ Former v. Never Smk (at FDR < 0.05)
library(readxl)
formerSmkEWAS <- read_excel(path = paste0(phenoDir,"joehanes_etal_supplemental_tables.xlsx"), 
                            sheet = 6, skip = 2)

smkCpG <- formerSmkEWAS |> 
  select(`Probe ID`) |> 
  rename(ID = `Probe ID`) |> 
  data.frame()


## Keep CpGs associated w/ Former Smk and having Strong IVs ##########

mrdocResF <- mrdocRes |> 
  filter(
    ## Associated w/ Former Smk
    cpg %in% smkCpG$ID,
    ## strong IV
    Fstat > 10,
    ## No missing SE
    rowSums(is.na(mrdocRes)) == 0,
  )


##  Inflation Factor using {bacon} ##############################

library(bacon)

#### Smk --> DNAm ###################

## Setting the seed for rng will ensure reproducibility
set.seed(2023) 
bc_g1 <- bacon(teststatistics = NULL, 
               effectsizes = mrdocResF$g1_hat, standarderrors = mrdocResF$g1_SE) 


#### DNAm --> Smk ###################

set.seed(2023) 
bc_g2 <- bacon(teststatistics = NULL, 
               effectsizes = mrdocResF$g2_hat, standarderrors = mrdocResF$g2_SE, na.exclude = T) 



## Make QQ plots #############

## g1
jpeg(file = paste0(plotDir,"g1_former_MRDoC2_qqPlot.jpeg"), 
     width = 6, height = 6, units = "in", res = 300)
GWASTools::qqPlot(mrdocResF$g1_p, col="steelblue4", 
                  main=paste("MR-DoC2: Former Smoking to DNA Methylation")) 
dev.off()

## g2
jpeg(file = paste0(plotDir,"g2_former_MRDoC2_qqPlot.jpeg"), 
     width = 6, height = 6, units = "in", res = 300)
GWASTools::qqPlot(mrdocResF$g2_p, col="indianred4", 
                  main=paste("MR-DoC2: DNA Methylation to Former Smoking")) 
dev.off()



## Significance ############

nTests <- nrow(mrdocResF)
bonfp <- 0.05/nTests
bonfZ <- qnorm(p = 1-bonfp/2)


#### FDR #####################

library(qvalue)
mrdocResF$g1_qval <- qvalue(mrdocResF$g1_p)$qvalue
mrdocResF$g2_qval <- qvalue(mrdocResF$g2_p)$qvalue

mrdocResF <- mrdocResF |> 
  arrange(desc(abs(g2_Z)))


## Annotation ##########

load(paste0(annoDir,"manifest.RData"))
manifest$cpg <- manifest$IlmnID

manifest    <- manifest[mrdocResF$cpg, 
                        c("cpg","CHR", "MAPINFO", "UCSC_RefGene_Name","UCSC_RefGene_Group","Genome_Build")]

load(paste0(annoDir,"cpgInfo_08112016.RData"))  ## Info on the nearest gene, if the CG is intergenic
info$cpg <- rownames(info)

nearestgene <- info[mrdocResF$cpg, c("cpg","SYMBOL")]

mrdocOut <- mrdocResF |> 
  left_join(manifest, by = "cpg") |> 
  left_join(nearestgene, by = "cpg")
dim(mrdocOut)

rownames(mrdocOut) <- mrdocOut$cpg


## save all results #################

save(mrdocOut, file = paste0(sumDir,"mrdoc2_formerSmk_dnam_Annotated.RData"))

## Subset
mrdocSave <- mrdocOut |> 
  select(cpg, CHR, MAPINFO, UCSC_RefGene_Name, UCSC_RefGene_Group, SYMBOL, 
         g1_hat, g1_SE, g1_Z, g1_p, g1_qval, 
         g2_hat, g2_SE, g2_Z, g2_p, g2_qval )

dim(mrdocSave)

save(mrdocSave, file = paste0(sumDir,"mrdoc2_formerSmk_dnam_results.RData"))

