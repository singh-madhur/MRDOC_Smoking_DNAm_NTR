## Process the Results for MRDoC DNAm --> Former Smk
## Madhur Singh

# clear the workspace
rm(list = ls()); gc()

# Libraries
library(data.table)
setDTthreads(16)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

#Directories
baseDir <- "/base_dir/" 
outDir  <- paste0(baseDir,"out/mrdoc_dnam_formerSmk/")
annoDir <- paste0(baseDir,"data/annotation/")
plotDir <- paste0(outDir,"plots/")
sumDir  <- paste0(outDir,"summary/")
phenoDir <- paste0(baseDir,"data/phenos/")


## Load results #########

batches <- 1:15

for (i in batches) {
  
  outFile   <- paste0(outDir,"mrdoc_DNAm_to_formerSmk_BATCH",i,".dat")
  mrdocOut  <- fread(outFile, header = T, data.table = F)
  print(paste0(i,": ",nrow(mrdocOut)))
  
  if(i == 1) { 
    mrdocRes <- mrdocOut 
    rm(mrdocOut); gc()
  } else {
    mrdocRes <- rbind(mrdocRes, mrdocOut)
    rm(mrdocOut); gc()
  }
}

## Examine ##################
dim(mrdocRes)
rownames(mrdocRes) <- mrdocRes$cpg

mrdocRes |> count(statusCode)
mrdocRes |> count(statusCode_wRE)
## Missing SEs
table(rowSums(is.na(mrdocRes)) == 0)


## Check IV R2 and Fstat #############

## IV R2 File
ivR2 <- fread(paste0(baseDir,"data/PRS_mQTL/Rsq/MRDoC_CpG_IV_with_maxR2.dat"), data.table = F)

mrdocRes <- mrdocRes |> 
  left_join(ivR2, by = c("cpg" = "CpG"))

table( rowSums( is.na(mrdocRes[ mrdocRes$Fstat > 10,]) ) == 0 )


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

table(is.element(smkCpG$ID, mrdocRes$cpg))
100*prop.table(table(is.element(smkCpG$ID, mrdocRes$cpg)))


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


# Inflation Factor using {bacon} ##############################

library(bacon)
library(qvalue)

## Setting the seed for rng will ensure reproducibility
set.seed(2023) 
bc <- bacon(teststatistics = NULL, 
            effectsizes = mrdocResF$g1_hat, standarderrors = mrdocResF$g1_SE) 

## wRE: Estimate genomic inflation using {bacon}
set.seed(2023) 
bc_wRE <- bacon(teststatistics = NULL, 
                effectsizes = mrdocResF$g1_wRE_hat, standarderrors = mrdocResF$g1_wRE_SE, na.exclude = T) 


## Make QQ plots #############

## g1_wPleio
jpeg(file = paste0(plotDir,"g1_wPleio_qqPlot.jpeg"), 
     width = 6, height = 6, units = "in", res = 300)
GWASTools::qqPlot(mrdocResF$g1_p, col="indianred4", 
                  main=paste("MR-DoC1\nDNA Methylation to Former Smoking with Pleiotropic Path")) 
dev.off()

## g1_wRE
jpeg(file = paste0(plotDir,"g1_wRE_qqPlot.jpeg"), 
     width = 6, height = 6, units = "in", res = 300)
GWASTools::qqPlot(mrdocResF$g1_wRE_p, col="indianred4", 
                  main=paste("MR-DoC1\nDNA Methylation to Former Smoking with rE")) 
dev.off()



## Significance ###########

nTests <- nrow(mrdocResF)
bonfp  <- 0.05/nTests
bonfZ  <- qnorm(p = 1-bonfp/2)

table(mrdocResF$g1_p < bonfp)
table(mrdocResF$g1_wRE_p < bonfp)

Sig_wPleio <- mrdocResF |> 
  filter(g1_p < bonfp)
Sig_wRE <- mrdocResF |> 
  filter(g1_wRE_p < bonfp)

#### FDR ############################

mrdocResF$qval         <- qvalue::qvalue(mrdocResF$g1_p)$qvalue
mrdocResF$qval_wRE     <- qvalue::qvalue(mrdocResF$g1_wRE_p)$qvalue

## save an indicator variable
mrdocResF <- mrdocResF |> 
  mutate( 
    bonfSig_wPleio = ifelse( g1_p < bonfp,          TRUE, FALSE),
    bonfSig_wRE    = ifelse( g1_wRE_p < bonfp,      TRUE, FALSE),
    fdrSig_wPleio  = ifelse( qval < 0.05,         TRUE, FALSE),
    fdrSig_wRE     = ifelse( qval_wRE < 0.05,     TRUE, FALSE),
  )


## Annotation ##########

## Load the Annotation Files
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

rownames(mrdocOut) <- mrdocOut$cpg

# Save results ####################

save(mrdocOut, file = paste0(sumDir,"mrdoc_dnam_formerSmk_Annotated.RData"))

## Subset
mrdocSave <- mrdocOut |> 
  select(cpg, CHR, MAPINFO, UCSC_RefGene_Name, UCSC_RefGene_Group, SYMBOL, 
         bonfSig_wPleio, bonfSig_wRE, fdrSig_wPleio, fdrSig_wRE,
         g1_hat, g1_SE, g1_Z, g1_p, qval, 
         g1_wRE_hat, g1_wRE_Z, g1_wRE_SE, g1_wRE_p, qval_wRE )

save(mrdocSave, file = paste0(sumDir,"mrdoc_dnam_formerSmk_dnam.RData"))

