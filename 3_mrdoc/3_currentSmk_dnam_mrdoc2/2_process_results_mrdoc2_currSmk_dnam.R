## Process the Results for MRDoC-2 Current Smk <--> DNAm
## Madhur Singh, 28-Feb-2023. Updated 7-Apr-2023

# clear the workspace
rm(list = ls()); gc()

# Libraries
library(data.table)
setDTthreads(16)
library(qvalue) ## To get the FDR

#Directories
baseDir <- "/base_dir/" 
outDir  <- paste0(baseDir,"out/mrdoc2_currSmk_dnam/")
annoDir <- paste0(baseDir,"data/annotation/")
plotDir <- paste0(outDir,"plots/")
sumDir  <- paste0(outDir,"summary/")


## Load results #########

batches <- 1:15

for (i in batches) {
  
  outFile   <- paste0(outDir,"mrdoc2_currentSmk_DNAm_BATCH",i,".dat")
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

## Examine ############
library(dplyr)
library(ggplot2)
library(forcats)

dim(mrdocRes)
rownames(mrdocRes) <- mrdocRes$cpg

table(mrdocRes$statusCode, useNA = "ifany")

## Missing SEs
table(rowSums(is.na(mrdocRes)) == 0)

## Filter low R2 #############

## IV R2 File
ivR2 <- fread(paste0(baseDir,"data/PRS_mQTL/Rsq/MRDoC_CpG_IV_with_maxR2.dat"), data.table = F)

mrdocRes <- mrdocRes |> 
  left_join(ivR2, by = c("cpg" = "CpG"))

## Fstat
table(mrdocRes$Fstat > 10)

## R2 in GEE
table(mrdocRes$R2percent > 0.5)


## Filter IV Fstat > 10 ###########################
mrdocResF <- mrdocRes |> 
  filter(Fstat > 10,
         rowSums(is.na(mrdocRes)) == 0)
dim(mrdocResF)



##  Inflation Factor using {bacon} ##############################

library(bacon)

#### Smk --> DNAm ###################

## Bacon generates a normal distribution to compare the test statistic with
## Setting the seed for rng will ensure reproducibility
set.seed(2023) 
bc_g1 <- bacon(teststatistics = NULL, 
               effectsizes = mrdocResF$g1_hat, standarderrors = mrdocResF$g1_SE) 


#### DNAm --> Smk ###################

set.seed(2023) 
bc_g2 <- bacon(teststatistics = NULL, 
               effectsizes = mrdocResF$g2_hat, standarderrors = mrdocResF$g2_SE, na.exclude = T) 

bc_g2


## Make QQ plots #############

## g1
jpeg(file = paste0(plotDir,"g1_MRDoC2_qqPlot.jpeg"), 
     width = 6, height = 6, units = "in", res = 300)
GWASTools::qqPlot(mrdocResF$g1_p, col="steelblue4", 
                  main=paste("MR-DoC2: Current Smoking to DNA Methylation")) 
dev.off()

## g2
jpeg(file = paste0(plotDir,"g2_MRDoC2_qqPlot.jpeg"), 
     width = 6, height = 6, units = "in", res = 300)
GWASTools::qqPlot(mrdocResF$g2_p, col="indianred4", 
                  main=paste("MR-DoC2: DNA Methylation to Current Smoking")) 
dev.off()


## Significance ############

nTests <- nrow(mrdocResF)
bonfp <- 0.05/nTests
bonfZ <- qnorm(p = 1-bonfp/2)

table(mrdocResF$g1_p < bonfp)
table(mrdocResF$g2_p < bonfp)

## FDR
mrdocResF$g1_qval <- qvalue(mrdocResF$g1_p)$qvalue
mrdocResF$g2_qval <- qvalue(mrdocResF$g2_p)$qvalue

table(mrdocResF$g1_qval < 0.05)
table(mrdocResF$g2_qval < 0.05)

table(mrdocResF$g1_qval < 0.05, mrdocResF$g2_qval < 0.05)


mrdocResF <- mrdocResF |> 
  arrange(desc(abs(g2_Z)))


## Annotation ##########

## Load the Annotation Files
load(paste0(annoDir,"manifest.RData"))
dim(manifest)
names(manifest)
manifest$cpg <- manifest$IlmnID


manifest    <- manifest[mrdocResF$cpg, 
                        c("cpg","CHR", "MAPINFO", "UCSC_RefGene_Name","UCSC_RefGene_Group","Genome_Build")]

load(paste0(annoDir,"cpgInfo_08112016.RData"))  ## Info on the nearest gene, if the CG is intergenic
info$cpg <- rownames(info)

nearestgene <- info[mrdocResF$cpg, c("cpg","SYMBOL")]

mrdocOut <- mrdocResF |> 
  left_join(manifest, by = "cpg") |> 
  left_join(nearestgene, by = "cpg")

rm(manifest,nearestgene); gc()

rownames(mrdocOut) <- mrdocOut$cpg


## Power to detect DNAm --> Smk depends on the strength of the mQTL score
cor.test(mrdocOut$by2_hat, mrdocOut$g2_SE)
cor.test(mrdocOut$by2_Z, mrdocOut$g2_SE)


## save all results #################

dim(mrdocOut)
save(mrdocOut, file = paste0(sumDir,"mrdoc2_currSmk_dnam_Annotated.RData"))

## Subset
mrdocSave <- mrdocOut |> 
  select(cpg, CHR, MAPINFO, SYMBOL, 
         g1_hat, g1_SE, g1_Z, g1_p, g1_qval, 
         g2_hat, g2_SE, g2_Z, g2_p, g2_qval )

save(mrdocSave, file = paste0(sumDir,"mrdoc2_currSmk_dnam_results.RData"))
