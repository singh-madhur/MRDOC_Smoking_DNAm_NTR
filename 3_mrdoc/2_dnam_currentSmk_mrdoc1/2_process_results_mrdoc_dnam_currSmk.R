## Process the Results for MRDoC DNAm --> Current Smk
## Madhur Singh

# clear the workspace
rm(list = ls()); gc()

# Libraries
library(data.table)
setDTthreads(16)

#Directories
baseDir <- "/base_dir/" 
outDir  <- paste0(baseDir,"out/mrdoc_dnam_currSmk/")
annoDir <- paste0(baseDir,"data/annotation/")
phenoDir<- paste0(baseDir, "data/phenos/")
plotDir <- paste0(outDir,"plots/")
sumDir  <- paste0(outDir,"summary/")

## Load results #########

batches <- 1:15

for (i in batches) {
  
  outFile   <- paste0(outDir,"mrdoc_DNAm_to_currentSmk_BATCH",i,".dat")
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


## Examine ############
library(dplyr)
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

## R2 in GEE
table(mrdocRes$R2percent > 0.5)
# FALSE  TRUE 
#  1850 11090  

table(mrdocRes$Fstat > 10)
# FALSE  TRUE 
#  1816 11124

## No missing S.E. among CpGs with mQTL Fstatistic >10
table( rowSums( is.na(mrdocRes[ mrdocRes$Fstat > 10,]) ) == 0 )


## Filter weak IVs ##########
mrdocResF <- mrdocRes |> 
  filter(Fstat > 10,
         rowSums(is.na(mrdocRes)) == 0)

## Compare the EWAS association size for CpGs retained for MR vs. those dropped ------------------------

## Read in the EWAS results from Joehanes et al
ewas <- readxl::read_xlsx( paste0(phenoDir,"joehanes_etal_supplemental_tables.xlsx"), sheet = 3, skip = 2  )

range(ewas$FDR)
# 3.250128e-42 4.998151e-02 (CpGs with FDR < 0.05)

## CpGs retained for smoking --> DNAm
mrdocRes_smk <- fread( paste0(baseDir,"out/mrdoc_currSmk_dnam/summary/mrdoc_currentSmk_DNAm_SmkRelatedCGs.dta"), data.table = F)
dim(mrdocRes_smk)

ewas_keep <- ewas |> 
  filter(`Probe ID` %in% mrdocRes_smk$cpg)

## Which of these are being assessed for DNAm --> Smk
length(mrdocResF$cpg)

ewas_keep <- ewas_keep |> 
  mutate(w_mQTL = ifelse(`Probe ID` %in% mrdocResF$cpg, 
                         "mQTL+", "mQTL-"),
         abs_Z = abs(`Z-value`),
         log10q = -1 * log10(FDR))

## plot
violin <- ggplot(ewas_keep, aes(w_mQTL, log10q #, fill = w_mQTL
                                )) +
  geom_jitter(aes(color = w_mQTL), shape = 1, alpha = 0.5, height = 0, width = 0.1) +
  geom_violin(# alpha = 0.95, 
              draw_quantiles = 0.5, 
              fill = NA,
              scale = "count") +
  labs(title = "EWAS Meta-Analysis Association Statistics of Smoking-Associated CpGs\nWith and Without an mQTL Allelic Score with F-statistic >10",
       x = NULL,
       y = "-log10(FDR) from EWAS of Current vs. Never Smoking"
  ) +
  scale_fill_manual(values = c("indianred","steelblue")) +
  scale_color_manual(values = c("indianred","steelblue")) +
  theme_bw(10) +
  theme(plot.title.position = "plot", 
        axis.text.x = element_text(size = 10),
        legend.position = "none")

jpeg(paste0(plotDir,"cpgs_with_without_strong_mQTL_score.jpeg"), 
     width = 6, height = 6, units = "in", res = 300)
violin 

dev.off()


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
                  main=paste("MR-DoC1\nDNA Methylation to Current Smoking with Pleiotropic Path")) 
dev.off()

## g2
jpeg(file = paste0(plotDir,"g1_wRE_qqPlot.jpeg"), 
     width = 6, height = 6, units = "in", res = 300)
GWASTools::qqPlot(mrdocResF$g1_wRE_p, col="indianred4", 
                  main=paste("MR-DoC1\nDNA Methylation to Current Smoking with rE")) 
dev.off()



## Significance ###########

nTests <- nrow(mrdocResF)
bonfp  <- 0.05/nTests
bonfZ  <- qnorm(p = 1-bonfp/2)

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

mrdocResF <- mrdocResF |> 
  arrange(g1_Z)


## Annotation ##########

## Load the Annotation Files
load(paste0(annoDir,"manifest.RData"))

manifest$cpg <- manifest$IlmnID

manifest    <- manifest[mrdocResF$cpg, 
                        c("cpg","CHR", "MAPINFO", "UCSC_RefGene_Name","UCSC_RefGene_Group","Genome_Build")]
dim(manifest)

load(paste0(annoDir,"cpgInfo_08112016.RData"))  ## Info on the nearest gene, if the CG is intergenic
info$cpg <- rownames(info)

nearestgene <- info[mrdocResF$cpg, c("cpg","SYMBOL")]

mrdocOut <- mrdocResF |> 
  left_join(manifest, by = "cpg") |> 
  left_join(nearestgene, by = "cpg")
dim(mrdocOut)

rownames(mrdocOut) <- mrdocOut$cpg


# Save results ####################
dim(mrdocOut)
save(mrdocOut, file = paste0(sumDir,"mrdoc_dnam_currSmk_Annotated.RData"))

## Subset
mrdocSave <- mrdocOut |> 
  select(cpg, CHR, SYMBOL, 
         bonfSig_wPleio, bonfSig_wRE, fdrSig_wPleio, fdrSig_wRE,
         g1_hat, g1_SE, g1_Z, g1_p, qval, 
         g1_wRE_hat, g1_wRE_Z, g1_wRE_SE, g1_wRE_p, qval_wRE )

dim(mrdocSave)

save(mrdocSave, file = paste0(sumDir,"mrdoc_dnam_currSmk_dnam.RData"))

