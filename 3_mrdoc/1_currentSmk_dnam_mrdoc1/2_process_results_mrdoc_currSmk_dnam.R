## Process the Results for MRDoC Current Smk --> DNAm
## Madhur Singh

# clear the workspace
rm(list = ls()); gc()

# Libraries
library(data.table)
setDTthreads(16)
library(tidyverse)

#Directories
baseDir <- "/base_dir/" 
outDir  <- paste0(baseDir,"out/mrdoc_currSmk_dnam/")
annoDir <- paste0(baseDir,"data/annotation/")
plotDir <- paste0(outDir,"plots/")
sumDir  <- paste0(outDir,"summary/")
phenoDir <- paste0(baseDir,"data/phenos/")


## Load results #########

batches <- 1:15

for (i in batches) {
  
  outFile   <- paste0(outDir,"mrdoc_currentSmk_to_DNAm_BATCH",i,".dat")
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

dim(mrdocRes)
rownames(mrdocRes) <- mrdocRes$cpg

table(mrdocRes$statusCode, useNA = "ifany")

## Missing SEs
table(rowSums(is.na(mrdocRes)) == 0)
round( prop.table( table(rowSums(is.na(mrdocRes)) == 0) ), 3 )


## Filter XY chr and Probes with SNPs #############

## Annotation
load(paste0(annoDir,"SNP.GONL.CG.Dataframe.RData"))
load(paste0(annoDir,"GFanno.RData"))
ls()


#### Annotate probes containing a SNP in the CG #################################
class(ProbSNP)
length(ProbSNP)

length (which(rownames(mrdocRes) %in% ProbSNP))
length (which(!rownames(mrdocRes) %in% ProbSNP))

## Add SNPprobe variable
mrdocRes <- mrdocRes |> 
  mutate( SNPprobe = ifelse(cpg %in% ProbSNP, "SNP", "NoSNP") ) 

mrdocRes |> 
  count(SNPprobe)


#### Annotate probes on the X or Y chromosome ##################################
## Annotation of all 450k microarray
str(GFanno_XY_withsplit_df)
head(GFanno_XY_withsplit_df) # file with chromosome and position (position=start)
rownames(GFanno_XY_withsplit_df) <- GFanno_XY_withsplit_df$name

## Filter out CpG sites already dropped during QC
GFanno <-GFanno_XY_withsplit_df[ rownames(mrdocRes), ]
dim(GFanno)
table(GFanno$seqnames)

XYcpg <- GFanno |> 
  filter( seqnames %in% c("chrX","chrY") ) |> 
  select(name)
length( XYcpg$name )

## Add XYprobe variable
mrdocRes <- mrdocRes |> 
  mutate( XYprobe = ifelse(cpg %in% XYcpg$name, "XY", "autosomal") ) 

mrdocRes |> 
  count(XYprobe)

#### Remove probes containing a SNP / CpG sites on chr X,Y ##################

mrdocRes |> 
  count( SNPprobe, XYprobe )

mrdocResAutoNoSNP <- mrdocRes |> 
  filter( SNPprobe == "NoSNP" & XYprobe == "autosomal" )

dim(mrdocResAutoNoSNP)


## Filter Missing SE #########################

## Filter
mrdocResF <- mrdocResAutoNoSNP[ which(rowSums(is.na(mrdocResAutoNoSNP)) == 0), ]
dim(mrdocResF)


## Filter the site with too low SE
mrdocResFF <- mrdocResF |> 
  filter( g1_wRE_SE > 0.05,
          g1_SE < 2) 

nGWtests <- nrow(mrdocResFF)



# GW genomic inflation using {bacon} ##############################

# BiocManager::install("qvalue", update=F)
# BiocManager::install("bacon", update=F)
library(bacon)
library(qvalue)

## Bacon generates a normal distribution to compare the test statistic with
## Setting the seed for rng will ensure reproducibility
set.seed(2023) 
bc <- bacon(teststatistics = NULL, effectsizes = mrdocResFF$g1_hat, standarderrors = mrdocResFF$g1_SE) 
# This will do 2 things: Estimate the amount of bias and inflation of test statistics, 
# and adjust test stastics (betas, standard errors, p-values) for this bias and inflation.

## wRE: Estimate genomic inflation using {bacon}
set.seed(2023) 
bc_wRE <- bacon(teststatistics = NULL, effectsizes = mrdocResFF$g1_wRE_hat, standarderrors = mrdocResFF$g1_wRE_SE, na.exclude = T) 

# compute FDR q-value
mrdocResFF$qvalue         <- qvalue(mrdocResFF$g1_p)$qvalue
mrdocResFF$qvalue_wRE     <- qvalue(mrdocResFF$g1_wRE_p)$qvalue


#### Compare with Sun et al 2021 =============

prior <- c("cg14391737","cg09338374","cg02978227","cg16841366","cg12956751","cg13849276",
           "cg00475490","cg22222502","cg14580211","cg15212295","cg13258799")

mrdocResFF |> 
  filter(cpg %in% prior)

## Save ###############################

fwrite( mrdocResFF, file = paste0(sumDir,"mrdoc_currentSmk_DNAm_gw_filt.dat"), 
        sep = "\t", quote = F, row.names = F, col.names = T )



## Make QQ plots #############
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("GWASTools")

## g1_wPleio
jpeg(file = paste0(plotDir,"g1_wPleio_qqPlot.jpeg"), 
     width = 6, height = 6, units = "in", res = 300)
GWASTools::qqPlot(mrdocResFF$g1_p, col="steelblue4", 
                  main=paste("Epigenome-wide MR-DoC1\nCurrent Smoking to DNA Methylation with Pleiotropic Path")) 
dev.off()

## g1_wRE
jpeg(file = paste0(plotDir,"g1_wRE_qqPlot.jpeg"), 
     width = 6, height = 6, units = "in", res = 300)
GWASTools::qqPlot(mrdocResFF$g1_wRE_p, col="steelblue4", 
                  main=paste("Epigenome-wide MR-DoC1\nCurrent Smoking to DNA Methylation with rE")) 
dev.off()



## GW Significance ###########


## Bonferroni 
nGWtests <- nrow(mrdocResFF)
bonfp <- 0.05/nGWtests
bonfp   
bonfZ <- qnorm(p = 1-bonfp/2)
bonfZ   

## Save
gwSig_wPleio <- mrdocResFF |> 
  dplyr::filter( qvalue < 0.05 )
gwSig_wRE <- mrdocResFF |> 
  dplyr::filter( qvalue_wRE < 0.05 )

fwrite( gwSig_wPleio, file = paste0(sumDir,"gw_significant_with_pleio.dat"), 
        sep = "\t", quote = F, row.names = F, col.names = T )
fwrite( gwSig_wRE, file = paste0(sumDir,"gw_significant_with_rE.dat"), 
        sep = "\t", quote = F, row.names = F, col.names = T )




## Filter Smoking-related CpG sites ############################


smkCpG <- fread(paste0(phenoDir,"joehanes_etal_cpg_current_v_never_smoking.csv"), data.table = F,
                header = F, col.names = "ID")

## Number of tests = 16,940

## Filter & Save
mrdocRes_smk <- mrdocResFF |> 
  filter( cpg %in% smkCpG$ID ) 
dim(mrdocRes_smk)  # 16940   118

fwrite(mrdocRes_smk,  file = paste0(sumDir,"mrdoc_currentSmk_DNAm_SmkRelatedCGs.dta"), 
       sep = "\t", quote = F, row.names = F, col.names = T)


# Smk-related genomic inflation using {bacon} ##############################

## Setting the seed for rng will ensure reproducibility
set.seed(2023) 
bc <- bacon(teststatistics = NULL, 
            effectsizes = mrdocRes_smk$g1_hat, standarderrors = mrdocRes_smk$g1_SE) 

## wRE: Estimate genomic inflation using {bacon}
set.seed(2023) 
bc_wRE <- bacon(teststatistics = NULL, 
                effectsizes = mrdocRes_smk$g1_wRE_hat, standarderrors = mrdocRes_smk$g1_wRE_SE, na.exclude = T) 


## Make QQ plots #############

## g1_wPleio
jpeg(file = paste0(plotDir,"g1_smkRel_wPleio_qqPlot.jpeg"), 
     width = 6, height = 6, units = "in", res = 300)
GWASTools::qqPlot(mrdocRes_smk$g1_p, col="steelblue4", 
                  main=paste("MR-DoC1\nCurrent Smoking to DNA Methylation with Pleiotropic Path")) 
dev.off()

## g2
jpeg(file = paste0(plotDir,"g1_smkRel_wRE_qqPlot.jpeg"), 
     width = 6, height = 6, units = "in", res = 300)
GWASTools::qqPlot(mrdocRes_smk$g1_wRE_p, col="steelblue4", 
                  main=paste("MR-DoC1\nCurrent Smoking to DNA Methylation with rE")) 
dev.off()


## Smk-related Significance ###########

## Bonferroni 
nSmktests <- nrow(mrdocRes_smk)
bonfp <- 0.05/nSmktests
bonfZ <- qnorm(p = 1-bonfp/2)


#### FDR ############################

## qvalue in smk-related subset
mrdocRes_smk$smk_qvalue         <- qvalue::qvalue(mrdocRes_smk$g1_p)$qvalue
mrdocRes_smk$smk_qvalue_wRE     <- qvalue::qvalue(mrdocRes_smk$g1_wRE_p)$qvalue

## save an indicator variable
mrdocRes_smk <- mrdocRes_smk |> 
  mutate( 
    bonfSig_wPleio = ifelse( g1_p < bonfp,          TRUE, FALSE),
    bonfSig_wRE    = ifelse( g1_wRE_p < bonfp,      TRUE, FALSE),
    fdrSig_wPleio  = ifelse( smk_qvalue < 0.05,         TRUE, FALSE),
    fdrSig_wRE     = ifelse( smk_qvalue_wRE < 0.05,     TRUE, FALSE),
  )

mrdocRes_smk <- mrdocRes_smk |> 
  arrange(g1_Z)



## Annotation ##########

## Load the Annotation Files
load(paste0(annoDir,"manifest.RData"))
dim(manifest)
names(manifest)

## How many CpGs from Sun etal are on 450k (they used EPIC)
table(prior %in% manifest$IlmnID)

manifest$cpg <- manifest$IlmnID

manifest    <- manifest[mrdocRes_smk$cpg, 
                        c("cpg","CHR", "MAPINFO", "UCSC_RefGene_Name","UCSC_RefGene_Group","Genome_Build")]

load(paste0(annoDir,"cpgInfo_08112016.RData"))  ## Info on the nearest gene, if the CG is intergenic
info$cpg <- rownames(info)

nearestgene <- info[mrdocRes_smk$cpg, c("cpg","SYMBOL")]

## Merge

mrdocOut <- mrdocRes_smk |> 
  left_join(manifest, by = "cpg") |> 
  left_join(nearestgene, by = "cpg")
dim(mrdocOut)

rownames(mrdocOut) <- mrdocOut$cpg

# Save results ####################

save(mrdocOut, file = paste0(sumDir,"mrdoc_currSmk_dnam_Annotated.RData"))

## Subset
mrdocSave <- mrdocOut |> 
  select(cpg, CHR, MAPINFO, SYMBOL, 
         bonfSig_wPleio, bonfSig_wRE, fdrSig_wPleio, fdrSig_wRE,
         g1_hat, g1_SE, g1_Z, g1_p, smk_qvalue, 
         g1_wRE_hat, g1_wRE_Z, g1_wRE_SE, g1_wRE_p, smk_qvalue_wRE )

save(mrdocSave, file = paste0(sumDir,"mrdoc_currSmk_dnam_results.RData"))


