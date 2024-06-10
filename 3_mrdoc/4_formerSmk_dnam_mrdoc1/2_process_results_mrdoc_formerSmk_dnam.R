## Process the Results for MRDoC Former Smk --> DNAm
## Madhur Singh

## Filter to keep only the CpGs significantly associated with Former Smk in prior EWAS

# clear the workspace
rm(list = ls()); gc()

# Libraries
library(data.table)
setDTthreads(24)

#Directories
baseDir <- "/base_dir/" 
outDir  <- paste0(baseDir,"out/mrdoc_formerSmk_dnam/")
annoDir <- paste0(baseDir,"data/annotation/")
plotDir <- paste0(outDir,"plots/")
sumDir  <- paste0(outDir,"summary/")
phenoDir <- paste0(baseDir,"data/phenos/")

## Load results #########

batches <- 1:15

for (i in batches) {
  
  outFile   <- paste0(outDir,"mrdoc_formerSmk_to_DNAm_BATCH",i,".dat")
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
library(tidyr)
library(readxl)

dim(mrdocRes)
rownames(mrdocRes) <- mrdocRes$cpg

table(mrdocRes$statusCode, useNA = "ifany")

## Missing SEs
table(rowSums(is.na(mrdocRes)) == 0)



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


#### Annotate probes on the X or Y chromosome ##################################
## Annotation of all 450k microarray
str(GFanno_XY_withsplit_df)
head(GFanno_XY_withsplit_df) # file with chromosome and position (position=start)
rownames(GFanno_XY_withsplit_df) <- GFanno_XY_withsplit_df$name

## Filter out CpG sites already dropped during QC
GFanno <-GFanno_XY_withsplit_df[ rownames(mrdocRes), ]

XYcpg <- GFanno |> 
  filter( seqnames %in% c("chrX","chrY") ) |> 
  select(name)

## Add XYprobe variable
mrdocRes <- mrdocRes |> 
  mutate( XYprobe = ifelse(cpg %in% XYcpg$name, "XY", "autosomal") ) 



#### Remove probes containing a SNP / CpG sites on chr X,Y ##################

mrdocResAutoNoSNP <- mrdocRes |> 
  filter( SNPprobe == "NoSNP" & XYprobe == "autosomal" )



## Filter Missing SE #########################
mrdocResF <- mrdocResAutoNoSNP[ which(rowSums(is.na(mrdocResAutoNoSNP)) == 0), ]


## Filter the site with too low SE
mrdocResFF <- mrdocResF |> 
  filter( g1_wRE_SE > 0.05) 


# Estimate genomic inflation using {bacon} ##############################

library(bacon)
library(qvalue)

## Setting the seed for rng will ensure reproducibility
set.seed(2023) 
bc <- bacon(teststatistics = NULL, effectsizes = mrdocResFF$g1_hat, standarderrors = mrdocResFF$g1_SE) 

## wRE: Estimate genomic inflation using {bacon}
set.seed(2023) 
bc_wRE <- bacon(teststatistics = NULL, effectsizes = mrdocResFF$g1_wRE_hat, standarderrors = mrdocResFF$g1_wRE_SE, na.exclude = T) 

# compute FDR q-value
mrdocResFF$qvalue         <- qvalue(mrdocResFF$g1_p)$qvalue
mrdocResFF$qvalue_wRE     <- qvalue(mrdocResFF$g1_wRE_p)$qvalue


## Save ###############################

fwrite( mrdocResFF, file = paste0(sumDir,"mrdoc_FormerSmk_DNAm_gw_filt.dat"), 
        sep = "\t", quote = F, row.names = F, col.names = T )


## GW Significance ###########


## Bonferroni 
nGWtests <- nrow(mrdocResFF)
bonfp <- 0.05/nGWtests
bonfZ <- qnorm(p = 1-bonfp/2)

## Save
gwSig_wPleio <- mrdocResFF |> 
  filter( qvalue < 0.05 )
gwSig_wRE <- mrdocResFF |> 
  filter( qvalue_wRE < 0.05 )

fwrite( gwSig_wPleio, file = paste0(sumDir,"gw_significant_with_pleio.dat"), 
        sep = "\t", quote = F, row.names = F, col.names = T )
fwrite( gwSig_wRE, file = paste0(sumDir,"gw_significant_with_rE.dat"), 
        sep = "\t", quote = F, row.names = F, col.names = T )




## Filter Former Smoking related CpG sites ############################

## EWAS results from Joehanes et al
## CpGs associated w/ Former v. Never Smk (at FDR < 0.05)
formerSmkEWAS <- read_excel(path = paste0(phenoDir,"joehanes_etal_supplemental_tables.xlsx"), 
                            sheet = 6, skip = 2)

smkCpG <- formerSmkEWAS |> 
  select(`Probe ID`) |> 
  rename(ID = `Probe ID`) |> 
  data.frame()

str(smkCpG)

table(is.element(smkCpG$ID, mrdocResAutoNoSNP$cpg))
table(is.element(smkCpG$ID, mrdocResFF$cpg))
## None of these CpGs had any missing SEs
## Number of tests = 2,330


## Filter & Save
mrdocRes_smk <- mrdocResFF |> 
  filter( cpg %in% smkCpG$ID )

fwrite(mrdocRes_smk,  file = paste0(sumDir,"mrdoc_FormerSmk_DNAm_SmkRelatedCGs.dta"), 
       sep = "\t", quote = F, row.names = F, col.names = T)


##### Estimate genomic inflation using {bacon} ##############################

## Setting the seed for rng will ensure reproducibility
set.seed(2023) 
bc <- bacon(teststatistics = NULL, effectsizes = mrdocRes_smk$g1_hat, standarderrors = mrdocRes_smk$g1_SE) 

## wRE: Estimate genomic inflation using {bacon}
set.seed(2023) 
bc_wRE <- bacon(teststatistics = NULL, effectsizes = mrdocRes_smk$g1_wRE_hat, standarderrors = mrdocRes_smk$g1_wRE_SE, na.exclude = T) 



##### Former Smk-related Significance ###########

## Bonferroni 
nSmktests <- nrow(mrdocRes_smk)
bonfp <- 0.05/nSmktests
bonfZ <- qnorm(p = 1-bonfp/2)

table(mrdocRes_smk$g1_p < bonfp)
table(mrdocRes_smk$g1_wRE_p < bonfp)

## qvalue in smk-related subset
mrdocRes_smk$smk_qvalue         <- qvalue(mrdocRes_smk$g1_p)$qvalue
mrdocRes_smk$smk_qvalue_wRE     <- qvalue(mrdocRes_smk$g1_wRE_p)$qvalue

mrdocRes_smk <- mrdocRes_smk |> 
  arrange(g1_Z)


## Annotation ##########

## Load the Annotation Files
load(paste0(annoDir,"manifest.RData"))
manifest$cpg <- manifest$IlmnID

manifest    <- manifest[mrdocRes_smk$cpg, 
                        c("cpg","CHR", "MAPINFO", "UCSC_RefGene_Name","UCSC_RefGene_Group","Genome_Build")]

load(paste0(annoDir,"cpgInfo_08112016.RData"))  ## Info on the nearest gene, if the CG is intergenic
info$cpg <- rownames(info)

nearestgene <- info[mrdocRes_smk$cpg, c("cpg","SYMBOL")]

mrdocOut <- mrdocRes_smk |> 
  left_join(manifest, by = "cpg") |> 
  left_join(nearestgene, by = "cpg")
dim(mrdocOut)
# 2330   126

rownames(mrdocOut) <- mrdocOut$cpg


# Save results ####################
dim(mrdocOut)
save(mrdocOut, file = paste0(sumDir,"mrdoc_formerSmk_dnam_Annotated.RData"))

## selected columns
mrdocOut_formerSmk <- mrdocOut |> 
  select(cpg, CHR, MAPINFO, UCSC_RefGene_Name, UCSC_RefGene_Group, SYMBOL, 
         g1_hat, g1_SE, g1_Z, g1_p, smk_qvalue, 
         g1_wRE_hat, g1_wRE_SE, g1_wRE_Z, g1_wRE_p, smk_qvalue_wRE) |> 
  rename(g1_q = smk_qvalue,
         g1_wRE_q = smk_qvalue_wRE)


save(mrdocOut_formerSmk, file = paste0(sumDir,"mrdoc_formerSmk_dnam_results.RData"))

