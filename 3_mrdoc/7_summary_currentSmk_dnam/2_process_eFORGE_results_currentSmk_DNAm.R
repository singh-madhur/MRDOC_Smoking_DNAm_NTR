## Examine the DNAM --> Smk causal estimates for brain-related CpGs from eFORGE
## Madhur Singh

# clear the workspace
rm(list = ls()); gc()

# Libraries
library(data.table)
setDTthreads(16)
library(tidyverse)
library(ggrepel)

#Directories
baseDir     <- "/base_dir/" 
eforgeDir   <- paste0(baseDir,"out/eFORGE/")
gscanDir    <- paste0(baseDir,"data/GSCAN_Sum_Stats_with_NTR")
outDir      <- paste0(baseDir,"out/mrdoc_summary_currentSmk_dnam/")
plotDir     <- paste0(outDir,"plots/")


## Load results #########

nominalH3File <- paste0(eforgeDir,
                        "Effects_of_DNAm_on_Smoking_CpGs_Enriched_for_Fetal_Brain.450k.erc2-H3-all.chart.tsv.gz")

## Sites enriched for fetal brain 
nominalH3 <- fread(nominalH3File, header = T, data.table = F)
dim(nominalH3) # 195   9

nominalH3 <- nominalH3 |> 
  arrange(Qvalue)

brainH3K4me3 <- nominalH3 |> 
  filter(Cell == "E082 Fetal Brain Female",
         Datatype == "H3K4me3") |> 
  select(Probe) |> 
  separate_longer_delim(Probe, delim = ",")

nrow(brainH3K4me3)

write.table(brainH3K4me3, file = paste0(eforgeDir,"brain_H3K4Me3_CpGs.tsv"), 
            quote = F, col.names = F, row.names = F)


## Causal Estimates ######################

dnam2smk_g2 <- read.table(paste0(outDir,"suppl_dnam2smk_g2est.tsv"), header = T)
str(dnam2smk_g2)

dnam2smk_g1 <- read.table(paste0(outDir,"suppl_dnam2smk_g1est.tsv"), header = T)
str(dnam2smk_g1)

table(dnam2smk_g2$IlmnID %in% brainH3K4me3$Probe)
table(dnam2smk_g1$IlmnID %in% brainH3K4me3$Probe)

## Add indicator if brain-sp eFORGE
dnam2smk_g2 <- dnam2smk_g2 |> 
  mutate(brain_specific_eFORGE = ifelse(IlmnID %in% brainH3K4me3$Probe, 
                                        TRUE, FALSE))
dnam2smk_g1 <- dnam2smk_g1 |> 
  mutate(brain_specific_eFORGE = ifelse(IlmnID %in% brainH3K4me3$Probe, 
                                        TRUE, FALSE))

dnam2smk_g2 |> 
  count(brain_specific_eFORGE)
dnam2smk_g1 |> 
  count(brain_specific_eFORGE)

## Save the added col
write.table(dnam2smk_g2, 
            file = paste0(outDir,"suppl_dnam2smk_g2est.tsv"), quote = F, row.names = F, col.names = T)
write.table(dnam2smk_g1, 
            file = paste0(outDir,"suppl_dnam2smk_g1est.tsv"), quote = F, row.names = F, col.names = T)


## Smk --> DNAm
smk2dnam_g1 <- read.table(paste0(outDir,"suppl_smk2dnam_g1est.tsv"), header = T)
str(smk2dnam_g1)

smk2dnam_g2 <- read.table(paste0(outDir,"suppl_smk2dnam_g2est.tsv"), header = T)
str(smk2dnam_g2)


## Save in excel.
# https://ycphs.github.io/openxlsx/articles/Introduction.html
outExcel <- list(
  "Current_Smk_to_DNAm_g1" = smk2dnam_g1, 
  "Current_Smk_to_DNAm_Reverse_g2" = smk2dnam_g2,
  "DNAm_to_Current_Smk_g2" = dnam2smk_g2, 
  "DNAm_to_Current_Smk_Reverse_g1" = dnam2smk_g1
)
openxlsx::write.xlsx(outExcel, 
                     file = paste0(outDir,"suppl_mrdoc_smk_dnam.xlsx"))



## Save the gene ====================

## Save the nearest gene for gene-set enrichment

brainH3K4me3_nearestGene <- dnam2smk_g2 |> 
  filter(brain_specific_eFORGE,
         !is.na(NearestGene)) |> 
  select(NearestGene) |> 
  t()

write.table(brainH3K4me3_nearestGene, file = paste0(outDir,"brain_H3K4Me3_nearestGenes.csv"), 
            sep = ",", col.names = F, row.names = F, quote = F)


## Save the gene for gene-set enrichment
brainH3K4me3_gene <- dnam2smk_g2 |> 
  filter(brain_specific_eFORGE,
         !is.na(UCSC_RefGene_Name)) |> 
  select(UCSC_RefGene_Name) |> 
  t()

write.table(brainH3K4me3_gene, file = paste0(outDir,"brain_H3K4Me3_genes.csv"), 
            sep = ",", col.names = F, row.names = F, quote = F)



## Plot g2 Causal Estimates ########################################

h3k4meEstimates <- dnam2smk_g2 |>
  filter(brain_specific_eFORGE) |> 
  mutate(cpg = paste(IlmnID,NearestGene, sep="\n"),
         cpg = fct_reorder(as_factor(cpg), - mrdoc1_wRE.g2_hat)
  ) |> 
  select(cpg, CHR, MAPINFO, NearestGene,
         mrdoc1_wPleio.g2_hat, mrdoc1_wPleio.g2_SE, mrdoc1_wPleio.g2_Z, mrdoc1_wPleio.g2_p, 
         mrdoc1_wRE.g2_hat, mrdoc1_wRE.g2_SE, mrdoc1_wRE.g2_Z, mrdoc1_wRE.g2_p, 
         mrdoc2.g2_hat, mrdoc2.g2_SE, mrdoc2.g2_Z, mrdoc2.g2_p) |> 
  pivot_longer(cols = -c(cpg, CHR, MAPINFO, NearestGene), 
               names_to = c("model","par"), names_sep = ".g2_") |>
  pivot_wider(id_cols = c("cpg","CHR","MAPINFO","NearestGene","model"), 
              names_from = "par") |> 
  mutate(
    model = fct_relevel(model, "mrdoc2"),
    model = fct_recode(model, 
                       "MR-DoC2" = "mrdoc2",
                       "MR-DoC1 w/ Pleiotropic Path" = "mrdoc1_wPleio",
                       "MR-DoC1 w/ rE" = "mrdoc1_wRE") 
  )

dim(h3k4meEstimates)


## Plot
jpeg(paste0(plotDir,"mrdoc_dnam2smk_brain_specific.jpeg"), 
     width = 5, height = 6, units = "in", res = 300)

h3k4meEstimates |> 
  ggplot(aes(x = hat, y = cpg )) + 
  geom_point(aes(color = model, 
                 fill = model, 
                 shape = model),
             size = 1.25, position = position_dodge(width = 0.5) ) +
  geom_errorbarh(aes(color = model, 
                     xmax = hat + qnorm(p=0.975)*SE,
                     xmin = hat - qnorm(p=0.975)*SE),
                 height = 0, 
                 linewidth = 0.5,
                 position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, 
             color = "midnightblue", 
             linetype = "dashed", 
             linewidth = 0.25) +
  labs(title = "Estimated Effects of DNA Methylation on Current Smoking\nAt 17 CpGs Showing Enrichment for Functional Elements in the Brain",
       x = "Causal Estimate (S.D.)", y = NULL, alpha = NULL,
       color = "Model", fill = "Model", shape = "Model") +
  scale_color_manual(values = c("indianred3","thistle4","mistyrose3")) +
  scale_fill_manual(values = c("indianred3","thistle4","mistyrose3")) +
  scale_shape_manual(values = c(15, 19, 17)) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme_light(6) +
  theme(
    plot.title.position = "plot",
    axis.text.y = element_text(size = 6),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom" 
  ) 


dev.off()





## Overlap with Bidirectional Effects ######################

dnam2smk_g1 |> 
  filter(brain_specific_eFORGE,g1_nominal) |> 
  count(IlmnID,CHR,NearestGene)

dnam2smk_g1 |> 
  filter(brain_specific_eFORGE,g1_robust) |> 
  count(IlmnID,NearestGene)

h3k4meEstimates |> 
  filter(alfa==1) |> 
  count(cpg)


## B cell sites ===============

## Check the input for further enrichment analyses
chromFile <- paste0(eforgeDir,
                    "Sites_with_Consistent_Effects_of_DNAm_on_Smoking_in_All_3_Models.450k.erc2-chromatin15state-all.chart.tsv.gz")
chromRes <- fread(chromFile, header = T, data.table = F)

bCellList <- chromRes |> 
  filter(Cell == "E031 Primary B cells from cord blood",
         Datatype == "Enh") |> 
  select(Probe) |> 
  separate_wider_delim(cols = Probe, delim = ",", names = paste0("Probe",1:18)) |> 
  t() |> 
  data.frame()
colnames(bCellList) <- "cpg"

nrow(bCellList)
# 18

## Annotate CpGs
bCellList <- bCellList |> 
  left_join(select(dnam2smk_g2, IlmnID, CHR, MAPINFO, NearestGene), 
            by = c("cpg"="IlmnID")) 


## Overlap with cell count GWAS signal. Vuckovic et al 2020 
cellCountTab1 <- readxl::read_xlsx(paste0(outDir,"VuckovicEtAl_TableS3.xlsx"))
cellCountTab1 <- cellCountTab1 |> 
  janitor::clean_names() 

glimpse(cellCountTab1)

cellCountTab2 <- readxl::read_xlsx(paste0(outDir,"VuckovicEtAl_TableS4.xlsx"))
cellCountTab2 <- cellCountTab2 |> 
  janitor::clean_names() 

glimpse(cellCountTab2)

cellCountTab1 <- cellCountTab1 |> 
  select(associated_blood_index,associated_blood_index_class,
         gene_symbol_s_for_most_serious_consequence,rs_id_where_available) |> 
  rename(blood_index = associated_blood_index,
         blood_index_class = associated_blood_index_class,
         gene = gene_symbol_s_for_most_serious_consequence, 
         id = rs_id_where_available)
head(cellCountTab1)

cellCountTab2 <- cellCountTab2 |> 
  select(associated_blood_index,gene_symbol_for_most_serious_consequence,id) |> 
  rename(blood_index = associated_blood_index,
         gene = gene_symbol_for_most_serious_consequence)

head(cellCountTab2)

# rbind
cellCountGWAS <- cellCountTab1 |> 
  bind_rows(cellCountTab2)

lymphCountGWAS <- cellCountGWAS |> 
  filter(grepl("LYMPH",blood_index))


## overlap
bCellList |> 
  distinct(NearestGene, .keep_all = T) |> 
  count(NearestGene %in% cellCountGWAS$gene)
#   NearestGene %in% cellCountGWAS$gene  n
# 1                               FALSE 11
# 2                                TRUE  5


bCellList |> 
  distinct(NearestGene, .keep_all = T) |> 
  count(NearestGene %in% lymphCountGWAS$gene)
#   NearestGene %in% lymphCountGWAS$gene  n
# 1                                FALSE 15
# 2                                 TRUE  1


