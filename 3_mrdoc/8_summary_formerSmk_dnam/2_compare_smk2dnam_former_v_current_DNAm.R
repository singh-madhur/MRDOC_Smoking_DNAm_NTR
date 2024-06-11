## Compare Former Smk --> DNAm and Current Smk <--> DNAm
## Madhur Singh

# clear the workspace
rm(list = ls()); gc()

# Libraries
library(data.table)
setDTthreads(16)
library(tidyverse)
library(ggrepel)
library(ComplexUpset)

#Directories
baseDir     <- "/base_dir/" 
currSmkDir  <- paste0(baseDir,"out/mrdoc_summary_currentSmk_dnam/")
formSmkDir  <- paste0(baseDir,"out/mrdoc_summary_formerSmk_dnam/")
outDir      <- formSmkDir
plotDir     <- paste0(outDir,"plots/")
phenoDir    <- paste0(baseDir,"data/phenos/")


## Load results #########

## https://stackoverflow.com/a/25455968/18446763
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

## Current Smk
currSmkRes <- loadRData( paste0(currSmkDir,"mrdoc_currentSmk_DNAm_results.RData") )

## Former Smk
formSmkRes <- loadRData( paste0(formSmkDir,"mrdoc_formerSmk_DNAm_results.RData") )



## Rename Columns ##################


#### Current Smk ####################

colnames(currSmkRes)
colnames(currSmkRes)[-c(1:4)] <- paste0("current.",colnames(currSmkRes))[-c(1:4)]
str(currSmkRes)

currSmkRes <- currSmkRes |> 
  mutate(CHR = as.character(CHR),
         CHR = as.numeric(CHR))


#### Keep Former-Smoking-related CpGs ############################

## EWAS results from Joehanes et al
## CpGs associated w/ Former v. Never Smk (at FDR < 0.05)
library(readxl)
formerSmkEWAS <- read_excel(path = paste0(phenoDir,"joehanes_etal_supplemental_tables.xlsx"), 
                            sheet = 6, skip = 2)

dim(formerSmkEWAS)

smkCpG <- formerSmkEWAS |> 
  select(`Probe ID`) |> 
  rename(ID = `Probe ID`) |> 
  data.frame()

str(smkCpG)

## Filter
currSmkResF <- currSmkRes |> 
  filter(cpg %in% smkCpG$ID)

dim(currSmkResF)


#### Former Smk ####################

colnames(formSmkRes)
colnames(formSmkRes)[-c(1:6)] <- paste0("former.",colnames(formSmkRes)[-c(1:6)])

str(formSmkRes)

formSmkRes <- formSmkRes |> 
  mutate(CHR = as.character(CHR),
         CHR = as.numeric(CHR))


## Merge ###################################

table(is.element(currSmkResF$cpg, formSmkRes$cpg ))

mrdoc_smk_dnam <- formSmkRes |> 
  full_join(currSmkResF,
            by = c("cpg", "CHR", "MAPINFO", "SYMBOL")) 

str(mrdoc_smk_dnam)


## Save ##########################
save(mrdoc_smk_dnam, file = paste0(outDir,"mrdoc_current_and_former_smk_DNAm.RData"))



## Estimates at Sites with  Consistent Effects of Former Smk ###########################


#### FDR < 0.05 in all three models

nFDR <- mrdoc_smk_dnam |>
  filter( (former.mrdoc2.g1_fdrSig &
             former.mrdoc1_wPleio.g1_fdrSig &
             former.mrdoc1_wRE.g1_fdrSig)  ) |> 
  nrow()


mrdoc_formerSmk_DNAm_FDR <- mrdoc_smk_dnam |> 
  filter( (former.mrdoc2.g1_fdrSig &
             former.mrdoc1_wPleio.g1_fdrSig &
             former.mrdoc1_wRE.g1_fdrSig)  ) |> 
  mutate(cpg = paste(cpg,SYMBOL, sep = "\n"),
         cpg = fct_reorder(cpg, current.mrdoc1_wPleio.g1_hat)) |> 
  select(cpg, CHR, MAPINFO, SYMBOL,
         former.mrdoc2.g1_hat, former.mrdoc2.g1_SE, former.mrdoc2.g1_Z,
         former.mrdoc1_wPleio.g1_hat, former.mrdoc1_wPleio.g1_SE, former.mrdoc1_wPleio.g1_Z,
         former.mrdoc1_wRE.g1_hat, former.mrdoc1_wRE.g1_SE, former.mrdoc1_wRE.g1_Z,
         current.mrdoc2.g1_hat, current.mrdoc2.g1_SE, current.mrdoc2.g1_Z,
         current.mrdoc1_wPleio.g1_hat, current.mrdoc1_wPleio.g1_SE, current.mrdoc1_wPleio.g1_Z,
         current.mrdoc1_wRE.g1_hat, current.mrdoc1_wRE.g1_SE, current.mrdoc1_wRE.g1_Z) |> 
  pivot_longer(cols = -c(cpg, CHR, MAPINFO, SYMBOL), 
               names_to = c("model","par"), names_sep = ".g1_") |> 
  pivot_wider(id_cols = c("cpg","CHR","MAPINFO","SYMBOL","model"), 
              names_from = "par", names_prefix = "g1_") |> 
  separate_wider_delim(model, delim = ".", names = c("Exposure","Model")) |> 
  mutate(Model = fct_relevel(Model, "mrdoc2"),
         Model = fct_recode(Model, 
                            "MR-DoC2" = "mrdoc2",
                            "MR-DoC1 w/ Pleiotropic Path" = "mrdoc1_wPleio",
                            "MR-DoC1 w/ rE" = "mrdoc1_wRE"),
         Exposure = case_when(Exposure == "former" ~ "Former Smoking",
                              Exposure == "current" ~ "Current Smoking")) |> 
  group_by(cpg,Exposure) |> 
  arrange(desc(abs(g1_hat))) |> 
  mutate(top = row_number(),
         plotLabel = ifelse(top == 1 & Exposure == "Current Smoking", SYMBOL, NA) ) |> 
  ungroup() |> 
  arrange(cpg)

head(mrdoc_formerSmk_DNAm_FDR)

mrdoc_formerSmk_DNAm_FDR_plot <- mrdoc_formerSmk_DNAm_FDR |> 
  ggplot(aes(g1_hat, cpg, label=plotLabel)) + 
  geom_errorbarh(aes(color = Model, alpha = Exposure,
                     linetype = Exposure, linewidth = Exposure,
                     xmax = g1_hat + qnorm(p=0.975)*g1_SE,
                     xmin = g1_hat - qnorm(p=0.975)*g1_SE),
                 height = 0, 
                 position = position_dodge(width = 0.8)) +
  geom_point(aes(color = Model, shape = Exposure, alpha = Exposure),
             size = 2.5, position = position_dodge(width = 0.8) ) +
  geom_vline(xintercept = 0, color = "grey50", linetype = "dashed") +
  facet_wrap(~ cpg, scales = "free_y", ncol = 1, strip.position = "left") + 
  labs( x = "Causal Estimate (95% C.I.)", y = NULL, # color = NULL, fill = NULL,
        title = "CpGs with Putative Effects of Former Smoking on DNA Methylation",
        subtitle = "Sites with FDR < 0.05 in All Three Models") +  
  scale_color_manual(values = c("grey25","steelblue3","skyblue1")) +
  scale_fill_manual(values = c("grey25","steelblue3","skyblue1")) +
  scale_alpha_manual(values = c(1,1)) + 
  scale_shape_manual(values = c(0,15)) +
  scale_linetype_manual(values = c("twodash","solid")) +
  scale_linewidth_manual(values = c(0.5,1)) +
  theme_bw(10) +
  theme(plot.title = element_text(size = 10),
        plot.subtitle = element_text(size = 8),
        plot.title.position = "plot",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.position = "right" ) 


jpeg(paste0(plotDir,"mrdoc_former_v_current_Smk_dnam_estimates.jpeg"), 
     width = 6, height = 6, units = "in", res = 300)

mrdoc_formerSmk_DNAm_FDR_plot

dev.off()

####
