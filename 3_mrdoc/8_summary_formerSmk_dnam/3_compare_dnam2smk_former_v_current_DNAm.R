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

str(formerSmkEWAS)

smkCpG <- formerSmkEWAS |> 
  select(`Probe ID`) |> 
  rename(ID = `Probe ID`) |> 
  data.frame()

table(is.element(currSmkRes$cpg, smkCpG$ID))

## Filter
currSmkResF <- currSmkRes |> 
  filter(cpg %in% smkCpG$ID)


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

glimpse(mrdoc_smk_dnam)


## Save ##########################
save(mrdoc_smk_dnam, file = paste0(outDir,"mrdoc_current_and_former_smk_DNAm.RData"))



## Estimates at Sites with Consistent Effects of DNAm on Current Smk ###########################


#### FDR < 0.05 in all three models

nFDR <- mrdoc_smk_dnam |>
  filter( (current.mrdoc2.g2_fdrSig &
             current.mrdoc1_wPleio.g2_fdrSig &
             current.mrdoc1_wRE.g2_fdrSig)  ) |> 
  nrow()

mrdoc_DNAm_currentSmk_FDR <- mrdoc_smk_dnam |> 
  filter( (current.mrdoc2.g2_fdrSig &
             current.mrdoc1_wPleio.g2_fdrSig &
             current.mrdoc1_wRE.g2_fdrSig)  ) |> 
  mutate(cpg = paste(cpg,SYMBOL, sep = "\n"),
         cpg = fct_reorder(cpg, current.mrdoc1_wPleio.g2_hat)) |> 
  select(cpg, CHR, MAPINFO, SYMBOL,
         former.mrdoc2.g2_hat, former.mrdoc2.g2_SE, former.mrdoc2.g2_Z,
         former.mrdoc1_wPleio.g2_hat, former.mrdoc1_wPleio.g2_SE, former.mrdoc1_wPleio.g2_Z,
         former.mrdoc1_wRE.g2_hat, former.mrdoc1_wRE.g2_SE, former.mrdoc1_wRE.g2_Z,
         current.mrdoc2.g2_hat, current.mrdoc2.g2_SE, current.mrdoc2.g2_Z,
         current.mrdoc1_wPleio.g2_hat, current.mrdoc1_wPleio.g2_SE, current.mrdoc1_wPleio.g2_Z,
         current.mrdoc1_wRE.g2_hat, current.mrdoc1_wRE.g2_SE, current.mrdoc1_wRE.g2_Z) |> 
  pivot_longer(cols = -c(cpg, CHR, MAPINFO, SYMBOL), 
               names_to = c("model","par"), names_sep = ".g2_") |> 
  pivot_wider(id_cols = c("cpg","CHR","MAPINFO","SYMBOL","model"), 
              names_from = "par", names_prefix = "g2_") |> 
  separate_wider_delim(model, delim = ".", names = c("Outcome","Model")) |> 
  mutate(Model = fct_relevel(Model, "mrdoc2"),
         Model = fct_recode(Model, 
                            "MR-DoC2" = "mrdoc2",
                            "MR-DoC1 w/ Pleiotropic Path" = "mrdoc1_wPleio",
                            "MR-DoC1 w/ rE" = "mrdoc1_wRE"),
         Outcome = case_when(Outcome == "former" ~ "Former Smoking",
                              Outcome == "current" ~ "Current Smoking")) |> 
  group_by(cpg,Outcome) |> 
  arrange(desc(abs(g2_hat))) |> 
  mutate(top = row_number(),
         plotLabel = ifelse(top == 1 & Outcome == "Current Smoking", SYMBOL, NA) ) |> 
  ungroup() |> 
  arrange(cpg)

head(mrdoc_DNAm_currentSmk_FDR)

### Plot

mrdoc_DNAm_currentSmk_FDR_plot <- mrdoc_DNAm_currentSmk_FDR |> 
  ggplot(aes(g2_hat, cpg, label=plotLabel)) + 
  geom_errorbarh(aes(color = Model, alpha = Outcome,
                     linetype = Outcome,
                     xmax = g2_hat + qnorm(p=0.975)*g2_SE,
                     xmin = g2_hat - qnorm(p=0.975)*g2_SE),
                 height = 0, linewidth = 0.5,
                 position = position_dodge(width = 0.8)) +
  geom_point(aes(color = Model, shape = Outcome, alpha = Outcome),
             size = 3, position = position_dodge(width = 0.8) ) +
  geom_vline(xintercept = 0, color = "grey50", linetype = "dashed") +
  facet_wrap(~ cpg, scales = "free_y", ncol = 1, strip.position = "left") + 
  labs( x = "Causal Estimate (95% C.I.)", y = NULL, # color = NULL, fill = NULL,
        title = "Putative Effects of DNA Methylation on Current Smoking\nCompared to the Estimated Effects on Former Smoking") +
  scale_color_manual(values = c("indianred3","thistle4","mistyrose3")) +
  scale_fill_manual(values = c("indianred3","thistle4","mistyrose3")) +
  scale_alpha_manual(values = c(1,1)) + 
  scale_shape_manual(values = c(1,16)) +
  scale_linetype_manual(values = c("twodash","solid")) +
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


jpeg(paste0(plotDir,"mrdoc_dnam_current_v_former_Smk_estimates.jpeg"), 
     width = 7, height = 5, units = "in", res = 300)

mrdoc_DNAm_currentSmk_FDR_plot

dev.off()


####
