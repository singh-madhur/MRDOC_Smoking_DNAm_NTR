## Combine the Results for MRDoC Current Smk <--> DNAm
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
smk2dnamDir <- paste0(baseDir,"out/mrdoc_currSmk_dnam/summary/")
dnam2smkDir <- paste0(baseDir,"out/mrdoc_dnam_currSmk/summary/")
mrdoc2Dir   <- paste0(baseDir,"out/mrdoc2_currSmk_dnam/summary/")
outDir      <- paste0(baseDir,"out/mrdoc_summary_currentSmk_dnam/")
plotDir     <- paste0(outDir,"plots/")

## Load results #########

## https://stackoverflow.com/a/25455968/18446763
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

## MRDoC1 Smk --> DNAm
smk2dnamRes <- loadRData( paste0(smk2dnamDir,"mrdoc_currSmk_dnam_results.RData") )

## MRDoC1 DNAm --> Smk
dnam2smkRes <- loadRData( paste0(dnam2smkDir,"mrdoc_dnam_currSmk_dnam.RData") )

## MRDoC2
mrdoc2Res   <- loadRData( paste0(mrdoc2Dir,"mrdoc2_currSmk_dnam_results.RData") )


## Rename Columns ##################

#### SMK --> DNAm ####################

smk2dnamRes <- smk2dnamRes |> 
  rename( mrdoc1_wPleio.g1_bonfSig = "bonfSig_wPleio",
          mrdoc1_wRE.g1_bonfSig    = "bonfSig_wRE",
          mrdoc1_wPleio.g1_fdrSig  = "fdrSig_wPleio",
          mrdoc1_wRE.g1_fdrSig     = "fdrSig_wRE",
          mrdoc1_wPleio.g1_hat     = "g1_hat",
          mrdoc1_wPleio.g1_SE      = "g1_SE",
          mrdoc1_wPleio.g1_Z       = "g1_Z",
          mrdoc1_wPleio.g1_p       = "g1_p",
          mrdoc1_wPleio.g1_qval    = "smk_qvalue",
          mrdoc1_wRE.g1_hat        = "g1_wRE_hat",
          mrdoc1_wRE.g1_SE         = "g1_wRE_SE",
          mrdoc1_wRE.g1_Z          = "g1_wRE_Z",
          mrdoc1_wRE.g1_p          = "g1_wRE_p",
          mrdoc1_wRE.g1_qval       = "smk_qvalue_wRE" )

smk2dnamRes <- smk2dnamRes |> 
  mutate(CHR = as.character(CHR),
         CHR = as.numeric(CHR))


#### DNAm --> Smk ####################

dnam2smkRes <- dnam2smkRes |> 
  rename( mrdoc1_wPleio.g2_bonfSig = "bonfSig_wPleio",
          mrdoc1_wRE.g2_bonfSig    = "bonfSig_wRE",
          mrdoc1_wPleio.g2_fdrSig  = "fdrSig_wPleio",
          mrdoc1_wRE.g2_fdrSig     = "fdrSig_wRE",
          mrdoc1_wPleio.g2_hat     = "g1_hat",
          mrdoc1_wPleio.g2_SE      = "g1_SE",
          mrdoc1_wPleio.g2_Z       = "g1_Z",
          mrdoc1_wPleio.g2_p       = "g1_p",
          mrdoc1_wPleio.g2_qval    = "qval",
          mrdoc1_wRE.g2_hat        = "g1_wRE_hat",
          mrdoc1_wRE.g2_SE         = "g1_wRE_SE",
          mrdoc1_wRE.g2_Z          = "g1_wRE_Z",
          mrdoc1_wRE.g2_p          = "g1_wRE_p",
          mrdoc1_wRE.g2_qval       = "qval_wRE" )

dnam2smkRes <- dnam2smkRes |> 
  mutate(CHR = as.character(CHR),
         CHR = as.numeric(CHR))



#### Bidirectional ####################

nTests <- nrow(mrdoc2Res)
bonfp <- 0.05/nTests

mrdoc2Res <- mrdoc2Res |> 
  mutate( 
    mrdoc2.g1_bonfSig = ifelse( g1_p < bonfp,   TRUE, FALSE),
    mrdoc2.g2_bonfSig = ifelse( g2_p < bonfp,   TRUE, FALSE),
    mrdoc2.g1_fdrSig  = ifelse( g1_qval < 0.05, TRUE, FALSE),
    mrdoc2.g2_fdrSig  = ifelse( g2_qval < 0.05, TRUE, FALSE)
  ) |> 
  rename(
    mrdoc2.g1_hat     = "g1_hat",
    mrdoc2.g1_SE      = "g1_SE",
    mrdoc2.g1_Z       = "g1_Z",
    mrdoc2.g1_p       = "g1_p",
    mrdoc2.g1_qval    = "g1_qval",
    mrdoc2.g2_hat     = "g2_hat",
    mrdoc2.g2_SE      = "g2_SE",
    mrdoc2.g2_Z       = "g2_Z",
    mrdoc2.g2_p       = "g2_p",
    mrdoc2.g2_qval    = "g2_qval"
  )

mrdoc2Res <- mrdoc2Res |> 
  mutate(CHR = as.character(CHR),
         CHR = as.numeric(CHR))




## Merge ###################################

nrow(smk2dnamRes) # Smk-related CpGs post-QC
nrow(dnam2smkRes) # Smk-related CpGs post-QC AND mQTLs with Fstat > 10
nrow(mrdoc2Res)   # Smk-related CpGs post-QC AND mQTLs with Fstat > 10 AND No missing SE

mrdoc_currSmk_Res <- mrdoc2Res |> 
  full_join(smk2dnamRes) |> 
  full_join( dnam2smkRes ) 


## Save ##########################
save(mrdoc_currSmk_Res, file = paste0(outDir,"mrdoc_currentSmk_DNAm_results.RData"))


# Estimates of Top Sites #####################################


### Smk --> DNAm #####################

nTest2 <- mrdoc_currSmk_Res |> 
  filter(!is.na(mrdoc2.g1_bonfSig)) |> 
  nrow()
bonfp2 <- 0.05/nTest2
bonfZ2 <- qnorm(p = 1-bonfp2/2)

nTest1 <- mrdoc_currSmk_Res |> 
  filter(!is.na(mrdoc1_wPleio.g1_bonfSig)) |> 
  nrow()
bonfp1 <- 0.05/nTest1
bonfZ1 <- qnorm(p = 1-bonfp1/2)

mrdoc_currSmk_DNAm_Est <- mrdoc_currSmk_Res |> 
  filter( mrdoc2.g1_bonfSig == T ) |> 
  mutate(cpg = fct_reorder(as_factor(cpg), - mrdoc2.g1_Z)) |> 
  select(cpg, CHR, MAPINFO, SYMBOL,
         mrdoc2.g1_hat, mrdoc2.g1_SE, mrdoc2.g1_Z,
         mrdoc1_wPleio.g1_hat, mrdoc1_wPleio.g1_SE, mrdoc1_wPleio.g1_Z,
         mrdoc1_wRE.g1_hat, mrdoc1_wRE.g1_SE, mrdoc1_wRE.g1_Z) |> 
  pivot_longer(cols = -c(cpg, CHR, MAPINFO, SYMBOL), 
               names_to = c("model","par"), names_sep = ".g1_") |> 
  pivot_wider(id_cols = c("cpg","CHR","MAPINFO","SYMBOL","model"), 
              names_from = "par", names_prefix = "g1_") |> 
  mutate(model = fct_relevel(model, "mrdoc2"),
         model = fct_recode(model, 
                            "MR-DoC2" = "mrdoc2",
                            "MR-DoC1 w/ Pleiotropic Path" = "mrdoc1_wPleio",
                            "MR-DoC1 w/ rE" = "mrdoc1_wRE") ) |> 
  group_by(cpg) |> 
  arrange(desc(abs(g1_hat))) |> 
  mutate(top = row_number(),
         plotLabel = ifelse(top == 1, SYMBOL, NA) ) |> 
  ungroup() |> 
  arrange(cpg) 

head(mrdoc_currSmk_DNAm_Est)


#### FDR Intersection Estimates ####################

## How many are in the MHC region?
# chr6:28,477,797-33,448,354
# https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37 

mrdoc_currSmk_Res |>
  filter( (mrdoc2.g1_fdrSig &
             mrdoc1_wPleio.g1_fdrSig &
             mrdoc1_wRE.g1_fdrSig) ,
          CHR == 6,
          MAPINFO >= 28477797,
          MAPINFO <= 33448354) |>
  select(cpg,CHR,MAPINFO,SYMBOL)

## One half of the FDR sites with negative estimate
mrdoc_currSmk_DNAm_FDRneg1 <- mrdoc_currSmk_Res |>
  filter( (mrdoc2.g1_fdrSig &
             mrdoc1_wPleio.g1_fdrSig &
             mrdoc1_wRE.g1_fdrSig) ,
          mrdoc2.g1_hat < 0 ) |> 
  ## Arrange by MR-DoC1 Z score
  arrange( mrdoc1_wPleio.g1_Z ) |> 
  ## Keep the 1st Half
  slice(1:30) |> 
  ## Pivot longer
  mutate(cpg = fct_reorder(as_factor(cpg), - mrdoc2.g1_Z)) |> 
  select(cpg, CHR, SYMBOL,
         mrdoc2.g1_hat, mrdoc2.g1_SE, mrdoc2.g1_Z,
         mrdoc1_wPleio.g1_hat, mrdoc1_wPleio.g1_SE, mrdoc1_wPleio.g1_Z,
         mrdoc1_wRE.g1_hat, mrdoc1_wRE.g1_SE, mrdoc1_wRE.g1_Z) |> 
  pivot_longer(cols = -c(cpg, CHR, SYMBOL), 
               names_to = c("model","par"), names_sep = ".g1_") |> 
  pivot_wider(id_cols = c("cpg","CHR","SYMBOL","model"), 
              names_from = "par", names_prefix = "g1_") |> 
  mutate(model = fct_relevel(model, "mrdoc2"),
         model = fct_recode(model, 
                            "MR-DoC2" = "mrdoc2",
                            "MR-DoC1 w/ Pleiotropic Path" = "mrdoc1_wPleio",
                            "MR-DoC1 w/ rE" = "mrdoc1_wRE") ) 

mrdoc_currSmk_DNAm_FDR_plot1 <- mrdoc_currSmk_DNAm_FDRneg1 |> 
  ggplot(aes(g1_hat, cpg)) + 
  geom_point(aes(color = model),
             size = 1.25, position = position_dodge(width = 0.5) ) +
  geom_errorbarh(aes(color = model,
                     xmax = g1_hat + qnorm(p=0.975)*g1_SE,
                     xmin = g1_hat - qnorm(p=0.975)*g1_SE),
                 height = 0, linewidth = 0.6,
                 position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, color = "grey50", linetype = "dashed") +
  labs( x = "Causal Estimate (95% C.I.)", y = NULL, color = NULL, fill = NULL) +
  scale_color_manual(values = c("steelblue4","slategray","lightskyblue3")) +
  scale_fill_manual(values = c("steelblue4","slategray","lightskyblue3")) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 6, face = "bold"),
        axis.title = element_text(size = 9),
        legend.text = element_text(size =8) )


## Second half of the FDR sites with negative estimate
mrdoc_currSmk_DNAm_FDRneg2 <- mrdoc_currSmk_Res |>
  filter( (mrdoc2.g1_fdrSig &
             mrdoc1_wPleio.g1_fdrSig &
             mrdoc1_wRE.g1_fdrSig) ,
          mrdoc2.g1_hat < 0 ) |> 
  ## Arrange by MR-DoC1 Z score
  arrange( mrdoc1_wPleio.g1_Z ) |> 
  ## Keep the 1st Half
  slice(31:59) |> 
  ## Pivot longer
  mutate(cpg = fct_reorder(as_factor(cpg), - mrdoc2.g1_Z)) |> 
  select(cpg, CHR, SYMBOL,
         mrdoc2.g1_hat, mrdoc2.g1_SE, mrdoc2.g1_Z,
         mrdoc1_wPleio.g1_hat, mrdoc1_wPleio.g1_SE, mrdoc1_wPleio.g1_Z,
         mrdoc1_wRE.g1_hat, mrdoc1_wRE.g1_SE, mrdoc1_wRE.g1_Z) |> 
  pivot_longer(cols = -c(cpg, CHR, SYMBOL), 
               names_to = c("model","par"), names_sep = ".g1_") |> 
  pivot_wider(id_cols = c("cpg","CHR","SYMBOL","model"), 
              names_from = "par", names_prefix = "g1_") |> 
  mutate(model = fct_relevel(model, "mrdoc2"),
         model = fct_recode(model, 
                            "MR-DoC2" = "mrdoc2",
                            "MR-DoC1 w/ Pleiotropic Path" = "mrdoc1_wPleio",
                            "MR-DoC1 w/ rE" = "mrdoc1_wRE") ) 


mrdoc_currSmk_DNAm_FDR_plot2 <- mrdoc_currSmk_DNAm_FDRneg2 |> 
  ggplot(aes(g1_hat, cpg)) + 
  geom_point(aes(color = model),
             size = 1.25, position = position_dodge(width = 0.5) ) +
  geom_errorbarh(aes(color = model,
                     xmax = g1_hat + qnorm(p=0.975)*g1_SE,
                     xmin = g1_hat - qnorm(p=0.975)*g1_SE),
                 height = 0, linewidth = 0.6,
                 position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, color = "grey50", linetype = "dashed") +
  labs( x = "Causal Estimate (95% C.I.)", y = NULL, color = NULL, fill = NULL) +
  scale_color_manual(values = c("steelblue4","slategray","lightskyblue3")) +
  scale_fill_manual(values = c("steelblue4","slategray","lightskyblue3")) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 6, face = "bold"),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8) )


## Combine
nFDRneg <- mrdoc_currSmk_Res |>
  filter( mrdoc2.g1_fdrSig &
            mrdoc1_wPleio.g1_fdrSig &
            mrdoc1_wRE.g1_fdrSig &
            mrdoc2.g1_hat < 0) |> 
  nrow()

library(patchwork)
jpeg(paste0(plotDir,"mrdoc_currentSmk_dnam_estimates_suppl.jpeg"), 
     width = 9, height = 6, units = "in", res = 300)

patchFDR <- mrdoc_currSmk_DNAm_FDR_plot1 + mrdoc_currSmk_DNAm_FDR_plot2 +
  plot_layout(guides = "collect")

patchFDR +
  plot_annotation(
    title = "CpG Sites where Current Smoking Likely Causes Hypomethylation",
    caption = paste0(nFDRneg," CG Sites with FDR < 0.05 in All Three Models")
  ) & 
  theme( plot.title = element_text(size = 10, face = "bold"),
         legend.position = "bottom" )

dev.off()


#### FDR Intersection Z-score MR-DoC_w_rE ####################

## Plot the estimates from MR-DoC2
nTests <- nrow(mrdoc_currSmk_Res)
bonfp  <- 0.05/nTests
bonfZ  <- qnorm(p = 1-bonfp/2)

## Plot Label = Bonferroni Significant in >1 models
mrdoc_currSmk_DNAm_FDR_Est <- mrdoc_currSmk_Res |> 
  rowwise() |> 
  mutate(
    sumOverlap = sum(mrdoc2.g1_bonfSig, 
                     mrdoc1_wPleio.g1_bonfSig, 
                     mrdoc1_wRE.g1_bonfSig, 
                     na.rm = T),
    plotLabel = ifelse( sumOverlap > 1  ,
                        SYMBOL, NA)
  ) 


## Factor for alpha = if FDR < 0.05 in all 3 models

mrdoc_currSmk_DNAm_FDR_Est <- mrdoc_currSmk_DNAm_FDR_Est |> 
  mutate(alfa = ifelse( (mrdoc2.g1_fdrSig & 
                           mrdoc1_wPleio.g1_fdrSig &
                           mrdoc1_wRE.g1_fdrSig ) , 
                        "sig", "nonSig") ) 

## Some alfa values are assigned NA (due to missing MR-DoC2 results)
## recode these as nonSig

mrdoc_currSmk_DNAm_FDR_Est <- mrdoc_currSmk_DNAm_FDR_Est |> 
  mutate(alfa = ifelse(is.na(alfa), "nonSig", alfa ) ) 

nFDR <- mrdoc_currSmk_DNAm_FDR_Est |> 
  filter(alfa == "sig") |> 
  nrow()

nBonf <- mrdoc_currSmk_DNAm_FDR_Est |> 
  filter(sumOverlap > 1) |> 
  nrow()

## Hard-coding "Manhattan" plots with Z scores
## Based on https://r-graph-gallery.com/101_Manhattan_plot.html

## CHR as numeric 
class(mrdoc_currSmk_DNAm_FDR_Est$CHR)

## Compute the cumulative position of CpG
mrdocOutPos <-  mrdoc_currSmk_DNAm_FDR_Est |>
  filter(!is.na(mrdoc1_wRE.g1_Z)) |> 
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(MAPINFO,na.rm = T)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(-chr_len) %>% 
  
  # Add this info to the initial dataset
  left_join(mrdoc_currSmk_DNAm_FDR_Est, ., by=c("CHR"="CHR")) %>% 
  
  # Add a cumulative position of each SNP
  arrange(MAPINFO) %>%
  mutate( BPcum = MAPINFO + tot)

summary(mrdocOutPos$BPcum)

## X-axis
axisdf <- mrdocOutPos %>% 
  group_by(CHR) %>% 
  summarize(center=( max(BPcum, na.rm = T) + min(BPcum, na.rm = T) ) / 2 )

glimpse(axisdf)

## Gene Labels = drop repeated plotlabel
mrdocOutPos <- mrdocOutPos |> 
  group_by(plotLabel) |> 
  arrange( mrdoc1_wRE.g1_p) |> 
  mutate(nCPG = row_number(), 
         geneLabel = ifelse(nCPG < 2, plotLabel, NA)) |> 
  ungroup() 

## Plot with highlighted sites with FDR < 0.05 in both models
mrdocPlot <- mrdocOutPos |> 
  ggplot(aes(x=BPcum, y= mrdoc1_wRE.g1_Z, label=geneLabel) ) +
  
  # Show all points
  geom_point( aes( color = as.factor(CHR), 
                   # fill = after_scale(alpha(color, 0.75)),
                   shape = alfa,
                   alpha = alfa ,
                   size = alfa ) ) +
  geom_label_repel(data = subset(mrdocOutPos,  mrdoc1_wRE.g1_Z < 0),
                   seed = 1, 
                   size = 3,
                   # mapping = aes(fontface = "bold"),
                   min.segment.length = 0, na.rm = T, position = position_nudge_repel(y=-1),
                   ylim = c(NA, -6)) +
  geom_label_repel(data = subset(mrdocOutPos,  mrdoc1_wRE.g1_Z > 0),
                   seed = 1, 
                   size = 3,
                   # mapping = aes(fontface = "bold"),
                   min.segment.length = 0, na.rm = T, position = position_nudge_repel(y=1) ) +
  scale_color_manual(values = rep(c( "steelblue3", "grey50"), 22 )) +
  scale_size_manual(values = c(1.25, 2.5)) +
  scale_alpha_manual(values = c(0.25, 0.9)) +
  scale_shape_manual(values = c(1, 19)) +
  
  # custom X axis:
  scale_x_continuous( expand = c(0.03, 0),
                      label = axisdf$CHR, breaks= axisdf$center, name = "Chromosome" ) +
  scale_y_continuous( name = "Causal Effect (Z Score)" ) +
  coord_cartesian(clip = "off") +
  geom_hline(yintercept = 0, color = "grey30") +
  geom_hline(yintercept = bonfZ, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -bonfZ, linetype = "dashed", color = "grey50") +
  labs(title = "Putative Causal Effects of Current Smoking on DNAm",
       subtitle = "Causal Estimates from MR-DoC1 Model with rE",
       caption = paste0("Solid Points indicate the ", nFDR," sites significant at FDR < 0.05 in all three models\n",
                        "Labeled Points indicate the ", nBonf," sites significant after Bonferroni correction in more than one model")) +
  # Custom the theme:
  theme_bw(8) +
  theme( plot.title = element_text(size = 10),
         plot.caption = element_text(size = 7),
         legend.position="none" , 
         # text = element_text(face = "bold"),
         panel.grid.minor.x = element_blank(),
         panel.border = element_blank()
  )

jpeg(paste0(plotDir,"mrdoc_currentSmk_dnam_wRE_fdr_estimates.jpeg"), 
     width = 9, height = 6, units = "in", res = 300)
mrdocPlot
dev.off()



### Evidence for bidirection effects? ===========================================

## At the 64 sites, look at the estimates of DNAm --> Smk

mrdoc_currSmk_bidir <- mrdoc_currSmk_Res |>
  ## filter FDR < 0.05 in all three
  filter(mrdoc2.g1_fdrSig & 
           mrdoc1_wPleio.g1_fdrSig &
           mrdoc1_wRE.g1_fdrSig) |> 
  mutate(
    ## if DNAm --> Smk has p < 0.05 in all three
    alfa = ifelse( ( 
      ## all three negative and p < 0.05
      (mrdoc2.g2_Z < -qnorm(p=0.975) &
         mrdoc1_wPleio.g2_Z < -qnorm(p=0.975) & 
         mrdoc1_wRE.g2_Z < -qnorm(p=0.975)) | 
        ## OR all three positive and p < 0.05
        (mrdoc2.g2_Z > qnorm(p=0.975) &
           mrdoc1_wPleio.g2_Z > qnorm(p=0.975) & 
           mrdoc1_wRE.g2_Z > qnorm(p=0.975)) ), 
      1, 0 ),
    alfa = factor(alfa),
    cpg = fct_reorder(as_factor(cpg), - mrdoc1_wRE.g1_Z)
  ) |> 
  select(cpg, CHR, MAPINFO, SYMBOL, alfa, 
         mrdoc1_wPleio.g1_hat, mrdoc1_wPleio.g1_SE, mrdoc1_wPleio.g1_Z, mrdoc1_wPleio.g1_p, 
         mrdoc1_wRE.g1_hat, mrdoc1_wRE.g1_SE, mrdoc1_wRE.g1_Z, mrdoc1_wRE.g1_p, 
         mrdoc2.g1_hat, mrdoc2.g1_SE, mrdoc2.g1_Z, mrdoc2.g1_p,
         mrdoc1_wPleio.g2_hat, mrdoc1_wPleio.g2_SE, mrdoc1_wPleio.g2_Z, mrdoc1_wPleio.g2_p, 
         mrdoc1_wRE.g2_hat, mrdoc1_wRE.g2_SE, mrdoc1_wRE.g2_Z, mrdoc1_wRE.g2_p, 
         mrdoc2.g2_hat, mrdoc2.g2_SE, mrdoc2.g2_Z, mrdoc2.g2_p) |> 
  pivot_longer(cols = -c(cpg, CHR, MAPINFO, SYMBOL, alfa), 
               names_to = c("model","par"), names_sep = ".g") |> 
  pivot_wider(id_cols = c("cpg","CHR","MAPINFO","SYMBOL","model","alfa"), 
              names_from = "par") |> 
  pivot_longer(cols = -c(cpg, CHR, MAPINFO, SYMBOL, model, alfa), 
               names_to = c("doc","par"), names_sep = "_") |> 
  pivot_wider(id_cols = c("cpg","CHR","MAPINFO","SYMBOL","doc","model","alfa"), 
              names_from = "par") |> 
  mutate(
    doc = factor(doc),
    doc = fct_recode(doc,
                     "Effect of Smoking on DNAm" = "1",
                     "Effect of DNAm on Smoking" = "2"),
    model = fct_relevel(model, "mrdoc2"),
    model = fct_recode(model, 
                       "MR-DoC2" = "mrdoc2",
                       "MR-DoC1 w/ Pleiotropic Path" = "mrdoc1_wPleio",
                       "MR-DoC1 w/ rE" = "mrdoc1_wRE") 
  ) 



jpeg(paste0(plotDir,"mrdoc_compare_bidir_FDR_estimates.jpeg"), 
     width = 7, height = 8, units = "in", res = 300)

mrdoc_currSmk_bidir |> 
  ggplot(aes(x = hat, y = cpg)) + 
  geom_point(aes(color = model, 
                 fill = model, 
                 # alpha = alfa, 
                 shape = model),
             size = 0.75, 
             alpha = 0.9,
             position = position_dodge(width = 0.6) ) +
  geom_errorbarh(aes(color = model, 
                     # alpha = alfa,
                     xmax = hat + qnorm(p=0.975)*SE,
                     xmin = hat - qnorm(p=0.975)*SE),
                 height = 0, 
                 linewidth = 0.25,
                 position = position_dodge(width = 0.6)) +
  facet_wrap(~ doc) +
  geom_vline(xintercept = 0, 
             color = "midnightblue", 
             linetype = "dashed", 
             linewidth = 0.25) +
  labs(title = "Bidirectional Causal Estimates between Current Smoking and DNAm",
       subtitle = "At 64 CpGs where Current Smoking Likely Affects DNAm", 
       x = "Causal Estimate (95% C.I.)", y = NULL, 
       color = "Model", fill = "Model", shape = "Model") +
  scale_color_manual(values = c("steelblue4","lightskyblue3","slategray")) +
  scale_fill_manual(values = c("steelblue4","lightskyblue3","slategray")) +
  # scale_alpha_manual(values = c(0.25, 1)) +
  scale_shape_manual(values = c(15, 19, 17)) +
  guides(color = guide_legend(override.aes = list(size = 3))) + 
  theme_light(10) +
  theme(
    plot.title.position = "plot",
    axis.text.y = element_text(size = 6),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom" 
  ) 


dev.off()






## DNAm --> Smk #####################

nTest2 <- mrdoc_currSmk_Res |> 
  filter(!is.na(mrdoc2.g2_bonfSig)) |> 
  nrow()
bonfp2 <- 0.05/nTest2
bonfZ2 <- qnorm(p = 1-bonfp2/2)

nTest1 <- mrdoc_currSmk_Res |> 
  filter(!is.na(mrdoc1_wPleio.g2_bonfSig)) |> 
  nrow()
bonfp1 <- 0.05/nTest1
bonfZ1 <- qnorm(p = 1-bonfp1/2)


#### Estimate ####################


mrdoc_DNAm_currSmk_Est <- mrdoc_currSmk_Res |> 
  filter( mrdoc2.g2_bonfSig == T ) |> 
  mutate(cpg = fct_reorder(as_factor(cpg), - mrdoc2.g2_Z)) |> 
  select(cpg, CHR, SYMBOL,
         mrdoc2.g2_hat, mrdoc2.g2_SE, mrdoc2.g2_Z,
         mrdoc1_wPleio.g2_hat, mrdoc1_wPleio.g2_SE, mrdoc1_wPleio.g2_Z,
         mrdoc1_wRE.g2_hat, mrdoc1_wRE.g2_SE, mrdoc1_wRE.g2_Z) |> 
  pivot_longer(cols = -c(cpg, CHR, SYMBOL), 
               names_to = c("model","par"), names_sep = ".g2_") |> 
  pivot_wider(id_cols = c("cpg","CHR","SYMBOL","model"), 
              names_from = "par", names_prefix = "g2_") |> 
  mutate(model = fct_relevel(model, "mrdoc2"),
         model = fct_recode(model, 
                            "MR-DoC2" = "mrdoc2",
                            "MR-DoC1 w/ Pleiotropic Path" = "mrdoc1_wPleio",
                            "MR-DoC1 w/ rE" = "mrdoc1_wRE") ) |> 
  group_by(cpg) |> 
  arrange(desc(abs(g2_hat))) |> 
  mutate(top = row_number(),
         plotLabel = ifelse(top == 1, SYMBOL, NA) ) |> 
  ungroup() |> 
  arrange(cpg) 

head(mrdoc_DNAm_currSmk_Est)



jpeg(paste0(plotDir,"mrdoc_dnam_currentSmk_estimates.jpeg"), width = 6, height = 6, units = "in", res = 300)

mrdoc_DNAm_currSmk_Est |> 
  ggplot(aes(g2_hat, cpg, label = plotLabel)) + 
  geom_point(aes(color = model),
             size = 3, position = position_dodge(width = 0.5) ) +
  geom_errorbarh(aes(color = model,
                     xmax = g2_hat + qnorm(p=0.975)*g2_SE,
                     xmin = g2_hat - qnorm(p=0.975)*g2_SE),
                 height = 0, size = 1,
                 position = position_dodge(width = 0.5)) +
  geom_text_repel(data = subset(mrdoc_DNAm_currSmk_Est,  g2_Z < 0),
                  seed = 1, mapping = aes(fontface = "bold"), 
                  segment.size = 0, na.rm = T, 
                  direction = "x", position = position_nudge_repel(x=-0.9) 
                  , hjust = "left" , xlim = c(NA,-1.1)
  ) +
  geom_text_repel(data = subset(mrdoc_DNAm_currSmk_Est,  g2_Z > 0),
                  seed = 1, mapping = aes(fontface = "bold"), 
                  segment.size = 0, na.rm = T, 
                  direction = "x", position = position_nudge_repel(x=0.75) 
                  , hjust = "right" , xlim = c(1.2,NA)
  ) +
  geom_vline(xintercept = 0, color = "grey50") +
  labs(title = "Putative Causal Effects of DNAm on Current Smoking",
       subtitle = paste("Sites with Bonferroni-significant Estimate in Bidirectional MR-DoC2"), 
       caption = paste0(nTest1," CpG Sites Tested"),
       x = "Causal Estimate (95% C.I.)", y = NULL, color = NULL, fill = NULL) +
  coord_cartesian(xlim = c(-1.74,1.8)) +
  scale_color_manual(values = c("indianred3","thistle4","mistyrose3")) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        legend.position = "bottom")

dev.off()


#### FDR Intersection Estimate ####################

mrdoc_DNAm_currSmk_Est <- mrdoc_currSmk_Res |> 
  filter( mrdoc2.g2_fdrSig == T &
            mrdoc1_wPleio.g2_fdrSig == T & 
            mrdoc1_wRE.g2_fdrSig == T ) |> 
  mutate(cpg = fct_reorder(as_factor(cpg), - mrdoc2.g2_Z)) |> 
  select(cpg, CHR, SYMBOL,
         mrdoc2.g2_hat, mrdoc2.g2_SE, mrdoc2.g2_Z,
         mrdoc1_wPleio.g2_hat, mrdoc1_wPleio.g2_SE, mrdoc1_wPleio.g2_Z,
         mrdoc1_wRE.g2_hat, mrdoc1_wRE.g2_SE, mrdoc1_wRE.g2_Z) |> 
  pivot_longer(cols = -c(cpg, CHR, SYMBOL), 
               names_to = c("model","par"), names_sep = ".g2_") |> 
  pivot_wider(id_cols = c("cpg","CHR","SYMBOL","model"), 
              names_from = "par", names_prefix = "g2_") |> 
  mutate(model = fct_relevel(model, "mrdoc2"),
         model = fct_recode(model, 
                            "MR-DoC2" = "mrdoc2",
                            "MR-DoC1 w/ Pleiotropic Path" = "mrdoc1_wPleio",
                            "MR-DoC1 w/ rE" = "mrdoc1_wRE") ) |> 
  group_by(cpg) |> 
  arrange(desc(abs(g2_hat))) |> 
  mutate(top = row_number(),
         plotLabel = ifelse(top == 1, SYMBOL, NA) ) |> 
  ungroup() |> 
  arrange(cpg) 

jpeg(paste0(plotDir,"mrdoc_dnam_currentSmk_fdr_estimates.jpeg"), width = 6, height = 6, units = "in", res = 300)

mrdoc_DNAm_currSmk_Est |> 
  ggplot(aes(g2_hat, cpg, label = plotLabel)) + 
  geom_point(aes(color = model, shape = model),
             size = 3, position = position_dodge(width = 0.5) ) +
  geom_errorbarh(aes(color = model,
                     xmax = g2_hat + qnorm(p=0.975)*g2_SE,
                     xmin = g2_hat - qnorm(p=0.975)*g2_SE),
                 height = 0, linewidth = 1,
                 position = position_dodge(width = 0.5)) +
  geom_text_repel(data = subset(mrdoc_DNAm_currSmk_Est,  g2_Z > 0),
                  seed = 1, mapping = aes(fontface = "bold"), 
                  segment.size = 0, na.rm = T, 
                  direction = "x", #position = position_nudge_repel(x=0.75), 
                  hjust = "right" , xlim = c(1.2,NA)
  ) +
  geom_vline(xintercept = 0, color = "grey50", linetype = "dashed") +
  labs(title = "Putative Causal Effects of DNAm on Current Smoking",
       subtitle = paste("Sites with FDR < 0.05 in All Three Models"), 
       caption = paste0(nTest1," CpG Sites Tested"),
       x = "Causal Estimate (95% C.I.)", y = NULL, color = NULL, shape = NULL) +
  coord_cartesian(xlim = c(0,1.5)) +
  scale_color_manual(values = c("indianred3","thistle4","mistyrose3")) +
  scale_shape_manual(values = c(17, 19, 15)) +
  theme_bw() +
  theme(#plot.title = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold", size = 12),
    # axis.text.x = element_text(face = "bold"),
    # axis.title.x = element_text(face = "bold"),
    legend.position = "bottom")

dev.off()




### Nominal Intersection Estimate ===============================================


mrdoc_DNAm_currSmk_nominal_Est <- mrdoc_currSmk_Res |> 
  filter( 
    ## all three negative and p < 0.05
    (mrdoc2.g2_Z < -qnorm(p=0.975) &
       mrdoc1_wPleio.g2_Z < -qnorm(p=0.975) & 
       mrdoc1_wRE.g2_Z < -qnorm(p=0.975)) | 
      ## OR all three positive and p < 0.05
      (mrdoc2.g2_Z > qnorm(p=0.975) &
         mrdoc1_wPleio.g2_Z > qnorm(p=0.975) & 
         mrdoc1_wRE.g2_Z > qnorm(p=0.975)) 
  ) |> 
  mutate(cpg = fct_reorder(as_factor(cpg), - mrdoc2.g2_Z)) |> 
  select(cpg, CHR, SYMBOL,
         mrdoc2.g2_hat, mrdoc2.g2_SE, mrdoc2.g2_Z,
         mrdoc1_wPleio.g2_hat, mrdoc1_wPleio.g2_SE, mrdoc1_wPleio.g2_Z,
         mrdoc1_wRE.g2_hat, mrdoc1_wRE.g2_SE, mrdoc1_wRE.g2_Z) |> 
  pivot_longer(cols = -c(cpg, CHR, SYMBOL), 
               names_to = c("model","par"), names_sep = ".g2_") |> 
  pivot_wider(id_cols = c("cpg","CHR","SYMBOL","model"), 
              names_from = "par", names_prefix = "g2_") |> 
  mutate(model = fct_relevel(model, "mrdoc2"),
         model = fct_recode(model, 
                            "MR-DoC2" = "mrdoc2",
                            "MR-DoC1 w/ Pleiotropic Path" = "mrdoc1_wPleio",
                            "MR-DoC1 w/ rE" = "mrdoc1_wRE") ) |> 
  group_by(cpg) |> 
  arrange(desc(abs(g2_hat))) |> 
  mutate(top = row_number(),
         plotLabel = ifelse(top == 1, SYMBOL, NA) ) |> 
  ungroup() |> 
  arrange(cpg) 


## Split into two parts: positive and negative
mrdoc_DNAm_currSmk_nominal_pos <- mrdoc_DNAm_currSmk_nominal_Est |> 
  filter(g2_hat > 0)
mrdoc_DNAm_currSmk_nominal_pos |> 
  count(cpg) |> 
  nrow()

mrdoc_DNAm_currSmk_nominal_neg <- mrdoc_DNAm_currSmk_nominal_Est |> 
  filter(g2_hat < 0)
mrdoc_DNAm_currSmk_nominal_neg |> 
  count(cpg) |> 
  nrow()

mrdoc_nominal_dnam_pos <- mrdoc_DNAm_currSmk_nominal_pos |> 
  ggplot(aes(g2_hat, cpg
             #, label = plotLabel
  )) + 
  geom_point(aes(color = model),
             size = 1, position = position_dodge(width = 0.5) ) +
  geom_errorbarh(aes(color = model,
                     xmax = g2_hat + qnorm(p=0.975)*g2_SE,
                     xmin = g2_hat - qnorm(p=0.975)*g2_SE),
                 height = 0, size = 0.5,
                 position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, color = "grey50", linetype = "dashed") +
  labs(x = "Causal Estimate (95% C.I.)", y = NULL, color = NULL, fill = NULL) +
  scale_color_manual(values = c("indianred3","thistle4","mistyrose3")) +
  theme_light() +
  theme(axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        legend.position = "bottom")



mrdoc_nominal_dnam_neg <- mrdoc_DNAm_currSmk_nominal_neg |> 
  ggplot(aes(g2_hat, cpg
  )) + 
  geom_point(aes(color = model),
             size = 1, position = position_dodge(width = 0.5) ) +
  geom_errorbarh(aes(color = model,
                     xmax = g2_hat + qnorm(p=0.975)*g2_SE,
                     xmin = g2_hat - qnorm(p=0.975)*g2_SE),
                 height = 0, size = 0.5,
                 position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, color = "grey50", linetype = "dashed") +
  labs(x = "Causal Estimate (95% C.I.)", y = NULL, color = NULL, fill = NULL) +
  scale_color_manual(values = c("indianred3","thistle4","mistyrose3")) +
  scale_y_discrete(position = "right") +
  theme_light() +
  theme(axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        legend.position = "bottom")


## Patch
jpeg(paste0(plotDir,"mrdoc_dnam_currentSmk_nominal_estimates.jpeg"), 
     width = 6, height = 8, units = "in", res = 300)

mrdoc_nominal_dnam_pos + mrdoc_nominal_dnam_neg + 
  plot_layout(guides = 'collect')  + 
  plot_annotation(
    title = "CpGs with Consistent Estimates of the Effects of DNAm on Current Smoking"
  )  & 
  theme(title = element_text(size = 8),
        legend.position = "bottom")

dev.off()





### Evidence for bidirection effects? ===========================================

## At the 3 sites, look at the estimates of Smk --> DNAm
mrdoc_currSmk_bidir2 <- mrdoc_currSmk_Res |>
  ## filter FDR < 0.05 in all three
  filter(mrdoc2.g2_fdrSig & 
           mrdoc1_wPleio.g2_fdrSig &
           mrdoc1_wRE.g2_fdrSig) |> 
  mutate(
    alfa = ifelse( ( 
      ## all three negative and p < 0.05
      (mrdoc2.g1_Z < -qnorm(p=0.975) &
         mrdoc1_wPleio.g1_Z < -qnorm(p=0.975) & 
         mrdoc1_wRE.g1_Z < -qnorm(p=0.975)) | 
        ## OR all three positive and p < 0.05
        (mrdoc2.g1_Z > qnorm(p=0.975) &
           mrdoc1_wPleio.g1_Z > qnorm(p=0.975) & 
           mrdoc1_wRE.g1_Z > qnorm(p=0.975)) ), 
      1, 0 ),
    alfa = factor(alfa),
    cpg = fct_reorder(as_factor(cpg), - mrdoc1_wRE.g2_Z)
  ) |> 
  select(cpg, CHR, MAPINFO, SYMBOL, alfa, 
         mrdoc1_wPleio.g1_hat, mrdoc1_wPleio.g1_SE, mrdoc1_wPleio.g1_Z, mrdoc1_wPleio.g1_p, 
         mrdoc1_wRE.g1_hat, mrdoc1_wRE.g1_SE, mrdoc1_wRE.g1_Z, mrdoc1_wRE.g1_p, 
         mrdoc2.g1_hat, mrdoc2.g1_SE, mrdoc2.g1_Z, mrdoc2.g1_p,
         mrdoc1_wPleio.g2_hat, mrdoc1_wPleio.g2_SE, mrdoc1_wPleio.g2_Z, mrdoc1_wPleio.g2_p, 
         mrdoc1_wRE.g2_hat, mrdoc1_wRE.g2_SE, mrdoc1_wRE.g2_Z, mrdoc1_wRE.g2_p, 
         mrdoc2.g2_hat, mrdoc2.g2_SE, mrdoc2.g2_Z, mrdoc2.g2_p) |> 
  pivot_longer(cols = -c(cpg, CHR, MAPINFO, SYMBOL, alfa), 
               names_to = c("model","par"), names_sep = ".g") |> 
  pivot_wider(id_cols = c("cpg","CHR","MAPINFO","SYMBOL","model","alfa"), 
              names_from = "par") |> 
  pivot_longer(cols = -c(cpg, CHR, MAPINFO, SYMBOL, model, alfa), 
               names_to = c("doc","par"), names_sep = "_") |> 
  pivot_wider(id_cols = c("cpg","CHR","MAPINFO","SYMBOL","doc","model","alfa"), 
              names_from = "par") |> 
  mutate(
    doc = fct_rev(factor(doc)),
    doc = fct_recode(doc,
                     "Effect of Smoking on DNAm" = "1",
                     "Effect of DNAm on Smoking" = "2"),
    model = fct_relevel(model, "mrdoc2"),
    model = fct_recode(model, 
                       "MR-DoC2" = "mrdoc2",
                       "MR-DoC1 w/ Pleiotropic Path" = "mrdoc1_wPleio",
                       "MR-DoC1 w/ rE" = "mrdoc1_wRE") 
  ) 


jpeg(paste0(plotDir,"mrdoc_compare_bidir2_FDR_estimates.jpeg"), 
     width = 4, height = 3, units = "in", res = 300)

mrdoc_currSmk_bidir2 |> 
  mutate(cpg = paste(cpg,SYMBOL, sep="\n")) |> 
  ggplot(aes(x = hat, y = cpg)) + 
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
  facet_wrap(~ doc) +
  geom_vline(xintercept = 0, 
             color = "midnightblue", 
             linetype = "dashed", 
             linewidth = 0.25) +
  labs(title = "Bidirectional Causal Estimates between Current Smoking and DNAm",
       subtitle = "At 3 CpGs where DNAm Likely Affects Current Smoking", 
       x = "Causal Estimate (S.D.)", y = NULL, 
       color = "Model", fill = "Model", shape = "Model") +
  scale_color_manual(values = c("indianred3","thistle4","mistyrose3")) +
  scale_fill_manual(values = c("indianred3","thistle4","mistyrose3")) +
  # scale_alpha_manual(values = c(0.25, 1)) +
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





## Potential Bidirectional Effects ######################################


## Part 1
mrdoc_currSmk_bidir_sub1 <- mrdoc_currSmk_bidir |> 
  filter(alfa == "1") |> 
  ## Combine cpg with gene symbol for plot axis text
  mutate(cpg = paste(cpg,SYMBOL, sep="\n"))

mrdoc_currSmk_bidir_sub1$robust <- "Stronger Evidence\nof Effects of\nSmoking on DNAm"


mrdoc_currSmk_bidir_sub1 |> 
  count(doc)


## Part 2
mrdoc_currSmk_bidir_sub2 <- mrdoc_currSmk_bidir2 |> 
  filter(alfa == "1") |> 
  ## Combine cpg with gene symbol for plot axis text
  mutate(cpg = paste(cpg,SYMBOL, sep="\n"))

mrdoc_currSmk_bidir_sub2$robust <- "Stronger Evidence\nof Effects of\nDNAm on Smoking"

## Combine
mrdoc_currSmk_bidir_comb <- mrdoc_currSmk_bidir_sub1 |> 
  rbind(mrdoc_currSmk_bidir_sub2) |> 
  mutate(robust = fct_rev(robust))


## Try making each column separately and then patchwork

## Combine and re-split by Direction of Causation 
mrdoc_currSmk_bidir_smk <- mrdoc_currSmk_bidir_sub1 |>
  rbind(mrdoc_currSmk_bidir_sub2) |>
  filter(doc == "Effect of Smoking on DNAm") |> 
  mutate(robust = fct_rev(robust), 
         SYMBOL = fct_relevel(SYMBOL, "GNG7", after = Inf))

mrdoc_currSmk_bidir_dnam <- mrdoc_currSmk_bidir_sub1 |>
  rbind(mrdoc_currSmk_bidir_sub2) |>
  filter(doc == "Effect of DNAm on Smoking") |> 
  mutate(robust = fct_rev(robust), 
         SYMBOL = fct_relevel(SYMBOL, "GNG7", after = Inf))

## Smk --> DNAm
bidir_smk <- mrdoc_currSmk_bidir_smk |>
  ggplot(aes(x = hat, y = cpg)) + 
  geom_point(aes(shape = model, color = model ),
             size = 2, 
             position = position_dodge(width = 0.5) ) +
  geom_errorbarh(aes(color = model, 
                     xmax = hat + qnorm(p=0.975)*SE,
                     xmin = hat - qnorm(p=0.975)*SE),
                 height = 0, 
                 linewidth = 0.5,
                 position = position_dodge(width = 0.5)) +
  facet_grid(rows = vars(robust), cols = vars(doc), 
             scales = "free_y", space = "free") +
  geom_vline(xintercept = 0, 
             color = "midnightblue", 
             linetype = "dashed", 
             linewidth = 0.2) +
  coord_cartesian(xlim = c(-1.5, 1.1)) +
  labs(x = "Causal Estimate (95% C.I.)", y = NULL, 
       color = NULL, 
       shape = NULL) +
  scale_color_manual(values = c("steelblue4","slategray","lightskyblue3")) +
  scale_shape_manual(values = c(17, 19, 15)) +
  guides(color = guide_legend(override.aes = list(size = 1.5))) +
  theme_bw(6) +
  theme(
    axis.text.y = element_text(size = 6),
    strip.text.x = element_text(#face = "bold",
      size = 6),
    strip.text.y = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "pt"),
    legend.position = "bottom" ,
    # legend.position = c(0.82,0.86),
    # legend.title = element_text(size = 5, hjust = 0.5),
    legend.text = element_text(size = 5)
  ) 


## DNAm --> Smk
bidir_dnam <- mrdoc_currSmk_bidir_dnam |>
  ggplot(aes(x = hat, y = cpg)) + 
  geom_point(aes(shape = model, color = model ),
             size = 2, 
             position = position_dodge(width = 0.5) ) +
  geom_errorbarh(aes(color = model, 
                     xmax = hat + qnorm(p=0.975)*SE,
                     xmin = hat - qnorm(p=0.975)*SE),
                 height = 0, 
                 linewidth = 0.5,
                 position = position_dodge(width = 0.5)) +
  facet_grid(rows = vars(robust), cols = vars(doc), 
             scales = "free_y", space = "free") +
  geom_vline(xintercept = 0, 
             color = "midnightblue", 
             linetype = "dashed", 
             linewidth = 0.2) +
  coord_cartesian(xlim = c(-1.5, 1.1)) +
  labs(x = "Causal Estimate (95% C.I.)", y = NULL, 
       color = NULL, 
       shape = NULL) +
  scale_color_manual(values = c("indianred3","thistle4","mistyrose3")) +
  scale_shape_manual(values = c(17, 19, 15)) +
  guides(color = guide_legend(override.aes = list(size = 1.5))) +
  theme_bw(6) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(#face = "bold",
      size = 6),
    plot.margin = margin(0, 0, 0, 0, "pt"),
    legend.position = "bottom" ,
    # legend.position = c(0.82,0.86),
    # legend.title = element_text(size = 5, hjust = 0.5),
    legend.text = element_text(size = 5)
  ) 

## Combine
library(patchwork)
patchBidir <-  bidir_smk + bidir_dnam + 
  plot_layout(guides = 'keep')

jpeg(paste0(plotDir,"bidir_causal_estimates.jpeg"), 
     width = 7, height = 4, units = "in",
     res = 300)

patchBidir + 
  plot_annotation(
    title = "Suggestive Bidirectional Causal Effects between Current Smoking and Blood DNA Methylation"
  )  & 
  theme(title = element_text(size = 6))

dev.off()




## Suggestive/Nominal Bidir Effects ######################

mrdoc_currSmk_nominal <- mrdoc_currSmk_Res |>
  mutate(
    g1_nominal = ifelse( ( 
      ## all three negative and p < 0.05
      (mrdoc2.g1_Z < -qnorm(p=0.975) &
         mrdoc1_wPleio.g1_Z < -qnorm(p=0.975) & 
         mrdoc1_wRE.g1_Z < -qnorm(p=0.975)) | 
        ## OR all three positive and p < 0.05
        (mrdoc2.g1_Z > qnorm(p=0.975) &
           mrdoc1_wPleio.g1_Z > qnorm(p=0.975) & 
           mrdoc1_wRE.g1_Z > qnorm(p=0.975)) ), 
      TRUE, FALSE ),
    g2_nominal = ifelse( ( 
      ## all three negative and p < 0.05
      (mrdoc2.g2_Z < -qnorm(p=0.975) &
         mrdoc1_wPleio.g2_Z < -qnorm(p=0.975) & 
         mrdoc1_wRE.g2_Z < -qnorm(p=0.975)) | 
        ## OR all three positive and p < 0.05
        (mrdoc2.g2_Z > qnorm(p=0.975) &
           mrdoc1_wPleio.g2_Z > qnorm(p=0.975) & 
           mrdoc1_wRE.g2_Z > qnorm(p=0.975)) ), 
      TRUE, FALSE ), 
    robust_smk2dnam = ifelse(mrdoc2.g1_fdrSig & 
                               mrdoc1_wPleio.g1_fdrSig &
                               mrdoc1_wRE.g1_fdrSig, 
                             TRUE, FALSE),
    robust_dnam2smk = ifelse(mrdoc2.g2_fdrSig & 
                               mrdoc1_wPleio.g2_fdrSig &
                               mrdoc1_wRE.g2_fdrSig, 
                             TRUE, FALSE)) 


#### Save the CpG IDs ===============

dnam2smk_intersect <- mrdoc_currSmk_nominal |> 
  filter(g2_nominal) |> 
  select(cpg)
length(dnam2smk_intersect$cpg)

write.table(dnam2smk_intersect, file = paste0(outDir,"dnam_2_smk_mrdoc_nominal_intersect.tsv"), 
            col.names = F, row.names = F, quote = F)



## How many are in the MHC region?
# chr6:28,477,797-33,448,354
# https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37 

mrdoc_currSmk_nominal |> 
  filter( g2_nominal ,
          CHR == 6,
          MAPINFO >= 28477797,
          MAPINFO <= 33448354) |>
  select(cpg,CHR,MAPINFO,SYMBOL)
## None

smk2dnam_intersect <- mrdoc_currSmk_nominal |> 
  filter(g1_nominal) |> 
  select(cpg)
length(smk2dnam_intersect$cpg)

write.table(smk2dnam_intersect, file = paste0(outDir,"smk_2_dnam_mrdoc_nominal_intersect.tsv"), 
            col.names = F, row.names = F, quote = F)



## Subset without MHC
smk2dnam_intersect_noMHC <- mrdoc_currSmk_nominal |> 
  filter( ! (CHR == 6 &
               MAPINFO >= 28477797 &
               MAPINFO <= 33448354) ) |>
  filter(g1_nominal) |> 
  select(cpg) 

write.table(smk2dnam_intersect_noMHC, file = paste0(outDir,"smk_2_dnam_mrdoc_nominal_intersect_noMHC.tsv"), 
            col.names = F, row.names = F, quote = F)


## robust
smk2dnam_robust <- mrdoc_currSmk_Res |>
  filter( (mrdoc2.g1_fdrSig &
             mrdoc1_wPleio.g1_fdrSig &
             mrdoc1_wRE.g1_fdrSig) ) |>
  select(cpg)

write.table(smk2dnam_robust, file = paste0(outDir,"smk_2_dnam_mrdoc_fdr_intersect.tsv"), 
            col.names = F, row.names = F, quote = F)



#### Compare with Sun et al 2021 =============

prior <- c("cg14391737","cg09338374","cg02978227","cg16841366","cg12956751","cg13849276",
           "cg00475490","cg22222502","cg14580211","cg15212295","cg13258799")

mrdoc_currSmk_Res |>
  filter( (mrdoc2.g1_fdrSig &
             mrdoc1_wPleio.g1_fdrSig &
             mrdoc1_wRE.g1_fdrSig) ) |>
  filter(cpg %in% prior)

mrdoc_currSmk_nominal |> 
  filter(g1_nominal) |> 
  filter(cpg %in% prior)

mrdoc_currSmk_Res |>
  filter(cpg %in% prior)


#### Compare with Jamieson et al 2020 =============

prior2 <- c("cg10255761","cg19758448","cg21201401","cg15059804","cg15951188",
            "cg24033122","cg09099830","cg15233611","cg06382664")

mrdoc_currSmk_nominal |>
  filter(cpg %in% prior2) |> 
  select(c(cpg,CHR,MAPINFO,SYMBOL,contains("nominal"))) |> 
  arrange(CHR,MAPINFO)

mrdoc_currSmk_Res |>
  filter( (mrdoc2.g1_fdrSig &
             mrdoc1_wPleio.g1_fdrSig &
             mrdoc1_wRE.g1_fdrSig) ) |> 
  filter(cpg %in% prior2)

mrdoc_currSmk_nominal <- data.frame(mrdoc_currSmk_nominal)
rownames(mrdoc_currSmk_nominal) <- mrdoc_currSmk_nominal$cpg
mrdoc_currSmk_nominal |>
  filter(cpg %in% prior2) |> 
  arrange(CHR,MAPINFO) |> 
  select(cpg,CHR,MAPINFO,SYMBOL,contains("g2")) |> 
  select(-c(contains("SE"),contains("Z"),contains("Sig")))


#### Compare with van Dongen et al 2023 =============

twins <- c("cg05575921","cg21566642","cg05951221","cg01940273","cg13411554",
           "cg01901332","cg21161138","cg00336149","cg22132788","cg21188533",
           "cg09935388","cg25648203","cg19089201")

table( twins %in% mrdoc_currSmk_Res$cpg)
# TRUE 
#   13 

mrdoc_currSmk_Res |>
  filter( (mrdoc2.g1_fdrSig &
             mrdoc1_wPleio.g1_fdrSig &
             mrdoc1_wRE.g1_fdrSig) ) |>
  filter(cpg %in% twins) |> 
  select(cpg,CHR,MAPINFO,SYMBOL) |> 
  arrange(CHR)

mrdoc_currSmk_nominal |> 
  filter(g1_nominal) |> 
  filter(cpg %in% twins)

mrdoc_currSmk_nominal |> 
  filter(g2_nominal) |> 
  filter(cpg %in% twins)
# None


#### Save the Gene names ===============

dnam2smk_gene <- mrdoc_currSmk_nominal |> 
  filter(g2_nominal) |> 
  arrange(mrdoc1_wPleio.g1_p) |> 
  select(SYMBOL) |> 
  t()
str(dnam2smk_gene)

write.table(dnam2smk_gene, file = paste0(outDir,"dnam_2_smk_mrdoc_nominal_intersect_genes.csv"), 
            sep = ",", col.names = F, row.names = F, quote = F)


smk2dnam_gene <- mrdoc_currSmk_nominal |> 
  filter(g1_nominal) |> 
  arrange(mrdoc1_wPleio.g1_p) |> 
  select(SYMBOL) |> 
  t()
str(smk2dnam_gene)

smk2dnam_gene_noMHC <- mrdoc_currSmk_nominal |> 
  filter( ! (CHR == 6 &
               MAPINFO >= 28477797 &
               MAPINFO <= 33448354) ) |>
  filter(g1_nominal) |> 
  arrange(mrdoc1_wPleio.g1_p) |> 
  select(SYMBOL) |> 
  t()
str(smk2dnam_gene_noMHC)
# 525

write.table(smk2dnam_gene_noMHC, file = paste0(outDir,"smk_2_dnam_mrdoc_nominal_intersect_genes_noMHC.csv"), 
            sep = ",", col.names = F, row.names = F, quote = F)


## Plot 

mrdoc_nominal_bidir <- mrdoc_currSmk_nominal |> 
  filter(g1_nominal, g2_nominal) |> 
  mutate(  
    cpg = fct_reorder(cpg, desc(mrdoc2.g1_Z)),
    robust = ifelse(robust_smk2dnam, "Stronger Evidence of\nEffects of Smoking on DNAm",
                    ifelse(robust_dnam2smk, "Stronger\nEvidence of\nEffects of\nDNAm on\nSmoking",
                           "Suggestive Evidence of\nEffects in Both Directions")),
    robust = fct_relevel(robust, "Stronger\nEvidence of\nEffects of\nDNAm on\nSmoking", after = Inf)
  ) |> 
  select(cpg, CHR, MAPINFO, SYMBOL, robust,
         mrdoc1_wPleio.g1_hat, mrdoc1_wPleio.g1_SE, mrdoc1_wPleio.g1_Z, mrdoc1_wPleio.g1_p, 
         mrdoc1_wRE.g1_hat, mrdoc1_wRE.g1_SE, mrdoc1_wRE.g1_Z, mrdoc1_wRE.g1_p, 
         mrdoc2.g1_hat, mrdoc2.g1_SE, mrdoc2.g1_Z, mrdoc2.g1_p,
         mrdoc1_wPleio.g2_hat, mrdoc1_wPleio.g2_SE, mrdoc1_wPleio.g2_Z, mrdoc1_wPleio.g2_p, 
         mrdoc1_wRE.g2_hat, mrdoc1_wRE.g2_SE, mrdoc1_wRE.g2_Z, mrdoc1_wRE.g2_p, 
         mrdoc2.g2_hat, mrdoc2.g2_SE, mrdoc2.g2_Z, mrdoc2.g2_p) |> 
  pivot_longer(cols = -c(cpg, CHR, MAPINFO, SYMBOL, robust), 
               names_to = c("model","par"), names_sep = ".g") |> 
  pivot_wider(id_cols = c("cpg","CHR","MAPINFO","SYMBOL","model","robust"), 
              names_from = "par") |> 
  pivot_longer(cols = -c(cpg, CHR, MAPINFO, SYMBOL, model, robust), 
               names_to = c("doc","par"), names_sep = "_") |> 
  pivot_wider(id_cols = c("cpg","CHR","MAPINFO","SYMBOL","doc","model","robust"), 
              names_from = "par") |> 
  mutate( 
    doc = fct_rev(factor(doc)),
    doc = fct_recode(doc,
                     "Effect of Smoking on DNAm" = "1",
                     "Effect of DNAm on Smoking" = "2"),
    model = fct_relevel(model, "mrdoc2"),
    model = fct_recode(model, 
                       "MR-DoC2" = "mrdoc2",
                       "MR-DoC1 w/ Pleiotropic Path" = "mrdoc1_wPleio",
                       "MR-DoC1 w/ rE" = "mrdoc1_wRE"),
    
  ) 


## Combine and re-split by Direction of Causation 
mrdoc_nominal_bidir_smk <- mrdoc_nominal_bidir |>
  filter(doc == "Effect of Smoking on DNAm") 

mrdoc_nominal_bidir_dnam <- mrdoc_nominal_bidir |>
  filter(doc == "Effect of DNAm on Smoking") 



## Plot Smk --> DNAm
nominal_bidir_smk <- mrdoc_nominal_bidir_smk |>
  ggplot(aes(x = hat, y = cpg)) + 
  geom_point(aes(shape = model, color = model ),
             size = 1.25, 
             position = position_dodge(width = 0.5) ) +
  geom_errorbarh(aes(color = model, 
                     xmax = hat + qnorm(p=0.975)*SE,
                     xmin = hat - qnorm(p=0.975)*SE),
                 height = 0, 
                 linewidth = 0.5,
                 position = position_dodge(width = 0.5)) +
  facet_grid(rows = vars(robust), cols = vars(doc), 
             scales = "free_y", space = "free") +
  geom_vline(xintercept = 0, 
             color = "midnightblue", 
             linetype = "dashed", 
             linewidth = 0.2) +
  coord_cartesian(xlim = c(-1.5, 1.1)) +
  labs(x = "Causal Estimate (95% C.I.)", y = NULL, 
       color = NULL, 
       shape = NULL) +
  scale_color_manual(values = c("steelblue4","slategray","lightskyblue3")) +
  scale_shape_manual(values = c(17, 19, 15)) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme_bw(6) +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text.x = element_text(size = 8),
    strip.text.y = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "pt"),
    legend.position = "bottom" ,
    legend.text = element_text(size = 6)
  ) 


## Plot DNAm --> Smk
nominal_bidir_dnam <- mrdoc_nominal_bidir_dnam |>
  ggplot(aes(x = hat, y = cpg)) + 
  geom_point(aes(shape = model, color = model ),
             size = 1.25, 
             position = position_dodge(width = 0.5) ) +
  geom_errorbarh(aes(color = model, 
                     xmax = hat + qnorm(p=0.975)*SE,
                     xmin = hat - qnorm(p=0.975)*SE),
                 height = 0, 
                 linewidth = 0.5,
                 position = position_dodge(width = 0.5)) +
  facet_grid(rows = vars(robust), cols = vars(doc), 
             scales = "free_y", space = "free") +
  geom_vline(xintercept = 0, 
             color = "midnightblue", 
             linetype = "dashed", 
             linewidth = 0.2) +
  coord_cartesian(xlim = c(-1.5, 1.1)) +
  labs(x = "Causal Estimate (95% C.I.)", y = NULL, 
       color = NULL, 
       shape = NULL) +
  scale_color_manual(values = c("indianred3","thistle4","mistyrose3")) +
  scale_shape_manual(values = c(17, 19, 15)) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme_bw(6) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text.x = element_text(size = 8),
    strip.text.y = element_text(size = 6),
    plot.margin = margin(0, 0, 0, 0, "pt"),
    legend.position = "bottom" ,
    legend.text = element_text(size = 6)
  ) 


## Combine
library(patchwork)
patchBidir <-  nominal_bidir_smk + nominal_bidir_dnam + 
  plot_layout(guides = 'keep')

nBidir <- mrdoc_currSmk_nominal |> 
  filter(g1_nominal, g2_nominal) |> 
  nrow()

jpeg(paste0(plotDir,"nominal_bidir_causal_estimates.jpeg"), 
     width = 8, height = 9, units = "in",
     res = 300)

patchBidir + 
  plot_annotation(
    title = paste("Suggestive Bidirectional Causal Effects between Smoking and DNAm at",
                  nBidir,"CpGs")
  )  & 
  theme(title = element_text(size = 8), 
        legend.position = "bottom")

dev.off()




## Save supplementary annonated tables (nominally significant) ==================


## Add annotation from Illumina manifest
annoDir <- paste0(baseDir,"data/annotation/")
## Load the Annotation Files
load(paste0(annoDir,"manifest.RData"))

manifest <- manifest[mrdoc_currSmk_Res$cpg, 
                     c("IlmnID","CHR","MAPINFO","Genome_Build","UCSC_RefGene_Name",
                       "UCSC_RefGene_Group","UCSC_CpG_Islands_Name","Relation_to_UCSC_CpG_Island")]

## Convert missing values to NA
manifestT <- manifest |> 
  mutate(IlmnID = as.character(IlmnID),
         CHR = as.character(CHR),
         UCSC_RefGene_Name = as.character(UCSC_RefGene_Name),
         UCSC_RefGene_Group = as.character(UCSC_RefGene_Group),
         UCSC_CpG_Islands_Name = as.character(UCSC_CpG_Islands_Name),
         Relation_to_UCSC_CpG_Island = as.character(Relation_to_UCSC_CpG_Island))

manifestT <- manifestT |> 
  mutate(IlmnID = ifelse(IlmnID == "", NA, IlmnID),
         CHR = ifelse(CHR == "", NA, CHR),
         UCSC_RefGene_Name = ifelse(UCSC_RefGene_Name == "", NA, UCSC_RefGene_Name),
         UCSC_RefGene_Group = ifelse(UCSC_RefGene_Group == "", NA, UCSC_RefGene_Group),
         UCSC_CpG_Islands_Name = ifelse(UCSC_CpG_Islands_Name == "", NA, UCSC_CpG_Islands_Name),
         Relation_to_UCSC_CpG_Island = ifelse(Relation_to_UCSC_CpG_Island == "", NA, Relation_to_UCSC_CpG_Island) )


load(paste0(annoDir,"cpgInfo_08112016.RData"))  ## Info on the nearest gene, if the CG is intergenic
info$IlmnID <- rownames(info)
info <- info |> 
  rename(NearestGene = SYMBOL)

nearestgene <- info |> 
  filter(IlmnID %in% manifestT$IlmnID) |> 
  select(IlmnID,NearestGene)

str(nearestgene)

mrdoc_currSmk_nominal <- mrdoc_currSmk_nominal |> 
  rename(IlmnID = cpg) |> 
  select(-c(CHR,MAPINFO,SYMBOL)) |> 
  left_join(manifestT, by = "IlmnID") |> 
  left_join(nearestgene, by = "IlmnID")

# rm(manifest,nearestgene); gc()


## Add indicator variable for MHC region
mrdoc_currSmk_nominal <- mrdoc_currSmk_nominal |> 
  mutate(MHC = ifelse( (CHR == 6 &
                          MAPINFO >= 28477797 &
                          MAPINFO <= 33448354),
                       TRUE, FALSE)
  )

mrdoc_currSmk_nominal |> 
  count(MHC)

## Make subsets
mrdoc_smk2DNAm_nominal_g1 <- mrdoc_currSmk_nominal |> 
  filter(g1_nominal) |> 
  rename(g1_robust = robust_smk2dnam,
         g2_robust = robust_dnam2smk) |> 
  select(-c(contains("g2"),g1_nominal)) |> 
  as.data.frame() |> 
  relocate(starts_with("mrdoc"), .after = last_col()) |> 
  relocate(contains("Sig"), .after = last_col()) |> 
  relocate(contains("robust"), .after = last_col())

mrdoc_smk2DNAm_nominal_g2 <- mrdoc_currSmk_nominal |> 
  filter(g1_nominal) |> 
  rename(g1_robust = robust_smk2dnam,
         g2_robust = robust_dnam2smk) |> 
  select(-c(contains("g1"))) |> 
  as.data.frame() |> 
  relocate(starts_with("mrdoc"), .after = last_col()) |> 
  relocate(contains("Sig"), .after = last_col()) |> 
  relocate(contains("nominal"), .after = last_col()) |> 
  relocate(contains("robust"), .after = last_col())

mrdoc_DNAm2smk_nominal_g1 <- mrdoc_currSmk_nominal |> 
  filter(g2_nominal) |> 
  rename(g1_robust = robust_smk2dnam,
         g2_robust = robust_dnam2smk) |> 
  select(-c(contains("g2"))) |> 
  as.data.frame() |> 
  relocate(starts_with("mrdoc"), .after = last_col()) |> 
  relocate(contains("Sig"), .after = last_col()) |> 
  relocate(contains("nominal"), .after = last_col()) |> 
  relocate(contains("robust"), .after = last_col())

mrdoc_DNAm2smk_nominal_g2 <- mrdoc_currSmk_nominal |> 
  filter(g2_nominal) |> 
  rename(g1_robust = robust_smk2dnam,
         g2_robust = robust_dnam2smk) |> 
  select(-c(contains("g1"),g2_nominal)) |> 
  as.data.frame() |> 
  relocate(starts_with("mrdoc"), .after = last_col()) |> 
  relocate(contains("Sig"), .after = last_col()) |> 
  relocate(contains("nominal"), .after = last_col()) |> 
  relocate(contains("robust"), .after = last_col())

str(mrdoc_smk2DNAm_nominal_g1)
str(mrdoc_smk2DNAm_nominal_g2)
str(mrdoc_DNAm2smk_nominal_g1)
str(mrdoc_DNAm2smk_nominal_g2)


## Save in excel.
# https://ycphs.github.io/openxlsx/articles/Introduction.html
outExcel <- list(
  "Current_Smk_to_DNAm_g1" = mrdoc_smk2DNAm_nominal_g1, 
  "Current_Smk_to_DNAm_Reverse_g2" = mrdoc_smk2DNAm_nominal_g2,
  "DNAm_to_Current_Smk_g2" = mrdoc_DNAm2smk_nominal_g2, 
  "DNAm_to_Current_Smk_Reverse_g1" = mrdoc_DNAm2smk_nominal_g1
)
openxlsx::write.xlsx(outExcel, 
                     file = paste0(outDir,"suppl_mrdoc_smk_dnam.xlsx"))

## Individual tables
write.table(mrdoc_smk2DNAm_nominal_g1, 
            file = paste0(outDir,"suppl_smk2dnam_g1est.tsv"), quote = F, row.names = F, col.names = T)
write.table(mrdoc_smk2DNAm_nominal_g2, 
            file = paste0(outDir,"suppl_smk2dnam_g2est.tsv"), quote = F, row.names = F, col.names = T)
write.table(mrdoc_DNAm2smk_nominal_g2, 
            file = paste0(outDir,"suppl_dnam2smk_g2est.tsv"), quote = F, row.names = F, col.names = T)
write.table(mrdoc_DNAm2smk_nominal_g1, 
            file = paste0(outDir,"suppl_dnam2smk_g1est.tsv"), quote = F, row.names = F, col.names = T)


## Fixer's Exact Test: Relationship with MHC ==========================

## nominal g1
mhc_by_nominal <- table(mrdoc_currSmk_nominal$MHC,mrdoc_currSmk_nominal$g1_nominal)
fisher.test(mhc_by_nominal)


## robust g1
mhc_by_robust <- table(mrdoc_currSmk_nominal$MHC,mrdoc_currSmk_nominal$robust_smk2dnam)
mhc_by_robust
#       FALSE  TRUE
# FALSE 16278    63
# TRUE    564     1

fisher.test(mhc_by_robust)


## Overlap with cell count GWAS signal. Vuckovic et al 2020 ==========================

cellCountTab1 <- readxl::read_xlsx(paste0(outDir,"VuckovicEtAl_TableS3.xlsx"))
cellCountTab1 <- cellCountTab1 |> 
  janitor::clean_names() 

cellCountTab2 <- readxl::read_xlsx(paste0(outDir,"VuckovicEtAl_TableS4.xlsx"))
cellCountTab2 <- cellCountTab2 |> 
  janitor::clean_names() 


cellCountTab1 <- cellCountTab1 |> 
  select(associated_blood_index,associated_blood_index_class,gene_symbol_s_for_most_serious_consequence,rs_id_where_available) |> 
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

mrdoc_currSmk_nominal |> 
  filter(g1_nominal) |> 
  count(SYMBOL %in% cellCountGWAS$gene)

mrdoc_currSmk_nominal |> 
  filter(g2_nominal) |> 
  count(SYMBOL %in% cellCountGWAS$gene)

mrdoc_currSmk_nominal <- mrdoc_currSmk_nominal |> 
  mutate(cellCountHit = ifelse(SYMBOL %in% cellCountGWAS$gene, TRUE, FALSE),
         lymphCountHit = ifelse(SYMBOL %in% lymphCountGWAS$gene, TRUE, FALSE))


mrdoc_currSmk_nominal |> 
  distinct(SYMBOL, .keep_all = T) |> 
  count(cellCountHit)

mrdoc_currSmk_nominal |> 
  distinct(SYMBOL, .keep_all = T) |> 
  count(lymphCountHit)

cellCount_by_g1 <- mrdoc_currSmk_nominal |> 
  distinct(SYMBOL, .keep_all = T) |> 
  filter(!is.na(g1_nominal)) 

cellCount_by_g1 |> 
  filter(g1_nominal) |> 
  nrow()


cellCount_by_g1 |> 
  count(g1_nominal,cellCountHit)

fisher.test(table(cellCount_by_g1$g1_nominal,cellCount_by_g1$cellCountHit), alternative = "g")

cellCount_by_g1 |> 
  count(g1_nominal,lymphCountHit)

cellCount_by_g2 <- mrdoc_currSmk_nominal |> 
  distinct(SYMBOL, .keep_all = T) |> 
  filter(!is.na(g2_nominal))

cellCount_by_g2 |> 
  filter(g2_nominal) |> 
  nrow()

cellCount_by_g2 |> 
  count(g2_nominal,cellCountHit)

fisher.test(table(cellCount_by_g2$g2_nominal,cellCount_by_g2$cellCountHit), alternative = "g")

cellCount_by_g2 |> 
  count(g2_nominal,lymphCountHit)

fisher.test(table(cellCount_by_g2$g2_nominal,cellCount_by_g2$lymphCountHit), alternative = "g")



## UpSet Plots #####################################

### Smk --> DNAm FDR UpSet Plot ######################

mrdoc_currSmk_DNAm_FDR <- mrdoc_currSmk_Res |> 
  select( cpg, CHR, SYMBOL, 
          mrdoc2.g1_fdrSig, mrdoc1_wPleio.g1_fdrSig, mrdoc1_wRE.g1_fdrSig ) |> 
  filter( mrdoc2.g1_fdrSig == T | mrdoc1_wPleio.g1_fdrSig == T | mrdoc1_wRE.g1_fdrSig == T  )


mrdoc_currSmk_DNAm_FDR <- mrdoc_currSmk_DNAm_FDR |> 
  rename( 
    "MR-DoC2" = 'mrdoc2.g1_fdrSig', 
    "MR-DoC1_wPleio" = "mrdoc1_wPleio.g1_fdrSig", 
    "MR-DoC1_wRE" = "mrdoc1_wRE.g1_fdrSig" 
  )

models <- colnames(mrdoc_currSmk_DNAm_FDR)[-(1:3)]


#### Plot

smk2dnam_upset <- upset(
  ## Data with binary vars with group membership
  data = mrdoc_currSmk_DNAm_FDR, 
  ## Columns with group membership
  intersect = models,  
  ## the label shown below the intersection matrix
  name="MR-DoC Model", 
  min_size = 0,
  width_ratio = 0.125 ,
  ## sort
  # sort_intersections_by = c('degree', 'cardinality'),
  # sort_intersections = 'descending',
  sort_intersections=FALSE,
  intersections=list(
    c( "MR-DoC1_wPleio","MR-DoC1_wRE","MR-DoC2" ),
    c( "MR-DoC1_wPleio","MR-DoC2" ),
    c( "MR-DoC1_wRE","MR-DoC1_wPleio" ),
    c( "MR-DoC1_wRE","MR-DoC2" ),
    "MR-DoC2",
    "MR-DoC1_wRE",
    "MR-DoC1_wPleio"
  ),
  ## display set size
  set_sizes=(
    upset_set_size()
    + geom_text(aes(label=after_stat(count)), 
                size=4, hjust=1.1, stat='count')
    # you can also add annotations on top of bars:
    + expand_limits(y=2600)
    + theme(axis.text.x = element_blank(),
            axis.title.x = element_text(size = 8))
  ), 
  queries=list(
    ## Color of the set size
    upset_query(
      set = "MR-DoC2",  fill = 'steelblue3'
    ),
    upset_query(
      set = "MR-DoC1_wPleio", fill = 'steelblue3'
    ),
    upset_query(
      set = "MR-DoC1_wRE", fill = 'steelblue3'
    ),
    ## Color of intersection of g1 
    upset_query(
      intersect = c('MR-DoC1_wPleio','MR-DoC1_wRE','MR-DoC2'), 
      color = 'steelblue3', 
      fill='steelblue3',
      only_components = c('intersections_matrix', 'Intersection size')
    )
  )
) +
  theme_minimal(9) + 
  labs(title = "Putative Causal Effects of Current Smoking on DNA Methylation\nFDR < 0.05",
       y = "Model",
       x = "Intersection of Significant Causal Estimates Across Models") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10))



jpeg(paste0(plotDir,"mrdoc_unidir_currentSmk_dnam_UpSet.jpeg"), width = 8, height = 5, units = "in", res = 300)
smk2dnam_upset
dev.off()


### DNAm --> Smk FDR UpSet Plot ######################

mrdoc_DNAm_currSmk_FDR <- mrdoc_currSmk_Res |> 
  select( cpg, CHR, SYMBOL, 
          mrdoc2.g2_fdrSig, mrdoc1_wPleio.g2_fdrSig, mrdoc1_wRE.g2_fdrSig ) |> 
  filter( mrdoc2.g2_fdrSig == T | mrdoc1_wPleio.g2_fdrSig == T | mrdoc1_wRE.g2_fdrSig == T  )


mrdoc_DNAm_currSmk_FDR <- mrdoc_DNAm_currSmk_FDR |> 
  rename( 
    "MR-DoC2" = 'mrdoc2.g2_fdrSig', 
    "MR-DoC1_wPleio" = "mrdoc1_wPleio.g2_fdrSig", 
    "MR-DoC1_wRE" = "mrdoc1_wRE.g2_fdrSig" 
  )

models <- colnames(mrdoc_DNAm_currSmk_FDR)[-(1:3)]


#### Plot

dnam2smk_upset <- upset(
  ## Data with binary vars with group membership
  data = mrdoc_DNAm_currSmk_FDR, 
  ## Columns with group membership
  intersect = models,  
  ## the label shown below the intersection matrix
  name="MR-DoC Model", 
  min_size = 0,
  width_ratio = 0.125 ,
  ## sort
  # sort_intersections_by = c('degree', 'cardinality'),
  # sort_intersections = 'descending',
  sort_intersections=FALSE,
  intersections=list(
    c( "MR-DoC1_wPleio","MR-DoC1_wRE","MR-DoC2" ),
    c( "MR-DoC1_wPleio","MR-DoC2" ),
    c( "MR-DoC1_wRE","MR-DoC1_wPleio" ),
    c( "MR-DoC1_wRE","MR-DoC2" ),
    "MR-DoC1_wRE",
    "MR-DoC2",
    "MR-DoC1_wPleio"
  ),
  ## display set size
  set_sizes=(
    upset_set_size()
    + geom_text(aes(label=after_stat(count)), hjust=1.1, stat='count')
    # you can also add annotations on top of bars:
    + expand_limits(y=2700)
    + theme(axis.text.x = element_blank(),
            axis.title.x = element_text(size = 8))
  ), 
  queries=list(
    ## Color of the set size
    upset_query(
      set = "MR-DoC2",  fill = 'indianred3'
    ),
    upset_query(
      set = "MR-DoC1_wPleio", fill = 'indianred3'
    ),
    upset_query(
      set = "MR-DoC1_wRE", fill = 'indianred3'
    ),
    ## Color of intersection of g1 
    upset_query(
      intersect = c('MR-DoC1_wPleio','MR-DoC1_wRE','MR-DoC2'), color = 'indianred3', fill='indianred3',
      only_components = c('intersections_matrix', 'Intersection size')
    )
  )
) + 
  theme_minimal(9) + 
  labs(title = "Putative Causal Effects of DNA Methylation on Current Smoking\nFDR < 0.05",
       y = "Model",
       x = "Intersection of Significant Causal Estimates Across Models") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10))



jpeg(paste0(plotDir,"mrdoc_unidir_dnam_currentSmk_UpSet.jpeg"), width = 8, height = 5, units = "in", res = 300)

dnam2smk_upset

dev.off()




### Smk --> DNAm Bonferroni UpSet Plot ######################


mrdoc_currSmk_DNAm_Bonf <- mrdoc_currSmk_Res |> 
  select( cpg, CHR, SYMBOL, 
          mrdoc2.g1_bonfSig, mrdoc1_wPleio.g1_bonfSig, mrdoc1_wRE.g1_bonfSig ) |> 
  filter( mrdoc2.g1_bonfSig == T   | mrdoc1_wPleio.g1_bonfSig == T | mrdoc1_wRE.g1_bonfSig == T )



mrdoc_currSmk_DNAm_Bonf <- mrdoc_currSmk_DNAm_Bonf |> 
  rename( 
    "MR-DoC2"        = 'mrdoc2.g1_bonfSig', 
    "MR-DoC1_wPleio" = "mrdoc1_wPleio.g1_bonfSig", 
    "MR-DoC1_wRE"    = "mrdoc1_wRE.g1_bonfSig"
  )

models <- colnames(mrdoc_currSmk_DNAm_Bonf)[-(1:3)]

#### Plot

jpeg(paste0(plotDir,"mrdoc_unidir_currentSmk_dnam_Bonf_UpSet.jpeg"), width = 9, height = 6, units = "in", res = 300)
upset(
  ## Data with binary vars with group membership
  data = mrdoc_currSmk_DNAm_Bonf, 
  ## Columns with group membership
  intersect = models,  
  ## the label shown below the intersection matrix
  name="MR-DoC Model", 
  min_size = 0,
  width_ratio = 0.125 ,
  ## sort
  # sort_intersections_by = c('degree', 'cardinality'),
  # sort_intersections = 'descending',
  sort_intersections=FALSE,
  intersections=list(
    c( "MR-DoC1_wPleio","MR-DoC1_wRE","MR-DoC2" ),
    c( "MR-DoC1_wRE","MR-DoC2" ),
    c( "MR-DoC1_wPleio","MR-DoC2" ),
    c( "MR-DoC1_wRE","MR-DoC1_wPleio" ),
    "MR-DoC2",
    "MR-DoC1_wRE",
    "MR-DoC1_wPleio"
  ),
  
  ## display set size
  set_sizes=(
    upset_set_size()
    + geom_text(aes(label=after_stat(count)), hjust=1.1, stat='count')
    # you can also add annotations on top of bars:
    + expand_limits(y=160)
    + theme(axis.text.x=element_text(angle=90))
  ), 
  queries=list(
    ## Color of the set size
    upset_query(
      set = "MR-DoC2",  fill = 'steelblue3'
    ),
    upset_query(
      set = "MR-DoC1_wPleio", fill = 'steelblue3'
    ),
    upset_query(
      set = "MR-DoC1_wRE", fill = 'steelblue3'
    ),
    ## Color of intersection
    upset_query(
      intersect = c('MR-DoC2','MR-DoC1_wRE','MR-DoC1_wPleio'), color = 'steelblue3', fill='steelblue3',
      only_components = c('intersections_matrix', 'Intersection size')
    ) 
  )
) +
  ggtitle("Putative Causal Effects of Current Smoking on DNA Methylation\nStatistically Significant after Bonferroni Correction") 

dev.off()



### DNAm --> Smk Bonferroni UpSet Plot ######################

mrdoc_DNAm_currSmk_Bonf <- mrdoc_currSmk_Res |> 
  select( cpg, CHR, SYMBOL, 
          mrdoc2.g2_bonfSig, mrdoc1_wPleio.g2_bonfSig, mrdoc1_wRE.g2_bonfSig ) |> 
  filter( mrdoc2.g2_bonfSig == T   | mrdoc1_wPleio.g2_bonfSig == T | mrdoc1_wRE.g2_bonfSig == T )



mrdoc_DNAm_currSmk_Bonf <- mrdoc_DNAm_currSmk_Bonf |> 
  rename( 
    "MR-DoC2"        = 'mrdoc2.g2_bonfSig', 
    "MR-DoC1_wPleio" = "mrdoc1_wPleio.g2_bonfSig", 
    "MR-DoC1_wRE"    = "mrdoc1_wRE.g2_bonfSig"
  )

models <- colnames(mrdoc_DNAm_currSmk_Bonf)[-(1:3)]

#### Plot

jpeg(paste0(plotDir,"mrdoc_unidir_dnam_currentSmk_Bonf_UpSet.jpeg"), width = 9, height = 6, units = "in", res = 300)

upset(
  ## Data with binary vars with group membership
  data = mrdoc_DNAm_currSmk_Bonf, 
  ## Columns with group membership
  intersect = models,  
  ## the label shown below the intersection matrix
  name="MR-DoC Model", 
  min_size = 0,
  width_ratio = 0.125 ,
  ## sort
  # sort_intersections_by = c('degree', 'cardinality'),
  # sort_intersections = 'descending',
  sort_intersections=FALSE,
  intersections=list(
    c( "MR-DoC1_wRE","MR-DoC2" ),
    "MR-DoC2",
    "MR-DoC1_wRE",
    "MR-DoC1_wPleio"
  ),
  ## display set size
  set_sizes=(
    upset_set_size()
    + geom_text(aes(label=after_stat(count)), hjust=1.1, stat='count')
    # you can also add annotations on top of bars:
    + expand_limits(y=100)
    + theme(axis.text.x = element_blank(),
            axis.title.x = element_text(size = 8))
  ), 
  queries=list(
    ## Color of the set size
    upset_query(
      set = "MR-DoC2",  fill = 'indianred3'
    ),
    upset_query(
      set = "MR-DoC1_wPleio", fill = 'indianred3'
    ),
    upset_query(
      set = "MR-DoC1_wRE", fill = 'indianred3'
    ),
    ## Color of intersection of g1 
    upset_query(
      intersect = c('MR-DoC1_wRE','MR-DoC2'), color = 'indianred3', fill='indianred3',
      only_components = c('intersections_matrix', 'Intersection size')
    )
  )
) + 
  theme_minimal(9) + 
  labs(title = "Putative Causal Effects of DNA Methylation on Current Smoking\nStatistically Significant after Bonferroni Correction",
       y = "Model",
       x = "Intersection of Significant Causal Estimates Across Models") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10))


dev.off()


