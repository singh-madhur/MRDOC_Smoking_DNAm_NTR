## Combine the Results for MRDoC Former Smk <--> DNAm
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
baseDir     <- "/data/msingh/MRDoC_Smk_DNAm/" 
smk2dnamDir <- paste0(baseDir,"out/mrdoc_formerSmk_dnam/summary/")
dnam2smkDir <- paste0(baseDir,"out/mrdoc_dnam_formerSmk/summary/")
mrdoc2Dir   <- paste0(baseDir,"out/mrdoc2_formerSmk_dnam/summary/")
outDir      <- paste0(baseDir,"out/mrdoc_summary_formerSmk_dnam/")
plotDir     <- paste0(outDir,"plots/")


## Load results #########

## https://stackoverflow.com/a/25455968/18446763
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

## MRDoC1 Smk --> DNAm
smk2dnamRes <- loadRData( paste0(smk2dnamDir,"mrdoc_formerSmk_dnam_results.RData") )

## MRDoC1 DNAm --> Smk
dnam2smkRes <- loadRData( paste0(dnam2smkDir,"mrdoc_dnam_formerSmk_dnam.RData") )

## MRDoC2
mrdoc2Res   <- loadRData( paste0(mrdoc2Dir,"mrdoc2_formerSmk_dnam_results.RData") )


## Rename Columns ##################

#### SMK --> DNAm ####################

dim(smk2dnamRes)
bonfp_smk2dnam <- 0.05 / nrow(smk2dnamRes)

smk2dnamRes <- smk2dnamRes |> 
  mutate( 
    mrdoc1_wPleio.g1_bonfSig = ifelse( g1_p     < bonfp_smk2dnam, TRUE, FALSE),
    mrdoc1_wRE.g1_bonfSig    = ifelse( g1_wRE_p < bonfp_smk2dnam, TRUE, FALSE),
    mrdoc1_wPleio.g1_fdrSig  = ifelse( g1_q     < 0.05,  TRUE, FALSE),
    mrdoc1_wRE.g1_fdrSig     = ifelse( g1_wRE_q < 0.05,  TRUE, FALSE)
  ) |> 
  rename( mrdoc1_wPleio.g1_hat     = "g1_hat",
          mrdoc1_wPleio.g1_SE      = "g1_SE",
          mrdoc1_wPleio.g1_Z       = "g1_Z",
          mrdoc1_wPleio.g1_p       = "g1_p",
          mrdoc1_wPleio.g1_qval    = "g1_q",
          mrdoc1_wRE.g1_hat        = "g1_wRE_hat",
          mrdoc1_wRE.g1_SE         = "g1_wRE_SE",
          mrdoc1_wRE.g1_Z          = "g1_wRE_Z",
          mrdoc1_wRE.g1_p          = "g1_wRE_p",
          mrdoc1_wRE.g1_qval       = "g1_wRE_q" )

smk2dnamRes <- smk2dnamRes |> 
  mutate(CHR = as.character(CHR),
         CHR = as.numeric(CHR))


#### DNAm --> Smk ####################

dim(dnam2smkRes)
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

dim(mrdoc2Res)
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

nrow(smk2dnamRes) # Former-Smk-related CpGs post-QC
nrow(dnam2smkRes) # Former-Smk-related CpGs post-QC AND mQTLs with Fstat > 10
nrow(mrdoc2Res)   # Former-Smk-related CpGs post-QC AND mQTLs with Fstat > 10

table(is.element(smk2dnamRes$cpg, dnam2smkRes$cpg ))
prop.table( table(is.element(smk2dnamRes$cpg, dnam2smkRes$cpg )) )

mrdoc_formerSmk_Res <- mrdoc2Res |> 
  full_join(smk2dnamRes) |> 
  full_join( dnam2smkRes ) 

dim(mrdoc_formerSmk_Res)


## Save ##########################
save(mrdoc_formerSmk_Res, file = paste0(outDir,"mrdoc_formerSmk_DNAm_results.RData"))


# Estimates of Top Sites #####################################


### Smk --> DNAm #####################

nTest2 <- mrdoc_formerSmk_Res |> 
  filter(!is.na(mrdoc2.g1_bonfSig)) |> 
  nrow()
bonfp2 <- 0.05/nTest2
bonfZ2 <- qnorm(p = 1-bonfp2/2)

nTest1 <- mrdoc_formerSmk_Res |> 
  filter(!is.na(mrdoc1_wPleio.g1_bonfSig)) |> 
  nrow()
bonfp1 <- 0.05/nTest1
bonfZ1 <- qnorm(p = 1-bonfp1/2)


#### FDR Intersection Estimates ####################

nFDR <- mrdoc_formerSmk_Res |>
  filter( mrdoc2.g1_fdrSig &
            mrdoc1_wPleio.g1_fdrSig &
            mrdoc1_wRE.g1_fdrSig ) |> 
  nrow()

mrdoc_formerSmk_Res |>
  filter( mrdoc2.g1_fdrSig &
            mrdoc1_wPleio.g1_fdrSig &
            mrdoc1_wRE.g1_fdrSig ) |> 
  select(cpg,SYMBOL,mrdoc2.g1_Z,mrdoc1_wPleio.g1_Z,mrdoc1_wRE.g1_Z) |> 
  arrange(mrdoc2.g1_Z)

## How many are in the MHC region?
# chr6:28,477,797-33,448,354
# https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37 

mrdoc_formerSmk_Res |>
  filter( (mrdoc2.g1_fdrSig &
             mrdoc1_wPleio.g1_fdrSig &
             mrdoc1_wRE.g1_fdrSig) ,
          CHR == 6,
          MAPINFO >= 28477797,
          MAPINFO <= 33448354) |>
  select(cpg,CHR,MAPINFO,SYMBOL)

## Plot the estimates
mrdoc_formerSmk_DNAm_FDRneg <- mrdoc_formerSmk_Res |> 
  filter( (mrdoc2.g1_fdrSig &
             mrdoc1_wPleio.g1_fdrSig &
             mrdoc1_wRE.g1_fdrSig)) |> 
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

nrow(mrdoc_formerSmk_DNAm_FDRneg)

mrdoc_formerSmk_DNAm_FDR_plot <- mrdoc_formerSmk_DNAm_FDRneg |> 
  ggplot(aes(g1_hat, cpg, label = plotLabel)) + 
  geom_point(aes(color = model),
             size = 2.5, position = position_dodge(width = 0.5) ) +
  geom_errorbarh(aes(color = model,
                     xmax = g1_hat + qnorm(p=0.975)*g1_SE,
                     xmin = g1_hat - qnorm(p=0.975)*g1_SE),
                 height = 0, linewidth = 1,
                 position = position_dodge(width = 0.5)) +
  geom_text_repel(seed = 1, mapping = aes(fontface = "bold"), 
                  segment.size = 0, na.rm = T, 
                  direction = "x", #position = position_nudge_repel(x=-0.9), 
                  hjust = "left" , xlim = c(NA,-1.75)
  ) +
  geom_vline(xintercept = 0, color = "grey50", linetype = "dashed") +
  coord_cartesian(xlim = c(-2.1,0)) +
  labs( x = "Causal Estimate (95% C.I.)", y = NULL, color = NULL, fill = NULL,
        title = "CpG Sites where Former Smoking Likely Causes Hypomethylation",
        caption = paste0(nFDR," CpG Sites with FDR < 0.05 in All Three Models")) +
  scale_color_manual(values = c("steelblue4","slategray","lightskyblue3")) +
  scale_fill_manual(values = c("steelblue4","slategray","lightskyblue3")) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold",
                                  size = 10),
        axis.text.y = element_text(face = "bold"),
        # axis.title = element_text(size = 9),
        # legend.text = element_text(size =8),
        legend.position = "bottom" ) 


jpeg(paste0(plotDir,"mrdoc_formerSmk_dnam_estimates_suppl.jpeg"), 
     width = 6, height = 6, units = "in", res = 300)
mrdoc_formerSmk_DNAm_FDR_plot
dev.off()


#### FDR Intersection Z-score MR-DoC_w_rE ####################

## Plot the estimates from MR-DoC2
nTests <- nrow(mrdoc_formerSmk_Res)
bonfp  <- 0.05/nTests
bonfZ  <- qnorm(p = 1-bonfp/2)

## Label = Bonferroni Significant in MR-DoC2
mrdoc_formerSmk_DNAm_FDR_Est <- mrdoc_formerSmk_Res |>
  mutate(
    plotLabel = ifelse( mrdoc2.g1_bonfSig  ,
                        SYMBOL, NA)
  )

mrdoc_formerSmk_DNAm_FDR_Est |>
  filter(!is.na(plotLabel)) |> 
  count(plotLabel) |> 
  arrange(-n) 

## Factor for alpha = if FDR < 0.05 in all 3 models
mrdoc_formerSmk_DNAm_FDR_Est <- mrdoc_formerSmk_DNAm_FDR_Est |> 
  mutate(alfa = ifelse( (mrdoc2.g1_fdrSig & 
                           mrdoc1_wPleio.g1_fdrSig &
                           mrdoc1_wRE.g1_fdrSig ) , 
                        "sig", "nonSig") ) 

## Some alfa values are assigned NA (due to missing MR-DoC2 results)
## recode these as nonSig

mrdoc_formerSmk_DNAm_FDR_Est <- mrdoc_formerSmk_DNAm_FDR_Est |> 
  mutate(alfa = ifelse(is.na(alfa), "nonSig", alfa ) ) 

nFDR <- mrdoc_formerSmk_DNAm_FDR_Est |> 
  filter(alfa == "sig") |> 
  nrow()
nFDR

nBonf <- mrdoc_formerSmk_DNAm_FDR_Est |> 
  filter( mrdoc2.g1_bonfSig ) |>
  nrow()
nBonf

## Hard-coding "Manhattan" plots with Z scores
## Based on https://r-graph-gallery.com/101_Manhattan_plot.html

## CHR as numeric 
class(mrdoc_formerSmk_DNAm_FDR_Est$CHR)

## Compute the cumulative position of CpG
mrdocOutPos <-  mrdoc_formerSmk_DNAm_FDR_Est |>
  filter(!is.na(mrdoc1_wRE.g1_Z)) |> 
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(MAPINFO,na.rm = T)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(-chr_len) %>% 
  
  # Add this info to the initial dataset
  left_join(mrdoc_formerSmk_DNAm_FDR_Est, ., by=c("CHR"="CHR")) %>% 
  
  # Add a cumulative position of each SNP
  arrange(MAPINFO) %>%
  mutate( BPcum = MAPINFO + tot)

## X-axis
axisdf <- mrdocOutPos %>% 
  group_by(CHR) %>% 
  summarize(center=( max(BPcum, na.rm = T) + min(BPcum, na.rm = T) ) / 2 )

## Gene Labels = drop repeated plotlabel
mrdocOutPos |> 
  group_by(plotLabel) |> 
  arrange( mrdoc1_wRE.g1_p) |> 
  mutate(nCPG = row_number(), 
         geneLabel = ifelse(nCPG < 2, plotLabel, NA)) |> 
  ungroup() |> 
  filter(!is.na(plotLabel)) |>
  count(CHR, MAPINFO, plotLabel, geneLabel,  mrdoc1_wRE.g1_p,  mrdoc1_wRE.g1_Z, nCPG) 

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
                   shape = alfa,
                   alpha = alfa ,
                   size = alfa ) ) +
  geom_label_repel(data = subset(mrdocOutPos,  mrdoc1_wRE.g1_Z < 0),
                   seed = 1, 
                   size = 3,
                   min.segment.length = 0, na.rm = T, position = position_nudge_repel(y=-1),
                   ylim = c(NA, -6)) +
  geom_label_repel(data = subset(mrdocOutPos,  mrdoc1_wRE.g1_Z > 0),
                   seed = 1, 
                   size = 3,
                   min.segment.length = 0, na.rm = T, position = position_nudge_repel(y=1) ) +
  scale_color_manual(values = rep(c( "steelblue3", "grey50"), 22 )) +
  scale_size_manual(values = c(1.25, 2.5)) +
  scale_alpha_manual(values = c(0.5, 1)) +
  scale_shape_manual(values = c(1, 19)) +
  
  # custom X axis:
  scale_x_continuous( expand = c(0.03, 0),
                      label = axisdf$CHR, breaks= axisdf$center, name = "Chromosome" ) +
  scale_y_continuous( name = "Causal Effect (Z Score)" ) +
  coord_cartesian(clip = "off") +
  geom_hline(yintercept = 0, color = "grey30") +
  geom_hline(yintercept = bonfZ, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -bonfZ, linetype = "dashed", color = "grey50") +
  labs(title = "Putative Causal Effects of Former Smoking on DNAm",
       subtitle = "Causal Estimates from MR-DoC1 Model with rE",
       caption = paste0(nFDR," sites significant at FDR < 0.05 in all three models [Solid Points]\nThe ",
                        nBonf," sites were also significant after Bonferroni correction in all three models")) +
  # Custom the theme:
  theme_bw(8) +
  theme( plot.title = element_text(size = 10),
         plot.caption = element_text(size = 8),
         legend.position="none" , 
         # text = element_text(face = "bold"),
         panel.grid.minor.x = element_blank(),
         panel.border = element_blank()
  )

jpeg(paste0(plotDir,"mrdoc_formerSmk_dnam_wRE_fdr_estimates.jpeg"), 
     width = 9, height = 6, units = "in", res = 300)
mrdocPlot
dev.off()




## DNAm --> Smk #####################

nTest2 <- mrdoc_formerSmk_Res |> 
  filter(!is.na(mrdoc2.g2_bonfSig)) |> 
  nrow()
bonfp2 <- 0.05/nTest2
bonfZ2 <- qnorm(p = 1-bonfp2/2)

nTest1 <- mrdoc_formerSmk_Res |> 
  filter(!is.na(mrdoc1_wPleio.g2_bonfSig)) |> 
  nrow()
bonfp1 <- 0.05/nTest1
bonfZ1 <- qnorm(p = 1-bonfp1/2)


#### FDR Intersection  ####################

mrdoc_formerSmk_Res |> 
  filter( mrdoc2.g2_fdrSig &
            mrdoc1_wPleio.g2_fdrSig & 
            mrdoc1_wRE.g2_fdrSig ) |> 
  nrow()
## 0


## Nominal Intersection Estimate ===============================================

mrdoc_DNAm_formerSmk_nominal_Est <- mrdoc_formerSmk_Res |> 
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
                            "MR-DoC1 w/ Pleiotropy" = "mrdoc1_wPleio",
                            "MR-DoC1 w/ rE" = "mrdoc1_wRE") ) |> 
  group_by(cpg) |> 
  arrange(desc(abs(g2_hat))) |> 
  mutate(top = row_number(),
         plotLabel = ifelse(top == 1, SYMBOL, NA) ) |> 
  ungroup() |> 
  arrange(cpg) 

mrdoc_DNAm_formerSmk_nominal_Est |> 
  count(cpg) |> 
  nrow()

## Split into two parts: positive and negative
mrdoc_DNAm_formerSmk_nominal_pos <- mrdoc_DNAm_formerSmk_nominal_Est |> 
  filter(g2_hat > 0)
nrow(mrdoc_DNAm_formerSmk_nominal_pos)

mrdoc_DNAm_formerSmk_nominal_neg <- mrdoc_DNAm_formerSmk_nominal_Est |> 
  filter(g2_hat < 0)
nrow(mrdoc_DNAm_formerSmk_nominal_neg)

## Plot
jpeg(paste0(plotDir,"mrdoc_dnam_formerSmk_nominal_estimates.jpeg"), 
     width = 8, height = 6, units = "in", res = 300)

mrdoc_DNAm_formerSmk_nominal_Est |> 
  mutate(cpg = paste(cpg,SYMBOL, sep = "\n")) |> 
  ggplot(aes(g2_hat, cpg )) + 
  geom_point(aes(color = model, shape = model),
             size = 2.5, position = position_dodge(width = 0.5) ) +
  geom_errorbarh(aes(color = model, 
                     xmax = g2_hat + qnorm(p=0.975)*g2_SE,
                     xmin = g2_hat - qnorm(p=0.975)*g2_SE),
                 height = 0, linewidth = 1,
                 position = position_dodge(width = 0.5)) +
geom_vline(xintercept = 0, color = "grey50", linetype = "dashed") +
  labs(title = "Putative Causal Effects of DNAm on Former Smoking",
    subtitle = paste("Sites with p < 0.05 in All Three Models"), 
    x = "Causal Estimate (95% C.I.)", y = NULL, 
    color = NULL, fill = NULL, shape = NULL) +
  # coord_cartesian(xlim = c(0,1.5)) +
  scale_color_manual(values = c("indianred3","thistle4","mistyrose3")) +
  theme_light() +
  theme(plot.title = element_text(size = 10),
    axis.text.y = element_text(),
    axis.text.x = element_text(size = 8),
    axis.title.x = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.position = "bottom")

dev.off()



mrdoc_formerSmk_nominal <- mrdoc_formerSmk_Res |>
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

mrdoc_formerSmk_nominal |> 
  count(g1_nominal, g2_nominal)

mrdoc_formerSmk_nominal |> 
  count(robust_smk2dnam, robust_dnam2smk)

mrdoc_formerSmk_nominal |> 
  filter(g1_nominal) |> 
  select(cpg, CHR, SYMBOL, 
         mrdoc2.g1_Z, 
         mrdoc2.g2_Z, 
         robust_smk2dnam ) |> 
  arrange(mrdoc2.g1_Z)

mrdoc_formerSmk_nominal |> 
  filter(g2_nominal) |> 
  select(cpg, CHR, SYMBOL, 
         mrdoc2.g1_Z, 
         mrdoc2.g2_Z, 
         robust_dnam2smk ) |> 
  arrange(mrdoc2.g2_Z)

mrdoc_formerSmk_nominal |> 
  filter(g1_nominal, g2_nominal) |> 
  select(cpg, CHR, SYMBOL, 
         mrdoc2.g1_Z, 
         mrdoc2.g2_Z, 
         robust_smk2dnam, robust_dnam2smk ) |> 
  arrange(mrdoc2.g1_Z)



## Save supplementary annonated tables (nominally significant) ==================

## Add annotation from Illumina manifest
annoDir <- paste0(baseDir,"data/annotation/")
## Load the Annotation Files
load(paste0(annoDir,"manifest.RData"))
head(rownames(manifest))

manifest <- manifest[mrdoc_formerSmk_Res$cpg, 
                     c("IlmnID","CHR","MAPINFO","Genome_Build","UCSC_RefGene_Name",
                       "UCSC_RefGene_Group","UCSC_CpG_Islands_Name","Relation_to_UCSC_CpG_Island")]
dim(manifest)
# 2330    8

## Convert missing values to NA
manifestT <- manifest |> 
  mutate(IlmnID = as.character(IlmnID),
         CHR = as.character(CHR),
         UCSC_RefGene_Name = as.character(UCSC_RefGene_Name),
         UCSC_RefGene_Group = as.character(UCSC_RefGene_Group),
         UCSC_CpG_Islands_Name = as.character(UCSC_CpG_Islands_Name),
         Relation_to_UCSC_CpG_Island = as.character(Relation_to_UCSC_CpG_Island))

str(manifestT)

manifestT <- manifestT |> 
  mutate(IlmnID = ifelse(IlmnID == "", NA, IlmnID),
         CHR = ifelse(CHR == "", NA, CHR),
         UCSC_RefGene_Name = ifelse(UCSC_RefGene_Name == "", NA, UCSC_RefGene_Name),
         UCSC_RefGene_Group = ifelse(UCSC_RefGene_Group == "", NA, UCSC_RefGene_Group),
         UCSC_CpG_Islands_Name = ifelse(UCSC_CpG_Islands_Name == "", NA, UCSC_CpG_Islands_Name),
         Relation_to_UCSC_CpG_Island = ifelse(Relation_to_UCSC_CpG_Island == "", NA, Relation_to_UCSC_CpG_Island) )

str(manifestT)

load(paste0(annoDir,"cpgInfo_08112016.RData"))  ## Info on the nearest gene, if the CG is intergenic
dim(info)
info$IlmnID <- rownames(info)
info <- info |> 
  rename(NearestGene = SYMBOL)

str(info)

nearestgene <- info |> 
  filter(IlmnID %in% manifestT$IlmnID) |> 
  select(IlmnID,NearestGene)


mrdoc_formerSmk_nominal <- mrdoc_formerSmk_nominal |> 
  rename(IlmnID = cpg) |> 
  select(-c(CHR,MAPINFO,SYMBOL)) |> 
  left_join(manifestT, by = "IlmnID") |> 
  left_join(nearestgene, by = "IlmnID")

dim(mrdoc_formerSmk_nominal)

# rm(manifest,nearestgene); gc()


## Make subsets
mrdoc_smk2DNAm_nominal_g1 <- mrdoc_formerSmk_nominal |> 
  filter(g1_nominal) |> 
  rename(g1_robust = robust_smk2dnam,
         g2_robust = robust_dnam2smk) |> 
  select(-c(contains("g2"),g1_nominal)) |> 
  as.data.frame() |> 
  relocate(starts_with("mrdoc"), .after = last_col()) |> 
  relocate(contains("Sig"), .after = last_col()) |> 
  relocate(contains("robust"), .after = last_col())

mrdoc_smk2DNAm_nominal_g2 <- mrdoc_formerSmk_nominal |> 
  filter(g1_nominal) |> 
  rename(g1_robust = robust_smk2dnam,
         g2_robust = robust_dnam2smk) |> 
  select(-c(contains("g1"))) |> 
  as.data.frame() |> 
  relocate(starts_with("mrdoc"), .after = last_col()) |> 
  relocate(contains("Sig"), .after = last_col()) |> 
  relocate(contains("nominal"), .after = last_col()) |> 
  relocate(contains("robust"), .after = last_col())

mrdoc_DNAm2smk_nominal_g1 <- mrdoc_formerSmk_nominal |> 
  filter(g2_nominal) |> 
  rename(g1_robust = robust_smk2dnam,
         g2_robust = robust_dnam2smk) |> 
  select(-c(contains("g2"))) |> 
  as.data.frame() |> 
  relocate(starts_with("mrdoc"), .after = last_col()) |> 
  relocate(contains("Sig"), .after = last_col()) |> 
  relocate(contains("nominal"), .after = last_col()) |> 
  relocate(contains("robust"), .after = last_col())

mrdoc_DNAm2smk_nominal_g2 <- mrdoc_formerSmk_nominal |> 
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
  "Former_Smk_to_DNAm_g1" = mrdoc_smk2DNAm_nominal_g1, 
  "Former_Smk_to_DNAm_Reverse_g2" = mrdoc_smk2DNAm_nominal_g2,
  "DNAm_to_Former_Smk_g2" = mrdoc_DNAm2smk_nominal_g2, 
  "DNAm_to_Former_Smk_Reverse_g1" = mrdoc_DNAm2smk_nominal_g1
)
openxlsx::write.xlsx(outExcel, 
                     file = paste0(outDir,"suppl_mrdoc_formerSmk_dnam.xlsx"))



write.table(mrdoc_smk2DNAm_nominal_g1, 
            file = paste0(outDir,"suppl_formerSmk2dnam_g1est.tsv"), quote = F, row.names = F, col.names = T)
write.table(mrdoc_smk2DNAm_nominal_g2, 
            file = paste0(outDir,"suppl_formerSmk2dnam_g2est.tsv"), quote = F, row.names = F, col.names = T)
write.table(mrdoc_DNAm2smk_nominal_g2, 
            file = paste0(outDir,"suppl_dnam2formerSmk_g2est.tsv"), quote = F, row.names = F, col.names = T)
write.table(mrdoc_DNAm2smk_nominal_g1, 
            file = paste0(outDir,"suppl_dnam2formerSmk_g1est.tsv"), quote = F, row.names = F, col.names = T)



#### Save the CpG IDs ===============

smk2dnam_intersect <- mrdoc_formerSmk_nominal |> 
  filter(g1_nominal) |> 
  select(IlmnID)
length(smk2dnam_intersect$IlmnID)

write.table(smk2dnam_intersect, file = paste0(outDir,"formerSmk_2_dnam_mrdoc_nominal_intersect.tsv"), 
            col.names = F, row.names = F, quote = F)



#### Save the Gene names ===============

smk2dnam_gene <- mrdoc_formerSmk_nominal |> 
  filter(g1_nominal) |> 
  arrange(mrdoc1_wPleio.g1_p) |> 
  select(NearestGene) |> 
  t()

write.table(smk2dnam_gene, file = paste0(outDir,"formerSmk_2_dnam_mrdoc_nominal_intersect_genes.csv"), 
            sep = ",", col.names = F, row.names = F, quote = F)



## UpSet Plots #####################################

### Smk --> DNAm FDR UpSet Plot ######################

mrdoc_formerSmk_DNAm_FDR <- mrdoc_formerSmk_Res |> 
  select( cpg, CHR, SYMBOL, 
          mrdoc2.g1_fdrSig, mrdoc1_wPleio.g1_fdrSig, mrdoc1_wRE.g1_fdrSig ) |> 
  filter( mrdoc2.g1_fdrSig == T | mrdoc1_wPleio.g1_fdrSig == T | mrdoc1_wRE.g1_fdrSig == T  )


mrdoc_formerSmk_DNAm_FDR <- mrdoc_formerSmk_DNAm_FDR |> 
  rename( 
    "MR-DoC2" = 'mrdoc2.g1_fdrSig', 
    "MR-DoC1_wPleio" = "mrdoc1_wPleio.g1_fdrSig", 
    "MR-DoC1_wRE" = "mrdoc1_wRE.g1_fdrSig" 
  )

models <- colnames(mrdoc_formerSmk_DNAm_FDR)[-(1:3)]

#### Plot

jpeg(paste0(plotDir,"mrdoc_unidir_formerSmk_dnam_UpSet.jpeg"), width = 8, height = 5, units = "in", res = 300)

upset(
  ## Data with binary vars with group membership
  data = mrdoc_formerSmk_DNAm_FDR, 
  ## Columns with group membership
  intersect = models,  
  ## the label shown below the intersection matrix
  name="MR-DoC Model", 
  min_size = 0,
  width_ratio = 0.125 ,
  ## sort
  sort_intersections=FALSE,
  intersections=list(
    c( "MR-DoC1_wPleio","MR-DoC1_wRE","MR-DoC2" ),
    "MR-DoC1_wPleio",
    "MR-DoC1_wRE"
  ),
  ## display set size
  set_sizes=(
    upset_set_size()
    + geom_text(aes(label=after_stat(count)), 
                size=4, hjust=1.1, stat='count')
    # you can also add annotations on top of bars:
    + expand_limits(y=550)
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
  labs(title = "Putative Causal Effects of Former Smoking on DNA Methylation\nFDR < 0.05",
       y = "Model",
       x = "Intersection of Significant Causal Estimates Across Models") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10))

dev.off()



### DNAm --> Smk FDR UpSet Plot ######################

mrdoc_DNAm_formerSmk_FDR <- mrdoc_formerSmk_Res |> 
  select( cpg, CHR, SYMBOL, 
          mrdoc2.g2_fdrSig, mrdoc1_wPleio.g2_fdrSig, mrdoc1_wRE.g2_fdrSig ) |> 
  filter( mrdoc2.g2_fdrSig == T | mrdoc1_wPleio.g2_fdrSig == T | mrdoc1_wRE.g2_fdrSig == T  )


mrdoc_DNAm_formerSmk_FDR <- mrdoc_DNAm_formerSmk_FDR |> 
  rename( 
    "MR-DoC2" = 'mrdoc2.g2_fdrSig', 
    "MR-DoC1_wPleio" = "mrdoc1_wPleio.g2_fdrSig", 
    "MR-DoC1_wRE" = "mrdoc1_wRE.g2_fdrSig" 
  )

models <- colnames(mrdoc_DNAm_formerSmk_FDR)[-(1:3)]


#### Plot

jpeg(paste0(plotDir,"mrdoc_unidir_dnam_formerSmk_UpSet.jpeg"), width = 8, height = 5, units = "in", res = 300)

upset(
  ## Data with binary vars with group membership
  data = mrdoc_DNAm_formerSmk_FDR, 
  ## Columns with group membership
  intersect = models,  
  ## the label shown below the intersection matrix
  name="MR-DoC Model", 
  min_size = 0,
  width_ratio = 0.125 ,
  ## sort
  sort_intersections=FALSE,
  intersections=list(
    c( "MR-DoC1_wRE","MR-DoC2" ),
    "MR-DoC1_wPleio"
  ),
  ## display set size
  set_sizes=(
    upset_set_size()
    + geom_text(aes(label=after_stat(count)), hjust=1.1, stat='count')
    # you can also add annotations on top of bars:
    + expand_limits(y=10)
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
  labs(title = "Putative Causal Effects of DNA Methylation on Former Smoking\nFDR < 0.05",
       y = "Model",
       x = "Intersection of Significant Causal Estimates Across Models") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10))
dev.off()

####

