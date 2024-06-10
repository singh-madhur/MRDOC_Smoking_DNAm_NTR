## MR-DOC2 model for bidirectional effects between Former Smoking and DNAm 
## BATCH1
## Madhur Singh

## DNAm is residualized for age, sex, cell counts, array row, and sample plate
## PRS is residualized for PCs and geno platform
## Include age and sex as covariates for Smk

# clear the workspace
rm(list = ls()); gc()

# Libraries
library(data.table)
setDTthreads(8)
Sys.setenv(OMP_NUM_THREADS=8)
library(OpenMx)
mxOption( NULL, "Default optimizer", "CSOLNP" )

## Libraries
baseDir      <- "/base_dir/"
phenoDir     <- paste0(baseDir,"data/phenos/")
mqtlResidDir <- paste0(baseDir,"data/mQTL_resid/")
mqtlDir      <- paste0(baseDir,"data/PRS_mQTL/top_score/")
dnamDir      <- paste0(baseDir,"data/DNAm_resid/")
godmcDir     <- paste0(baseDir,"data/GoDMC/batches/")
outDir       <- paste0(baseDir,"out/mrdoc2_formerSmk_dnam/")
outFile      <- paste0(outDir,"mrdoc2_formerSmk_DNAm_BATCH1.dat")

## Load the processed dataset with Smk, age, sex, and residualized PRS (long format)
twin_dat_long <- fread(paste0(phenoDir,"twin_dat_long_mrdoc_resid.csv"), header = T)
dim(twin_dat_long)

twin_dat_long$fisnumber <- as.numeric(twin_dat_long$fisnumber)

## Load the CpG IDs by batch
CpGs <- fread(paste0(godmcDir,"BATCH1_IDs.dat"), header = F)
colnames(CpGs) <- "ID"
head(CpGs)
nCpG <- length(CpGs$ID)
nCpG

mqtlFilesRaw <- list.files(mqtlResidDir)
mqtlFiles    <- mqtlFilesRaw[grep("_godmc_resid_prs.dat",mqtlFilesRaw)]

## File with GEE R-sq and estimate per mQTL
## Some mQTL have a negative estimate --> Need to multiply the score with -1
mqtlR2       <- fread(paste0(mqtlDir,"MRDoC_CpG_IV_with_maxR2.dat"), header = T, data.table = F)
table(mqtlR2$negEst)
# FALSE  TRUE 
# 12615   325 

negMQTL_IDs <- mqtlR2[mqtlR2$Estimate<0,]$CpG
length(negMQTL_IDs)
# 325

# intersection
length(intersect(CpGs$ID, negMQTL_IDs))

## Load the residualized beta values
load(file = paste0(dnamDir,"residualized_beta_blood450k_BATCH1.RData"))
dim(resBeta)  

cpgDat <- as.data.frame(resBeta)
cpgDat <- cpgDat[, colnames(cpgDat) %in% CpGs$ID ]
dim(cpgDat)

rm(resBeta); gc()


# MR-DOC Model #################################################################

#Starting Values
x=   1    # x is the standard deviation of the PGS for x
y=   1    # y is the standard deviation of the PGS for y
ax= sqrt(.7)
cx= sqrt(.1)
ex= sqrt(.2)    # a2 + c2 + e2 = 1
ay= sqrt(.7) 
cy= sqrt(.1) 
ey= sqrt(.2)    # a2 + c2 + e2 = 1
ra= .10    # the larger this is - when unmodeled the more it will inflate g1
rc= .10    # the larger this is - when unmodeled the more it will inflate g1  
re= .10    # the larger this is - when unmodeled the more it will inflate g1 
rIV= .10
bx1= .1    # larger values give larger genetic instrument strength
by1= .0    # non-zero values are a violation of the non-pleiotropy assumption
bx2= .0    # non-zero values are a violation of the non-pleiotropy assumption
by2= .1    # larger values give larger genetic instrument strength
g1=  .10   # x->y causal effect size 
g2=  .10   # y->x would mean that there is also reverse causation 

## Means and threshold
frMV      <- c(F,T,T,T)                # free status for variables (Smk,CpG,PRS,mQTL)
svMe      <- c(0,0,0,0)                # start value for means
svTh      <- .1                        # start value for thresholds

# Select Binary Exposure (smoking) Variable
expos       <- "formerSmk"            
nExpos      <- 2                     # number of smoking variables
expVar      <- paste(expos,c(1,2),sep="_")

# Select Continuous Outcome (CpG DNAm) Variable
out       <- 'cpg'              # list of variables names

# Select the genetic IV
IVx     <- 'resPRS'
IVy     <- 'resMQTL'

# Total phenotypes
vars    <- c(expos,out,IVx,IVy)
nv      <- length(vars)
ntv     <- nv*2

vnames  <- paste(vars, c(rep(1,nv),rep(2,nv)), sep = "_")

# Select Definition variables/Covariates for Analysis

covar  <- c("age", "female")
nCov   <- length(covar)
indVar <- paste(covar,c(rep(1,nCov),rep(2,nCov)),sep="_")

selVars <- c(vnames, indVar)



#### DEFINE MATRICES FOR OPENMX ################################################

# We first define the openmx matrices L, B and Psy using mxMatrix 
# mxMatrix requires various input matrices which we define here: 
# mxfree_ : which parameters should be estimated in the matrices BETA and PSY 
# mxlabels_: parameter names
# mxvalues_: parameter starting values or fixed values


ny=8
ne=20


# Dimnames for Beta and Psi matrices
matDim <- c(vnames,                                            # 1:8
            paste(c("A","C","E"),expos,"1",sep = "_"),         # 9:11
            paste(c("A","C","E"),out,"1",sep = "_"),           # 12:14
            paste(c("A","C","E"),expos,"2",sep = "_"),         # 15:17  
            paste(c("A","C","E"),out,"2",sep = "_"))           # 18:20


##### MATRIX BETA ##############################################################

mxfreeBE=matrix(F,ne,ne,dimnames=list(matDim,matDim))

# [1] "formerSmk_1"   "cpg_1"          "resPRS_1"      
# [4] "resMQTL_1"      "formerSmk_2"   "cpg_2"         
# [7] "resPRS_2"       "resMQTL_2"      "A_formerSmk_1"
# [10] "C_formerSmk_1" "E_formerSmk_1" "A_cpg_1"       
# [13] "C_cpg_1"        "E_cpg_1"        "A_formerSmk_2"
# [16] "C_formerSmk_2" "E_formerSmk_2" "A_cpg_2"       
# [19] "C_cpg_2"        "E_cpg_2"    

# Twin 1
mxfreeBE[1,2]=T #'g2' 
mxfreeBE[1,3]=T	#'bx1'
# mxfreeBE[1,4]=F	#'bx2'  # IV2 --> X Fixed at 0
mxfreeBE[2,1]=T	#'g1'
# mxfreeBE[2,3]=F	#'by1'  # IV1 --> Y Fixed at 0
mxfreeBE[2,4]=T	#'by2'

# Twin 2
mxfreeBE[5,6]=T #'g2' 
mxfreeBE[5,7]=T	#'bx1'
# mxfreeBE[5,8]=F	#'bx2'  # IV2 --> X Fixed at 0
mxfreeBE[6,5]=T	#'g1'
# mxfreeBE[6,7]=F	#'by1'  # IV1 --> Y Fixed at 0
mxfreeBE[6,8]=T	#'by2'
#
mxfreeBE[1, 9]=T	#'aSmk'
mxfreeBE[1,10]=T  #'cSmk'
mxfreeBE[1,11]=T	#'eSmk'
mxfreeBE[2,12]=T	#'aCpG'
mxfreeBE[2,13]=T	#'cCpG' 
mxfreeBE[2,14]=T	#'eCpG'
mxfreeBE[5,15]=T	#'aSmk'
mxfreeBE[5,16]=T 	#'cSmk'
mxfreeBE[5,17]=T	#'eSmk'
mxfreeBE[6,18]=T	#'aCpG'
mxfreeBE[6,19]=T	#'cCpG'
mxfreeBE[6,20]=T	#'eCpG'


mxlabelsBE=matrix(NA,ne,ne,dimnames=list(matDim,matDim))
mxlabelsBE[1,2]='g2'
mxlabelsBE[1,3]='bx1'
mxlabelsBE[1,4]='bx2'
mxlabelsBE[2,1]='g1'
mxlabelsBE[2,3]='by1'
mxlabelsBE[2,4]='by2'
#
mxlabelsBE[5,6]='g2'
mxlabelsBE[5,7]='bx1'
mxlabelsBE[5,8]='bx2'
mxlabelsBE[6,5]='g1'
mxlabelsBE[6,7]='by1'
mxlabelsBE[6,8]='by2'
#
mxlabelsBE[1, 9]='aSmk'
mxlabelsBE[1,10]='cSmk'
mxlabelsBE[1,11]='eSmk'
mxlabelsBE[2,12]='aCpG'
mxlabelsBE[2,13]='cCpG'
mxlabelsBE[2,14]='eCpG'
mxlabelsBE[5,15]='aSmk'
mxlabelsBE[5,16]='cSmk'
mxlabelsBE[5,17]='eSmk'
mxlabelsBE[6,18]='aCpG'
mxlabelsBE[6,19]='cCpG'
mxlabelsBE[6,20]='eCpG'

# VALUES

mxvaluesBE=matrix(0,ne,ne,dimnames=list(matDim,matDim))
mxvaluesBE[1,2]=g2
mxvaluesBE[1,3]=bx1
mxvaluesBE[1,4]=bx2
mxvaluesBE[2,1]=g1
mxvaluesBE[2,3]=by1
mxvaluesBE[2,4]=by2
#
mxvaluesBE[5,6]=g2
mxvaluesBE[5,7]=bx1
mxvaluesBE[5,8]=bx2
mxvaluesBE[6,5]=g1
mxvaluesBE[6,7]=by1
mxvaluesBE[6,8]=by2
#
mxvaluesBE[1, 9]=ax 
mxvaluesBE[1,10]=cx 
mxvaluesBE[1,11]=ex 
mxvaluesBE[2,12]=ay 
mxvaluesBE[2,13]=cy  
mxvaluesBE[2,14]=ey
mxvaluesBE[5,15]=ax 
mxvaluesBE[5,16]=cx 
mxvaluesBE[5,17]=ex 
mxvaluesBE[6,18]=ay 
mxvaluesBE[6,19]=cy 
mxvaluesBE[6,20]=ey 




##### MATRIX PSY ################################################################
# [1] "formerSmk_1"   "cpg_1"          "resPRS_1"      
# [4] "resMQTL_1"      "formerSmk_2"   "cpg_2"         
# [7] "resPRS_2"       "resMQTL_2"      "A_formerSmk_1"
# [10] "C_formerSmk_1" "E_formerSmk_1" "A_cpg_1"       
# [13] "C_cpg_1"        "E_cpg_1"        "A_formerSmk_2"
# [16] "C_formerSmk_2" "E_formerSmk_2" "A_cpg_2"       
# [19] "C_cpg_2"        "E_cpg_2"  
mxfreePSY=matrix(F,ne,ne,dimnames= list(matDim,matDim)) 
mxfreePSY[12, 9]=T # 'ra'
mxfreePSY[15,12]=T # 'ra'
mxfreePSY[18, 9]=T # 'ra'
mxfreePSY[18,15]=T # 'ra'
#
mxfreePSY[13,10]=T # 'rc'
mxfreePSY[19,16]=T # 'rc'
mxfreePSY[16,13]=T # 'rc'
mxfreePSY[19,10]=T # 'rc'
#
mxfreePSY[14,11]=T # re  
mxfreePSY[20,17]=T # re   
#
mxfreePSY[ 4, 3]=T # rIV
mxfreePSY[ 8, 7]=T # rIV 
mxfreePSY[ 8, 3]=T # rIV 
mxfreePSY[ 7, 4]=T # rIV 


# symmetric fill of upper half of symm matrix
for (i in 1:ne) { for (j in 1:i) { 
  mxfreePSY[j,i]=mxfreePSY[i,j]}}


mxlabelsPSY=matrix(NA,ne,ne,dimnames= list(matDim,matDim)) 
mxlabelsPSY[12, 9]= 'ra'
mxlabelsPSY[15,12]= 'ra' 
mxlabelsPSY[18, 9]= 'ra' 
mxlabelsPSY[18,15]= 'ra' 
mxlabelsPSY[13,10]= 'rc' 
mxlabelsPSY[16,13]= 'rc'
mxlabelsPSY[19,10]= 'rc' 
mxlabelsPSY[19,16]= 'rc' 
mxlabelsPSY[14,11]= 're' 
mxlabelsPSY[20,17]= 're' 
mxlabelsPSY[ 4, 3]= 'rIV'
mxlabelsPSY[ 8, 7]= 'rIV' 
mxlabelsPSY[ 8, 3]= 'rIV' 
mxlabelsPSY[ 7, 4]= 'rIV' 

# symmetric fill of upper half of symm matrix
for (i in 1:ne) { for (j in 1:i) { 
  mxlabelsPSY[j,i]=mxlabelsPSY[i,j]}}

# VALUES
# [1] "formerSmk_1"   "cpg_1"          "resPRS_1"      
# [4] "resMQTL_1"      "formerSmk_2"   "cpg_2"         
# [7] "resPRS_2"       "resMQTL_2"      "A_formerSmk_1"
# [10] "C_formerSmk_1" "E_formerSmk_1" "A_cpg_1"       
# [13] "C_cpg_1"        "E_cpg_1"        "A_formerSmk_2"
# [16] "C_formerSmk_2" "E_formerSmk_2" "A_cpg_2"       
# [19] "C_cpg_2"        "E_cpg_2" 
mxvaluesPSY=matrix(0,ne,ne,dimnames=list(matDim,matDim)) 
diag(mxvaluesPSY)=c(0,0,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
# fixed
mxvaluesPSY[15, 9]=1   # r(Ax_twin1, Ax_twin2)
mxvaluesPSY[16,10]=1   # r(Cx_twin1, Cx_twin2)
mxvaluesPSY[17,11]=0   # r(Ex_twin1, Ex_twin2)
mxvaluesPSY[18,12]=1   # r(Ay_twin1, Ay_twin2)
mxvaluesPSY[19,13]=1   # r(Cy_twin1, Cy_twin2)
mxvaluesPSY[20,14]=0   # r(Ey_twin1, Ry_twin2)
mxvaluesPSY[ 7, 3]=1   # correlation pgs twin1 - twin2
mxvaluesPSY[ 8, 4]=1   # correlation mQTL twin1 - twin2
# free
mxvaluesPSY[12, 9]= ra
mxvaluesPSY[15,12]= ra 
mxvaluesPSY[18, 9]= ra 
mxvaluesPSY[18,15]= ra 
mxvaluesPSY[13,10]= rc 
mxvaluesPSY[16,13]= rc
mxvaluesPSY[19,10]= rc 
mxvaluesPSY[19,16]= rc 
mxvaluesPSY[14,11]= re 
mxvaluesPSY[20,17]= re 
mxvaluesPSY[ 4, 3]= rIV
mxvaluesPSY[ 8, 7]= rIV 
mxvaluesPSY[ 8, 3]= rIV 
mxvaluesPSY[ 7, 4]= rIV 

# symmetric fill of upper half of symm matrix
for (i in 1:ne) { for (j in 1:i) { 
  mxvaluesPSY[j,i]=mxvaluesPSY[i,j]}}

# values for DZ ; we will later multiply the MZ-PSY matrix with this one to get DZ-PSY2, 
# hence the initial '1' everywhere
mxvaluesPSY2=matrix(1,ne,ne,dimnames=list(matDim,matDim)) 
mxvaluesPSY2[15, 9]= 0.5  # r(Ax_twin1, Ax_twin2)
mxvaluesPSY2[15,12]= 0.5  # 'ra' (X_t2, Y_t1)
mxvaluesPSY2[18, 9]= 0.5  # 'ra' (X_t1, Y_t2)
mxvaluesPSY2[18,12]= 0.5  # r(Ay_twin1, Ay_twin2)
mxvaluesPSY2[ 7, 3]= 0.5  # correlation pgs twin1 - twin2
mxvaluesPSY2[ 8, 4]= 0.5  # correlation mQTL twin1 - twin2
mxvaluesPSY2[ 8, 3]= 0.5  # r(PGS_t1, mQTL_t2)
mxvaluesPSY2[ 7, 4]= 0.5  # r(PGS_t2, mQTL_t1)

for (i in 1:ne) { for (j in 1:i) { 
  mxvaluesPSY2[j,i]=mxvaluesPSY2[i,j]}}



##### MATRIX LAMBDA ############################################################## 

mxfreeLY=matrix(F,ny,ne,dimnames=list(matDim[1:ny],matDim))
mxfreeLY[3,3]=mxfreeLY[7,7]=T  # IVx std dev
mxfreeLY[4,4]=mxfreeLY[8,8]=T  # IVx std dev


mxlabelsLY=matrix(NA,ny,ne,dimnames=list(matDim[1:ny],matDim))
mxlabelsLY[3,3]= mxlabelsLY[7,7]='SDprs' # pgs variance
mxlabelsLY[4,4]= mxlabelsLY[8,8]='SDmqtl'

# VALUES

mxvaluesLY=matrix(0,ny,ne,dimnames=list(matDim[1:ny],matDim))
mxvaluesLY[1:ny,1:ny]=diag(ny)
mxvaluesLY[3,3]=mxvaluesLY[7,7]=x  # IVx std dev
mxvaluesLY[4,4]=mxvaluesLY[8,8]=y  # IVx std dev


# What is the role of the mysterious matrix mXFilt in
# mxMatrix(type='Full',nrow=5,ncol=6,free=F,values=mxFilt, name='F')?
# This simply removes the row and column relating to the second twin's pgs in the MZ group. 
# We remove this because it is correlated 1 with the first twin's pgs in the MZ group,
# causing the MZ covariance matrix to be singular. 
# If the covariance matrix is singular, OpenMx cannot estimate the maximum likelihood parameters.  
ny1=ny-2
#
mxFilt=matrix(c(
  1,0,0,0,0,0,0,0,
  0,1,0,0,0,0,0,0,
  0,0,1,0,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,0,1,0,0,0,
  0,0,0,0,0,1,0,0),ny1,ny,byrow=T,
  dimnames=list(vnames[1:ny1],vnames) )



##### MEANS and COVARIATES ######################################################

# # Matrix to hold definition variables for the regression
# Twin 1
dCov1 <- mxMatrix( type = "Full", nrow = nCov, ncol = 1, free = FALSE,
                   labels = paste("data.",paste(covar,c(rep(1,nCov)),sep="_"),sep=""), 
                   name = "dCov1" )
# Twin 2
dCov2 <- mxMatrix( type = "Full", nrow = nCov, ncol = 1, free = FALSE,
                   labels = paste("data.",paste(covar,c(rep(2,nCov)),sep="_"),sep=""), 
                   name = "dCov2" )

# Matrix to hold betas for the covariates
bCov <- mxMatrix( type = "Full", nrow = nv, ncol = nCov, 
                  free = c(rep(T,nCov),                      # Smk  - regressed on covariates
                           rep(F,nCov),                      # CpG  - not regressed - already residualized
                           rep(F,nCov),                      # PRS  - not regressed on any covariate (exogenous IV)
                           rep(F,nCov)),                     # mQTL - not regressed on any covariate (exogenous IV)
                  values = c(rep(0,nCov*nv)), byrow = T,
                  labels = c(paste("bSmk",covar,sep="_"),
                             rep(NA,nCov),
                             rep(NA,nCov),
                             rep(NA,nCov)), 
                  dimnames = list(vars,covar),
                  name = "bCov" )


# Create Algebra for expected Mean & Threshold Matrices
meanT1     <- mxMatrix( type="Full", nrow=1, ncol=nv, free=frMV, values=svMe, 
                        labels=paste("mean", vars, sep = "_"), dimnames=list("mean",vars), 
                        name="meanT1" )
# In MZ twins, prs_twin1 == prs_twin2
# So, Twin2 in MZ pairs does not have the prs/mQTL variables
meanT2MZ     <- mxMatrix( type="Full", nrow=1, ncol=2, free=frMV[1:2], values=svMe[1:2],  
                          labels=paste("mean", c(expos,out), sep = "_"), dimnames=list("mean",vars[1:2]), 
                          name="meanT2MZ" )
# In DZ twins, prs_twin1 != prs_twin2
meanT2DZ     <- mxMatrix( type="Full", nrow=1, ncol=nv, free=frMV, values=svMe,  
                          labels=paste("mean", vars, sep = "_"), dimnames=list("mean",vars), 
                          name="meanT2DZ" )

expMeanMZ <- mxAlgebra(expression = cbind(meanT1 + t(bCov%*%dCov1),
                                          meanT2MZ + t(bCov[1:2,] %*%dCov2)
) ,  name="expMeanMZ" )   
expMeanDZ <- mxAlgebra(expression = cbind(meanT1 + t(bCov%*%dCov1),
                                          meanT2DZ + t(bCov%*%dCov2)
) ,  name="expMeanDZ" )   

thresh    <- mxMatrix( type="Full", nrow=1, ncol=2, free=TRUE, values=svTh, 
                       labels=c("thrSmk","thrSmk"), name="thresh" )


## Computational settings
funML    <-    mxFitFunctionML()
plan     <-    omxDefaultComputePlan()


## Matrices to store results #####################################################

mrdocEst  <- matrix(NA_real_, 1, 19)
mrdocSE   <- matrix(NA_real_, 1, 19)
mrdocZ    <- matrix(NA_real_, 1, 19)
mrdocPval <- matrix(NA_real_, 1, 4)  # only computing for g1, g2, bx1, by2

namesGoF <- c("statusCode","nStats","nObs","nPars","nConst","df","minus2LL",
              "AICpar","AICsample","BICpar","BICsample","wallTime")
nGoF     <- length(namesGoF)
mrdocFit  <- matrix(NA_real_, 1, nGoF)


## START LOOP #######################################################################

t1 <- Sys.time()

for (i in 1:nCpG) {
  
  # Extract the CpG column
  cpgT = colnames(cpgDat)[i]
  cpgT
  
  smk_cpg_i <- as.data.frame(cbind(IdatName = rownames(cpgDat),
                                   cpg = cpgDat[, colnames(cpgDat) %in% cpgT ]))
  smk_cpg_i$cpg <- as.numeric(smk_cpg_i$cpg)
  smk_cpg_i$cpg <- as.numeric(scale(smk_cpg_i$cpg))
  
  ## Load the mQTL
  mqtlT <- mqtlFiles[grep(cpgT,mqtlFiles)]
  mQTL <- fread(paste0(mqtlResidDir,mqtlT), header = T)
  mQTL$fisnumber <- as.numeric(mQTL$fisnumber)
  
  ## If the mQTL had a negative association in GEE, multiply it by -1
  if(cpgT %in% negMQTL_IDs) {
    mQTL$resMQTL <- -1 * mQTL$resMQTL
  }
  
  ## Merge
  twin_dat_mQTL <- merge(twin_dat_long, mQTL, by = "fisnumber", all.x = T)
  
  # Bind with DNAm data
  twin_dat_mrdoc <- merge(twin_dat_mQTL, smk_cpg_i, by = "IdatName", all.x = T)
  
  ## Pivot to wide format
  twin_dat <- reshape(
    data = twin_dat_mrdoc,  direction = "wide",
    v.names = c("fisnumber", "IdatName", "age", "female", "resPRS", "resMQTL",
                "everSmk", "currentSmk", "formerSmk", "smkCess",
                "ncigday", "nyearstop", "cpg", "Platform_2" , "Platform_3",
                "PC1_sc", "PC2_sc", "PC3_sc", "PC4_sc", "PC5_sc", 
                "PC6_sc","PC7_sc","PC8_sc","PC9_sc","PC10_sc",
                "Sample_Plate", "Array_rownum", "Neut_Perc", "Mono_Perc", "Eos_Perc" ), 
    idvar = c("FamilyNumber","zygo"),
    timevar = "Twin_nr", sep = "_"
  )
  
  
  
  ##### MR-DOC #####################################################################
  
  # Select Data for Analysis
  mzData <- subset(twin_dat, zygo==1, selVars)
  dzData <- subset(twin_dat, zygo==2, selVars)
  # Transform binary Smk pheno into an ordered factor
  mzDataF   <- cbind(formerSmk_1 = mxFactor( x=mzData[,1], levels=c(0:1)), mzData[,2:4],
                     formerSmk_2 = mxFactor( x=mzData[,5], levels=c(0:1)), mzData[,-c(1:5)]) 
  dzDataF   <- cbind(formerSmk_1 = mxFactor( x=dzData[,1], levels=c(0:1)), dzData[,2:4],
                     formerSmk_2 = mxFactor( x=dzData[,5], levels=c(0:1)), dzData[,-c(1:5)]) 
 
  # Matrices
  Matmodel=mxModel('TwRM',
                   mxMatrix(type='Full',nrow=ny1,ncol=ny,free=F,values=mxFilt, name='Filt'),
                   # #
                   mxMatrix(type='Full',nrow=ne,ncol=ne, free=mxfreeBE,value=mxvaluesBE,labels=mxlabelsBE,
                            ubound=10,lbound=-10,dimnames=list(matDim,matDim),name='BE'),
                   #
                   mxMatrix(type='Iden',nrow=ne,ncol=ne,name='I'),
                   mxMatrix(type='Full',nrow=ny,ncol=ne, free=mxfreeLY,value=mxvaluesLY,labels=mxlabelsLY,
                            ubound=10,lbound=-10,dimnames=list(vnames,matDim),name='LY'),
                   mxMatrix(type='Symm',nrow=ne,ncol=ne, free=mxfreePSY,value=mxvaluesPSY,labels=mxlabelsPSY,
                            ubound=10,lbound=-10,dimnames=list(matDim,matDim),name='PS'),
                   mxMatrix(type='Symm',nrow=ne,ncol=ne, free=F,value=mxvaluesPSY2,labels=NA, 
                            dimnames=list(matDim,matDim),name='PS2'),
                   # Constrain Variance of Binary Smk Variable
                   # Create Algebra for expected Variance
                   mxAlgebra( expression= (BE[1,9])^2 + (BE[1,10])^2 + (BE[1,11])^2, name="VSmk" ),
                   # Constrain Variance
                   mxConstraint( expression=diag2vec(VSmk)==1, name="Var1Smk" ),
                   
                   mxAlgebra(expression=solve(I-BE),name='iBE'),
                   #
                   mxAlgebra(expression=Filt%*%LY%*%iBE%*%PS%*%t(iBE)%*%t(LY)%*%t(Filt),name='Smz'),
                   mxAlgebra(expression=LY%*%iBE%*%(PS*PS2)%*%t(iBE)%*%t(LY),name='Sdz')   
                   
  )
  
  
  #
  # a model the data, the fit function (MZ)
  MZmodel=mxModel("MZ",
                  meanT1, meanT2MZ, expMeanMZ, thresh,
                  mxData( observed=mzDataF, type="raw", verbose = 0L),   # the data minus pgst2
                  # covariates
                  dCov1, dCov2, bCov,
                  mxExpectationNormal(covariance="TwRM.Smz",means="expMeanMZ",
                                      thresholds="thresh", threshnames=expVar,
                                      dimnames=vnames[1:ny1]),
                  # uncomment one of the 2 fit functions
                  funML,
                  # funWLS,
                  plan
  )
  
  # a model the data, the fit function (DZ)
  DZmodel=mxModel("DZ",
                  meanT1, meanT2DZ, expMeanDZ, thresh,
                  mxData(observed=dzDataF, type="raw", verbose = 0L),
                  # covariates
                  dCov1, dCov2, bCov,
                  mxExpectationNormal(covariance="TwRM.Sdz",means="expMeanDZ",
                                      thresholds="thresh", threshnames=expVar,
                                      dimnames = vnames), 
                  # uncomment one of the 2 fit functions
                  funML,
                  # funWLS,
                  plan
  )
  
  mrdocModel <-  mxModel( "mrdoc", Matmodel, MZmodel, DZmodel,
                          mxFitFunctionMultigroup(c("MZ","DZ")) )
  
  mrdocModel <- omxSetParameters(mrdocModel, labels = c("cSmk","rc","cCpG"), free = F, values = 0)
  
  # fit the model 
  mrdocOut <- mxTryHardOrdinal(mrdocModel, exhaustive = F, OKstatuscodes = 0)
  sumMRDOC <- summary(mrdocOut)
  
  ## check for any missing SE 
  mrdocStdErr <- as.vector(mrdocOut$output$standardErrors)
  if( sum(is.na(mrdocStdErr)) > 0 ) {
    ## Exhaustive mxTryHard
    mrdocOut <- mxTryHardOrdinal(mrdocOut, exhaustive = T, OKstatuscodes = 0)
    sumMRDOC <- summary(mrdocOut)
  }
  
  
  ### SAVE RESULTS ##################################################################
  ## Save 
  mrdocCoef   <- as.vector(coef(mrdocOut))
  mrdocStdErr <- as.vector(mrdocOut$output$standardErrors)
  mrdocZedd   <- as.vector( mrdocCoef / mrdocStdErr )
  mrdocP      <- 2*( 1-pnorm( abs(mrdocZedd[1:4]) ) )
  mrdocGoF    <- as.vector(c(mrdocOut$output$status$code,
                             sumMRDOC$observedStatistics,
                             sumMRDOC$numObs,
                             sumMRDOC$estimatedParameters,
                             sumMRDOC$constraints,
                             sumMRDOC$degreesOfFreedom,
                             sumMRDOC$Minus2LogLikelihood,
                             sumMRDOC$informationCriteria["AIC:","par"],
                             sumMRDOC$informationCriteria["AIC:","sample"],
                             sumMRDOC$informationCriteria["BIC:","par"],
                             sumMRDOC$informationCriteria["BIC:","sample"],
                             mrdocOut$output$wallTime)) 
  
  
  if(i==1) {
    mrdocEst[1,]  <- mrdocCoef;      rownames(mrdocEst) = colnames(cpgDat)[i]
    mrdocSE[1,]   <- mrdocStdErr;    rownames(mrdocSE) = colnames(cpgDat)[i]
    mrdocZ[1,]    <- mrdocZedd;      rownames(mrdocZ) = colnames(cpgDat)[i]
    mrdocPval[1,] <- mrdocP;         rownames(mrdocPval) = colnames(cpgDat)[i]
    mrdocFit[1,]  <- mrdocGoF;       rownames(mrdocFit) = colnames(cpgDat)[i]
  } else {
    mrdocEst   <- rbind(mrdocEst, mrdocCoef);    rownames(mrdocEst)[i] = colnames(cpgDat)[i]
    mrdocSE    <- rbind(mrdocSE, mrdocStdErr);   rownames(mrdocSE)[i] = colnames(cpgDat)[i]
    mrdocZ     <- rbind(mrdocZ, mrdocZedd);      rownames(mrdocZ)[i] = colnames(cpgDat)[i]
    mrdocPval  <- rbind(mrdocPval, mrdocP);      rownames(mrdocPval)[i] = colnames(cpgDat)[i]
    mrdocFit   <- rbind(mrdocFit, mrdocGoF);     rownames(mrdocFit)[i] = colnames(cpgDat)[i]
  }
  
  ## colnames
  colnames(mrdocEst)  = paste( names(coef(mrdocOut)), "hat", sep = "_" )
  colnames(mrdocSE)   = paste( names(coef(mrdocOut)), "SE", sep = "_" )
  colnames(mrdocZ)    = paste( names(coef(mrdocOut)), "Z", sep = "_" )
  colnames(mrdocPval) = paste( names(coef(mrdocOut))[1:4], "p", sep = "_" )
  colnames(mrdocFit)  = namesGoF
  ## cbind the matrices
  mrdocRes <- as.data.frame(cbind(mrdocEst,mrdocSE,mrdocZ,mrdocPval,mrdocFit))
  ## Fisher's transformation of rA and rE
  mrdocRes$ra_fisherZ <- DescTools::FisherZ(mrdocRes$ra_hat)
  mrdocRes$re_fisherZ <- DescTools::FisherZ(mrdocRes$re_hat)
  ## cpg ID
  mrdocRes$cpg = colnames(cpgDat)[1:i]
  
  ## Save
  fwrite(mrdocRes, file = outFile, sep = "\t", row.names = T, col.names = T, quote = F)
  #
  if( i%%100 == 0 ) { print(paste("Ran Cpg",i,"of",nCpG)) }
  
}

t2 <- Sys.time()
print(paste("Time cost =",(t2-t1)))

## END LOOP #####################################################################

