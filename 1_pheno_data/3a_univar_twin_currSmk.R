## Madhur Singh

## clear workspace
rm(list = ls()); gc()

## Directories #####################

baseDir   <- "/base_dir/"
phenoDir  <- paste0(baseDir,"data/phenos/")
scriptDir <- paste0(baseDir,"scripts/")

# Libraries
library(data.table)
setDTthreads(16)
library(dplyr)
library(tidyr)
library(polycor)      # for polyserial correlation of binary outcome

Sys.setenv(OMP_NUM_THREADS=16)
library(OpenMx)
source(paste0(scriptDir,"miFunctions.R"))
mxOption( NULL, "Default optimizer", "CSOLNP" )


## LOAD TWIN COVARIATE DATA ################################

## Load the processed dataset with Smk, age, sex, and residualized PRS (long format)
twin_dat_long <- fread(paste0(phenoDir,"twin_dat_long_mrdoc_resid.csv"), header = T)
twin_dat_long$fisnumber <- as.numeric(twin_dat_long$fisnumber)

str(twin_dat_long)

dim(twin_dat_long[!is.na(twin_dat_long$fisnumber),])

## Pivot to wide format
twin_dat <- twin_dat_long |>
  select(FamilyNumber, zygo, Twin_nr, fisnumber, age, female, currentSmk) |>
  pivot_wider(id_cols = c(FamilyNumber, zygo), names_from = Twin_nr,
              values_from = c(fisnumber, age, female, currentSmk),
              names_sep = "_")

dim(twin_dat)


## OpenMX model ----------------------------------------------------------------

# Select Variables for Analysis
vars      <- 'currentSmk'                 # list of variables names
nv        <- length(vars)                 # number of variables
ntv       <- nv*2                         # number of total variables
depVar    <- paste(vars,c(rep(1,nv),rep(2,nv)),sep="_")

# Select Covariates for Analysis
covar  <- c("age", "female") 
nCov   <- length(covar)
indVar <- paste(covar,c(rep(1,nCov),rep(2,nCov)),sep="_")

selVars <- c(depVar,indVar)

# Select Data for Analysis
mzData <- subset(twin_dat, zygo==1, selVars)
dzData <- subset(twin_dat, zygo==2, selVars)

mzDataF   <- cbind(mxFactor( x=mzData[,depVar], levels=c(0:1)), mzData[,indVar] ) 
dzDataF   <- cbind(mxFactor( x=dzData[,depVar], levels=c(0:1)), dzData[,indVar] ) 

dim(mzDataF); head(mzDataF)
dim(dzDataF); head(dzDataF)

# Generate Descriptive Statistics
hetcor(mzDataF$currentSmk_1, mzDataF$currentSmk_2)$cor
hetcor(dzDataF$currentSmk_1, dzDataF$currentSmk_2)$cor

# Set Starting Values
svBe      <- 0.01                      # start value for regressions
svTh      <- .15                       # start value for thresholds
svCor     <- .75                       # start value for correlations
lbCor     <- -0.999                    # lower bounds for correlations
ubCor     <- 0.999                     # upper bounds for correlations

# Create Data Objects for Multiple Groups

dataMZ <- mxData( observed=mzDataF, type="raw" )
dataDZ <- mxData( observed=dzDataF, type="raw" )

# Matrix to hold definition variables for the regression
dCov <- mxMatrix( type = "Full", nrow = nCov, ncol = 2, free = FALSE, 
                  labels = paste("data.",indVar,sep=""), name = "dCov" )

# Matrix to hold betas for the covariates
bCov <- mxMatrix( type = "Full", nrow = 1, ncol = nCov, free = TRUE, values = svBe,
                  labels = paste("b",covar,sep = "_"), name = "bCov" )

# Create Algebra for expected Mean & Threshold Matrices
meanG     <- mxMatrix( type="Zero", nrow=1, ncol=ntv, name="meanG" )
expMeanMZ <- mxAlgebra(expression = (meanG + bCov%*%dCov ) ,  name="expMeanMZ" )
expMeanDZ <- mxAlgebra(expression = (meanG + bCov%*%dCov ) ,  name="expMeanDZ" )
thresh    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svTh, labels=c("thr","thr"), name="thresh" )

# Create Algebra for expected Correlation Matrices
corMZ     <- mxMatrix( type="Stand", nrow=ntv, ncol=ntv, free=TRUE, values=svCor, 
                       lbound=lbCor, ubound=ubCor, labels="rMZ", name="corMZ" )
corDZ     <- mxMatrix( type="Stand", nrow=ntv, ncol=ntv, free=TRUE, values=svCor, 
                       lbound=lbCor, ubound=ubCor, labels="rDZ", name="corDZ" )

# Create Expectation Objects for Multiple Groups
expMZ     <- mxExpectationNormal( covariance="corMZ", means="expMeanMZ", dimnames=depVar, thresholds="thresh" )
expDZ     <- mxExpectationNormal( covariance="corDZ", means="expMeanDZ", dimnames=depVar, thresholds="thresh" )
funML     <- mxFitFunctionML()

# Create Model Objects for Multiple Groups
pars      <- c(meanG, bCov)
defs      <- dCov

modelMZ   <- mxModel( "MZ", pars, defs, expMeanMZ, corMZ, thresh, dataMZ, expMZ, funML, name="MZ" )  
modelDZ   <- mxModel( "DZ", pars, defs, expMeanDZ, corDZ, thresh, dataDZ, expDZ, funML, name="DZ" )  
multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )

# Create Confidence Interval Objects
ciCor     <- mxCI( c('MZ.corMZ','DZ.corDZ') )
ciThre    <- mxCI( c('MZ.threMZ','DZ.threDZ') )

# Build Saturated Model with Confidence Intervals
modelSAT  <- mxModel( "currSmkSAT", pars, modelMZ, modelDZ, multi, ciCor, ciThre )

# RUN MODEL

# Run Saturated Model
fitSAT    <- mxRun( modelSAT, intervals=F )
summary(fitSAT)


## ACE Model ####################################################################

# Set Starting Values
svPa      <- .8                        # start value for A 
svPc      <- .1                        # start value for C
svPe      <- .1                        # start value for E

# Create Matrices for Variance Components
covA      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VA11", name="VA" ) 
covC      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPc, label="VC11", name="VC" )
covE      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPe, label="VE11", name="VE" )

# Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covP      <- mxAlgebra( expression= VA+VC+VE, name="V" )
covMZ     <- mxAlgebra( expression= VA+VC, name="cMZ" )
covDZ     <- mxAlgebra( expression= 0.5%x%VA+ VC, name="cDZ" )
expCovMZ  <- mxAlgebra( expression= rbind( cbind(V, cMZ), cbind(t(cMZ), V)), name="expCovMZ" )
expCovDZ  <- mxAlgebra( expression= rbind( cbind(V, cDZ), cbind(t(cDZ), V)), name="expCovDZ" )

# Constrain Variance of Binary Variables
var1     <- mxConstraint( expression=diag2vec(V)==1, name="Var1" )

# Create Expectation Objects for Multiple Groups
expMZ     <- mxExpectationNormal( covariance="expCovMZ", means="expMeanMZ", dimnames=depVar, thresholds="thresh" )
expDZ     <- mxExpectationNormal( covariance="expCovDZ", means="expMeanDZ", dimnames=depVar, thresholds="thresh" )
funML     <- mxFitFunctionML()

# Create Model Objects for Multiple Groups
pars      <- list( bCov, meanG, thresh, covA, covC, covE, covP )
defs      <- list( dCov )
modelMZ   <- mxModel( pars, defs, expMeanMZ, covMZ, expCovMZ, dataMZ, expMZ, funML, name="MZ" )
modelDZ   <- mxModel( pars, defs, expMeanDZ, covDZ, expCovDZ, dataDZ, expDZ, funML, name="DZ" )
multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )

# Create Algebra for Unstandardized and Standardized Variance Components
rowUS     <- rep('US',nv)
colUS     <- rep(c('VA','VC','VE','SA','SC','SE'),each=nv)
estUS     <- mxAlgebra( expression=cbind(VA,VC,VE,VA/V,VC/V,VE/V), name="US", dimnames=list(rowUS,colUS) )

# Create Confidence Interval Objects
ciACE     <- mxCI( "US[1,]" )

# Build Model with Confidence Intervals
modelACE  <- mxModel( "currSmkACE", pars, var1, modelMZ, modelDZ, multi, estUS, ciACE )

# ----------------------------------------------------------------------------------------------------------------------
# RUN MODEL

# Run ACE Model
fitACE    <- mxTryHardOrdinal( modelACE, intervals=F )
sumACE    <- summary( fitACE )
sumACE


# Run AE model
modelAE   <- mxModel( fitACE, name="currSmkAE" )
modelAE   <- omxSetParameters( modelAE, labels="VC11", free=FALSE, values=0 )
fitAE     <- mxRun( modelAE, intervals=T )

# Run CE model
modelCE   <- mxModel( fitACE, name="currSmkCE" )
modelCE   <- omxSetParameters( modelCE, labels="VA11", free=FALSE, values=0 )
fitCE     <- mxRun( modelCE )

# Run E model
modelE    <- mxModel( fitAE, name="currSmkE" )
modelE    <- omxSetParameters( modelE, labels="VA11", free=FALSE, values=0 )
fitE      <- mxRun( modelE)

mxCompare(fitACE, list(fitAE,fitCE,fitE))
round(rbind(fitACE$US$result,fitAE$US$result),4)
round(fitAE$output$confidenceIntervals, 3)
