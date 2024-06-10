## MR-DOC model for the effect of Former Smk on DNAm 
## BATCH1
## Madhur Singh, 15-April-2023

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
mxOption( key='Number of Threads', value=8 )

## Libraries
baseDir  <- "/data/msingh/MRDoC_Smk_DNAm/"
phenoDir <- paste0(baseDir,"data/phenos/")
dnamDir  <- paste0(baseDir,"data/DNAm_resid/")
outDir   <- paste0(baseDir,"out/mrdoc_formerSmk_dnam/")
outFile  <- paste0(outDir,"mrdoc_formerSmk_to_DNAm_BATCH1.dat")

## Data #####################

## Load the processed dataset with Smk, age, sex, and residualized PRS (long format)
twin_dat_long <- fread(paste0(phenoDir,"twin_dat_long_mrdoc_resid.csv"), header = T)
twin_dat_long$fisnumber <- as.numeric(twin_dat_long$fisnumber)


## Residualized DNAm 
load(file = paste0(dnamDir,"residualized_beta_blood450k_BATCH1.RData"))
dim(resBeta) 

cpgDat <- as.data.frame(resBeta)
cpgBeta <- colnames(cpgDat)
nCpG <- length(cpgBeta)

rm(resBeta); gc()


# MR-DOC Model #################################################################

## We'll fit AE models for both Smoking and DNAm
## However, the initial script is with ACE models for easier adaptation to other traits or the ACE model
## The C components and rC are fixed at 0 before model-fitting

#Starting Values
x=   1    # x is the standard deviation of the PGS for x
ax= sqrt(.7)
cx= sqrt(.1)
ex= sqrt(.2)    # ax+cx+ex=1
ay= sqrt(.7)
cy= sqrt(.1)
ey= sqrt(.2)    # ay+cy+ey=1
ra= .50   # the larger this is - when unmodeled the more it will inflate g1
rc= .0    # the larger this is - when unmodeled the more it will inflate g1  
re= .0    # the larger this is - when unmodeled the more it will inflate g1 
b1= .1    # larger values give larger genetic instrument strength
b2= .0    # non-zero values are a violation of the non-pleiotropy assumption
g1= -.10  # true causal effect size 
g2= .00   # set to zero ; y->x would mean that there is also reverse causation 

## Means and threshold
frMV      <- c(F,T,T)                  # free status for variables (Smk,CpG,PRS)
svMe      <- c(0,0,0)                  # start value for means
svTh      <- .1                        # start value for thresholds

# Select Binary Exposure (smoking) Variable
expos       <- "formerSmk"            
nExpos      <- 2                     # number of smoking variables
expVar      <- paste(expos,c(1,2),sep="_")

# Select Continuous Outcome (CpG DNAm) Variable
out       <- 'cpg'              # list of variables names
nOut      <- 2                         # number of total variables
outVar    <- paste(out,c(1,2),sep="_")

# Select the genetic IV
iv        <- 'resPRS'
nIV       <- 2
instVars  <- paste(iv,c(1,2),sep="_")

# Total phenotypes
vars    <- c(expos,out,iv)
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
#
ne=18 # restate
ny=6  # restate 
#

# Dimnames for Beta and Psi matrices
matDim <- c(vnames,                                            # 1:6
            paste(c("A","C","E"),expos,"1",sep = "_"),         # 7:9
            paste(c("A","C","E"),out,"1",sep = "_"),           # 10:12
            paste(c("A","C","E"),expos,"2",sep = "_"),         # 13:15  
            paste(c("A","C","E"),out,"2",sep = "_"))           # 16:18


##### MATRIX BETA ##############################################################

mxfreeBE=matrix(F,ne,ne)
# Vars: Smk1 CpG1 PRS1 Smk2 CpG2 PRS2

# Twin 1
mxfreeBE[1,2]=F # 'g2' not freely estimated but set to zero
mxfreeBE[1,3]=T	#'b1'
mxfreeBE[2,1]=T	#'g1'
mxfreeBE[2,3]=T	#'b2'

# Twin 2
mxfreeBE[4, 5]=F  # 'g2' not freely estimated but set to zero
mxfreeBE[4, 6]=T	#'b1'
mxfreeBE[5, 4]=T	#'g1'
mxfreeBE[5, 6]=T	#'b2' 
#
mxfreeBE[1, 7]=T	#'aSmk'
mxfreeBE[1, 8]=T  #'cSmk'
mxfreeBE[1, 9]=T	#'eSmk'
mxfreeBE[2,10]=T	#'aCpG'
mxfreeBE[2,11]=T	#'cCpG' 
mxfreeBE[2,12]=T	#'eCpG'
mxfreeBE[4,13]=T	#'aSmk'
mxfreeBE[4,14]=T 	#'cSmk'
mxfreeBE[4,15]=T	#'eSmk'
mxfreeBE[5,16]=T	#'aCpG'
mxfreeBE[5,17]=T	#'cCpG'
mxfreeBE[5,18]=T	#'eCpG'


mxlabelsBE=matrix(NA,ne,ne)
mxlabelsBE[1,2]='g2'
mxlabelsBE[1,3]='b1'
mxlabelsBE[2,1]='g1'
mxlabelsBE[2,3]='b2'
#
mxlabelsBE[4,5]='g2'
mxlabelsBE[4,6]='b1'
mxlabelsBE[5,4]='g1'
mxlabelsBE[5,6]='b2'
#
mxlabelsBE[1, 7]='aSmk'
mxlabelsBE[1, 8]='cSmk'
mxlabelsBE[1, 9]='eSmk'
mxlabelsBE[2,10]='aCpG'
mxlabelsBE[2,11]='cCpG'
mxlabelsBE[2,12]='eCpG'
mxlabelsBE[4,13]='aSmk'
mxlabelsBE[4,14]='cSmk'
mxlabelsBE[4,15]='eSmk'
mxlabelsBE[5,16]='aCpG'
mxlabelsBE[5,17]='cCpG'
mxlabelsBE[5,18]='eCpG'

# VALUES

mxvaluesBE=matrix(0,ne,ne)
mxvaluesBE[1,2]=g2 
mxvaluesBE[1,3]=b1 
mxvaluesBE[2,1]=g1 
mxvaluesBE[2,3]=b2 
#
mxvaluesBE[4,5]=g2   
mxvaluesBE[4,6]=b1 
mxvaluesBE[5,4]=g1  
mxvaluesBE[5,6]=b2 
#
mxvaluesBE[1, 7]=ax 
mxvaluesBE[1, 8]=cx 
mxvaluesBE[1, 9]=ex 
mxvaluesBE[2,10]=ay 
mxvaluesBE[2,11]=cy  
mxvaluesBE[2,12]=ey
mxvaluesBE[4,13]=ax 
mxvaluesBE[4,14]=cx 
mxvaluesBE[4,15]=ex 
mxvaluesBE[5,16]=ay 
mxvaluesBE[5,17]=cy 
mxvaluesBE[5,18]=ey 




##### MATRIX PSY ################################################################

mxfreePSY=matrix(F,ne,ne) 
mxfreePSY[10, 7]=T # 'ra'
mxfreePSY[13,10]=T # 'ra'
mxfreePSY[16, 7]=T # 'ra'
mxfreePSY[16,13]=T # 'ra'
#
mxfreePSY[11, 8]=T # 'rc'
mxfreePSY[17,14]=T # 'rc'
mxfreePSY[14,11]=T # 'rc'
mxfreePSY[17, 8]=T # 'rc'
#
mxfreePSY[12, 9]=F # re  NOT free  
mxfreePSY[18,15]=F # re  NOT free  

# symmetric fill of upper half of symm matrix
for (i in 1:ne) { for (j in 1:i) { 
  mxfreePSY[j,i]=mxfreePSY[i,j]}}


mxlabelsPSY=matrix(NA,ne,ne) 
mxlabelsPSY[10, 7]='ra'
mxlabelsPSY[13,10]='ra'
mxlabelsPSY[16, 7]='ra'
mxlabelsPSY[16,13]='ra'
mxlabelsPSY[11, 8]='rc'
mxlabelsPSY[14,11]='rc'
mxlabelsPSY[17,14]='rc'
mxlabelsPSY[17, 8]='rc'
mxlabelsPSY[12, 9]='re'  
mxlabelsPSY[18,15]='re' 

# symmetric fill of upper half of symm matrix
for (i in 1:ne) { for (j in 1:i) { 
  mxlabelsPSY[j,i]=mxlabelsPSY[i,j]}}

# VALUES

mxvaluesPSY=matrix(0,ne,ne) 
diag(mxvaluesPSY)=c(0,0,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1)
# fixed
mxvaluesPSY[13, 7]=1 
mxvaluesPSY[14, 8]=1 
mxvaluesPSY[15, 9]=0
mxvaluesPSY[16,10]=1 
mxvaluesPSY[17,11]=1 
mxvaluesPSY[18,12]=0
mxvaluesPSY[ 6, 3]=1   # correlation pgs twin1 - twin2
# free
mxvaluesPSY[10, 7]= ra 
mxvaluesPSY[13,10]= ra 
mxvaluesPSY[16, 7]= ra 
mxvaluesPSY[16,13]= ra 
mxvaluesPSY[11, 8]= rc 
mxvaluesPSY[14,11]= rc
mxvaluesPSY[17, 8]= rc 
mxvaluesPSY[17,14]= rc 
# fixed
mxvaluesPSY[12, 9]= re 
mxvaluesPSY[18,15]= re 

# symmetric fill of upper half of symm matrix
for (i in 1:ne) { for (j in 1:i) { 
  mxvaluesPSY[j,i]=mxvaluesPSY[i,j]}}

# values for DZ ; we will later multiply the MZ-PSY matrix with this one to get DZ-PSY2, 
# hence the initial '1' everywhere
mxvaluesPSY2=matrix(1,ne,ne) 
mxvaluesPSY2[13, 7]= 0.5 
mxvaluesPSY2[13,10]= 0.5  
mxvaluesPSY2[16, 7]= 0.5  # 'ra'
mxvaluesPSY2[16,10]= 0.5  # 'ra'
mxvaluesPSY2[ 6, 3]= 0.5  # correlation pgs twin1 - twin2

for (i in 1:ne) { for (j in 1:i) { 
  mxvaluesPSY2[j,i]=mxvaluesPSY2[i,j]}}



##### MATRIX LAMBDA ############################################################## 

mxfreeLY=matrix(F,ny,ne)
mxfreeLY[3,3]=mxfreeLY[6,6]=T  # pgs variance


mxlabelsLY=matrix(NA,ny,ne)
mxlabelsLY[3,3]= mxlabelsLY[6,6]='x' # pgs variance

# VALUES

mxvaluesLY=matrix(0,ny,ne)
mxvaluesLY[1:ny,1:ny]=diag(ny)
mxvaluesLY[3,3]=mxvaluesLY[6,6]=x # use of a standardized PGS would lead to x = 1


# What is the role of the mysterious matrix mXFilt in
# mxMatrix(type='Full',nrow=5,ncol=6,free=F,values=mxFilt, name='F')?
# This simply removes the row and column relating to the second twin's pgs in the MZ group. 
# We remove this because it is correlated 1 with the first twin's pgs in the MZ group,
# causing the MZ covariance matrix to be singular. 
# If the covariance matrix is singular, OpenMx cannot estimate the maximum likelihood parameters.  
ny1=ny-1
mxFilt=matrix(c(
  1,0,0,0,0,0,
  0,1,0,0,0,0,
  0,0,1,0,0,0,
  0,0,0,1,0,0,
  0,0,0,0,1,0),ny1,ny,byrow=T)


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
# row 1 = pheno 1 (smk)
# row 2 = pheno 2 (cpg)
# row 3 = prs (betas fixed at zero)
bCov <- mxMatrix( type = "Full", nrow = 3, ncol = nCov, 
                  free = c(rep(T,nCov),                      # Smk  - regressed on covariates
                           rep(F,nCov),                      # CpG  - not regressed 
                           rep(F,nCov)),                     # PRS  - not regressed on any covariate (exogenous IV)
                  values = c(rep(0,nCov*3)), byrow = T,
                  labels = c(paste("bSmk",covar,sep="_"),
                             rep(NA,nCov),
                             rep(NA,nCov)), 
                  dimnames = list(vars,covar),
                  name = "bCov" )


# Create Algebra for expected Mean & Threshold Matrices
meanT1     <- mxMatrix( type="Full", nrow=1, ncol=nv, free=frMV, values=svMe, 
                        labels=c("mean_smk","mean_dnam","mean_prs"), dimnames=list("mean",vars), 
                        name="meanT1" )
# In MZ twins, prs_twin1 == prs_twin2
# So, Twin2 in MZ pairs does not have the prs variable
meanT2MZ     <- mxMatrix( type="Full", nrow=1, ncol=(nv-1), free=frMV[1:2], values=svMe[1:2],  
                          labels=c("mean_smk","mean_dnam"), dimnames=list("mean",vars[1:2]), 
                          name="meanT2MZ" )
# In DZ twins, prs_twin1 != prs_twin2
meanT2DZ     <- mxMatrix( type="Full", nrow=1, ncol=nv, free=frMV, values=svMe,  
                          labels=c("mean_smk","mean_dnam","mean_prs"), dimnames=list("mean",vars), 
                          name="meanT2DZ" )

expMeanMZ <- mxAlgebra(expression = cbind(meanT1 + t(bCov%*%dCov1),
                                          meanT2MZ + t(bCov[1:2,] %*%dCov2)
) ,  name="expMeanMZ" )   
expMeanDZ <- mxAlgebra(expression = cbind(meanT1 + t(bCov%*%dCov1),
                                          meanT2DZ + t(bCov%*%dCov2)
) ,  name="expMeanDZ" )   

thresh    <- mxMatrix( type="Full", nrow=1, ncol=2, free=TRUE, values=svTh, labels=c("thr","thr"), name="thresh" )


## Computational settings
funML    <-    mxFitFunctionML()
plan     <-    omxDefaultComputePlan()

## Matrices to store results
mrdocEst  <- matrix(NA_real_, 1, 28)
mrdocSE   <- matrix(NA_real_, 1, 28)
mrdocZ    <- matrix(NA_real_, 1, 28)
mrdocPval <- matrix(NA_real_, 1, 5)  # only computing for g1, b1, b2

namesGoF <- c("statusCode","nStats","nObs","nPars","nConst","df","minus2LL",
              "AICpar","AICsample","BICpar","BICsample","wallTime")
nGoF     <- length(namesGoF)
mrdocFit  <- matrix(NA_real_, 1, 2*nGoF)



## START LOOP #######################################################################

t1 <- Sys.time()

for (i in 1:nCpG) {
  
  # Extract the CpG column
  cpgT <- colnames(cpgDat)[i]
  cpgT
  
  #### Read in the CpG Data ##########################################################
  
  smk_cpg_i <- as.data.frame(cbind(IdatName = rownames(cpgDat),
                                   cpg = cpgDat[, colnames(cpgDat) %in% cpgT ]))
  smk_cpg_i$cpg <- as.numeric(smk_cpg_i$cpg)
  smk_cpg_i$cpg <- as.numeric(scale(smk_cpg_i$cpg))
  
  # Bind with DNAm data
  twin_dat_mrdoc <- merge(twin_dat_long, smk_cpg_i, by = "IdatName", all.x = T)
  
  twin_dat_mrdoc$FamilyNumber <- as.numeric(twin_dat_mrdoc$FamilyNumber)
  twin_dat_mrdoc$Twin_nr <- as.numeric(twin_dat_mrdoc$Twin_nr)
  twin_dat_mrdoc <- twin_dat_mrdoc[order(twin_dat_mrdoc$FamilyNumber,
                                         twin_dat_mrdoc$Twin_nr),]
  
  ## Pivot to wide format
  twin_dat <- reshape(
    data = twin_dat_mrdoc,  direction = "wide",
    v.names = c("fisnumber", "IdatName", "age", "female", "resPRS", 
                "everSmk", "currentSmk", "formerSmk", "smkCess",
                "ncigday", "nyearstop", "cpg", "Platform_2" , "Platform_3",
                "PC1_sc", "PC2_sc", "PC3_sc", "PC4_sc", "PC5_sc", 
                "PC6_sc","PC7_sc","PC8_sc","PC9_sc","PC10_sc",
                "Sample_Plate", "Array_rownum", "Neut_Perc", "Mono_Perc", "Eos_Perc" ), 
    idvar = c("FamilyNumber","zygo"),
    timevar = "Twin_nr", sep = "_"
  )
  
  
  #### MR-DOC #####################################################################
  twin_dat <- data.frame(twin_dat)
  # Select Data for Analysis
  mzData <- subset(twin_dat, zygo==1, selVars)
  dzData <- subset(twin_dat, zygo==2, selVars)
  # Transform binary Smk pheno into an ordered factor
  mzDataF   <- cbind(formerSmk_1 = mxFactor( x=mzData[,1], levels=c(0:1)), mzData[,2:3],
                     formerSmk_2 = mxFactor( x=mzData[,4], levels=c(0:1)), mzData[,-c(1:4)]) 
  # head(mzDataF)
  dzDataF   <- cbind(formerSmk_1 = mxFactor( x=dzData[,1], levels=c(0:1)), dzData[,2:3],
                     formerSmk_2 = mxFactor( x=dzData[,4], levels=c(0:1)), dzData[,-c(1:4)]) 
  # head(dzDataF)
  
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
                   mxAlgebra( expression= (BE[1,7])^2 + (BE[1,8])^2 + (BE[1,9])^2, name="VSmk" ),
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
  
  ## Constrain C components to be 0
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
  
  ## Save output
  mrdocCoef   <- as.vector(coef(mrdocOut))
  mrdocStdErr <- as.vector(mrdocOut$output$standardErrors)
  mrdocZedd   <- as.vector( mrdocCoef / mrdocStdErr )
  mrdocP      <- 2*( 1-pnorm( abs(mrdocZedd[1:3]) ) )
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
  
  
  #### Fit MR-DoC with free rE ##################################################
  
  mrdoc_rE_Model <- omxSetParameters(mrdocOut, labels = c("b2"), free = F, values = 0)
  mrdoc_rE_Model <- omxSetParameters(mrdoc_rE_Model, labels = c("re"), free = T, values = .1)
  
  # fit the model 
  mrdoc_rE_Out <- mxTryHardOrdinal(mrdoc_rE_Model, exhaustive = F, OKstatuscodes = 0)
  sumMRDOC_rE  <- summary(mrdoc_rE_Out)
  
  ## check for any missing SE 
  mrdocStdErr_re <- as.vector(mrdoc_rE_Out$output$standardErrors)
  if( sum(is.na(mrdocStdErr_re)) > 0 ) {
    ## Exhaustive mxTryHard
    mrdoc_rE_Out <- mxTryHardOrdinal(mrdoc_rE_Out, exhaustive = T, OKstatuscodes = 0)
    sumMRDOC_rE  <- summary(mrdoc_rE_Out)
  }
  
  ## Save output
  mrdocCoef_re   <- as.vector(coef(mrdoc_rE_Out))
  mrdocStdErr_re <- as.vector(mrdoc_rE_Out$output$standardErrors)
  mrdocZedd_re   <- as.vector( mrdocCoef_re / mrdocStdErr_re )
  mrdocP_re      <- 2*( 1-pnorm( abs(mrdocZedd_re[1:2]) ) )
  mrdocGoF_re    <- as.vector(c(mrdoc_rE_Out$output$status$code,
                                sumMRDOC_rE$observedStatistics,
                                sumMRDOC_rE$numObs,
                                sumMRDOC_rE$estimatedParameters,
                                sumMRDOC_rE$constraints,
                                sumMRDOC_rE$degreesOfFreedom,
                                sumMRDOC_rE$Minus2LogLikelihood,
                                sumMRDOC_rE$informationCriteria["AIC:","par"],
                                sumMRDOC_rE$informationCriteria["AIC:","sample"],
                                sumMRDOC_rE$informationCriteria["BIC:","par"],
                                sumMRDOC_rE$informationCriteria["BIC:","sample"],
                                mrdoc_rE_Out$output$wallTime)) 
  
  
  
  ### SAVE RESULTS ##################################################################
  
  if(i==1) {
    
    mrdocEst[1,]  <- c(mrdocCoef,mrdocCoef_re);      rownames(mrdocEst) = colnames(cpgDat)[i]
    mrdocSE[1,]   <- c(mrdocStdErr,mrdocStdErr_re);  rownames(mrdocSE) = colnames(cpgDat)[i]
    mrdocZ[1,]    <- c(mrdocZedd,mrdocZedd_re);      rownames(mrdocZ) = colnames(cpgDat)[i]
    mrdocPval[1,] <- c(mrdocP,mrdocP_re);            rownames(mrdocPval) = colnames(cpgDat)[i]
    mrdocFit[1,]  <- c(mrdocGoF,mrdocGoF_re);        rownames(mrdocFit) = colnames(cpgDat)[i]
    
  } else {
    
    mrdocEst   <- rbind(mrdocEst, c(mrdocCoef,mrdocCoef_re));     rownames(mrdocEst)[i] = colnames(cpgDat)[i]
    mrdocSE    <- rbind(mrdocSE, c(mrdocStdErr,mrdocStdErr_re));  rownames(mrdocSE)[i] = colnames(cpgDat)[i]
    mrdocZ     <- rbind(mrdocZ, c(mrdocZedd,mrdocZedd_re));       rownames(mrdocZ)[i] = colnames(cpgDat)[i]
    mrdocPval  <- rbind(mrdocPval, c(mrdocP,mrdocP_re));          rownames(mrdocPval)[i] = colnames(cpgDat)[i]
    mrdocFit   <- rbind(mrdocFit, c(mrdocGoF,mrdocGoF_re));       rownames(mrdocFit)[i] = colnames(cpgDat)[i]
    
  }
  
  ## colnames
  colnames(mrdocEst)  = paste( c(names(coef(mrdocOut)), 
                                 paste0(names(coef(mrdoc_rE_Out)),"_wRE")), "hat", sep = "_" )
  colnames(mrdocSE)   = paste( c(names(coef(mrdocOut)), 
                                 paste0(names(coef(mrdoc_rE_Out)),"_wRE")), "SE", sep = "_" )
  colnames(mrdocZ)    = paste( c(names(coef(mrdocOut)), 
                                 paste0(names(coef(mrdoc_rE_Out)),"_wRE")), "Z", sep = "_" )
  colnames(mrdocPval) = paste( c(names(coef(mrdocOut))[1:3], 
                                 paste0(names(coef(mrdoc_rE_Out)[1:2]),"_wRE")), "p", sep = "_" )
  colnames(mrdocFit)  = c(namesGoF, paste0(namesGoF,"_wRE"))
  
  ## cbind the matrices
  mrdocRes <- as.data.frame(cbind(mrdocEst,mrdocSE,mrdocZ,mrdocPval,mrdocFit))
  ## cpg ID
  mrdocRes$cpg = colnames(cpgDat)[1:i]
  
  ## Save
  fwrite(mrdocRes, file = outFile, sep = "\t", col.names = T, row.names = F, quote = F)
  #
  if( i%%100 == 0 ) { print(paste("Ran Cpg",i,"of",nCpG)) }
  if( i == nCpG ) { print(paste("Ran Cpg",i,"of",nCpG)) }
  
}

t2 <- Sys.time()
print(paste("Time cost =",(t2-t1)))

## END LOOP #####################################################################
## Save

