## Constructing the allelic score from mQTLs
## Madhur Singh

rm(list = ls()); gc()

## Libraries
library(data.table)
setDTthreads(12)
library(dplyr)
library(tidyr)

baseDir   <- "/base_dir/" 
dataDir   <- paste0(baseDir,"data/GoDMC/")
godmcDir  <- paste0(baseDir,"data/GoDMC/batches/")
phenoDir  <- paste0(baseDir,"data/phenos/")
cisDir    <- paste0(dataDir, "cpg_mqtls/cis_mqtls/")
qcGenoDor <- paste0(baseDir, "out/geno/")
prsDir    <- paste0(baseDir, "out/prs_mqtls/")
scoreDir  <- paste0(baseDir, "out/prs_mqtls/scores/")
tempdir   <- paste0(baseDir,"scripts/temp/")
logDir    <- paste0(baseDir,"logs/")

plinkDir  <- "/base_dir/Util/plink1.9/"

## Files
refFile   <- paste0(qcGenoDor,"1KGP_GoDMC_SNPs")
bedFile   <- paste0(qcGenoDor,"MRG16_PRSQCed_GoDMC_matched")
pvalRange <- paste0(prsDir, "range_list_2")
## This has only a single p-val threshold
## Because cohort-level thresholding prior to the meta-analysis GWAS, 
## the p-values are already "kind of" thresholded
## In test runs with a few CpGs, using different p-val thresholds in C&T did not make much difference

## CpGs with cis-mQTLs
cpgIDs <- fread(paste0(phenoDir,"godmc_cis_cpg_ids.dta"),
                header = F, sep = "\t")
nCpG <- length(cpgIDs$V1)
nCpG

snp_per_prs <- data.frame(matrix(NA_real_, nrow = 1, ncol = 3))
colnames(snp_per_prs) <- c("CpG","clumpLDr2","nSNP")

# Start CpG Loop ##################################################################

for (jj in 1:nCpG) {     
  
  cpgT     <- cpgIDs$V1[jj]
  
  cpgDir   <- paste0(prsDir,cpgT,"/")
  sumStats <- paste0(cisDir,cpgT,"_GoDMC_cis_sumStats_noBIOS_QC.dat")
  
  system(paste0("mkdir ",cpgDir))
  
  r2Range <- c(0.5, 0.1)
  nLDr2   <- length(r2Range)
  
  ## Clumping ######################################################################
  
  for (ii in 1:nLDr2) {
    
    r2cutoff <- r2Range[ii]      
    outFile  <- paste0(cpgDir, cpgT, "_godmc_cis_prs_r2_",10*r2cutoff)
    
    dateCheck <- gsub( pattern="-", replacement="", x=Sys.Date() )
    jobName   <- paste0("plink_clump_cis_r2_",10*r2cutoff,"_",dateCheck)
    shellName <- paste0(cpgDir,jobName,".sh",sep="")
    logFile   <- paste0(cpgDir,jobName)
    
    plinkCom <- paste0(
      plinkDir,"plink \\
--bfile ", refFile, " \\
--clump-p1 1 \\
--clump-kb 250 \\
--clump-r2 ", r2cutoff, " \\
--clump ", sumStats," \\
--clump-snp-field kgpSNP \\
--clump-field pval \\
--threads 12 \\
--out ", outFile
    )
    
    write(plinkCom, file = shellName)
    system(paste0("chmod +x ",shellName))
    system(paste0("nohup ",shellName," > ",logFile,".log 2> ",logFile,".err &"))
    Sys.sleep(2)
    
    print(paste(ii,"of",nLDr2,"clumping job submitted for", cpgT))
    
  }  
  
  ## END Clumping Loop
  
  
  # Scoring #######################################################################
  
  for (ii in 1:nLDr2) {
    
    r2cutoff <- r2Range[ii]  
    outFile  <- paste0(cpgDir, cpgT, "_godmc_cis_prs_r2_",10*r2cutoff)
    outFile1  <- paste0(scoreDir, cpgT, "_godmc_cis_prs_LDr2_",10*r2cutoff,"_p")
    
    ## get the clumped KGP SNPid
    ## Col 3 in the *.clumped file
    awkCom <- paste0("awk 'NR!=1{print $3}' ", outFile, ".clumped >  ", outFile,".snp")
    system(awkCom)
    
    ## get the corresponding NTR SNPid
    gwasSum   <- fread(sumStats, data.table = F)
    clumpSNP  <- fread(paste0(outFile,".snp"), data.table = F, header = F)
    
    ## For the PRS, we need to specify the col numbers in this order:
    ## 2 SNP ID
    ## 9 Effect Allele
    ## 5 Effect Size
    
    ## filter Sum stats for the clumped SNPs
    clumpSum  <- gwasSum |> 
      filter(kgpSNP %in% clumpSNP$V1)
    
    ## Save the clumped NTR SNPid
    snpID  <- clumpSum |> 
      select(SNP)
    
    fwrite(snpID, file = paste0(outFile,".valid.snp"), sep = "\t",
           quote = F, row.names = F, col.names = F )
    
    ## Save the SNPid + pval
    fwrite(clumpSum[,c("SNP","pval")], file = paste0(outFile,".snp.pval"), sep = "\t",
           quote = F, row.names = F, col.names = T )
    
    ## Save the number of SNPs in the PRS
    nSNP <- clumpSum |> 
      filter(pval < 0.05) |> 
      nrow()
    
    if(ii*jj == 1) {
      snp_per_prs[ii*jj,] <- c(cpgT,r2cutoff,nSNP)
    } else {
      snp_per_prs <- rbind( snp_per_prs, c(cpgT,r2cutoff,nSNP) )
    }

    ## PLINK 
    dateCheck <- gsub( pattern="-", replacement="", x=Sys.Date() )
    jobName   <- paste0("plink_cis_score_LDr2_",10*r2cutoff,"_",dateCheck)
    shellName <- paste0(cpgDir,jobName,".sh",sep="")
    logFile   <- paste0(cpgDir,jobName)
    snpPval   <- paste0(outFile,".snp.pval")
    validSNP  <- paste0(outFile,".valid.snp")
    
    ## Check the min P-val
    minP <- min(clumpSum$pval)
    
    ## Apply p-val threshold only if the min p-val is <0.05
    ## Else Keep all the SNPs that we have in the Sum Stats
    ## The derived PRS will likely not perform very well in the latter case
    ## But we can still compute it and get the R2
    
    if(minP < 0.05) {
      
      plinkCom <- paste0(
        plinkDir,"plink \\
--bfile ", bedFile, " \\
--score ", sumStats," 2 9 5 header \\
--q-score-range ", pvalRange," ", snpPval," \\
--extract ", validSNP," \\
--threads 12 \\
--out ", outFile1
      )
      
    } else {
      
      ## Plink command without p-val threshold
      plinkCom <- paste0(
        plinkDir,"plink \\
--bfile ", bedFile, " \\
--score ", sumStats," 2 9 5 header \\
--extract ", validSNP," \\
--threads 12 \\
--out ", outFile1
      )
      
    }
  
    write(plinkCom, file = shellName)
    system(paste0("chmod +x ",shellName))
    system(paste0("nohup ",shellName," > ",logFile,".log 2> ",logFile,".err &"))
    
    Sys.sleep(1)
    
    print(paste(ii,"of",nLDr2,"scoring job submitted for",cpgT))
    
  }  
  
  ## END scoring Loop
  
  
  # Top mQTL score ##################################################################
  
  topMQTL <- gwasSum |> 
    arrange(pval) |> 
    slice(1)
  
  topMQTLid   <- topMQTL$SNP
  topMQTLpval <- topMQTL |> select(SNP, pval)
  
  ## Save the SNP/mQTL ID
  write(topMQTLid, file = paste0(cpgDir, cpgT, "_godmc_top_cis_mQTL.snp"))
  ## Save the SNPid + pval
  fwrite(topMQTLpval, file = paste0(cpgDir, cpgT, "_godmc_top_cis_mQTL.snp.pval"), sep = "\t",
         quote = F, row.names = F, col.names = T )
  
  
  outFile2  <- paste0(scoreDir, cpgT, "_godmc_top_cis_mQTL_score")
  
  dateCheck <- gsub( pattern="-", replacement="", x=Sys.Date() )
  jobName   <- paste0("plink_score_top_cis_mQTL_",dateCheck)
  shellName <- paste0(cpgDir,jobName,".sh",sep="")
  logFile   <- paste0(cpgDir,jobName)
  mqtlSNP   <- paste0(cpgDir, cpgT, "_godmc_top_cis_mQTL.snp")
  
  plinkCom <- paste0(
    plinkDir,"plink \\
--bfile ", bedFile, " \\
--score ", sumStats," 2 9 5 header \\
--extract ", mqtlSNP," \\
--threads 12 \\
--out ", outFile2
  )
  
  write(plinkCom, file = shellName)
  system(paste0("chmod +x ",shellName))
  system(paste0("nohup ",shellName," > ",logFile,".log 2> ",logFile,".err &"))
  
  print(paste("PRS calculation jobs submitted for cpg", jj,"of",nCpG))
  
  Sys.sleep(1)
  
  ## Save the SNP_per_PRS file
  if(jj == nCpG) {
    snp_per_prs$nSNP <- as.numeric(snp_per_prs$nSNP)
    fwrite(snp_per_prs, file = paste0(prsDir,"n_snps_p_0.05_by_clumpLDr2.dta"), 
           sep = "\t", quote = F, col.names = T, row.names = F)
  }
 
}

########### END ##################

