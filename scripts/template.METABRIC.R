#! /usr/bin/R

##script to run survival analysis on METABRIC
##NB needs to be run over CWL with METABRIC.cwl
##therefore need to change variables in METABRIC.yml
##also runs in a Docker which should have library(survival) installed
#################
## Preliminary ##
#################
library('survival')
options(scipen=999)

##take command line
args <- commandArgs(trailingOnly = TRUE)

##load METABRIC data, assign to vars
load(args[1])

##NB these are sorted, all sample IDs are same (MB-xxxx)
mbexpr <- METABRIC_df
clin <- Complete_METABRIC_Clinical_Features_Data

##NB these are identical, so use DSS only!
##from https://www.synapse.org/#!Synapse:syn1688369/wiki/27311)
dss <- Complete_METABRIC_Clinical_Survival_Data_DSS
#os <- Complete_METABRIC_Clinical_Survival_Data_OS

##userargs
userargs <- read.table(args[2],header=T)

##genes of interest, matchargs and covaraites:
gois <- strsplit(as.vector(userargs[,1]),",")[[1]]
matchargs <- strsplit(as.vector(userargs[,2]),",")[[1]]
covars <- strsplit(as.vector(userargs[,3]),",")[[1]]

##user args
##iterate over gois
for (goi in 1:length(gois)){

  #############################
  ## Survival Analysis Setup ##
  #############################
  goi <- gois[goi]

  ##iterate over all userargs 2..n
  ##set clin to equal that
  clin1 <- clin

  for(x in 2:length(colnames(userargs))){
    for (y in 1:length(matchargs)){
      argsmatch <- strsplit(matchargs[1],";")[[1]]
      clin1 <- clin1[do.call(argsmatch[2],list(clin1[,argsmatch[1]],argsmatch[3])),]
    }
  }

  ##define what variables to use as covariates
  clin2 <- clin1[,colnames(clin1) %in% covars]

  ##make expression set for goi
  ilmn_probe_goi <- as.vector(ilmn_anno[grep(goi,ilmn_anno[,4]),1])
  expruse <- mbexpr[rownames(mbexpr) %in% ilmn_probe_goi,
                    colnames(mbexpr) %in% rownames(clin2)]

  ##possible that multiple probes per gene
  ##take most highly expressed in absence of other indicators
  expruse <- expruse[grep(max(rowSums(expruse)),rowSums(expruse)),]

  ##similarly define DSS
  survuse <- dss[rownames(dss) %in% rownames(clin2),]

  ###############################
  ## Survival Analysis Running ##
  ###############################
  uniMultiList <- survival_UM_continuous(survuse,expruse)
  write.table(lapply(uniMultiList[[1]][7],round,4),
              paste0(goi,".coxPHmultivar.tab"),
              quote=F,row=T,col=T,sep="\t")

  pdf(paste0("METABRIC","_",types[x],"_",goi,"-median_survival.pdf"))
  plot(survFit, lty = 1:4, xlab = "Years To Event", ylab = "Proportion Surviving",
  lwd = 2,col=c("forestgreen","red"),main=paste0("METABRIC ", types[x]," ",goi," (",ilmn_probe_goi,") median survival"))
  legend(x = 6, y = 1, legend = c(paste0(goi," below median (n=",svd$n[1],")"),paste0(goi," above median (n=",svd$n[2],")")), lty = c(1:4), lwd = 2, cex = 0.8,col=c("forestgreen","red"))
  text(8, 0.7, pval, pos = 4)
  dev.off()

}
