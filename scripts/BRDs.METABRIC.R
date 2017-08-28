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

##remove NA
clinnoNA <- clin2[rownames(clin2) %in% names(unlist(lapply(apply(!is.na(clin2),1,unique),
       function(x){if(length(x)==1){return(x)}}))),]

##make expression set for goi
ilmn_probe_goi <- as.vector(ilmn_anno[grep(goi,ilmn_anno[,4]),1])
expruse <- mbexpr[rownames(mbexpr) %in% ilmn_probe_goi,
          colnames(mbexpr) %in% rownames(clinnoNA)]

##possible that multiple probes per gene
##take most highly expressed in absence of other indicators
expruse <- expruse[grep(max(rowSums(expruse)),rowSums(expruse)),]

##similarly define DSS
survuse <- dss[rownames(dss) %in% rownames(clinnoNA),]

###############################
## Survival Analysis Running ##
###############################
uniMultiList <- survival_UM_continuous(survuse, expruse, clinnoNA)

##write out table
write.table(summary(uniMultiList[[2]])[7],
	paste0(goi,".coxPHmultivar.tab"),
    	quote=F,row=T,col=NA,sep="\t")

##KM plot on uniMultiList[[1]], univar
pdf(paste0(goi,".METABRIC.medianExpr_survival.pdf"),onefile=F)

plot(uniMultiList[[3]],
lty = 1:4,
xlab = "Days To Event",
ylab = "Proportional Survival",
lwd = 2,
col=c("forestgreen","red"),
main=paste0("METABRIC ",goi," (",rownames(expruse),") median survival"))

legend(x="topright",
legend = c(paste0(goi," below median (n=", uniMultiList[[3]]$n[1],")"),
         paste0(goi," above median (n=", uniMultiList[[3]]$n[2],")")),
	 lty = c(1:4), lwd = 2, cex = 0.8,
	 col=c("forestgreen","red"))

text(4250, 0.85, strsplit(capture.output(uniMultiList[[1]])[7]," ")[[1]][8], pos = 4)

dev.off()
}
