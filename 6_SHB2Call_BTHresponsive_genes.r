#####################
#README
#The script takes minimally filtered logTPM file, constructs LMM for each gene
#and performs permutation to identify genes whose expression variation across samples can be significantly attributed to BTH treatment
  #per-gene permutation threshold is established, do a *two-sided* test
  #rankings of test statistic is transformed into p-value, and FDR-adjusted
  #the separately-ran permutation script is included below in a separate section
#####################
#Config
library(RColorBrewer)
library(scales)
library(vioplot)
library(lme4)
wd<-"/ebio/abt6_projects7/SHB2/data/1_analysis/6_BTHdiffgenes/"
outbasename<-"run197_RowFilter_BTHlmm_"
triocol<-colorRampPalette(brewer.pal(8,"Set2"))
metafile<-"Run197_rmCMRTPHBS_rm2MSeqtkDupOL_corr.meta.txt"
logTPMfile<-"run197_BTHdiffgene_rowfiltered_logTPM.txt"
additivegenefile<-"run197_RowFilter_addpredict_slp0.4Sigma0.65_geneID.txt"
######################
#Data input
setwd(wd)
meta<-read.table(metafile,header=T,sep="\t",row.names=1)
logTPM<-read.table(logTPMfile,header=T,sep="\t",row.names=1)
AdditiveGenes<-read.table(additivegenefile,header=F,sep="\t")[,1]

meta$Treatment<-factor(meta$Treatment,levels=c("Mock","BTH"))
meta$IsHybrid<-factor(meta$IsHybrid,levels=c("Inbred","Hybrid"))
#######################
#use lmm and oermutation to identify genes whose expression can be significantly attributed to BTH treatment
#Ran once
Treatlmm<-data.frame()
for(g in 1:nrow(logTPM)){
  gene<-as.numeric(logTPM[g,])
  fit<-lmer(gene~Treatment+IsHybrid+Treatment*IsHybrid+(1|PlantBatch),data=meta,REML=F)
  returnvec<-c(fixef(fit),summary(fit)$coefficients[,3], summary(fit)$sigma)
  Treatlmm<-rbind(Treatlmm,returnvec)
}

names(Treatlmm)<-c("Mock","BTH","Hybrid","BTHxHybrid","t-Mock","t-BTH","t-Hybrid","t-BTHxHybrid","Sigma")
rownames(Treatlmm)<-rownames(logTPM)

write.table(Treatlmm,paste0(outbasename,"SimpleModel_FixEffect.txt"),quote=F,sep="\t")

#################################################
#################################################
##PERMUTATION SCRIPT, TO BE RUN INDEPENDENTLY
#each single script perform 100X permutation, 
#with a shell script to automate 100X iteration of such process, for a total of 10k permutation
#################
args<-commandArgs(trailingOnly = T)
idx<-args[1]
#perform permutation of scalelogTPM to test for BTH treatment effect
library(lme4)
setwd("/ebio/abt6_projects7/SHB2/data/1_analysis/6_BTHdiffgenes/")
meta<-read.table("Run197_rmCMRTPHBS_rm2MSeqtkDupOL_corr.meta.txt",header=T,sep="\t",row.names=1)
logTPM<-read.table("run197_BTHdiffgene_rowfiltered_logTPM.txt",header=T,sep="\t",row.names=1)
meta$Treatment<-factor(meta$Treatment,levels=c("Mock","BTH"))
meta$IsHybrid<-factor(meta$IsHybrid,levels=c("Inbred","Hybrid"))
outbasename<-"Run197_BTHlmm_perm_"
print(dim(logTPM))
print(head(meta))

function_permTreatlmm<-function(x){
  #scramble per gene for its association with meta data
  permgene<-sample(as.numeric(x),length(x),replace = F)
  fit<-lmer(permgene~Treatment+IsHybrid+Treatment*IsHybrid+(1|PlantBatch),data=meta,REML=F)
  tfixed<-summary(fit)$coefficients[,3] #get the t-statistic for each fixed effect term
  return(tfixed) #instead of returning one value, this will be four values, for mock, BTH, mock*hybrid, and BTH*hybrid
}

#permutation
print("starting permutation, round:")
permEffBTH<-matrix(nrow=nrow(logTPM),ncol=400)
#timer<-c()
for (i in 1:100){
  startT<-Sys.time()
  print (i)
  permTreatlmm<-apply(logTPM,1,function_permTreatlmm)
  permTreatlmm<-t(permTreatlmm)
  permEffBTH[,(4*i-3):(4*i)]<-permTreatlmm
  endT<-Sys.time()
  timer<-c(timer,(endT-startT))
}

permEffBTH<-as.data.frame(permEffBTH)
rownames(permEffBTH)<-rownames(logTPM)
write.table(permEffBTH,file=paste0(outbasename,idx,"_nullt-stats.txt"), quote=F,sep="\t")
#########################################################################
#########################################################################

#integrating permutation result and FDR
##loop through all 100 perm files, and record the rank.minus1 and rank.plus1 for each gene in each file 

BTHrank<-matrix(nrow=nrow(Treatlmm),ncol=200)
permrnd=10000
for(i in 1:100){
  print (paste0("i=",i))
  permranks<-matrix(nrow=nrow(Treatlmm),ncol=2)
  perm<-paste0("Run197_BTHlmm_perm_",i,"_nullt-stats.txt")
  print(perm)
  con<-file(perm)
  open(con)
  geneperm<-readLines(con,n=1) #get rid of the header
  for(j in 1:nrow(Treatlmm)){
    geneBTH<-Treatlmm[j,]$t.BTH #t-stat for BTH term
    geneperm<-readLines(con,n=1)
    geneperm<-unlist(strsplit(geneperm,"\t"))
    genename<-geneperm[1]
    geneperm<-as.numeric(geneperm[2:length(geneperm)])
    #this permutation result include t-stats for Mock, BTH, MockxHybrid, and BTHxHybrid. 
    geneperm<-geneperm[seq(2,length(geneperm),4)]#For our current purpose, only t-stats for BTH is useful (aka every 2nd column in groups of 4)
    
    #remove NAs from permutation, but count them as more extreme than the observed value
    sumNA<-sum(is.na(geneperm))
    if(sumNA!=0){print(c(j,genename,sumNA))}
    geneperm<-geneperm[!is.na(geneperm)]
    rank.minus1<-ifelse(geneBTH<0, sum(geneperm<geneBTH)+sumNA,sum(geneperm<geneBTH)) #add sumNA into the extreme rank
    rank.plus1<-ifelse(geneBTH>0,sum(geneperm>geneBTH)+sumNA,sum(geneperm>geneBTH))
    if(j%%1000==1){print(c(j,genename,rank.minus1,rank.plus1),sep=" ")}
    permranks[j,1]<-rank.minus1
    permranks[j,2]<-rank.plus1
  }
  close(con)
  BTHrank[,2*i-1]<-permranks[,1]
  BTHrank[,2*i]<-permranks[,2]
}

##calculate the overall plus and minus rank, and subsequently, p-value
ttl.rank.minus1<-apply(BTHrank[,seq(1,ncol(BTHrank),2)],1,sum)
ttl.rank.plus1<-apply(BTHrank[,seq(2,ncol(BTHrank),2)],1,sum)
##two-sided P-value
pvalBTH<-sapply(1:nrow(Treatlmm),function(x){ifelse(Treatlmm$t.BTH[x]>0,2*(ttl.rank.plus1[x]+1)/(permrnd+1),2*(ttl.rank.minus1[x]+1)/(permrnd+1))})
pvalBTH<-ifelse(pvalBTH>1,1,pvalBTH)

######################
#REMOVE ADDITIVE GENES
Treatlmm[,"IsAdditive"]<-ifelse(rownames(Treatlmm)%in%AdditiveGenes,1,0) #additive =1, not additive =0 (slp0.4,sigma 0.66)
Treatlmm_rmadd<-Treatlmm[Treatlmm$IsAdditive==0,]
Treatlmm_rmadd[,"padjBTH_rmadd"]<-p.adjust(Treatlmm_rmadd$Perm10k_pvalBTH,method="BH")
write.table(Treatlmm_rmadd,paste0(outbasename,"rmadd_SimpleModel_FixEffect.txt"),quote=F,sep="\t")
