#############################
## READ ME
#This code takes logTPM file, 
#1) remove duplicate samples
#2) calculate all possible combinations of mid-parent value (MPV)
#3) perform Kolmogorov-Smirnov test on log2 (TPM+1) and log2 (MPV+1), then output test values
############################
#Config
library(scales)
wd<-"/ebio/abt6_projects9/SigHeterosis_Batch1/data/1_analysis/4_DEanalysis/run125_rmCMRTPHBS_seqtk"
metaname<-"run125_RSEM_rmCMRTPHBS_rmColOL3Mdup_seqtk28_Shootmetainfo.txt"
dataname<-"run125_RSEM_rmCMRTPHBS_rmColOL3Mdup_seqtk28_ShootTPM.txt"
outbasename<-"run125_rmCMRTPHBS3.5dup_seqtk28_Shoot"
drop=80 #threshold of max dropout samples
varcoefcut=0.2 #threshold of min coefficient of variance
lowexpcut=0.5 #threshold of min average expression (trimmed mean on log scale)
trimperc=0.20 #the fraction trimmed when calculating trimmed mean

##############
setwd(wd)
meta<-read.table(metaname,header=T,sep="\t",row.names=1)
TPM<-read.table(dataname,header=T,sep="\t",row.names=1)

#remove S97 from the files, leaving only the seqtk subsampled one
maxdepth<-which(meta$ttlMapped.rsem==max(meta$ttlMapped.rsem))
meta<-meta[-maxdepth,]
TPM<-TPM[,-maxdepth]

TPMdrop<-TPM[apply(TPM,1,function(x){sum(x==0)})<=drop,]
TPMdrop<-TPMdrop[apply(TPMdrop,1,function(x){!any(x==Inf)}),]
TPMvar<-TPMdrop[apply(TPMdrop,1,function(x){(sd(x)/mean(x))>=varcoefcut}),]
TPMvarlog<-as.data.frame(apply(TPMvar,c(1,2),function(x){log2(x+1)}))#transform into log2(TPM+1)
trimmean<-function(x,perc){
  xrank<-sort(x)
  extremepercent<-floor(length(x)*perc)
  trimmedx<-xrank[extremepercent:(length(x)-extremepercent)]
  trimmedmean<-mean(trimmedx)
  return(trimmedmean)
}
TPMvarlog<-TPMvarlog[apply(TPMvarlog,1,trimmean,perc=trimperc)>=lowexpcut,]
print(dim(TPMvarlog))
print(dim(meta))
TPMvar<-TPMvar[rownames(TPMvar)%in%rownames(TPMvarlog),] #lazy work-around, index genes on untransformed scale that fits the filtering criteria


## calculate pairwise MPV
TPM2use=TPMvar
inbredidx<-which(meta$IsHybrid=="Inbred")
inbredcomb<-t(combn(inbredidx,2))
calcMPV<-function(x,TPMtable){
  parent1<-x[1]
  parent2<-x[2]
  par1TPM<-TPMtable[,parent1]
  par2TPM<-TPMtable[,parent2]
  #print(c(colnames(TPMtable)[parent1],colnames(TPMtable)[parent2]))
  pairMPV<-(par1TPM+par2TPM)/2
  return(pairMPV)
}
allpossibleMPV<-as.data.frame(apply(inbredcomb,1,calcMPV,TPMtable=TPM2use))
#testpossibleMPV<-apply(inbredcomb[1:4,],1,calcMPV,TPMtable=TPMtable[1:10,])
rownames(allpossibleMPV)<-rownames(TPM2use)
colnames(allpossibleMPV)<-paste0("Comb",inbredcomb[,1],"_",inbredcomb[,2])

## calculate Kolmogorov-Smirnov statistic
#Compare distribution of each gene's expression in hybrids vs. that of all possible combination of Mid-Parental Value (MPV).
#Test if the hybrids are drawn from the same distribution as the distribution of MPVs. 
#Significant deviation from MPV indicates some form of dominance in expression pattern

function_ks<-function(x,Hybgroup,MPVgroup){
  hyb<-as.numeric(log2(Hybgroup[x,]+1)) #log2 tranform here
  mpv<-as.numeric(log2(MPVgroup[x,]+1))
  kspergene<-ks.test(hyb,mpv,alternative="two.sided")  
  ksstat<-kspergene$statistic
  kspval<-kspergene$p.value
  returnvec<-c(ksstat,kspval)
  #if(x%%100==1){print(x)}
  print(x)
  return(returnvec)
}

HvMPV_KS<-t(sapply(1:nrow(TPM2use),function_ks,Hybgroup=TPM2use[,which(meta$IsHybrid=="Hybrid")], MPVgroup=allpossibleMPV))
TPM2use[,"HvMPV_KS.D"]<-HvMPV_KS[,1]
TPM2use[,"HvMPV_KS.p"]<-HvMPV_KS[,2]
TPM2use[,"HvMPV_KS.fdr"]<-p.adjust(HvMPV_KS[,2],"BH")
TPMoutname<-paste0(outbasename,"_TPM_HvMPVtest.txt")
write.table(TPM2use,file=file.path(wd,TPMoutname),quote=F,sep="\t")
q0.05<-sum(TPM2use$HvMPV_KS.fdr<0.05) #number of non-additive genes with q<0.05
q0.01<-sum(TPM2use$HvMPV_KS.fdr<0.01) #number of non-additive genes with q<0.01
q0.001<-sum(TPM2use$HvMPV_KS.fdr<0.001) #number of non-additive genes with q<0.001
#Kolmogorov-Smirnov statistic appended to the log TPM file.

#Output plots:
#Rank genes by Kolmogorov-Smirnov stat and plot overall distribution along the rank

pdfname<-paste0(outbasename,"_HvMPV_example_distribution.pdf")
Inbredcol<-alpha("purple4",0.6)
Hybridcol<-alpha("turquoise4",0.6)
Hybgroup<-TPM2use[,which(meta$IsHybrid=="Hybrid")]
MPVgroup<-allpossibleMPV

#generate a list of coordinate to plot Hybrid vs MPV distribution
lastthousand<-ifelse(round(nrow(TPM2use),digit=-3)<=nrow(TPMvar),round(nrow(TPM2use),digit=-3),(round(nrow(TPM2use),digit=-3)-1000))#the last full thousand for the table
lasthundred<-ifelse(round(nrow(TPM2use),digit=-2)<=nrow(TPMvar),round(nrow(TPM2use),digit=-2),(round(nrow(TPM2use),digit=-2)-100))
lastten<-ifelse(round(nrow(TPM2use),digit=-1)<=nrow(TPMvar),round(nrow(TPM2use),digit=-1),(round(nrow(TPM2use),digit=-1)-10))

#plotsamples<-c(1:9,seq(10,490,10),seq(500,8500,500),seq(9000,(lastthousand-1),1000),seq(lastthousand,(lasthundred-1),100),seq(lasthundred,(lastten-1),10),lastten:nrow(TPM2use))
plotsamples<-c(1:9,seq(10,490,10),500,seq(1000,(lastthousand-1),1000),seq(lastthousand,(lasthundred-1),100),seq(lasthundred,(lastten-1),10),lastten:nrow(TPM2use))

KSfdr<-p.adjust(HvMPV_KS[,2],"BH")
orderbyKS<-order(KSfdr)
plot_orderbyKS<-orderbyKS[plotsamples]

pdf(file.path(wd,pdfname),width=30,height=36)
par(mfrow=c(6,5))
for (i in 1:length(plotsamples)){
   cat(i)
   whichrow<-plot_orderbyKS[i]
   genename<-rownames(TPM2use)[whichrow]
   KSP<-KSfdr[whichrow]
   KSD<-HvMPV_KS[whichrow,1]
   hyb<-as.numeric(log2(Hybgroup[whichrow,]+1))
   mpv<-as.numeric(log2(allpossibleMPV[whichrow,]+1))
   densmpv<-density(mpv)
   denshyb<-density(hyb)
   hist(mpv,breaks=100,col=Inbredcol,border=Inbredcol,freq=F,ylim=c(0,max(densmpv$y,denshyb$y)*2.5),xlim=c(min(densmpv$x,denshyb$x),max(densmpv$x,denshyb$x)),main=paste0(genename,":k.s.stat=", round(KSD,digits=4),", k.s.fdr=", round(KSP,digits=4),", rank=",plotsamples[i]))
   hist(hyb,breaks=100,col=Hybridcol,border=Hybridcol,freq=F,add=T)
   lines(densmpv,lwd=2,col="purple4")
   lines(denshyb,lwd=2,col="turquoise4")
}

dev.off()






