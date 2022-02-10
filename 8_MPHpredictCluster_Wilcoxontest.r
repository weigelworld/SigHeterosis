##################
#README
#the script takes spline-based cluster information (61 clusters), 
#and performs Wilcoxon signed rank sum test (non parametric, paired) on the ends of the fitted cluster spline,
#to test whether the clusters are correctly assigned to the super cluster 
#(i.e. whether MPH values in the groups at either end of the spline are significantly different)

#test strategy:
#calculate mean gene expression and size for each trio*treatment across all reps
#find the two ends of expression values (or the median for the quad test), and perform paired test for each cluster (pairing are the genes)
#Use the 10% from each tail of the spline (and the 10% at the median value) 
##################
#Configuration
library(vioplot)
library(scales)
wd<-"/ebio/abt6_projects7/SHB2/data/1_analysis/8_MPHsplineRanksumTest/"
logTPMfile<-"run197_BTHdiffgene_rowfiltered_logTPM.txt"
metafile<-"Run197_rmCMRTPHBS_rm2MSeqtkDupOL_corr.meta.txt"
MPHidxfile<-"run197_BTHdiffgene_MetaRowIdx_forpairedMPHcalc_manualcorr.txt"
clusterinfo<-"/ebio/abt6_projects7/SHB2/data/1_analysis/6_BTHdiffgenes/run197_RowFilter_BTHlmm_rmadd_fdr0.001_rmflat_MPHpredsplinek2-200clustID.txt"
##################

setwd(wd)
logTPM<-read.table(logTPMfile,header=T,sep="\t",row.names=1)
meta<-read.table(metafile,header=T,sep="\t",row.names=1)
MPHidx<-read.table(MPHidxfile,header=T,sep="\t")
clusterfile<-read.table(clusterinfo,header=T,sep="\t",row.names = 1)


meta$Treatment<-factor(meta$Treatment,levels=c("Mock","BTH"))
meta$IsHybrid<-factor(meta$IsHybrid,levels=c("Inbred","Hybrid"))
MPHidxM<-which(MPHidx$Treatment=="Mock") # these are just row numbers, used for indexing/sorting/spline fitting in the subsequent step
MPHidxB<-which(MPHidx$Treatment=="BTH")
sublogTPM<-logTPM[rownames(logTPM)%in%rownames(clusterfile),] #subset expression values to only the genes who participated in clustering
scalelogTPM<-t(scale(t(sublogTPM)))

#calculate MPH for each gene in each trio
calc_trioxpnMPH<-function(idx,gene){
  fxp<-gene[as.numeric(idx[7])]
  mxp<-gene[as.numeric(idx[8])]
  hxp<-gene[as.numeric(idx[9])]
  mphxp<-(hxp-0.5*(fxp+mxp))
}
xpnMPH<-t(apply(scalelogTPM,1,function(x){apply(MPHidx,1,calc_trioxpnMPH,gene=x)}))

clustID<-clusterfile$K61

#index genes according to supercluster information
clustergenelist<-list()
clustergenelist$"flat-flat"<-c(24,28,30,32,40,43,53,59)
clustergenelist$"flat-pos"<-c(16,23,45)
clustergenelist$"flat-neg"<-c(7,22,41,58,18,19,26,38,54,56)
clustergenelist$"pos-flat"<-c(4,6,34,44,50,55)
clustergenelist$"pos-pos"<-c(3,20,46,48,49,52)
clustergenelist$"pos-neg"<-c(35,51)
clustergenelist$"neg-flat"<-c(8,13,42,47)
clustergenelist$"neg-pos"<-c(10)
clustergenelist$"neg-neg"<-c(5,11,21,29,31,33,36,37,39,60)
clustergenelist$"quad-flat"<-c(9,17,25,27,57,61)
clustergenelist$"quad-pos"<-c(1,12,15)
clustergenelist$"quad-neg"<-c(2,14)



#acquire average lines of individual genes in each cluster
trioszmean<-tapply(MPHidx$MPHmm2,paste0(MPHidx$Treatment,".",MPHidx$TrioID),mean)
Mockszmean<-trioszmean[grep("Mock",names(trioszmean))]
BTHszmean<-trioszmean[grep("BTH",names(trioszmean))]

#dataframe to record statistics
Wilcox.output<-data.frame()

#######
#1. flat-flat
########
pdf("run197_BTHdiffgene_MPHpredict_1.flat-flat_windowvioplot.pdf",height=10,width=20)
par(mfrow=c(2,4))
for(i in clustergenelist[[1]]){
  print(i)
  genelist<-xpnMPH[clustID==i,]
  #calculate trio*treatment mean
  trioxpnmean<-t(apply(genelist,1,function(x){tapply(x,paste0(MPHidx$Treatment,".",MPHidx$TrioID),mean)})) #rows are the genes within the chosen cluster, cols are mean genlist of the trio
  #separate mock and BTH
  Mockxpnmean<-trioxpnmean[,grep("Mock",colnames(trioxpnmean))] #the mean genlist for each mock trio and each gene
  BTHxpnmean<-trioxpnmean[,grep("BTH",colnames(trioxpnmean))] 
  
  #acquire average value in an extreme window (data points=4)
  #Mock
  minxM<-apply(Mockxpnmean,1,function(x){mean(sort(x)[1:4])}) #min expression 
  maxxM<-apply(Mockxpnmean,1,function(x){mean(sort(x)[(length(x)-3):length(x)])}) #max expression
  minyM<-apply(Mockxpnmean,1,function(x){mean(Mockszmean[order(x)[1:4]])}) # size corresponding to min expression
  maxyM<-apply(Mockxpnmean,1,function(x){mean(Mockszmean[order(x)[(length(x)-3):length(x)]])}) # size corresponding to max expression
  #BTH
  minxB<-apply(BTHxpnmean,1,function(x){mean(sort(x)[1:4])})
  maxxB<-apply(BTHxpnmean,1,function(x){mean(sort(x)[(length(x)-3):length(x)])})
  minyB<-apply(BTHxpnmean,1,function(x){mean(BTHszmean[order(x)[1:4]])})
  maxyB<-apply(BTHxpnmean,1,function(x){mean(BTHszmean[order(x)[(length(x)-3):length(x)]])})
  
  #plot
  vioplot(minyM, minyB, maxyM,maxyB, 
          col=alpha(rep(c("springgreen4","darkorange2"),2),0.4),at=c(1,2,4,5),
          xlim=c(0,6),ylim=c(-20,100),
          names=c("Mock.min","BTH.min","Mock.max","BTH.max"),las=2,
          main=paste0("SplineCluster",i," flat-flat n=",length(minyM)),ylab="Rosette size MPH")
  points(jitter(rep(1,length(minyM)),amount=0.1),minyM,pch=20,col="springgreen4")
  points(jitter(rep(2,length(minyB)),amount=0.1),minyB,pch=20,col="darkorange2")
  points(jitter(rep(4,length(maxyM)),amount=0.1),maxyM,pch=20,col="springgreen4")
  points(jitter(rep(5,length(maxyB)),amount=0.1),maxyB,pch=20,col="darkorange2")
  
  #wilcoxon test
  Mocktest<-wilcox.test(minyM, maxyM, paired = TRUE, alternative = "two.sided")
  BTHtest<-wilcox.test(minyB, maxyB, paired = TRUE, alternative = "two.sided")
  statdf<-data.frame("Cluster"=rep(i,2),"Mode"=rep("flat-flat",2),"Treatment"=c("Mock","BTH"),
                     "ContrastA"=rep("Min.xpn.sz",2),"ContrastB"=rep("Max.xpn.sz",2),"Sided"=rep("two",2),
                     "Stat"=c(Mocktest$statistic,BTHtest$statistic),"p.value"=c(Mocktest$p.value,BTHtest$p.value))
  Wilcox.output<-rbind(Wilcox.output,statdf)
}
dev.off()
########
#2.flat-pos
########
pdf("run197_BTHdiffgene_MPHpredict_2.flat-pos_windowvioplot.pdf",height=10,width=10)
par(mfrow=c(2,2))
for(i in clustergenelist[[2]]){
  print(i)
  genelist<-xpnMPH[clustID==i,]
  #calculate trio*treatment mean
  trioxpnmean<-t(apply(genelist,1,function(x){tapply(x,paste0(MPHidx$Treatment,".",MPHidx$TrioID),mean)})) #rows are the genes within the chosen cluster, cols are mean genlist of the trio
  #separate mock and BTH
  Mockxpnmean<-trioxpnmean[,grep("Mock",colnames(trioxpnmean))] #the mean genlist for each mock trio and each gene
  BTHxpnmean<-trioxpnmean[,grep("BTH",colnames(trioxpnmean))] 
  
  #acquire average value in an extreme window (data points=4)
  #Mock
  minxM<-apply(Mockxpnmean,1,function(x){mean(sort(x)[1:4])}) #min expression 
  maxxM<-apply(Mockxpnmean,1,function(x){mean(sort(x)[(length(x)-3):length(x)])}) #max expression
  minyM<-apply(Mockxpnmean,1,function(x){mean(Mockszmean[order(x)[1:4]])}) # size corresponding to min expression
  maxyM<-apply(Mockxpnmean,1,function(x){mean(Mockszmean[order(x)[(length(x)-3):length(x)]])}) # size corresponding to max expression
  #BTH
  minxB<-apply(BTHxpnmean,1,function(x){mean(sort(x)[1:4])})
  maxxB<-apply(BTHxpnmean,1,function(x){mean(sort(x)[(length(x)-3):length(x)])})
  minyB<-apply(BTHxpnmean,1,function(x){mean(BTHszmean[order(x)[1:4]])})
  maxyB<-apply(BTHxpnmean,1,function(x){mean(BTHszmean[order(x)[(length(x)-3):length(x)]])})
  
  #plot
  vioplot(minyM, minyB, maxyM,maxyB, 
          col=alpha(rep(c("springgreen4","darkorange2"),2),0.4),at=c(1,2,4,5),
          xlim=c(0,6),ylim=c(-20,100),
          names=c("Mock.min","BTH.min","Mock.max","BTH.max"),las=2,
          main=paste0("SplineCluster",i," flat-pos n=",length(minyM)),ylab="Rosette size MPH")
  points(jitter(rep(1,length(minyM)),amount=0.1),minyM,pch=20,col="springgreen4")
  points(jitter(rep(2,length(minyB)),amount=0.1),minyB,pch=20,col="darkorange2")
  points(jitter(rep(4,length(maxyM)),amount=0.1),maxyM,pch=20,col="springgreen4")
  points(jitter(rep(5,length(maxyB)),amount=0.1),maxyB,pch=20,col="darkorange2")
  
  #wilcoxon test
  Mocktest<-wilcox.test(minyM, maxyM, paired = TRUE, alternative = "two.sided")
  BTHtest<-wilcox.test(minyB, maxyB, paired = TRUE, alternative = "less")
  statdf<-data.frame("Cluster"=rep(i,2),"Mode"=rep("flat-pos",2),"Treatment"=c("Mock","BTH"),
                     "ContrastA"=rep("Min.xpn.sz",2),"ContrastB"=rep("Max.xpn.sz",2),"Sided"=c("two","less"),
                     "Stat"=c(Mocktest$statistic,BTHtest$statistic),"p.value"=c(Mocktest$p.value,BTHtest$p.value))
  Wilcox.output<-rbind(Wilcox.output,statdf)
}
dev.off()
########
#3.flat-neg
########
pdf("run197_BTHdiffgene_MPHpredict_3.flat-neg_windowvioplot.pdf",height=10,width=25)
par(mfrow=c(2,5))
for(i in clustergenelist[[3]]){
  print(i)
  genelist<-xpnMPH[clustID==i,]
  #calculate trio*treatment mean
  trioxpnmean<-t(apply(genelist,1,function(x){tapply(x,paste0(MPHidx$Treatment,".",MPHidx$TrioID),mean)})) #rows are the genes within the chosen cluster, cols are mean genlist of the trio
  #separate mock and BTH
  Mockxpnmean<-trioxpnmean[,grep("Mock",colnames(trioxpnmean))] #the mean genlist for each mock trio and each gene
  BTHxpnmean<-trioxpnmean[,grep("BTH",colnames(trioxpnmean))] 
  
  #acquire average value in an extreme window (data points=4)
  #Mock
  minxM<-apply(Mockxpnmean,1,function(x){mean(sort(x)[1:4])}) #min expression 
  maxxM<-apply(Mockxpnmean,1,function(x){mean(sort(x)[(length(x)-3):length(x)])}) #max expression
  minyM<-apply(Mockxpnmean,1,function(x){mean(Mockszmean[order(x)[1:4]])}) # size corresponding to min expression
  maxyM<-apply(Mockxpnmean,1,function(x){mean(Mockszmean[order(x)[(length(x)-3):length(x)]])}) # size corresponding to max expression
  #BTH
  minxB<-apply(BTHxpnmean,1,function(x){mean(sort(x)[1:4])})
  maxxB<-apply(BTHxpnmean,1,function(x){mean(sort(x)[(length(x)-3):length(x)])})
  minyB<-apply(BTHxpnmean,1,function(x){mean(BTHszmean[order(x)[1:4]])})
  maxyB<-apply(BTHxpnmean,1,function(x){mean(BTHszmean[order(x)[(length(x)-3):length(x)]])})
  
  #plot
  vioplot(minyM, minyB, maxyM,maxyB, 
          col=alpha(rep(c("springgreen4","darkorange2"),2),0.4),at=c(1,2,4,5),
          xlim=c(0,6),ylim=c(-20,100),
          names=c("Mock.min","BTH.min","Mock.max","BTH.max"),las=2,
          main=paste0("SplineCluster",i," flat-neg n=",length(minyM)),ylab="Rosette size MPH")
  points(jitter(rep(1,length(minyM)),amount=0.1),minyM,pch=20,col="springgreen4")
  points(jitter(rep(2,length(minyB)),amount=0.1),minyB,pch=20,col="darkorange2")
  points(jitter(rep(4,length(maxyM)),amount=0.1),maxyM,pch=20,col="springgreen4")
  points(jitter(rep(5,length(maxyB)),amount=0.1),maxyB,pch=20,col="darkorange2")
  
  #wilcoxon test
  Mocktest<-wilcox.test(minyM, maxyM, paired = TRUE, alternative = "two.sided")
  BTHtest<-wilcox.test(minyB, maxyB, paired = TRUE, alternative = "greater")
  statdf<-data.frame("Cluster"=rep(i,2),"Mode"=rep("flat-neg",2),"Treatment"=c("Mock","BTH"),
                     "ContrastA"=rep("Min.xpn.sz",2),"ContrastB"=rep("Max.xpn.sz",2),"Sided"=c("two","greater"),
                     "Stat"=c(Mocktest$statistic,BTHtest$statistic),"p.value"=c(Mocktest$p.value,BTHtest$p.value))
  Wilcox.output<-rbind(Wilcox.output,statdf)
}
dev.off()
########
#4.pos-flat
########
pdf("run197_BTHdiffgene_MPHpredict_4.pos-flat_windowvioplot.pdf",height=10,width=15)
par(mfrow=c(2,3))
for(i in clustergenelist[[4]]){
  print(i)
  genelist<-xpnMPH[clustID==i,]
  #calculate trio*treatment mean
  trioxpnmean<-t(apply(genelist,1,function(x){tapply(x,paste0(MPHidx$Treatment,".",MPHidx$TrioID),mean)})) #rows are the genes within the chosen cluster, cols are mean genlist of the trio
  #separate mock and BTH
  Mockxpnmean<-trioxpnmean[,grep("Mock",colnames(trioxpnmean))] #the mean genlist for each mock trio and each gene
  BTHxpnmean<-trioxpnmean[,grep("BTH",colnames(trioxpnmean))] 
  
  #acquire average value in an extreme window (data points=4)
  #Mock
  minxM<-apply(Mockxpnmean,1,function(x){mean(sort(x)[1:4])}) #min expression 
  maxxM<-apply(Mockxpnmean,1,function(x){mean(sort(x)[(length(x)-3):length(x)])}) #max expression
  minyM<-apply(Mockxpnmean,1,function(x){mean(Mockszmean[order(x)[1:4]])}) # size corresponding to min expression
  maxyM<-apply(Mockxpnmean,1,function(x){mean(Mockszmean[order(x)[(length(x)-3):length(x)]])}) # size corresponding to max expression
  #BTH
  minxB<-apply(BTHxpnmean,1,function(x){mean(sort(x)[1:4])})
  maxxB<-apply(BTHxpnmean,1,function(x){mean(sort(x)[(length(x)-3):length(x)])})
  minyB<-apply(BTHxpnmean,1,function(x){mean(BTHszmean[order(x)[1:4]])})
  maxyB<-apply(BTHxpnmean,1,function(x){mean(BTHszmean[order(x)[(length(x)-3):length(x)]])})
  
  #plot
  vioplot(minyM, minyB, maxyM,maxyB, 
          col=alpha(rep(c("springgreen4","darkorange2"),2),0.4),at=c(1,2,4,5),
          xlim=c(0,6),ylim=c(-20,100),
          names=c("Mock.min","BTH.min","Mock.max","BTH.max"),las=2,
          main=paste0("SplineCluster",i," pos-flat n=",length(minyM)),ylab="Rosette size MPH")
  points(jitter(rep(1,length(minyM)),amount=0.1),minyM,pch=20,col="springgreen4")
  points(jitter(rep(2,length(minyB)),amount=0.1),minyB,pch=20,col="darkorange2")
  points(jitter(rep(4,length(maxyM)),amount=0.1),maxyM,pch=20,col="springgreen4")
  points(jitter(rep(5,length(maxyB)),amount=0.1),maxyB,pch=20,col="darkorange2")
  
  #wilcoxon test
  Mocktest<-wilcox.test(minyM, maxyM, paired = TRUE, alternative = "less")
  BTHtest<-wilcox.test(minyB, maxyB, paired = TRUE, alternative = "two.sided")
  statdf<-data.frame("Cluster"=rep(i,2),"Mode"=rep("pos-flat",2),"Treatment"=c("Mock","BTH"),
                     "ContrastA"=rep("Min.xpn.sz",2),"ContrastB"=rep("Max.xpn.sz",2),"Sided"=c("less","two"),
                     "Stat"=c(Mocktest$statistic,BTHtest$statistic),"p.value"=c(Mocktest$p.value,BTHtest$p.value))
  Wilcox.output<-rbind(Wilcox.output,statdf)
}
dev.off()
########
#5.pos-pos
########
pdf("run197_BTHdiffgene_MPHpredict_5.pos-pos_windowvioplot.pdf",height=10,width=15)
par(mfrow=c(2,3))
for(i in clustergenelist[[5]]){
  print(i)
  genelist<-xpnMPH[clustID==i,]
  #calculate trio*treatment mean
  trioxpnmean<-t(apply(genelist,1,function(x){tapply(x,paste0(MPHidx$Treatment,".",MPHidx$TrioID),mean)})) #rows are the genes within the chosen cluster, cols are mean genlist of the trio
  #separate mock and BTH
  Mockxpnmean<-trioxpnmean[,grep("Mock",colnames(trioxpnmean))] #the mean genlist for each mock trio and each gene
  BTHxpnmean<-trioxpnmean[,grep("BTH",colnames(trioxpnmean))] 

  #acquire average value in an extreme window (data points=4)
  #Mock
  minxM<-apply(Mockxpnmean,1,function(x){mean(sort(x)[1:4])}) #min expression 
  maxxM<-apply(Mockxpnmean,1,function(x){mean(sort(x)[(length(x)-3):length(x)])}) #max expression
  minyM<-apply(Mockxpnmean,1,function(x){mean(Mockszmean[order(x)[1:4]])}) # size corresponding to min expression
  maxyM<-apply(Mockxpnmean,1,function(x){mean(Mockszmean[order(x)[(length(x)-3):length(x)]])}) # size corresponding to max expression
  #BTH
  minxB<-apply(BTHxpnmean,1,function(x){mean(sort(x)[1:4])})
  maxxB<-apply(BTHxpnmean,1,function(x){mean(sort(x)[(length(x)-3):length(x)])})
  minyB<-apply(BTHxpnmean,1,function(x){mean(BTHszmean[order(x)[1:4]])})
  maxyB<-apply(BTHxpnmean,1,function(x){mean(BTHszmean[order(x)[(length(x)-3):length(x)]])})
  
  #plot
  vioplot(minyM, minyB, maxyM,maxyB, 
          col=alpha(rep(c("springgreen4","darkorange2"),2),0.4),at=c(1,2,4,5),
          xlim=c(0,6),ylim=c(-20,100),
          names=c("Mock.min","BTH.min","Mock.max","BTH.max"),las=2,
          main=paste0("SplineCluster",i," pos-pos n=",length(minyM)),ylab="Rosette size MPH")
  points(jitter(rep(1,length(minyM)),amount=0.1),minyM,pch=20,col="springgreen4")
  points(jitter(rep(2,length(minyB)),amount=0.1),minyB,pch=20,col="darkorange2")
  points(jitter(rep(4,length(maxyM)),amount=0.1),maxyM,pch=20,col="springgreen4")
  points(jitter(rep(5,length(maxyB)),amount=0.1),maxyB,pch=20,col="darkorange2")
  
  #wilcoxon test
  Mocktest<-wilcox.test(minyM, maxyM, paired = TRUE, alternative = "less")
  BTHtest<-wilcox.test(minyB, maxyB, paired = TRUE, alternative = "less")
  statdf<-data.frame("Cluster"=rep(i,2),"Mode"=rep("pos-flat",2),"Treatment"=c("Mock","BTH"),
                     "ContrastA"=rep("Min.xpn.sz",2),"ContrastB"=rep("Max.xpn.sz",2),"Sided"=c("less","less"),
                     "Stat"=c(Mocktest$statistic,BTHtest$statistic),"p.value"=c(Mocktest$p.value,BTHtest$p.value))
  Wilcox.output<-rbind(Wilcox.output,statdf)
}
dev.off()
########
#6.pos-neg
########
pdf("run197_BTHdiffgene_MPHpredict_6.pos-neg_windowvioplot.pdf",height=5,width=10)
par(mfrow=c(1,2))
for(i in clustergenelist[[6]]){
  print(i)
  genelist<-xpnMPH[clustID==i,]
  #calculate trio*treatment mean
  trioxpnmean<-t(apply(genelist,1,function(x){tapply(x,paste0(MPHidx$Treatment,".",MPHidx$TrioID),mean)})) #rows are the genes within the chosen cluster, cols are mean genlist of the trio
  #separate mock and BTH
  Mockxpnmean<-trioxpnmean[,grep("Mock",colnames(trioxpnmean))] #the mean genlist for each mock trio and each gene
  BTHxpnmean<-trioxpnmean[,grep("BTH",colnames(trioxpnmean))] 

  #acquire average value in an extreme window (data points=4)
  #Mock
  minxM<-apply(Mockxpnmean,1,function(x){mean(sort(x)[1:4])}) #min expression 
  maxxM<-apply(Mockxpnmean,1,function(x){mean(sort(x)[(length(x)-3):length(x)])}) #max expression
  minyM<-apply(Mockxpnmean,1,function(x){mean(Mockszmean[order(x)[1:4]])}) # size corresponding to min expression
  maxyM<-apply(Mockxpnmean,1,function(x){mean(Mockszmean[order(x)[(length(x)-3):length(x)]])}) # size corresponding to max expression
  #BTH
  minxB<-apply(BTHxpnmean,1,function(x){mean(sort(x)[1:4])})
  maxxB<-apply(BTHxpnmean,1,function(x){mean(sort(x)[(length(x)-3):length(x)])})
  minyB<-apply(BTHxpnmean,1,function(x){mean(BTHszmean[order(x)[1:4]])})
  maxyB<-apply(BTHxpnmean,1,function(x){mean(BTHszmean[order(x)[(length(x)-3):length(x)]])})
  
  #plot
  vioplot(minyM, minyB, maxyM,maxyB, 
          col=alpha(rep(c("springgreen4","darkorange2"),2),0.4),at=c(1,2,4,5),
          xlim=c(0,6),ylim=c(-20,100),
          names=c("Mock.min","BTH.min","Mock.max","BTH.max"),las=2,
          main=paste0("SplineCluster",i," pos-neg n=",length(minyM)),ylab="Rosette size MPH")
  points(jitter(rep(1,length(minyM)),amount=0.1),minyM,pch=20,col="springgreen4")
  points(jitter(rep(2,length(minyB)),amount=0.1),minyB,pch=20,col="darkorange2")
  points(jitter(rep(4,length(maxyM)),amount=0.1),maxyM,pch=20,col="springgreen4")
  points(jitter(rep(5,length(maxyB)),amount=0.1),maxyB,pch=20,col="darkorange2")
  
  #wilcoxon test
  Mocktest<-wilcox.test(minyM, maxyM, paired = TRUE, alternative = "less")
  BTHtest<-wilcox.test(minyB, maxyB, paired = TRUE, alternative = "greater")
  statdf<-data.frame("Cluster"=rep(i,2),"Mode"=rep("pos-neg",2),"Treatment"=c("Mock","BTH"),
                     "ContrastA"=rep("Min.xpn.sz",2),"ContrastB"=rep("Max.xpn.sz",2),"Sided"=c("less","greater"),
                     "Stat"=c(Mocktest$statistic,BTHtest$statistic),"p.value"=c(Mocktest$p.value,BTHtest$p.value))
  Wilcox.output<-rbind(Wilcox.output,statdf)
}
dev.off()
########
#7. neg-flat
########
pdf("run197_BTHdiffgene_MPHpredict_7.neg-flat_windowvioplot.pdf",height=10,width=10)
par(mfrow=c(2,2))
for(i in clustergenelist[[7]]){
  print(i)
  genelist<-xpnMPH[clustID==i,]
  #calculate trio*treatment mean
  trioxpnmean<-t(apply(genelist,1,function(x){tapply(x,paste0(MPHidx$Treatment,".",MPHidx$TrioID),mean)})) #rows are the genes within the chosen cluster, cols are mean genlist of the trio
  #separate mock and BTH
  Mockxpnmean<-trioxpnmean[,grep("Mock",colnames(trioxpnmean))] #the mean genlist for each mock trio and each gene
  BTHxpnmean<-trioxpnmean[,grep("BTH",colnames(trioxpnmean))] 
  #acquire average value in an extreme window (data points=4)
  #Mock
  minxM<-apply(Mockxpnmean,1,function(x){mean(sort(x)[1:4])}) #min expression 
  maxxM<-apply(Mockxpnmean,1,function(x){mean(sort(x)[(length(x)-3):length(x)])}) #max expression
  minyM<-apply(Mockxpnmean,1,function(x){mean(Mockszmean[order(x)[1:4]])}) # size corresponding to min expression
  maxyM<-apply(Mockxpnmean,1,function(x){mean(Mockszmean[order(x)[(length(x)-3):length(x)]])}) # size corresponding to max expression
  #BTH
  minxB<-apply(BTHxpnmean,1,function(x){mean(sort(x)[1:4])})
  maxxB<-apply(BTHxpnmean,1,function(x){mean(sort(x)[(length(x)-3):length(x)])})
  minyB<-apply(BTHxpnmean,1,function(x){mean(BTHszmean[order(x)[1:4]])})
  maxyB<-apply(BTHxpnmean,1,function(x){mean(BTHszmean[order(x)[(length(x)-3):length(x)]])})
  
  #plot
  vioplot(minyM, minyB, maxyM,maxyB, 
          col=alpha(rep(c("springgreen4","darkorange2"),2),0.4),at=c(1,2,4,5),
          xlim=c(0,6),ylim=c(-20,100),
          names=c("Mock.min","BTH.min","Mock.max","BTH.max"),las=2,
          main=paste0("SplineCluster",i," neg-flat n=",length(minyM)),ylab="Rosette size MPH")
  points(jitter(rep(1,length(minyM)),amount=0.1),minyM,pch=20,col="springgreen4")
  points(jitter(rep(2,length(minyB)),amount=0.1),minyB,pch=20,col="darkorange2")
  points(jitter(rep(4,length(maxyM)),amount=0.1),maxyM,pch=20,col="springgreen4")
  points(jitter(rep(5,length(maxyB)),amount=0.1),maxyB,pch=20,col="darkorange2")
  
  #wilcoxon test
  Mocktest<-wilcox.test(minyM, maxyM, paired = TRUE, alternative = "greater")
  BTHtest<-wilcox.test(minyB, maxyB, paired = TRUE, alternative = "two.sided")
  statdf<-data.frame("Cluster"=rep(i,2),"Mode"=rep("neg-flat",2),"Treatment"=c("Mock","BTH"),
                     "ContrastA"=rep("Min.xpn.sz",2),"ContrastB"=rep("Max.xpn.sz",2),"Sided"=c("greater","two"),
                     "Stat"=c(Mocktest$statistic,BTHtest$statistic),"p.value"=c(Mocktest$p.value,BTHtest$p.value))
  Wilcox.output<-rbind(Wilcox.output,statdf)
}
dev.off()
########
#8.neg-pos
########
pdf("run197_BTHdiffgene_MPHpredict_8.neg-pos_windowvioplot.pdf",height=5,width=5)
#par(mfrow=c(2,2))
for(i in clustergenelist[[8]]){
  print(i)
  genelist<-xpnMPH[clustID==i,]
  #calculate trio*treatment mean
  trioxpnmean<-t(apply(genelist,1,function(x){tapply(x,paste0(MPHidx$Treatment,".",MPHidx$TrioID),mean)})) #rows are the genes within the chosen cluster, cols are mean genlist of the trio
  #separate mock and BTH
  Mockxpnmean<-trioxpnmean[,grep("Mock",colnames(trioxpnmean))] #the mean genlist for each mock trio and each gene
  BTHxpnmean<-trioxpnmean[,grep("BTH",colnames(trioxpnmean))] 
 
  #acquire average value in an extreme window (data points=4)
  #Mock
  minxM<-apply(Mockxpnmean,1,function(x){mean(sort(x)[1:4])}) #min expression 
  maxxM<-apply(Mockxpnmean,1,function(x){mean(sort(x)[(length(x)-3):length(x)])}) #max expression
  minyM<-apply(Mockxpnmean,1,function(x){mean(Mockszmean[order(x)[1:4]])}) # size corresponding to min expression
  maxyM<-apply(Mockxpnmean,1,function(x){mean(Mockszmean[order(x)[(length(x)-3):length(x)]])}) # size corresponding to max expression
  #BTH
  minxB<-apply(BTHxpnmean,1,function(x){mean(sort(x)[1:4])})
  maxxB<-apply(BTHxpnmean,1,function(x){mean(sort(x)[(length(x)-3):length(x)])})
  minyB<-apply(BTHxpnmean,1,function(x){mean(BTHszmean[order(x)[1:4]])})
  maxyB<-apply(BTHxpnmean,1,function(x){mean(BTHszmean[order(x)[(length(x)-3):length(x)]])})
  
  #plot
  vioplot(minyM, minyB, maxyM,maxyB, 
          col=alpha(rep(c("springgreen4","darkorange2"),2),0.4),at=c(1,2,4,5),
          xlim=c(0,6),ylim=c(-20,100),
          names=c("Mock.min","BTH.min","Mock.max","BTH.max"),las=2,
          main=paste0("SplineCluster",i," neg-pos n=",length(minyM)),ylab="Rosette size MPH")
  points(jitter(rep(1,length(minyM)),amount=0.1),minyM,pch=20,col="springgreen4")
  points(jitter(rep(2,length(minyB)),amount=0.1),minyB,pch=20,col="darkorange2")
  points(jitter(rep(4,length(maxyM)),amount=0.1),maxyM,pch=20,col="springgreen4")
  points(jitter(rep(5,length(maxyB)),amount=0.1),maxyB,pch=20,col="darkorange2")
  
  #wilcoxon test
  Mocktest<-wilcox.test(minyM, maxyM, paired = TRUE, alternative = "greater")
  BTHtest<-wilcox.test(minyB, maxyB, paired = TRUE, alternative = "less")
  statdf<-data.frame("Cluster"=rep(i,2),"Mode"=rep("neg-pos",2),"Treatment"=c("Mock","BTH"),
                     "ContrastA"=rep("Min.xpn.sz",2),"ContrastB"=rep("Max.xpn.sz",2),"Sided"=c("greater","less"),
                     "Stat"=c(Mocktest$statistic,BTHtest$statistic),"p.value"=c(Mocktest$p.value,BTHtest$p.value))
  Wilcox.output<-rbind(Wilcox.output,statdf)
}
dev.off()
########
#9.neg-neg
########
pdf("run197_BTHdiffgene_MPHpredict_9.neg-neg_windowvioplot.pdf",height=10,width=25)
par(mfrow=c(2,5))
for(i in clustergenelist[[9]]){
  print(i)
  genelist<-xpnMPH[clustID==i,]
  #calculate trio*treatment mean
  trioxpnmean<-t(apply(genelist,1,function(x){tapply(x,paste0(MPHidx$Treatment,".",MPHidx$TrioID),mean)})) #rows are the genes within the chosen cluster, cols are mean genlist of the trio
  #separate mock and BTH
  Mockxpnmean<-trioxpnmean[,grep("Mock",colnames(trioxpnmean))] #the mean genlist for each mock trio and each gene
  BTHxpnmean<-trioxpnmean[,grep("BTH",colnames(trioxpnmean))] 
 
  #acquire average value in an extreme window (data points=4)
  #Mock
  minxM<-apply(Mockxpnmean,1,function(x){mean(sort(x)[1:4])}) #min expression 
  maxxM<-apply(Mockxpnmean,1,function(x){mean(sort(x)[(length(x)-3):length(x)])}) #max expression
  minyM<-apply(Mockxpnmean,1,function(x){mean(Mockszmean[order(x)[1:4]])}) # size corresponding to min expression
  maxyM<-apply(Mockxpnmean,1,function(x){mean(Mockszmean[order(x)[(length(x)-3):length(x)]])}) # size corresponding to max expression
  #BTH
  minxB<-apply(BTHxpnmean,1,function(x){mean(sort(x)[1:4])})
  maxxB<-apply(BTHxpnmean,1,function(x){mean(sort(x)[(length(x)-3):length(x)])})
  minyB<-apply(BTHxpnmean,1,function(x){mean(BTHszmean[order(x)[1:4]])})
  maxyB<-apply(BTHxpnmean,1,function(x){mean(BTHszmean[order(x)[(length(x)-3):length(x)]])})
  
  #plot
  vioplot(minyM, minyB, maxyM,maxyB, 
          col=alpha(rep(c("springgreen4","darkorange2"),2),0.4),at=c(1,2,4,5),
          xlim=c(0,6),ylim=c(-20,100),
          names=c("Mock.min","BTH.min","Mock.max","BTH.max"),las=2,
          main=paste0("SplineCluster",i," neg-neg n=",length(minyM)),ylab="Rosette size MPH")
  points(jitter(rep(1,length(minyM)),amount=0.1),minyM,pch=20,col="springgreen4")
  points(jitter(rep(2,length(minyB)),amount=0.1),minyB,pch=20,col="darkorange2")
  points(jitter(rep(4,length(maxyM)),amount=0.1),maxyM,pch=20,col="springgreen4")
  points(jitter(rep(5,length(maxyB)),amount=0.1),maxyB,pch=20,col="darkorange2")
  
  #wilcoxon test
  Mocktest<-wilcox.test(minyM, maxyM, paired = TRUE, alternative = "greater")
  BTHtest<-wilcox.test(minyB, maxyB, paired = TRUE, alternative = "greater")
  statdf<-data.frame("Cluster"=rep(i,2),"Mode"=rep("neg-neg",2),"Treatment"=c("Mock","BTH"),
                     "ContrastA"=rep("Min.xpn.sz",2),"ContrastB"=rep("Max.xpn.sz",2),"Sided"=c("greater","greater"),
                     "Stat"=c(Mocktest$statistic,BTHtest$statistic),"p.value"=c(Mocktest$p.value,BTHtest$p.value))
  Wilcox.output<-rbind(Wilcox.output,statdf)
}
dev.off()
########
#10.quad-flat
########
pdf("run197_BTHdiffgene_MPHpredict_10.quad-flat_windowvioplot.pdf",height=10,width=15)
par(mfrow=c(2,3))
for(i in clustergenelist[[10]]){
  print(i)
  genelist<-xpnMPH[clustID==i,]
  #calculate trio*treatment mean
  trioxpnmean<-t(apply(genelist,1,function(x){tapply(x,paste0(MPHidx$Treatment,".",MPHidx$TrioID),mean)})) #rows are the genes within the chosen cluster, cols are mean genlist of the trio
  #separate mock and BTH
  Mockxpnmean<-trioxpnmean[,grep("Mock",colnames(trioxpnmean))] #the mean genlist for each mock trio and each gene
  BTHxpnmean<-trioxpnmean[,grep("BTH",colnames(trioxpnmean))] 
 
  #acquire average value in an extreme window (data points=4)
  #Mock
  minxM<-apply(Mockxpnmean,1,function(x){mean(sort(x)[1:4])}) #min expression 
  maxxM<-apply(Mockxpnmean,1,function(x){mean(sort(x)[(length(x)-3):length(x)])}) #max expression
  medxM<-apply(Mockxpnmean,1,function(x){mean(sort(x)[15:18])})
  minyM<-apply(Mockxpnmean,1,function(x){mean(Mockszmean[order(x)[1:4]])}) # size corresponding to min expression
  maxyM<-apply(Mockxpnmean,1,function(x){mean(Mockszmean[order(x)[(length(x)-3):length(x)]])}) # size corresponding to max expression
  medyM<-apply(Mockxpnmean,1,function(x){mean(Mockszmean[order(x)[15:18]])}) # size corresponding to median expression
  #BTH
  minxB<-apply(BTHxpnmean,1,function(x){mean(sort(x)[1:4])})
  maxxB<-apply(BTHxpnmean,1,function(x){mean(sort(x)[(length(x)-3):length(x)])})
  minyB<-apply(BTHxpnmean,1,function(x){mean(BTHszmean[order(x)[1:4]])})
  maxyB<-apply(BTHxpnmean,1,function(x){mean(BTHszmean[order(x)[(length(x)-3):length(x)]])})
  
  #plot
  vioplot(minyM,minyB, medyM, maxyM,maxyB, 
          col=alpha(c("springgreen4","darkorange2","springgreen4","springgreen4","darkorange2"),0.4),at=c(1,2,3,4,5),
          xlim=c(0,6),ylim=c(-20,100),
          names=c("Mock.min","BTH.min","Mock.median","Mock.max","BTH.max"),las=2,
          main=paste0("SplineCluster",i," quad-flat n=",length(minyM)),ylab="Rosette size MPH")
  points(jitter(rep(1,length(minyM)),amount=0.1),minyM,pch=20,col="springgreen4")
  points(jitter(rep(2,length(minyB)),amount=0.1),minyB,pch=20,col="darkorange2")
  points(jitter(rep(3,length(medyM)),amount=0.1),medyM,pch=20,col="springgreen4")
  points(jitter(rep(4,length(maxyM)),amount=0.1),maxyM,pch=20,col="springgreen4")
  points(jitter(rep(5,length(maxyB)),amount=0.1),maxyB,pch=20,col="darkorange2")
  
  #wilcoxon test
  Mocktest1<-wilcox.test(minyM, medyM, paired = TRUE, alternative = "less")
  Mocktest2<-wilcox.test(maxyM, medyM, paired = TRUE, alternative = "less")
  BTHtest<-wilcox.test(minyB, maxyB, paired = TRUE, alternative = "two.sided")
  statdf<-data.frame("Cluster"=rep(i,3),"Mode"=rep("quad-flat",3),"Treatment"=c("Mock","Mock","BTH"),
                     "ContrastA"=c("Min.xpn.sz","Max.xpn.sz","Min.xpn.sz"),"ContrastB"=c("Median.xpn.sz","Median.xpn.sz","Max.xpn.sz"),"Sided"=c("less","less","two"),
                     "Stat"=c(Mocktest1$statistic,Mocktest2$statistic,BTHtest$statistic),"p.value"=c(Mocktest1$p.value,Mocktest2$p.value,BTHtest$p.value))
  Wilcox.output<-rbind(Wilcox.output,statdf)
}
dev.off()
########
#11. quad-pos
########
pdf("run197_BTHdiffgene_MPHpredict_11.quad-pos_windowvioplot.pdf",height=10,width=10)
par(mfrow=c(2,2))
for(i in clustergenelist[[11]]){
  print(i)
  genelist<-xpnMPH[clustID==i,]
  #calculate trio*treatment mean
  trioxpnmean<-t(apply(genelist,1,function(x){tapply(x,paste0(MPHidx$Treatment,".",MPHidx$TrioID),mean)})) #rows are the genes within the chosen cluster, cols are mean genlist of the trio
  #separate mock and BTH
  Mockxpnmean<-trioxpnmean[,grep("Mock",colnames(trioxpnmean))] #the mean genlist for each mock trio and each gene
  BTHxpnmean<-trioxpnmean[,grep("BTH",colnames(trioxpnmean))] 
  
  #acquire average value in an extreme window (data points=4)
  #Mock
  minxM<-apply(Mockxpnmean,1,function(x){mean(sort(x)[1:4])}) #min expression 
  maxxM<-apply(Mockxpnmean,1,function(x){mean(sort(x)[(length(x)-3):length(x)])}) #max expression
  medxM<-apply(Mockxpnmean,1,function(x){mean(sort(x)[15:18])})
  minyM<-apply(Mockxpnmean,1,function(x){mean(Mockszmean[order(x)[1:4]])}) # size corresponding to min expression
  maxyM<-apply(Mockxpnmean,1,function(x){mean(Mockszmean[order(x)[(length(x)-3):length(x)]])}) # size corresponding to max expression
  medyM<-apply(Mockxpnmean,1,function(x){mean(Mockszmean[order(x)[15:18]])}) # size corresponding to median expression
  #BTH
  minxB<-apply(BTHxpnmean,1,function(x){mean(sort(x)[1:4])})
  maxxB<-apply(BTHxpnmean,1,function(x){mean(sort(x)[(length(x)-3):length(x)])})
  minyB<-apply(BTHxpnmean,1,function(x){mean(BTHszmean[order(x)[1:4]])})
  maxyB<-apply(BTHxpnmean,1,function(x){mean(BTHszmean[order(x)[(length(x)-3):length(x)]])})
  
  #plot
  vioplot(minyM,minyB, medyM, maxyM,maxyB, 
          col=alpha(c("springgreen4","darkorange2","springgreen4","springgreen4","darkorange2"),0.4),at=c(1,2,3,4,5),
          xlim=c(0,6),ylim=c(-20,100),
          names=c("Mock.min","BTH.min","Mock.median","Mock.max","BTH.max"),las=2,
          main=paste0("SplineCluster",i," quad-pos n=",length(minyM)),ylab="Rosette size MPH")
  points(jitter(rep(1,length(minyM)),amount=0.1),minyM,pch=20,col="springgreen4")
  points(jitter(rep(2,length(minyB)),amount=0.1),minyB,pch=20,col="darkorange2")
  points(jitter(rep(3,length(medyM)),amount=0.1),medyM,pch=20,col="springgreen4")
  points(jitter(rep(4,length(maxyM)),amount=0.1),maxyM,pch=20,col="springgreen4")
  points(jitter(rep(5,length(maxyB)),amount=0.1),maxyB,pch=20,col="darkorange2")
  
  #wilcoxon test
  Mocktest1<-wilcox.test(minyM, medyM, paired = TRUE, alternative = "less")
  Mocktest2<-wilcox.test(maxyM, medyM, paired = TRUE, alternative = "less")
  BTHtest<-wilcox.test(minyB, maxyB, paired = TRUE, alternative = "less")
  statdf<-data.frame("Cluster"=rep(i,3),"Mode"=rep("quad-pos",3),"Treatment"=c("Mock","Mock","BTH"),
                     "ContrastA"=c("Min.xpn.sz","Max.xpn.sz","Min.xpn.sz"),"ContrastB"=c("Median.xpn.sz","Median.xpn.sz","Max.xpn.sz"),"Sided"=c("less","less","less"),
                     "Stat"=c(Mocktest1$statistic,Mocktest2$statistic,BTHtest$statistic),"p.value"=c(Mocktest1$p.value,Mocktest2$p.value,BTHtest$p.value))
  Wilcox.output<-rbind(Wilcox.output,statdf)
}
dev.off()
########
#12. quad-neg
########
pdf("run197_BTHdiffgene_MPHpredict_12.quad-neg_windowvioplot.pdf",height=5,width=10)
par(mfrow=c(1,2))
for(i in clustergenelist[[12]]){
  print(i)
  genelist<-xpnMPH[clustID==i,]
  #calculate trio*treatment mean
  trioxpnmean<-t(apply(genelist,1,function(x){tapply(x,paste0(MPHidx$Treatment,".",MPHidx$TrioID),mean)})) #rows are the genes within the chosen cluster, cols are mean genlist of the trio
  #separate mock and BTH
  Mockxpnmean<-trioxpnmean[,grep("Mock",colnames(trioxpnmean))] #the mean genlist for each mock trio and each gene
  BTHxpnmean<-trioxpnmean[,grep("BTH",colnames(trioxpnmean))] 
 
  #acquire average value in an extreme window (data points=4)
  #Mock
  minxM<-apply(Mockxpnmean,1,function(x){mean(sort(x)[1:4])}) #min expression 
  maxxM<-apply(Mockxpnmean,1,function(x){mean(sort(x)[(length(x)-3):length(x)])}) #max expression
  medxM<-apply(Mockxpnmean,1,function(x){mean(sort(x)[15:18])})
  minyM<-apply(Mockxpnmean,1,function(x){mean(Mockszmean[order(x)[1:4]])}) # size corresponding to min expression
  maxyM<-apply(Mockxpnmean,1,function(x){mean(Mockszmean[order(x)[(length(x)-3):length(x)]])}) # size corresponding to max expression
  medyM<-apply(Mockxpnmean,1,function(x){mean(Mockszmean[order(x)[15:18]])}) # size corresponding to median expression
  #BTH
  minxB<-apply(BTHxpnmean,1,function(x){mean(sort(x)[1:4])})
  maxxB<-apply(BTHxpnmean,1,function(x){mean(sort(x)[(length(x)-3):length(x)])})
  minyB<-apply(BTHxpnmean,1,function(x){mean(BTHszmean[order(x)[1:4]])})
  maxyB<-apply(BTHxpnmean,1,function(x){mean(BTHszmean[order(x)[(length(x)-3):length(x)]])})
  
  #plot
  vioplot(minyM,minyB, medyM, maxyM,maxyB, 
          col=alpha(c("springgreen4","darkorange2","springgreen4","springgreen4","darkorange2"),0.4),at=c(1,2,3,4,5),
          xlim=c(0,6),ylim=c(-20,100),
          names=c("Mock.min","BTH.min","Mock.median","Mock.max","BTH.max"),las=2,
          main=paste0("SplineCluster",i," quad-neg n=",length(minyM)),ylab="Rosette size MPH")
  points(jitter(rep(1,length(minyM)),amount=0.1),minyM,pch=20,col="springgreen4")
  points(jitter(rep(2,length(minyB)),amount=0.1),minyB,pch=20,col="darkorange2")
  points(jitter(rep(3,length(medyM)),amount=0.1),medyM,pch=20,col="springgreen4")
  points(jitter(rep(4,length(maxyM)),amount=0.1),maxyM,pch=20,col="springgreen4")
  points(jitter(rep(5,length(maxyB)),amount=0.1),maxyB,pch=20,col="darkorange2")
  
  #wilcoxon test
  Mocktest1<-wilcox.test(minyM, medyM, paired = TRUE, alternative = "less")
  Mocktest2<-wilcox.test(maxyM, medyM, paired = TRUE, alternative = "less")
  BTHtest<-wilcox.test(minyB, maxyB, paired = TRUE, alternative = "greater")
  statdf<-data.frame("Cluster"=rep(i,3),"Mode"=rep("quad-neg",3),"Treatment"=c("Mock","Mock","BTH"),
                     "ContrastA"=c("Min.xpn.sz","Max.xpn.sz","Min.xpn.sz"),"ContrastB"=c("Median.xpn.sz","Median.xpn.sz","Max.xpn.sz"),"Sided"=c("less","less","greater"),
                     "Stat"=c(Mocktest1$statistic,Mocktest2$statistic,BTHtest$statistic),"p.value"=c(Mocktest1$p.value,Mocktest2$p.value,BTHtest$p.value))
  Wilcox.output<-rbind(Wilcox.output,statdf)
}
dev.off()

########
#correct for multiple hypothesis testing, and output result table
Wilcox.output[,"bonferroni"]<-p.adjust(Wilcox.output$p.value,"bonferroni")
write.table(Wilcox.output,"run197_BTHdiffgene_MPHpredict_extremewindowmean_Wilcoxon.txt",quote=F,sep="\t",row.names=F)
