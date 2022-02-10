############
#README
#The script takes per-accession rare allele count, index to positive- and negative clusters, and associate rare allele count with expression ranking
#(expression rank~ RA count relationship,is there on average more RA for a (given) more extreme expression rank than a more moderate one?)
#only genes with non-zero SNPs among accessions were kept

#index relevant clusters

### ++
#    K3       Positive        1.48E-03        Positive        6.12E-07
#    K20
#    K46      Positive        3.12E-05        Positive        1.81E-09
#    K48      Positive        4.88E-10        Positive        5.36E-14  
#    K49      Positive        8.14E-12        Positive        5.27E-12
#    K52
### --
#    K5       Negative        8.07E-12        Negative        3.21E-13  
#    K11      Negative        5.36E-06        Negative        3.63E-06
#    K21      Negative        6.39E-17        Negative        5.90E-18   
#    K29      Negative        5.75E-13        Negative        1.18E-14  
#    K31      Negative        2.41E-13        Negative        2.28E-17  
#    K33      Negative        3.87E-15        Negative        1.38E-11  
#    K36      Negative        1.19E-10        Negative        1.95E-10
#    K37      Negative        4.95E-06        Negative        4.20E-03
#    K39      Negative        1.39E-15        Negative        5.33E-07
#    K60      Negative        4.63E-04        Negative        2.65E-16

##############
#Config and data input
library(vioplot)
library(scales)
setwd("/ebio/abt6_projects7/SHB2/data/1_analysis/10_Subgenome_Genetdist")
Acc.rareallele<-read.table("run197_BTHdiffgene_1kbupstreamONLY_MAF0.005_SHB2acc_pergeneRareallelecount.txt",header=T,sep="\t",row.names=1)
seqmeta<-read.table("/ebio/abt6_projects7/SHB2/data/1_analysis/6_BTHdiffgenes/Run197_rmCMRTPHBS_rm2MSeqtkDupOL_corr.meta.txt",header=T,sep="\t",row.names=1)
#get the cluster assignment of genes
clusttab<-read.table("/ebio/abt6_projects7/SHB2/data/1_analysis/6_BTHdiffgenes/run197_RowFilter_BTHlmm_rmadd_fdr0.001_rmflat_MPHpredsplinek2-200clustID.txt",header=T,sep="\t",row.names=1)
clust61<-clusttab$K61
PPclust<-c(3,20,46,48,49,52) #positive clusters
NNclust<-c(5,11,21,29,31,33,36,37,39,60) #negative clusters
#import the expression value
ScaledTPM<-read.table("/ebio/abt6_projects7/SHB2/data/1_analysis/6_BTHdiffgenes/run197_BTHdiffgene_rowfiltered_scaledlogTPM.txt",header=T,sep="\t",row.names=1)
#input phenotype info
trioMPH<-read.table("/ebio/abt6_projects7/SHB2/data/1_analysis/6_BTHdiffgenes/run197_BTHdiffgene_MetaRowIdx_forpairedMPHcalc_manualcorr.txt",
                    header=T,sep="\t")


#shorten the gene list
ScaledTPM<-ScaledTPM[rownames(ScaledTPM)%in% rownames(clusttab),] 
#index the columns with accession samples (Mock + BTH)
meta.acc<-subset(seqmeta,IsHybrid=="Inbred")
ScaledTPM<-ScaledTPM[,which(seqmeta$IsHybrid=="Inbred")]

#get the genes from each cluster
###PPclusters
PPclustRAmock<-list()
PPclustRAbth<-list()
i<-1
for(p in PPclust){
  print(paste0("clust",p))
  genes<-rownames(clusttab)[clust61==p]
  print(paste0("n=",length(genes)))
  
  #calculate average expression rank for each gene across replicates, in Mock and BTH,
  #and use the rank to sort the RA table
  clustRAmock<-matrix(nrow=length(genes),ncol=64)
  clustRAbth<-matrix(nrow=length(genes),ncol=64)
  j<-1
  for(g in genes){
    xpn<-as.numeric(ScaledTPM[rownames(ScaledTPM)==g,])
    if(!any(rownames(Acc.rareallele)==g)){
      print(paste0(g," has no upstream rare allele"))
      next
    }
    RA<-Acc.rareallele[rownames(Acc.rareallele)==g,]
    meanxpn<-tapply(xpn,paste0(meta.acc$Genotype,"_",meta.acc$Treatment),mean)
    mockxpn<-meanxpn[grepl("Mock",names(meanxpn))]
    BTHxpn<-meanxpn[grepl("BTH",names(meanxpn))]
    mockrank<-as.character(sapply(names(sort(rank(mockxpn))),function(x){substr(x,1,nchar(x)-5)})) #rank value from low to high
    BTHrank<-as.character(sapply(names(sort(rank(BTHxpn))),function(x){substr(x,1,nchar(x)-4)}))
    RAmock<-RA[as.numeric(sapply(mockrank,function(x){which(names(RA)==x)}))]
    RAbth<-RA[sapply(BTHrank,function(x){which(names(RA)==x)})]
    clustRAmock[j,]<-as.numeric(RAmock)
    clustRAbth[j,]<-as.numeric(RAbth)
    j<-(j+1)
  }
  PPclustRAmock[[i]]<-na.omit(clustRAmock)
  PPclustRAbth[[i]]<-na.omit(clustRAbth)
  i<-(i+1)
}
############
# [1] "clust3"
# [1] "n=115"
# [1] "AT2G03510 has no upstream rare allele"
# [1] "AT2G07180 has no upstream rare allele"
# [1] "AT2G41450 has no upstream rare allele"
# [1] "AT5G39660 has no upstream rare allele"
# [1] "AT5G41480 has no upstream rare allele"
# [1] "AT5G64370 has no upstream rare allele"
# [1] "clust20"
# [1] "n=32"
# [1] "clust46"
# [1] "n=59"
# [1] "AT1G14880 has no upstream rare allele"
# [1] "AT1G64150 has no upstream rare allele"
# [1] "AT4G02530 has no upstream rare allele"
# [1] "AT5G45930 has no upstream rare allele"
# [1] "AT5G59670 has no upstream rare allele"
# [1] "clust48"
# [1] "n=88"
# [1] "AT2G04690 has no upstream rare allele"
# [1] "AT3G27350 has no upstream rare allele"
# [1] "AT4G02920 has no upstream rare allele"
# [1] "AT4G08930 has no upstream rare allele"
# [1] "AT4G23130 has no upstream rare allele"
# [1] "AT5G39785 has no upstream rare allele"
# [1] "AT5G44580 has no upstream rare allele"
# [1] "clust49"
# [1] "n=93"
# [1] "AT1G09932 has no upstream rare allele"
# [1] "AT2G17705 has no upstream rare allele"
# [1] "AT3G07195 has no upstream rare allele"
# [1] "AT3G45290 has no upstream rare allele"
# [1] "AT4G04700 has no upstream rare allele"
# [1] "AT4G05060 has no upstream rare allele"
# [1] "AT4G23150 has no upstream rare allele"
# [1] "AT5G08120 has no upstream rare allele"
# [1] "AT5G36700 has no upstream rare allele"
# [1] "AT5G36790 has no upstream rare allele"
# [1] "clust52"
# [1] "n=112"
# [1] "AT3G15900 has no upstream rare allele"
# [1] "AT4G04740 has no upstream rare allele"
# [1] "AT5G24230 has no upstream rare allele"
# [1] "AT5G44578 has no upstream rare allele"
############

###NNclusters
NNclustRAmock<-list()
NNclustRAbth<-list()
i<-1
for(n in NNclust){
  print(paste0("clust",n))
  genes<-rownames(clusttab)[clust61==n]
  print(paste0("n=",length(genes)))
 
  #calculate average expression rank for each gene across replicates, in Mock and BTH,
  #and use the rank to sort the RA table
  clustRAmock<-matrix(nrow=length(genes),ncol=64)
  clustRAbth<-matrix(nrow=length(genes),ncol=64)
  j<-1
  for(g in genes){
    xpn<-as.numeric(ScaledTPM[rownames(ScaledTPM)==g,])
    if(!any(rownames(Acc.rareallele)==g)){
      print(paste0(g," has no upstream rare allele"))
      next
    }
    RA<-Acc.rareallele[rownames(Acc.rareallele)==g,]
    meanxpn<-tapply(xpn,paste0(meta.acc$Genotype,"_",meta.acc$Treatment),mean)
    mockxpn<-meanxpn[grepl("Mock",names(meanxpn))]
    BTHxpn<-meanxpn[grepl("BTH",names(meanxpn))]
    mockrank<-as.character(sapply(names(sort(rank(mockxpn))),function(x){substr(x,1,nchar(x)-5)})) #rank value from low to high
    BTHrank<-as.character(sapply(names(sort(rank(BTHxpn))),function(x){substr(x,1,nchar(x)-4)}))
    RAmock<-RA[as.numeric(sapply(mockrank,function(x){which(names(RA)==x)}))]
    RAbth<-RA[sapply(BTHrank,function(x){which(names(RA)==x)})]
    clustRAmock[j,]<-as.numeric(RAmock)
    clustRAbth[j,]<-as.numeric(RAbth)
    j<-(j+1)
  }
  NNclustRAmock[[i]]<-na.omit(clustRAmock)
  NNclustRAbth[[i]]<-na.omit(clustRAbth)
  i<-(i+1)
}

#####
# [1] "clust5"
# [1] "n=82"
# [1] "AT1G60890 has no upstream rare allele"
# [1] "AT2G17870 has no upstream rare allele"
# [1] "AT3G54950 has no upstream rare allele"
# [1] "AT5G09460 has no upstream rare allele"
# [1] "AT5G09461 has no upstream rare allele"
# [1] "AT5G09462 has no upstream rare allele"
# [1] "AT5G09463 has no upstream rare allele"
# [1] "AT5G41580 has no upstream rare allele"
# [1] "clust11"
# [1] "n=39"
# [1] "AT3G16450 has no upstream rare allele"
# [1] "AT4G09030 has no upstream rare allele"
# [1] "clust21"
# [1] "n=110"
# [1] "AT2G35110 has no upstream rare allele"
# [1] "AT3G51770 has no upstream rare allele"
# [1] "AT4G17140 has no upstream rare allele"
# [1] "AT5G07690 has no upstream rare allele"
# [1] "AT5G64340 has no upstream rare allele"
# [1] "AT5G64341 has no upstream rare allele"
# [1] "AT5G64342 has no upstream rare allele"
# [1] "AT5G64343 has no upstream rare allele"
# [1] "clust29"
# [1] "n=90"
# [1] "AT1G71010 has no upstream rare allele"
# [1] "AT2G03090 has no upstream rare allele"
# [1] "AT2G33620 has no upstream rare allele"
# [1] "AT5G43560 has no upstream rare allele"
# [1] "clust31"
# [1] "n=125"
# [1] "AT1G43580 has no upstream rare allele"
# [1] "AT1G44910 has no upstream rare allele"
# [1] "AT1G50280 has no upstream rare allele"
# [1] "AT1G52080 has no upstream rare allele"
# [1] "AT1G70560 has no upstream rare allele"
# [1] "AT2G21840 has no upstream rare allele"
# [1] "AT3G23590 has no upstream rare allele"
# [1] "AT4G00820 has no upstream rare allele"
# [1] "AT4G02710 has no upstream rare allele"
# [1] "AT5G40830 has no upstream rare allele"
# [1] "clust33"
# [1] "n=102"
# [1] "AT1G36070 has no upstream rare allele"
# [1] "AT1G69450 has no upstream rare allele"
# [1] "AT1G73360 has no upstream rare allele"
# [1] "AT3G28740 has no upstream rare allele"
# [1] "AT3G29390 has no upstream rare allele"
# [1] "AT3G30300 has no upstream rare allele"
# [1] "AT3G61480 has no upstream rare allele"
# [1] "AT4G11800 has no upstream rare allele"
# [1] "AT5G52230 has no upstream rare allele"
# [1] "AT5G63950 has no upstream rare allele"
# [1] "clust36"
# [1] "n=79"
# [1] "AT1G12930 has no upstream rare allele"
# [1] "AT1G17990 has no upstream rare allele"
# [1] "AT1G47380 has no upstream rare allele"
# [1] "AT1G58100 has no upstream rare allele"
# [1] "AT3G19270 has no upstream rare allele"
# [1] "AT4G34588 has no upstream rare allele"
# [1] "AT4G34590 has no upstream rare allele"
# [1] "AT5G18270 has no upstream rare allele"
# [1] "AT5G35730 has no upstream rare allele"
# [1] "clust37"
# [1] "n=28"
# [1] "clust39"
# [1] "n=120"
# [1] "AT2G05760 has no upstream rare allele"
# [1] "AT2G25560 has no upstream rare allele"
# [1] "AT2G44280 has no upstream rare allele"
# [1] "AT3G53240 has no upstream rare allele"
# [1] "AT4G00755 has no upstream rare allele"
# [1] "clust60"
# [1] "n=105"
# [1] "AT3G44160 has no upstream rare allele"
# [1] "AT5G43670 has no upstream rare allele"
# [1] "AT5G49520 has no upstream rare allele"
#####

#associate the count of rare alleles to gene expression

#calculate mean rare alleles per rank across all cluster members
PPclust_meanRAmock<-lapply(PPclustRAmock,function(x){apply(x,2,mean)}) #one mean RA per cluster
PPclust_meanRAbth<-lapply(PPclustRAbth,function(x){apply(x,2,mean)})
NNclust_meanRAmock<-lapply(NNclustRAmock,function(x){apply(x,2,mean)})
NNclust_meanRAbth<-lapply(NNclustRAbth,function(x){apply(x,2,mean)})

#merge all clusters within PP or NN, plot the overall

unlsPPmock<-do.call("rbind",PPclustRAmock) #467 genes with rare alleles in the PP clusters
unlsPPbth<-do.call("rbind",PPclustRAbth)
unlsPPmock.mean<-apply(unlsPPmock,2,mean)
unlsPPbth.mean<-apply(unlsPPbth,2,mean)

unlsNNmock<-do.call("rbind",NNclustRAmock) #821 genes with RA in the NN clusters
unlsNNbth<-do.call("rbind",NNclustRAbth)
unlsNNmock.mean<-apply(unlsNNmock,2,mean)
unlsNNbth.mean<-apply(unlsNNbth,2,mean)

#use paired Wilcoxon test to see if low rank and high rank has different average RA count. 
# To avoid test statistic solely driven by outliers, 
# average each gene's rare allele count among the 10% lowest (7 out of 64) & 10% highest ranks (7 out of 64)
PPmocklow0.1mean<-apply(unlsPPmock[,1:7],1,mean)
PPmockhigh0.1mean<-apply(unlsPPmock[,(ncol(unlsPPmock)-6):ncol(unlsPPmock)],1,mean)
PPmockLvH<-wilcox.test(PPmocklow0.1mean,PPmockhigh0.1mean,paired=T,alternative = "greater")
# Wilcoxon signed rank test with continuity correction
# 
# data:  PPmocklow0.1mean and PPmockhigh0.1mean
# V = 36063, p-value = 0.02777
# alternative hypothesis: true location shift is greater than 0

PPbthlow0.1mean<-apply(unlsPPbth[,1:7],1,mean)
PPbthhigh0.1mean<-apply(unlsPPbth[,(ncol(unlsPPbth)-6):ncol(unlsPPbth)],1,mean)
PPbthLvH<-wilcox.test(PPbthlow0.1mean,PPbthhigh0.1mean,paired=F,alternative = "greater")
# Wilcoxon rank sum test with continuity correction
# 
# data:  PPbthlow0.1mean and PPbthhigh0.1mean
# W = 114556, p-value = 0.08721
#alternative hypothesis: true location shift is greater than 0

NNmocklow0.1mean<-apply(unlsNNmock[,1:7],1,mean)
NNmockhigh0.1mean<-apply(unlsNNmock[,(ncol(unlsNNmock)-6):ncol(unlsNNmock)],1,mean)
NNmockLvH<-wilcox.test(NNmocklow0.1mean,NNmockhigh0.1mean,paired=T,alternative = "less")
# Wilcoxon signed rank test with continuity correction
# 
# data:  NNmocklow0.1mean and NNmockhigh0.1mean
# V = 79558, p-value = 2.173e-06
# alternative hypothesis: true location shift is less than 0

NNbthlow0.1mean<-apply(unlsNNbth[,1:7],1,mean)
NNbthhigh0.1mean<-apply(unlsNNbth[,(ncol(unlsNNbth)-6):ncol(unlsNNbth)],1,mean)
NNbthLvH<-wilcox.test(NNbthlow0.1mean,NNbthhigh0.1mean,paired=F,alternative = "less")
# Wilcoxon rank sum test with continuity correction
# 
# data:  NNbthlow0.1mean and NNbthhigh0.1mean
# W = 289102, p-value = 1.761e-07
# alternative hypothesis: true location shift is less than 0

pdf("run197_MPHpredict_expressionrank_rareallele.1.pdf",width=15,height=15)
par(mfrow=c(2,2))
plot(1:64,unlsPPmock.mean,pch=20,col="springgreen4",ylim=c(0.2,1.2), xlim=c(0,65),
     main="Average of all positively-predicting genes",
     xlab="Expression Rank",ylab="Mean Rare-allele Count")
points(1:64,unlsPPbth.mean,pch=20,col="darkorange2")
lines(lowess(1:64,unlsPPmock.mean),lwd=3,col="springgreen4")
lines(lowess(1:64,unlsPPbth.mean),lwd=3,col="darkorange2")
legend(1,1.2,legend=c("Mock","BTH"),pch=20,lwd=2,col=c("springgreen4","darkorange2"),bty="n")
mtext(text=c("LOW","HIGH"),side=1,line=3,at=c(0,62))

plot(1:64,unlsNNmock.mean,pch=20,col="springgreen4",ylim=c(0.2,1.2), xlim=c(0,65),
     main="Average of all negatively-predicting genes",
     xlab="Expression Rank",ylab="Mean Rare-allele Count")
points(1:64,unlsNNbth.mean,pch=20,col="darkorange2")
lines(lowess(1:64,unlsNNmock.mean),lwd=3,col="springgreen4")
lines(lowess(1:64,unlsNNbth.mean),lwd=3,col="darkorange2")
legend(50,1.2,legend=c("Mock","BTH"),pch=20,lwd=2,col=c("springgreen4","darkorange2"),bty="n")
mtext(text=c("LOW","HIGH"),side=1,line=3,at=c(0,62))

#violin plot of low vs high- rank RAcount
#PP
vioplot(PPmocklow0.1mean,PPmockhigh0.1mean,PPbthlow0.1mean,PPbthhigh0.1mean,
        col=alpha(rep(c("springgreen4","darkorange2"),each=2),0.7),ylim=c(0,17),
        main="Positive-positive genes mean low- vs high- rank RA",las=2,
        names=c("M.L.","M.H.","B.L.","B.H."), ylab="Mean rare-allele count")
set.seed(97)
points(x=rep(1,length(PPmocklow0.1mean))+rnorm(length(PPmocklow0.1mean),0,0.1),
       y=PPmocklow0.1mean,pch=20,col=alpha("springgreen4",0.5))

points(x=rep(2,length(PPmockhigh0.1mean))+rnorm(length(PPmockhigh0.1mean),0,0.1),
       y=PPmockhigh0.1mean,pch=20,col=alpha("springgreen4",0.5))

points(x=rep(3,length(PPbthlow0.1mean))+rnorm(length(PPbthlow0.1mean),0,0.1),
       y=PPbthlow0.1mean,pch=20,col=alpha("darkorange2",0.5))

points(x=rep(4,length(PPbthhigh0.1mean))+rnorm(length(PPbthhigh0.1mean),0,0.1),
       y=PPbthhigh0.1mean,pch=20,col=alpha("darkorange2",0.5))

segments(1,14,1,16,lwd=3)
segments(2,8,2,16,lwd=3)
segments(2,16,1.7,16,lwd=3)
segments(1,16,1.3,16,lwd=3)
text(1.5,16,paste0("p=",round(PPmockLvH$p.value,digits = 3)),cex=0.7,srt=45,pos=1)

segments(3,12,3,16,lwd=3)
segments(4,5,4,16,lwd=3)
segments(3,16,3.3,16,lwd=3)
segments(4,16,3.7,16,lwd=3)
text(3.5,16,paste0("p=",round(PPbthLvH$p.value,digits = 3)),cex=0.7,srt=45,pos=1)

##NN
vioplot(NNmocklow0.1mean,NNmockhigh0.1mean,NNbthlow0.1mean,NNbthhigh0.1mean,
        col=alpha(rep(c("springgreen4","darkorange2"),each=2),0.7),ylim=c(0,17),
        main="Negative-negative genes mean low- vs high- rank RA",las=2,
        names=c("M.L.","M.H.","B.L.","B.H."), ylab="Mean rare-allele count")
set.seed(97)
points(x=rep(1,length(NNmocklow0.1mean))+rnorm(length(NNmocklow0.1mean),0,0.1),
       y=NNmocklow0.1mean,pch=20,col=alpha("springgreen4",0.5))

points(x=rep(2,length(NNmockhigh0.1mean))+rnorm(length(NNmockhigh0.1mean),0,0.1),
       y=NNmockhigh0.1mean,pch=20,col=alpha("springgreen4",0.5))

points(x=rep(3,length(NNbthlow0.1mean))+rnorm(length(NNbthlow0.1mean),0,0.1),
       y=NNbthlow0.1mean,pch=20,col=alpha("darkorange2",0.5))

points(x=rep(4,length(NNbthhigh0.1mean))+rnorm(length(NNbthhigh0.1mean),0,0.1),
       y=NNbthhigh0.1mean,pch=20,col=alpha("darkorange2",0.5))

segments(1,8,1,16,lwd=3)
segments(2,13,2,16,lwd=3)
segments(1,16,1.3,16,lwd=3)
segments(2,16,1.7,16,lwd=3)
text(1.5,16,paste0("p=",formatC(NNmockLvH$p.value,format="e",digits = 2)),cex=0.7,srt=45,pos=1)


segments(3,9,3,16,lwd=3)
segments(4,15,4,16,lwd=3)
segments(3,16,3.3,16,lwd=3)
segments(4,16,3.7,16,lwd=3)
text(3.5,16,paste0("p=",formatC(NNbthLvH$p.value,format="e",digits = 2)),cex=0.7,srt=45,pos=1)

dev.off()


###############################################################################
#ask whether F1 expression can be associated with parental rare-allele counts
#use mean rare alleles in parents as indicator

ScaledTPM<-read.table("/ebio/abt6_projects7/SHB2/data/1_analysis/6_BTHdiffgenes/run197_BTHdiffgene_rowfiltered_scaledlogTPM.txt",header=T,sep="\t",row.names=1) #read in expression value again
#shorten the gene list
ScaledTPM<-ScaledTPM[rownames(ScaledTPM)%in% rownames(clusttab),] 
#index the columns with accession samples (Mock + BTH)
meta.F1<-subset(seqmeta,IsHybrid=="Hybrid") #the meta file specified the
ScaledTPM<-ScaledTPM[,which(seqmeta$IsHybrid=="Hybrid")]

#get the genes from each cluster, rank their F1 expression (record the trio association)
#simulatneously get the RA counts from all parents (record the trio association)
#then calculate parental-RA vs. F1-xpn association based on the above-mentioned 5 criteria

###PPclusts
PPclustRAacc<-list()
PPclustF1rank<-list()

i<-1
for(p in PPclust){
  print(paste0("clust",p))
  genes<-rownames(clusttab)[clust61==p]
  print(paste0("n=",length(genes)))
  
  #calculate average expression rank for each gene across replicates, in Mock and BTH,
  #and extract the parental RA count
  clustRA<-matrix(nrow=length(genes),ncol=64)
  F1xpnRank<-matrix(nrow=length(genes),ncol=64) #first 32 columns for mock, and the remaining 32 columns for BTH
  
  j<-1
  for(g in genes){
    xpn<-as.numeric(ScaledTPM[rownames(ScaledTPM)==g,])
    if(!any(rownames(Acc.rareallele)==g)){
      print(paste0(g," has no upstream rare allele"))
      next
    }
    RA<-Acc.rareallele[rownames(Acc.rareallele)==g,]
    meanxpn<-tapply(xpn,paste0(meta.F1$Trio,"_",meta.F1$Treatment),mean)
    mockxpn<-meanxpn[grepl("Mock",names(meanxpn))]
    BTHxpn<-meanxpn[grepl("BTH",names(meanxpn))]
    mockrank<-as.character(sapply(names(sort(rank(mockxpn))),function(x){substr(x,1,nchar(x)-5)})) #rank value from low to high
    BTHrank<-as.character(sapply(names(sort(rank(BTHxpn))),function(x){substr(x,1,nchar(x)-4)}))
    clustRA[j,]<-as.numeric(RA)
    F1xpnRank[j,]<-c(paste0("M",mockrank),paste0("B",BTHrank))
    j<-(j+1)
  }
  
  PPclustRAacc[[i]]<-na.omit(clustRA)
  PPclustF1rank[[i]]<-na.omit(F1xpnRank)
  i<-(i+1)
}
PPclustRAacc<-do.call("rbind",PPclustRAacc)
PPclustF1rank<-do.call("rbind",PPclustF1rank) #1:32 are mock ranks for each gene, 33:64 are BTH ranks

# now we have a list of accession RAs and a list of F1 xpn ranks, both in corresponding order
acc2trio<-sapply(colnames(Acc.rareallele),function(x){unique(seqmeta$Trio[which(seqmeta$Genotype==x)])})

#####
# calculate for each gene: F1 xpn rank
trioRA<-sapply(acc2trio,function(x){substr(x,1,nchar(x)-1)})
meanRA<-t(apply(PPclustRAacc,1,function(x){sapply(unique(trioRA),function(y){mean(x[trioRA==y])})})) #467*32,each row is a gene, each col is a trio
#sort RA by rank for each gene
sortingrankMock<-c()
sortingrankBTH<-c()#each row is a gene in the PP cluster, each column is a given rank.
for(g in 1:nrow(PPclustRAacc)){ 
  if(g%%10==1){print(g)}
  RA<-PPclustRAacc[g,]
  F1rankM<-sapply(PPclustF1rank[g,1:32],function(x){substr(x,start=2,stop=nchar(x)-1)}) #F1ranking by Mock
  F1rankB<-sapply(PPclustF1rank[g,33:64], function(x){substr(x,start=2,stop=nchar(x)-1)})#F1 ranking by BTH
  RArankM<-as.numeric(sapply(F1rankM,function(x){which(unique(trioRA)==x)}))
  RArankB<-as.numeric(sapply(F1rankB,function(x){which(unique(trioRA)==x)}))
  sortingrankMock<-rbind(sortingrankMock,RArankM)
  sortingrankBTH<-rbind(sortingrankBTH,RArankB)
}

#rank each gene's RA by the table, calculate mean over each given rank
meanRA_rankMock<-t(sapply(1:nrow(PPclustRAacc),function(x){return(as.numeric(meanRA[x,sortingrankMock[x,]]))})) #each row is a gene, each column is a rank (low->high)
meanRA_rankBTH<-t(sapply(1:nrow(PPclustRAacc),function(x){return(as.numeric(meanRA[x,sortingrankBTH[x,]]))})) #each row is a gene, each column is a rank (low->high)
meanRA_meanrankMock<-apply(meanRA_rankMock,2,mean)
meanRA_meanrankBTH<-apply(meanRA_rankBTH,2,mean)

#use paired wilcoxon test to see if low rank and high rank has different average RA count. 
#average each gene's rare allele count among the 10% lowest (4 out of 32) & 10% highest ranks
PPmocklow0.1mean<-apply(meanRA_rankMock[,1:4],1,mean)
PPmockhigh0.1mean<-apply(meanRA_rankMock[,(ncol(meanRA_rankMock)-3):ncol(meanRA_rankMock)],1,mean)
PPmockLvH<-wilcox.test(PPmocklow0.1mean,PPmockhigh0.1mean,paired=T,alternative = "greater")
# Wilcoxon signed rank test with continuity correction
# 
# data:  PPmocklow0.1mean and PPmockhigh0.1mean
# V = 38772, p-value = 0.003697
# alternative hypothesis: true location shift is greater than 0

PPbthlow0.1mean<-apply(meanRA_rankBTH[,1:4],1,mean)
PPbthhigh0.1mean<-apply(meanRA_rankBTH[,(ncol(meanRA_rankBTH)-3):ncol(meanRA_rankBTH)],1,mean)
PPbthLvH<-wilcox.test(PPbthlow0.1mean,PPbthhigh0.1mean,paired=F,alternative = "greater")
# Wilcoxon rank sum test with continuity correction
# 
# data:  PPbthlow0.1mean and PPbthhigh0.1mean
# W = 112124, p-value = 0.2245
# alternative hypothesis: true location shift is greater than 0

########
#Negative clusts
NNclustRAacc<-list()
NNclustF1rank<-list()

i<-1
for(n in NNclust){
  print(paste0("clust",n))
  genes<-rownames(clusttab)[clust61==n]
  print(paste0("n=",length(genes)))
  
  #calculate average expression rank for each gene across replicates, in Mock and BTH,
  #and use the rank to sort the RA table
  clustRA<-matrix(nrow=length(genes),ncol=64)
  F1xpnRank<-matrix(nrow=length(genes),ncol=64) #first 32 columns for mock, and the remaining 32 columns for BTH
  
  j<-1
  for(g in genes){
    xpn<-as.numeric(ScaledTPM[rownames(ScaledTPM)==g,])
    if(!any(rownames(Acc.rareallele)==g)){
      print(paste0(g," has no upstream rare allele"))
      next
    }
    RA<-Acc.rareallele[rownames(Acc.rareallele)==g,]
    meanxpn<-tapply(xpn,paste0(meta.F1$Trio,"_",meta.F1$Treatment),mean)
    mockxpn<-meanxpn[grepl("Mock",names(meanxpn))]
    BTHxpn<-meanxpn[grepl("BTH",names(meanxpn))]
    mockrank<-as.character(sapply(names(sort(rank(mockxpn))),function(x){substr(x,1,nchar(x)-5)})) #rank value from low to high
    BTHrank<-as.character(sapply(names(sort(rank(BTHxpn))),function(x){substr(x,1,nchar(x)-4)}))
    clustRA[j,]<-as.numeric(RA)
    F1xpnRank[j,]<-c(paste0("M",mockrank),paste0("B",BTHrank))
    j<-(j+1)
  }
  NNclustRAacc[[i]]<-na.omit(clustRA)
  NNclustF1rank[[i]]<-na.omit(F1xpnRank)
  i<-(i+1)
}
NNclustRAacc<-do.call("rbind",NNclustRAacc)
NNclustF1rank<-do.call("rbind",NNclustF1rank)
#####
# calculate for each gene: F1 xpn rank
NNtrioRA<-sapply(acc2trio,function(x){substr(x,1,nchar(x)-1)})
NNmeanRA<-t(apply(NNclustRAacc,1,function(x){sapply(unique(trioRA),function(y){mean(x[trioRA==y])})}))

#sort RA by rank for each gene
NNsortingrankMock<-c()
NNsortingrankBTH<-c()#each row is a gene in the cluster, each column is a given rank.
for(g in 1:nrow(NNclustRAacc)){ 
  if(g%%10==1){print(g)}
  RA<-NNclustRAacc[g,]
  F1rankM<-sapply(NNclustF1rank[g,1:32],function(x){substr(x,start=2,stop=nchar(x)-1)}) #F1ranking by Mock
  F1rankB<-sapply(NNclustF1rank[g,33:64], function(x){substr(x,start=2,stop=nchar(x)-1)})#F1 ranking by BTH
  RArankM<-as.numeric(sapply(F1rankM,function(x){which(unique(trioRA)==x)}))
  RArankB<-as.numeric(sapply(F1rankB,function(x){which(unique(trioRA)==x)}))
  NNsortingrankMock<-rbind(NNsortingrankMock,RArankM)
  NNsortingrankBTH<-rbind(NNsortingrankBTH,RArankB)
}

#rank each gene's RA by the table, calculate mean over each given rank

NNmeanRA_rankMock<-t(sapply(1:nrow(NNclustRAacc),function(x){return(as.numeric(NNmeanRA[x,NNsortingrankMock[x,]]))})) #each row is a gene, each column is a rank (low->high)
NNmeanRA_rankBTH<-t(sapply(1:nrow(NNclustRAacc),function(x){return(as.numeric(NNmeanRA[x,NNsortingrankBTH[x,]]))})) #each row is a gene, each column is a rank (low->high)
NNmeanRA_meanrankMock<-apply(NNmeanRA_rankMock,2,mean)
NNmeanRA_meanrankBTH<-apply(NNmeanRA_rankBTH,2,mean)


NNmocklow0.1mean<-apply(NNmeanRA_rankMock[,1:4],1,mean)
NNmockhigh0.1mean<-apply(NNmeanRA_rankMock[,(ncol(NNmeanRA_rankMock)-3):ncol(NNmeanRA_rankMock)],1,mean)
NNmockLvH<-wilcox.test(NNmocklow0.1mean,NNmockhigh0.1mean,paired=T,alternative = "less")
# Wilcoxon signed rank test with continuity correction
# 
# data:  NNmocklow0.1mean and NNmockhigh0.1mean
# V = 70638, p-value = 2.46e-13
# alternative hypothesis: true location shift is less than 0

NNbthlow0.1mean<-apply(NNmeanRA_rankBTH[,1:4],1,mean)
NNbthhigh0.1mean<-apply(NNmeanRA_rankBTH[,(ncol(NNmeanRA_rankBTH)-3):ncol(NNmeanRA_rankBTH)],1,mean)
NNbthLvH<-wilcox.test(NNbthlow0.1mean,NNbthhigh0.1mean,paired=F,alternative = "less")
# Wilcoxon rank sum test with continuity correction
# 
# data:  NNbthlow0.1mean and NNbthhigh0.1mean
# W = 286838, p-value = 5.575e-08
# alternative hypothesis: true location shift is less than 0

###########
#plot
pdf("run197_MPHpredict_PosNegClust_F1expressionrank_rareallele.1.pdf",width=15,height=15)
par(mfcol=c(2,2))

#Positives

plot(1:32,meanRA_meanrankMock,pch=20,col="springgreen4",ylim=c(0.2,1.2),
     main="Mean parental rare alleles in positively-predicting genes",
     xlab="F1 Expression Rank",ylab="Mean parental Rare-allele Count",xpd=NA)
points(1:32,meanRA_meanrankBTH,pch=20,col="darkorange2")
lines(lowess(1:32,meanRA_meanrankMock),lwd=3,col="springgreen4")
lines(lowess(1:32,meanRA_meanrankBTH),lwd=3,col="darkorange2")
legend(1,1.2,legend=c("Mock","BTH"),pch=20,lwd=2,col=c("springgreen4","darkorange2"),bty="n")
mtext(text=c("LOW","HIGH"),side=1,line=3,at=c(0,31))

#violin plot of low vs high- rank RAcount
#PP
vioplot(PPmocklow0.1mean,PPmockhigh0.1mean,PPbthlow0.1mean,PPbthhigh0.1mean,
        col=alpha(rep(c("springgreen4","darkorange2"),each=2),0.7),ylim=c(0,13),
        main="Positive-positive genes mean low- vs high- rank RA",las=2,
        names=c("M.L.","M.H.","B.L.","B.H."), ylab="Mean rare-allele count")
set.seed(97)
points(x=rep(1,length(PPmocklow0.1mean))+rnorm(length(PPmocklow0.1mean),0,0.1),
       y=PPmocklow0.1mean,pch=20,col=alpha("springgreen4",0.5))

points(x=rep(2,length(PPmockhigh0.1mean))+rnorm(length(PPmockhigh0.1mean),0,0.1),
       y=PPmockhigh0.1mean,pch=20,col=alpha("springgreen4",0.5))

points(x=rep(3,length(PPbthlow0.1mean))+rnorm(length(PPbthlow0.1mean),0,0.1),
       y=PPbthlow0.1mean,pch=20,col=alpha("darkorange2",0.5))

points(x=rep(4,length(PPbthhigh0.1mean))+rnorm(length(PPbthhigh0.1mean),0,0.1),
       y=PPbthhigh0.1mean,pch=20,col=alpha("darkorange2",0.5))

segments(1,10,1,12,lwd=3)
segments(2,6,2,12,lwd=3)
segments(2,12,1.7,12,lwd=3)
segments(1,12,1.3,12,lwd=3)
text(1.5,12,paste0("p=",formatC(PPmockLvH$p.value,format="e",digits = 2)),cex=0.7,srt=45,pos=1)

segments(3,8,3,12,lwd=3)
segments(4,6,4,12,lwd=3)
segments(3,12,3.3,12,lwd=3)
segments(4,12,3.7,12,lwd=3)
text(3.5,12,paste0("p=",round(PPbthLvH$p.value,digits = 3)),cex=0.7,srt=45,pos=1)

#####

#Negatives
plot(1:32,NNmeanRA_meanrankMock,pch=20,col="springgreen4",ylim=c(0.2,1.2),
     main="Mean parental rare alleles in negatively-predicting genes",
     xlab="F1 Expression Rank",ylab="Mean parental Rare-allele Count",xpd=NA)
points(1:32,NNmeanRA_meanrankBTH,pch=20,col="darkorange2")
lines(lowess(1:32,NNmeanRA_meanrankMock),lwd=3,col="springgreen4")
lines(lowess(1:32,NNmeanRA_meanrankBTH),lwd=3,col="darkorange2")
legend(27,1.2,legend=c("Mock","BTH"),pch=20,lwd=2,col=c("springgreen4","darkorange2"),bty="n")
mtext(text=c("LOW","HIGH"),side=1,line=3,at=c(0,31))

vioplot(NNmocklow0.1mean,NNmockhigh0.1mean,NNbthlow0.1mean,NNbthhigh0.1mean,
        col=alpha(rep(c("springgreen4","darkorange2"),each=2),0.7),ylim=c(0,13),
        main="Negative-negative genes mean low- vs high- rank RA",las=2,
        names=c("M.L.","M.H.","B.L.","B.H."), ylab="Mean rare-allele count")
set.seed(97)
points(x=rep(1,length(NNmocklow0.1mean))+rnorm(length(NNmocklow0.1mean),0,0.1),
       y=NNmocklow0.1mean,pch=20,col=alpha("springgreen4",0.5))

points(x=rep(2,length(NNmockhigh0.1mean))+rnorm(length(NNmockhigh0.1mean),0,0.1),
       y=NNmockhigh0.1mean,pch=20,col=alpha("springgreen4",0.5))

points(x=rep(3,length(NNbthlow0.1mean))+rnorm(length(NNbthlow0.1mean),0,0.1),
       y=NNbthlow0.1mean,pch=20,col=alpha("darkorange2",0.5))

points(x=rep(4,length(NNbthhigh0.1mean))+rnorm(length(NNbthhigh0.1mean),0,0.1),
       y=NNbthhigh0.1mean,pch=20,col=alpha("darkorange2",0.5))

segments(1,7,1,12,lwd=3)
segments(2,10,2,12,lwd=3)
segments(1,12,1.3,12,lwd=3)
segments(2,12,1.7,12,lwd=3)
text(1.5,12,paste0("p=",formatC(NNmockLvH$p.value,format="e",digits = 2)),cex=0.7,srt=45,pos=1)


segments(3,5,3,12,lwd=3)
segments(4,8,4,12,lwd=3)
segments(3,12,3.3,12,lwd=3)
segments(4,12,3.7,12,lwd=3)
text(3.5,12,paste0("p=",formatC(NNbthLvH$p.value,format="e",digits = 2)),cex=0.7,srt=45,pos=1)

dev.off()
##########
#plot rare-allele vs MPH
#do "mean MPH per rank of F1 expression" exhibit a pattern agains rare-allele counts?
#get mean size MPH value per trio * treatment, 
#use absolute MPH

trioMPHmean<-tapply(trioMPH$MPHmm2,paste0(trioMPH$TrioID,trioMPH$Treatment),mean)
trioMockmean<-trioMPHmean[grepl("Mock",names(trioMPHmean))]
names(trioMockmean)<-as.character(sapply(names(trioMockmean),function(x){substr(x,1,nchar(x)-4)}))
trioBTHmean<-trioMPHmean[grepl("BTH",names(trioMPHmean))]
names(trioBTHmean)<-as.character(sapply(names(trioBTHmean),function(x){substr(x,1,nchar(x)-3)}))

####PPclust
#repetitively sort size by F1 expression rank
SizeMPH_PPrankMock<-c()
SizeMPH_PPrankBTH<-c()#each row is a gene in the cluster, each column is rank-sorted MPH value
for(g in 1:nrow(PPclustF1rank)){ 
  # if(g%%10==1){print(g)}
  F1rankM<-sapply(PPclustF1rank[g,1:32],function(x){substr(x,start=2,stop=nchar(x)-1)}) #F1ranking by Mock
  F1rankB<-sapply(PPclustF1rank[g,33:64], function(x){substr(x,start=2,stop=nchar(x)-1)})#F1 ranking by BTH
  
  MPHrankM<-trioMockmean[as.numeric(sapply(F1rankM,function(x){which(names(trioMockmean)==x)}))]
  MPHrankB<-trioBTHmean[as.numeric(sapply(F1rankB,function(x){which(names(trioBTHmean)==x)}))]
  
  SizeMPH_PPrankMock<-rbind(SizeMPH_PPrankMock,MPHrankM) #467genes*32 rank
  SizeMPH_PPrankBTH<-rbind(SizeMPH_PPrankBTH,MPHrankB)
}

#mean over rank
SizeMPH_PPrankMock<-apply(SizeMPH_PPrankMock,2,mean)
SizeMPH_PPrankBTH<-apply(SizeMPH_PPrankBTH,2,mean)

#simple pearson correlation
corPPMock<-cor.test(meanRA_meanrankMock,SizeMPH_PPrankMock)
# Pearson's product-moment correlation
# 
# data:  meanRA_meanrankMock and SizeMPH_PPrankMock
# t = -5.3767, df = 30, p-value = 8.041e-06
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.8432427 -0.4655542
# sample estimates:
#        cor 
# -0.7005274 
corPPbth<-cor.test(meanRA_meanrankBTH,SizeMPH_PPrankBTH)
# Pearson's product-moment correlation
# 
# data:  meanRA_meanrankBTH and SizeMPH_PPrankBTH
# t = -2.3946, df = 30, p-value = 0.02309
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.65744656 -0.06031725
# sample estimates:
#       cor 
# -0.400586 

###NNclust
#repetitively sort size by F1 expression rank
SizeMPH_NNrankMock<-c()
SizeMPH_NNrankBTH<-c()#each row is a gene in the cluster, each column is rank-sorted MPH value
for(g in 1:nrow(NNclustF1rank)){ 
  #if(g%%10==1){print(g)}
  F1rankM<-sapply(NNclustF1rank[g,1:32],function(x){substr(x,start=2,stop=nchar(x)-1)}) #F1ranking by Mock
  F1rankB<-sapply(NNclustF1rank[g,33:64], function(x){substr(x,start=2,stop=nchar(x)-1)})#F1 ranking by BTH
  MPHrankM<-trioMockmean[as.numeric(sapply(F1rankM,function(x){which(names(trioMockmean)==x)}))]
  MPHrankB<-trioBTHmean[as.numeric(sapply(F1rankB,function(x){which(names(trioBTHmean)==x)}))]
  SizeMPH_NNrankMock<-rbind(SizeMPH_NNrankMock,MPHrankM)
  SizeMPH_NNrankBTH<-rbind(SizeMPH_NNrankBTH,MPHrankB)
}
#mean over rank
SizeMPH_NNrankMock<-apply(SizeMPH_NNrankMock,2,mean)
SizeMPH_NNrankBTH<-apply(SizeMPH_NNrankBTH,2,mean)

corNNmock<-cor.test(NNmeanRA_meanrankMock,SizeMPH_NNrankMock)
# Pearson's product-moment correlation
# 
# data:  NNmeanRA_meanrankMock and SizeMPH_NNrankMock
# t = -7.7505, df = 30, p-value = 1.2e-08
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.9070550 -0.6542773
# sample estimates:
#        cor 
# -0.8166571 
corNNbth<-cor.test(NNmeanRA_meanrankBTH,SizeMPH_NNrankBTH)
# Pearson's product-moment correlation
# 
# data:  NNmeanRA_meanrankBTH and SizeMPH_NNrankBTH
# t = -7.9793, df = 30, p-value = 6.612e-09
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.9111948 -0.6677142
# sample estimates:
#        cor 
# -0.8244526 

#####
#plot RA rank vs Size rank
pdf("run197_MPHpredict_PosNegClust_F1sizeMPH_rareallele.1.pdf",width=15,height=8)
par(mfrow=c(1,2))
plot(meanRA_meanrankMock,SizeMPH_PPrankMock,pch=20,col="springgreen4",ylim=c(5,60),xlim=c(0.35,0.6),
     main="F1 MPH vs mean parental rare alleles in positively-predicting genes",
     ylab="Mean F1 size MPH (mm2)",xlab="Mean parental Rare-allele Count")
points(meanRA_meanrankBTH,SizeMPH_PPrankBTH,pch=20,col="darkorange2")
lines(lowess(meanRA_meanrankMock,SizeMPH_PPrankMock),lwd=3,col="springgreen4")
lines(lowess(meanRA_meanrankBTH,SizeMPH_PPrankBTH),lwd=3,col="darkorange2")
legend(0.55,60,legend=c("Mock","BTH"),pch=20,lwd=2,col=c("springgreen4","darkorange2"),bty="n")
text(0.35,57,paste0("R=",round(corPPMock$estimate,digits=3),
                    ", p=",formatC(corPPMock$p.value,format="e",digits=2)),
      col="springgreen4",pos=4)
text(0.35,10,paste0("R=",round(corPPbth$estimate,digits=3),
                    ", p=",formatC(corPPbth$p.value,format="e",digits=2)),
     col="darkorange2",pos=4)

plot(NNmeanRA_meanrankMock,SizeMPH_NNrankMock,pch=20,col="springgreen4",ylim=c(5,70),xlim=c(0.25,0.55),
     main="F1 MPH vs mean parental rare alleles in negatively-predicting genes",
     ylab="Mean F1 size MPH (mm2)",xlab="Mean parental Rare-allele Count")
points(NNmeanRA_meanrankBTH,SizeMPH_NNrankBTH,pch=20,col="darkorange2")
lines(lowess(NNmeanRA_meanrankMock,SizeMPH_NNrankMock),lwd=3,col="springgreen4")
lines(lowess(NNmeanRA_meanrankBTH,SizeMPH_NNrankBTH),lwd=3,col="darkorange2")
text(0.25,67,paste0("R=",round(corNNmock$estimate,digits=3),
                    ", p=",formatC(corNNmock$p.value,format="e",digits=2)),
     col="springgreen4",pos=4)
text(0.25,20,paste0("R=",round(corNNbth$estimate,digits=3),
                    ", p=",formatC(corNNbth$p.value,format="e",digits=2)),
     col="darkorange2",pos=4)

dev.off()
