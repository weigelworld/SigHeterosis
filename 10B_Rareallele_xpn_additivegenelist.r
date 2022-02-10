###########
#README
#The script perform same kind of analysis as 10A, but on common additive genes
###########
#Config and data input
setwd("/ebio/abt6_projects7/SHB2/data/1_analysis/10_Subgenome_Genetdist")
additive<-read.table("/ebio/abt6_projects7/SHB2/data/1_analysis/6_BTHdiffgenes/run125_run197_loosethreshold_additiveintersect.txt",header=F,sep="\t")[,1] #only read in the 300 over-lapping additive genes
Acc.rareallele<-read.table("run197_BTHdiffgene_1kbupstreamONLY_MAF0.005_SHB2acc_pergeneRareallelecount.txt",header=T,sep="\t",row.names=1)
seqmeta<-read.table("/ebio/abt6_projects7/SHB2/data/1_analysis/6_BTHdiffgenes/Run197_rmCMRTPHBS_rm2MSeqtkDupOL_corr.meta.txt",header=T,sep="\t",row.names=1)
trioMPH<-read.table("/ebio/abt6_projects7/SHB2/data/1_analysis/6_BTHdiffgenes/run197_BTHdiffgene_MetaRowIdx_forpairedMPHcalc_manualcorr.txt",
                    header=T,sep="\t")
###########

#import the expression value
ScaledTPM<-read.table("/ebio/abt6_projects7/SHB2/data/1_analysis/6_BTHdiffgenes/run197_BTHdiffgene_rowfiltered_scaledlogTPM.txt",header=T,sep="\t",row.names=1)
#shorten the gene list
ScaledTPM<-ScaledTPM[rownames(ScaledTPM)%in% additive,] 


#index the columns with accession samples (Mock + BTH)
meta.acc<-subset(seqmeta,IsHybrid=="Inbred")
meta.F1<-subset(seqmeta,IsHybrid=="Hybrid") #the meta file specified the

ScaledTPM.acc<-ScaledTPM[,which(seqmeta$IsHybrid=="Inbred")]
ScaledTPM.F1<-ScaledTPM[,which(seqmeta$IsHybrid=="Hybrid")]

#get accession expression value and RA of the additive genes
AddRAmock<-matrix(nrow=length(additive),ncol=64)
AddRAbth<-matrix(nrow=length(additive),ncol=64)
AddnoRA<-c()
j<-1
for(g in additive){
  xpn<-as.numeric(ScaledTPM.acc[rownames(ScaledTPM.acc)==g,])
  if(!any(rownames(Acc.rareallele)==g)){
    print(paste0(g," has no upstream rare allele"))
    AddnoRA<-c(AddnoRA,g)
    next
  }
  RA<-Acc.rareallele[rownames(Acc.rareallele)==g,]
  meanxpn<-tapply(xpn,paste0(meta.acc$Genotype,"_",meta.acc$Treatment),mean)
  mockxpn<-meanxpn[grepl("Mock",names(meanxpn))]
  BTHxpn<-meanxpn[grepl("BTH",names(meanxpn))]
  mockrank<-as.character(sapply(names(sort(rank(mockxpn))),function(x){substr(x,1,nchar(x)-5)})) #rank value from low to high
  BTHrank<-as.character(sapply(names(sort(rank(BTHxpn))),function(x){substr(x,1,nchar(x)-4)}))
  RAmock<-RA[as.numeric(sapply(mockrank,function(x){which(names(RA)==x)}))]
  RAbth<-RA[as.numeric(sapply(BTHrank,function(x){which(names(RA)==x)}))]
  AddRAmock[j,]<-as.numeric(RAmock)
  AddRAbth[j,]<-as.numeric(RAbth)
  j<-(j+1)
}
AddRAmock<-na.omit(AddRAmock) #211/300 genes have rare alleles
AddRAbth<-na.omit(AddRAbth)
Add_meanRAmock<-apply(AddRAmock,2,mean) #one mean RA per cluster
Add_meanRAbth<-apply(AddRAbth,2,mean)
cat(AddnoRA,file="run125_run197_intersectAddtivegenes_noRareAllele_geneID.txt",sep="\n")

#Wilcoxon test for accession RA-rank relationship
Addmocklow0.1mean<-apply(AddRAmock[,1:7],1,mean)
Addmockhigh0.1mean<-apply(AddRAmock[,(ncol(AddRAmock)-6):ncol(AddRAmock)],1,mean)
AddmockLvH<-wilcox.test(Addmocklow0.1mean,Addmockhigh0.1mean,paired=T,alternative = "greater")
# Wilcoxon signed rank test with continuity correction
# 
# data:  Addmocklow0.1mean and Addmockhigh0.1mean
# V = 10655, p-value = 5.158e-06
# alternative hypothesis: true location shift is greater than 0

Addbthlow0.1mean<-apply(AddRAbth[,1:7],1,mean)
Addbthhigh0.1mean<-apply(AddRAbth[,(ncol(AddRAbth)-6):ncol(AddRAbth)],1,mean)
AddbthLvH<-wilcox.test(Addbthlow0.1mean,Addbthhigh0.1mean,paired=T,alternative = "greater")
# Wilcoxon signed rank test with continuity correction
# 
# data:  Addbthlow0.1mean and Addbthhigh0.1mean
# V = 10536, p-value = 5.281e-06
# alternative hypothesis: true location shift is greater than 0

######
#plotting part1
pdf("run197_MPHpredict_300Additivegenes_expressionrank_rareallele.1.pdf",width=15,height=30)
par(mfrow=c(4,2))
plot(1:64,Add_meanRAmock,pch=20,col="springgreen4",ylim=c(0.3,1.8),
     main="Average of all additive genes in Accessions",
     xlab="Expression Rank of Accessions",ylab="Mean Rare-allele Count",xpd=NA)
points(1:64,Add_meanRAbth,pch=20,col="darkorange2")
lines(lowess(1:64,Add_meanRAmock),lwd=3,col="springgreen4")
lines(lowess(1:64,Add_meanRAbth),lwd=3,col="darkorange2")
legend(0,1.8,legend=c("Mock","BTH"),pch=20,lwd=2,col=c("springgreen4","darkorange2"),bty="n")
mtext(text=c("LOW","HIGH"),side=1,line=3,at=c(0,62))

vioplot(Addmocklow0.1mean,Addmockhigh0.1mean,Addbthlow0.1mean,Addbthhigh0.1mean,
        col=alpha(rep(c("springgreen4","darkorange2"),each=2),0.7),ylim=c(0,20),
        main="Additive genes mean low- vs high- rank RA",las=2,
        names=c("M.L.","M.H.","B.L.","B.H."), ylab="Mean rare-allele count")
set.seed(97)
points(x=rep(1,length(Addmocklow0.1mean))+rnorm(length(Addmocklow0.1mean),0,0.1),
       y=Addmocklow0.1mean,pch=20,col=alpha("springgreen4",0.5))

points(x=rep(2,length(Addmockhigh0.1mean))+rnorm(length(Addmockhigh0.1mean),0,0.1),
       y=Addmockhigh0.1mean,pch=20,col=alpha("springgreen4",0.5))

points(x=rep(3,length(Addbthlow0.1mean))+rnorm(length(Addbthlow0.1mean),0,0.1),
       y=Addbthlow0.1mean,pch=20,col=alpha("darkorange2",0.5))

points(x=rep(4,length(Addbthhigh0.1mean))+rnorm(length(Addbthhigh0.1mean),0,0.1),
       y=Addbthhigh0.1mean,pch=20,col=alpha("darkorange2",0.5))

segments(1,16,1,18,lwd=3)
segments(2,8,2,18,lwd=3)
segments(1,18,2,18,lwd=3)
text(1.5,19,paste0("p=",formatC(AddmockLvH$p.value,format="e",digits = 2)),cex=0.7)

segments(3,16,3,18,lwd=3)
segments(4,13,4,18,lwd=3)
segments(3,18,4,18,lwd=3)
text(3.5,19,paste0("p=",formatC(AddbthLvH$p.value,format="e",digits = 2)),cex=0.7)
#######

#model F1 expression with parental rare-allele counts
#use mean rare alleles from both parents

#calculate average expression rank for each gene across replicates, in Mock and BTH,
#and extract the parental RA count
AddRAacc<-matrix(nrow=length(additive),ncol=64)
AddF1xpnRank<-matrix(nrow=length(additive),ncol=64) #first 32 columns for mock, and the remaining 32 columns for BTH
j<-1
for(g in additive){
  xpn<-as.numeric(ScaledTPM.F1[rownames(ScaledTPM.F1)==g,])
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
  AddRAacc[j,]<-as.numeric(RA)
  AddF1xpnRank[j,]<-c(paste0("M",mockrank),paste0("B",BTHrank))
  j<-(j+1)
}
AddRAacc<-na.omit(AddRAacc)
AddF1xpnRank<-na.omit(AddF1xpnRank)

# now we have a list of accession RAs and a list of F1 xpn ranks, both in corresponding order
acc2trio<-sapply(colnames(Acc.rareallele),function(x){unique(seqmeta$Trio[which(seqmeta$Genotype==x)])})

# calculate for each gene-F1 xpn rank: sum of rare alleles in parents
trioRA<-sapply(acc2trio,function(x){substr(x,1,nchar(x)-1)})
meanRA<-t(apply(AddRAacc,1,function(x){sapply(unique(trioRA),function(y){mean(x[trioRA==y])})})) #each row is a PP gene, each column the sum RA of a trio

#sort RA by rank for each gene,note that the sum/mean/max/min/diffRA table right now all have the same column order (ie unique(trioRA)), and the F1 rank will also be shared
sortingrankMock<-c()
sortingrankBTH<-c()#each row is a gene in the PP cluster, each column is a given rank.
for(g in 1:nrow(AddRAacc)){ 
  #if(g%%10==1){print(g)}
  RA<-AddRAacc[g,]
  F1rankM<-sapply(AddF1xpnRank[g,1:32],function(x){substr(x,start=2,stop=nchar(x)-1)}) #F1ranking by Mock
  F1rankB<-sapply(AddF1xpnRank[g,33:64], function(x){substr(x,start=2,stop=nchar(x)-1)})#F1 ranking by BTH
  RArankM<-as.numeric(sapply(F1rankM,function(x){which(unique(trioRA)==x)}))
  RArankB<-as.numeric(sapply(F1rankB,function(x){which(unique(trioRA)==x)}))
  sortingrankMock<-rbind(sortingrankMock,RArankM)
  sortingrankBTH<-rbind(sortingrankBTH,RArankB)
}

#rank each gene's RA by the table, calculate mean over each given rank
meanRA_rankMock<-t(sapply(1:nrow(AddRAacc),function(x){return(as.numeric(meanRA[x,sortingrankMock[x,]]))})) #each row is a gene, each column is a rank (low->high)
meanRA_rankBTH<-t(sapply(1:nrow(AddRAacc),function(x){return(as.numeric(meanRA[x,sortingrankBTH[x,]]))})) #each row is a gene, each column is a rank (low->high)
meanRA_meanrankMock<-apply(meanRA_rankMock,2,mean)
meanRA_meanrankBTH<-apply(meanRA_rankBTH,2,mean)

Addmocklow0.1mean<-apply(meanRA_rankMock[,1:4],1,mean)
Addmockhigh0.1mean<-apply(meanRA_rankMock[,(ncol(meanRA_rankMock)-3):ncol(meanRA_rankMock)],1,mean)
AddmockLvH<-wilcox.test(Addmocklow0.1mean,Addmockhigh0.1mean,paired=T,alternative = "greater")
# Wilcoxon signed rank test with continuity correction
# 
# data:  Addmocklow0.1mean and Addmockhigh0.1mean
# V = 10166, p-value = 2.963e-05
# alternative hypothesis: true location shift is greater than 0

Addbthlow0.1mean<-apply(meanRA_rankBTH[,1:4],1,mean)
Addbthhigh0.1mean<-apply(meanRA_rankBTH[,(ncol(meanRA_rankBTH)-3):ncol(meanRA_rankBTH)],1,mean)
AddbthLvH<-wilcox.test(Addbthlow0.1mean,Addbthhigh0.1mean,paired=T,alternative = "greater")
# Wilcoxon signed rank test with continuity correction
# 
# data:  Addbthlow0.1mean and Addbthhigh0.1mean
# V = 10460, p-value = 0.0002553
# alternative hypothesis: true location shift is greater than 0

#######
#plotting part 2 (continued)
plot(1:32,meanRA_meanrankMock,pch=20,col="springgreen4",ylim=c(0.3,1.5),
     main="Mean parental rare alleles in additive genes",
     xlab="F1 Expression Rank",ylab="Mean parental Rare-allele Count")
points(1:32,meanRA_meanrankBTH,pch=20,col="darkorange2")
lines(lowess(1:32,meanRA_meanrankMock),lwd=3,col="springgreen4")
lines(lowess(1:32,meanRA_meanrankBTH),lwd=3,col="darkorange2")
legend(1,1.5,legend=c("Mock","BTH"),pch=20,lwd=2,col=c("springgreen4","darkorange2"),bty="n")
mtext(text=c("LOW","HIGH"),side=1,line=3,at=c(0,31))

vioplot(Addmocklow0.1mean,Addmockhigh0.1mean,Addbthlow0.1mean,Addbthhigh0.1mean,
        col=alpha(rep(c("springgreen4","darkorange2"),each=2),0.7),ylim=c(0,20),
        main="Additive genes mean low- vs high- rank RA",las=2,
        names=c("M.L.","M.H.","B.L.","B.H."), ylab="Mean rare-allele count")
set.seed(97)
points(x=rep(1,length(Addmocklow0.1mean))+rnorm(length(Addmocklow0.1mean),0,0.1),
       y=Addmocklow0.1mean,pch=20,col=alpha("springgreen4",0.5))

points(x=rep(2,length(Addmockhigh0.1mean))+rnorm(length(Addmockhigh0.1mean),0,0.1),
       y=Addmockhigh0.1mean,pch=20,col=alpha("springgreen4",0.5))

points(x=rep(3,length(Addbthlow0.1mean))+rnorm(length(Addbthlow0.1mean),0,0.1),
       y=Addbthlow0.1mean,pch=20,col=alpha("darkorange2",0.5))

points(x=rep(4,length(Addbthhigh0.1mean))+rnorm(length(Addbthhigh0.1mean),0,0.1),
       y=Addbthhigh0.1mean,pch=20,col=alpha("darkorange2",0.5))

segments(1,16,1,18,lwd=3)
segments(2,8,2,18,lwd=3)
segments(1,18,2,18,lwd=3)
text(1.5,19,paste0("p=",formatC(AddmockLvH$p.value,format="e",digits = 2)),cex=0.7)

segments(3,16,3,18,lwd=3)
segments(4,13,4,18,lwd=3)
segments(3,18,4,18,lwd=3)
text(3.5,19,paste0("p=",formatC(AddbthLvH$p.value,format="e",digits = 2)),cex=0.7)
#######

#rare-allele vs MPH
#do "mean MPH per rank of F1 expression" exhibit a pattern agains rare-allele counts?
trioMPHmean<-tapply(trioMPH$MPHmm2,paste0(trioMPH$TrioID,trioMPH$Treatment),mean)
trioMockmean<-trioMPHmean[grepl("Mock",names(trioMPHmean))]
names(trioMockmean)<-as.character(sapply(names(trioMockmean),function(x){substr(x,1,nchar(x)-4)}))
trioBTHmean<-trioMPHmean[grepl("BTH",names(trioMPHmean))]
names(trioBTHmean)<-as.character(sapply(names(trioBTHmean),function(x){substr(x,1,nchar(x)-3)}))

#repetitively sort size by F1 expression rank
SizeMPH_AddrankMock<-c()
SizeMPH_AddrankBTH<-c()#each row is a gene in the cluster, each column is rank-sorted MPH value
for(g in 1:nrow(AddF1xpnRank)){ 
  if(g%%10==1){print(g)}
  F1rankM<-sapply(AddF1xpnRank[g,1:32],function(x){substr(x,start=2,stop=nchar(x)-1)}) #F1ranking by Mock
  F1rankB<-sapply(AddF1xpnRank[g,33:64], function(x){substr(x,start=2,stop=nchar(x)-1)})#F1 ranking by BTH
  MPHrankM<-trioMockmean[as.numeric(sapply(F1rankM,function(x){which(names(trioMockmean)==x)}))]
  MPHrankB<-trioBTHmean[as.numeric(sapply(F1rankB,function(x){which(names(trioBTHmean)==x)}))]
  SizeMPH_AddrankMock<-rbind(SizeMPH_AddrankMock,MPHrankM)
  SizeMPH_AddrankBTH<-rbind(SizeMPH_AddrankBTH,MPHrankB)
}
#mean over rank
SizeMPH_AddrankMock<-apply(SizeMPH_AddrankMock,2,mean)
SizeMPH_AddrankBTH<-apply(SizeMPH_AddrankBTH,2,mean)

#simple pearson correlation
corAddMock<-cor.test(meanRA_meanrankMock,SizeMPH_AddrankMock)
# Pearson's product-moment correlation
# 
# data:  meanRA_meanrankMock and SizeMPH_AddrankMock
# t = -0.88823, df = 30, p-value = 0.3815
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.4818735  0.1997680
# sample estimates:
#        cor 
# -0.1600766 

corAddbth<-cor.test(meanRA_meanrankBTH,SizeMPH_AddrankBTH)
# Pearson's product-moment correlation
# 
# data:  meanRA_meanrankBTH and SizeMPH_AddrankBTH
# t = -1.9861, df = 30, p-value = 0.05622
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.61632877  0.00884879
# sample estimates:
#       cor 
# -0.340897 

#####
#plotting part 4 (continued)
plot(meanRA_meanrankMock,SizeMPH_AddrankMock,pch=20,col="springgreen4",ylim=c(5,60),xlim=c(0.3,1.2),
     main="F1 MPH vs mean parental rare alleles in additive genes",
     ylab="Mean F1 size MPH (mm2)",xlab="Mean parental Rare-allele Count")
points(meanRA_meanrankBTH,SizeMPH_AddrankBTH,pch=20,col="darkorange2")
lines(lowess(meanRA_meanrankMock,SizeMPH_AddrankMock),lwd=3,col="springgreen4")
lines(lowess(meanRA_meanrankBTH,SizeMPH_AddrankBTH),lwd=3,col="darkorange2")
legend(1,60,legend=c("Mock","BTH"),pch=20,lwd=2,col=c("springgreen4","darkorange2"),bty="n")

text(0.3,57,paste0("R=",round(corAddMock$estimate,digits=3),
                    ", p=",formatC(corAddMock$p.value,format="e",digits=2)),
     col="springgreen4",pos=4)
text(0.3,10,paste0("R=",round(corAddbth$estimate,digits=3),
                    ", p=",formatC(corAddbth$p.value,format="e",digits=2)),
     col="darkorange2",pos=4)


dev.off()
########
#generate FDR/bonferroni correction for all wilcoxon p-values
group<-paste(rep(c("Add","Pos","Neg"),each=4),
             rep(rep(c("Acc","F1"),each=2),3),
             rep(rep(c("M","B"),2),3),sep=".")

pval<-c(5.16e-6,5.28e-6,2.96e-5,2.55e-4,0.028,0.087,3.7e-3,0.225,2.17e-6,1.76e-7,2.46e-13,5.58e-8)
pvaldf<-data.frame(group,pval)
pvaldf[,"fdr"]<-formatC(p.adjust(pval,"BH"),format="e",digits=2)
pvaldf[,"bonferroni"]<-formatC(p.adjust(pval,"bonferroni"),format="e",digits=2)

