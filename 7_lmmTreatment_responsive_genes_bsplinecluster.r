##################
#README
#The script takes BTH treatment-responsive gnees (additive genes removed, fdr<0.001)
#and performs b-spline fitting of rosette size MPH (y-axis) against expression MPH of each gene (x-axis).
#The spline coefficients are then used for K-means clustering
##################
#Config
library(cluster)
library(clValid)
library(splines)
library(scales)

wd<-"/ebio/abt6_projects7/SHB2/data/1_analysis/6_BTHdiffgenes/"
outbasename<-"run197_RowFilter_BTHlmm_"
outbasename1<-paste0(outbasename,"rmadd_fdr0.001_")
metafile<-"Run197_rmCMRTPHBS_rm2MSeqtkDupOL_corr.meta.txt"
logTPMfile<-"run197_BTHdiffgene_rowfiltered_logTPM.txt"
Treatgenelist<-"run197_RowFilter_BTHlmm_rmadd_SimpleModel_FixEffect.txt"
MPHidxfile<-"run197_BTHdiffgene_MetaRowIdx_forpairedMPHcalc_manualcorr.txt"
###################
#data input
setwd(wd)
meta<-read.table(metafile,header=T,sep="\t",row.names=1)
meta$Treatment<-factor(meta$Treatment,levels=c("Mock","BTH"))
meta$IsHybrid<-factor(meta$IsHybrid,levels=c("Inbred","Hybrid"))

logTPM<-read.table(logTPMfile,header=T,sep="\t",row.names=1)
Treatlmm_rmadd<-read.table(Treatgenelist,header=T,sep="\t",row.names = 1)

MPHidx<-read.table(MPHidxfile,header=T,sep="\t")
MPHidxM<-which(MPHidx$Treatment=="Mock") # these are just row numbers, used for indexing/sorting/spline fitting in the subsequent step
MPHidxB<-which(MPHidx$Treatment=="BTH")
############################################
#Subset by fdr<0.001
BTHrmadd0.001<-Treatlmm_rmadd[Treatlmm_rmadd$padjBTH_rmadd<0.001,] #8102 genes
logTPM_rmadd0.001<-logTPM[rownames(logTPM)%in%rownames(BTHrmadd0.001),]

#calculate MPH for each of the 8K+ genes in each trio, then perform spline clustering with either all 8k+ genes
scaleBTHgenes<-t(scale(t(logTPM_rmadd0.001)))

calc_trioxpnMPH<-function(idx,gene){
  fxp<-gene[as.numeric(idx[7])]
  mxp<-gene[as.numeric(idx[8])]
  hxp<-gene[as.numeric(idx[9])]
  mphxp<-(hxp-0.5*(fxp+mxp))
}
xpnMPH<-t(apply(scaleBTHgenes,1,function(x){apply(MPHidx,1,calc_trioxpnMPH,gene=x)}))#each gene correspond to one row, within which is the MPH for each trio*treatment*rep

#########
#construct b-spline on quantiles of x (expression value), then cluster the coefficients
#note that spline should be fitted for mock and BTH separately
splinefitM<-list()
splinefitB<-list()
splinecoef<-data.frame()
for(i in 1:nrow(xpnMPH)){
  print(i)
  gene<-rownames(xpnMPH)[i]
  dfM<-data.frame("x"=xpnMPH[i,MPHidxM],y=MPHidx$MPHmm2[MPHidxM])
  dfB<-data.frame("x"=xpnMPH[i,MPHidxB],y=MPHidx$MPHmm2[MPHidxB])
  dfM<-dfM[order(dfM$x),]
  dfB<-dfB[order(dfB$x),]
  nsfitM<-ns(dfM$x,knots=quantile(dfM$x,probs=c(0.33,0.67)))
  nsfitB<-ns(dfB$x,knots=quantile(dfB$x,probs=c(0.33,0.67)))
  modelnsM<-lm(dfM$y~nsfitM)
  modelnsB<-lm(dfB$y~nsfitB)
  splinefitM[[i]]<-data.frame("x"=dfM$x,"spline"=fitted(modelnsM))
  splinefitB[[i]]<-data.frame("x"=dfB$x,"spline"=fitted(modelnsB))
  coef<-c(modelnsM$coefficients, modelnsB$coefficients)
  splinecoef<-rbind(splinecoef,coef)
}

rownames(splinecoef)<-rownames(xpnMPH)
colnames(splinecoef)<-paste0(rep(c("Intercept","nsfit1","nsfit2","nsfit3"),2),rep(c(".M",".B"),each=4))

#output spline coefficient
write.table(splinecoef,paste0(outbasename1,"splinecoefficient.txt"),quote=F,sep="\t")

#########
#clustering
maxk=100
genlist<-splinecoef
set.seed(109) 

#calculate elbow, silhouette, and dunn index
ds<-dist(genlist)
avg_elb.sil<-data.frame()
testclust<-c()
for(k in 2:maxk){
  print(k)
  km <- kmeans(genlist, centers = k, nstart=100,iter.max =3000000)
  testclust<-cbind(testclust,km$cluster)
  print("elbow")
  perc_within<-km$tot.withinss/km$totss
  print("silhouette")
  ss <- silhouette(km$cluster, ds)
  meanss<-mean(ss[, 3])
  print("Dunn")
  dunnIdx<-dunn(distance=ds,clusters=km$cluster)
  avg_elb.sil<-rbind(avg_elb.sil,c(perc_within,meanss,dunnIdx))
}

# compute gap statistic
 gap_stat <- clusGap(genlist, FUN = kmeans, nstart = 30, K.max = maxk, B = 50)
 print(gap_stat, method = "firstmax")

pdf(paste0(outbasename1,"MPHpredspline_clust_OptimalK.pdf"),width=10,height=15)
par(mfrow=c(3,1))
plot(2:maxk, type='l', avg_elb.sil[,1], main="Elbow plot",xlab='Number of clusters', ylab='Total within sum of squares', frame=FALSE) 
text(2:maxk,avg_elb.sil[,1],2:maxk,col="indianred4",cex=0.7,srt=45)
plot(2:maxk, type='l', avg_elb.sil[,2], main="Silhouette plot",xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE)
text(2:maxk,avg_elb.sil[,2],2:maxk,col="indianred4",cex=0.7,srt=45)
plot(2:maxk, type='l',avg_elb.sil[,3], main="Dunn Index",xlab="Number of clusters",ylab="Dunn Index",frame=F)
text(2:maxk,avg_elb.sil[,3],2:maxk,col="indianred4",cex=0.7,srt=45)
plot(1:maxk,gap_stat$Tab[,3],type="b",main="Gap statistic",xlab="Number of clusters",ylab="Gap statistic",frame=FALSE,pch=20)
arrows(x0=1:maxk,y0=(gap_stat$Tab[,3]-gap_stat$Tab[,4]),x1=1:maxk, y1=(gap_stat$Tab[,3]+gap_stat$Tab[,4]),code=3,angle=90,length=0.05)
dev.off()

colnames(testclust)<-paste0("K",2:maxk)
write.table(testclust,file=paste0(outbasename1,"MPHpredsplinek2-200clustID.txt"),quote=F,sep="\t")

################
#find strategy to determine optimal cluster numbers
#histogram cluster size per division
pdf(paste0(outbasename1,"2-200clustsize.pdf"),width=50,height=25)
par(mfrow=c(5,10))
for(i in 1:ncol(testclust)){
  k<-(i+1)
  clustsz<-table(testclust[,i])
  barplot(sort(clustsz),axes=F,main=paste0("K=",k),xlab="Cluster ID", ylab="No. cluster members",cex.main=2,cex.lab=2)
  axis(side=2,las=2)#one tick mark per 5% of genes
  for(quant in quantile(clustsz,probs=c(0.05,0.25,0.5,0.75,0.95))) {
    abline(h=quant,col="indianred4",lwd=3,lty=3)
  }
  abline(h=ceiling(nrow(testclust)*0.05), col="purple",lwd=3)
  text(0,ceiling(nrow(testclust)*0.05),"5%",pos=3)
  abline(h=ceiling(nrow(testclust)*0.025), col="purple",lwd=3)
  text(0,ceiling(nrow(testclust)*0.025),"2.5%",pos=3)
  abline(h=ceiling(nrow(testclust)*0.01), col="purple",lwd=3)
  text(0,ceiling(nrow(testclust)*0.01),"1%",pos=3)
}
dev.off()
#####
#plot results: cluster mean, individual cluster spline
#from the above result, picked local peaks of best clusters k=7,15,56,77

#prep for base spline model matrix. This is originated from values within the range of scaled expression MPH (-4,4)
rangeserie<-seq(-4,4,0.1)
nsrangeserie<-ns(rangeserie,knots=quantile(rangeserie,probs=c(0.33,0.67)))
rangeseriemat<-model.matrix(~nsrangeserie)  #use this to get fitted average spline value

K<-c(7,15,56,77)
for(k in K){
  print(k)
  clustID<-testclust[,(k-1)]
  print(table(clustID))
  
  #get the mean coefficient for each cluster
  ClustSplineCoef<-c()
  for(clust in 1:k){
    ClustSplineCoef<-rbind(ClustSplineCoef,apply(splinecoef[which(clustID==clust),],2,mean))
  }
  
  pdfname<-paste0(outbasename1,"bsplineK",k,"_MPHprediction.pdf")
  pdf(pdfname,width=60,height=35)
  par(mfrow=c(7,12)) #plot Mock and BTH side by side
  
  for(clust in 1:k){
    print(clust)
    print("Mock")
    #retrieve the average cluster coefficient
    meanCoefM<-ClustSplineCoef[clust,1:4]
    meanSplinefitM<-rangeseriemat%*%meanCoefM
    
    plot(seq(-4,4,length=10),seq(-50,200,length=10),xlim=c(-4,4),ylim=c(-80,220),type="n",
         xlab="expression MPH",ylab="rosette size MPH",main=paste0("Cluster", clust, " Mock n=", table(clustID)[clust]),
         cex.main=1.8,cex.lab=1.5)
    #first plot individual lines for each gene
    for(i in which(clustID==clust)){
      geneM<-splinefitM[[i]]
      lines(geneM$x,geneM$spline,col=alpha("springgreen4",0.15),lwd=0.3)
    }
    abline(h=0,lty=2,lwd=3)
    abline(v=0,lty=2,lwd=3)
    lines(rangeserie,meanSplinefitM,col="springgreen4",lwd=4)
    
    
    print("BTH")
    #retrieve the average cluster coefficient
    meanCoefB<-ClustSplineCoef[clust,5:8]
    meanSplinefitB<-rangeseriemat%*%meanCoefB
    
  #  plot(seq(-4,4,length=10),seq(-50,200,length=10),xlim=c(-4,4),ylim=c(-80,220),type="n",
  #       xlab="expression MPH",ylab="rosette size MPH",main=paste0("Cluster", clust, " BTH n=", table(clustID)[clust]))
    #first plot individual lines for each gene
    for(i in which(clustID==clust)){
      geneB<-splinefitB[[i]]
      lines(geneB$x,geneB$spline,col=alpha("darkorange2",0.15),lwd=0.3)
    }
    abline(h=0,lty=2,lwd=3)
    abline(v=0,lty=2,lwd=3)
    lines(rangeserie,meanSplinefitB,col="darkorange2",lwd=4)
    
  }
  
  dev.off()
  
}

################
#of course the more cluster there are, the more refined each cluster is. But the indi cluster size suffer
#use medium-level cluster numbers and make super cluster would be a good idea
#genes from the "flat" clusters are removed, then clustering was ran again on the remaining genes.