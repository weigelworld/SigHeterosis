################
#README
################
#Config
library(scales)
library(lme4)
library(RColorBrewer)
library(arm)
library(splines)
library(gplots)

wd<-"/ebio/abt6_projects9/SigHeterosis_Batch1/data/1_analysis/4_DEanalysis/run125_rmCMRTPHBS_seqtk"
metaname<-"run125_RSEM_rmCMRTPHBS_rmColOL3Mdup_seqtk28_Shootmetainfo.txt"
fulltriofile<-"run125_RSEM_rmCMRTPHBS_rmColOL3Mdup_seqtk28_ShootFullTrio_filtered.txt"
sortedclusterfile<-"run125_rmCMRTPHBS3MOLdup_seqtk28_Shootq0.001Kmeans_Optimalkmeans_sortedclustergeneID.txt"
env.path<-file.path(wd,"run125_rmCMRTPHBS3MOLdup_seqtk28_ShootHvMPV_q0.001Kmeans.elbsildunn.RData")
outbasename<-"run125_rmCMRTPHBS3MOLdup_seqtk28_Shootq0.001Kmeans"
myKpalette<-c("black",brewer.pal(12,"Paired"))
##################
setwd(wd)
meta<-read.table(metaname,header=T,sep="\t",row.names=1)
trio<-read.table(fulltriofile,header=T,sep="\t")
reassignclustdf<-read.table(sortedclusterfile,header=T,sep="\t",row.names=1)
load(env.path)

maxdepth<-which(meta$ttlMapped.rsem==max(meta$ttlMapped.rsem))
meta<-meta[-maxdepth,] #keep the seqtk sub-sampled entry
namelab<-meta$Genotype

genelist<-TPM[rownames(TPM)%in%rownames(testclust),]
genelist<-genelist[,1:nrow(meta)]
genelist<-apply(genelist,c(1,2),function(x){log2(x+1)})
colnames(genelist)<-namelab

reordermeta<-meta[order(meta$IsHybrid,meta$lastdaymm2),]
reordermeta[,"Size.z"]<-scale(reordermeta$lastdaymm2)
list1scale<-t(scale(t(genelist),center=T,scale=T))
HI_col_reorder<-ifelse(reordermeta$IsHybrid=="Inbred","purple4","turquoise4")

reassignclustdf<-read.table(sortedclusterfile,header=T,sep="\t",row.names=1)

#########fit spline models
set.seed(3295)
nsim=10000

sumSplineLMMlist<-list()
for(op in 1:ncol(reassignclustdf)){
  reassignclust<-reassignclustdf[,op]
  k<-optimalK[op]
  print(paste0("No.K group=",k))
  pdf(paste0(outbasename,k,"_perclustSplineLMM.pdf"),height=20,width=15)
  par(mfrow=c(4,1))
  list2scale<-list1scale[order(reassignclust),order(meta$IsHybrid,meta$lastdaymm2)]
  perclustmean<-apply(list2scale,2,function(x){tapply(x,sort(reassignclust),mean)})
  
  perKlist<-list()
  for(clust in 1:k){
    perclustlist<-list()
    
    print(paste0("Cluster ",clust))
    Kavg<-perclustmean[clust,]
    Kgenes<-list2scale[sort(reassignclust)==clust,]
    Kgenesvec<-c()
    for(i in 1:nrow(Kgenes)) Kgenesvec=c(Kgenesvec,Kgenes[i,])
    #build LMM data, z-transform the size
    LMMdf<-data.frame("Genotype"=rep(reordermeta$Genotype,nrow(Kgenes)),
                      "IsHybrid"=rep(reordermeta$IsHybrid,nrow(Kgenes)),
                      "LibBatch"=rep(reordermeta$LibBatch,nrow(Kgenes)),
                      "lastdaymm2"=rep(reordermeta$lastdaymm2,nrow(Kgenes)),
                      "Size.z"=rep(reordermeta$Size.z,nrow(Kgenes)),
                      "Gene"=as.factor(rep(1:nrow(Kgenes),each=nrow(reordermeta))),
                      "Genexpn"=Kgenesvec)
    ######fitting lmespline model#######
    print("modspline")
    modsp<-lmer(Genexpn~IsHybrid+ns(Size.z,knots=quantile(Size.z,probs=c(0.33,0.67)))+ns(Size.z,knots=quantile(Size.z,probs=c(0.33,0.67))):IsHybrid
                +(IsHybrid|LibBatch),data=LMMdf,REML=F,control = lmerControl(optimizer ="Nelder_Mead"))
    print("Bayesian")
    bsim.sp<-sim(modsp,n.sim=nsim)
    newdat.sp<-expand.grid(IsHybrid=levels(LMMdf$IsHybrid),lastdaymm2=seq(1,300,length=100))
    newdat.sp$Size.z<-(newdat$lastdaymm2-mean(LMMdf$lastdaymm2))/sd(LMMdf$lastdaymm2)
    Xmat.sp<-model.matrix(~IsHybrid+ns(Size.z,knots=quantile(Size.z,probs=c(0.33,0.67)))+ns(Size.z,knots=quantile(Size.z,probs=c(0.33,0.67))):IsHybrid,data=newdat.sp)
    fitmat.sp<-matrix(nrow=nsim,ncol=nrow(newdat.sp))
    for(j in 1:nsim) fitmat.sp[j,]<-Xmat.sp%*%bsim.sp@fixef[j,]
    colnames(fitmat.sp)<-paste0(rep(c("H_","I_"),100),round(newdat.sp$lastdaymm2,digits=2))
    newdat.sp$fitted<-Xmat.sp%*%fixef(modsp)
    newdat.sp$lower<-apply(fitmat.sp,2,quantile,prob=0.025)
    newdat.sp$upper<-apply(fitmat.sp,2,quantile,prob=0.975)
    
    perclustlist$splSummary<-summary(modsp)
    perclustlist$splFitmat<-fitmat.sp
    perclustlist$splCrI<-newdat.sp
    
    ###########plotting######
    print("plot")
    plot(reordermeta$lastdaymm2,perclustmean[clust,],col=HI_col_reorder,pch=ifelse(reordermeta$IsHybrid=="Hybrid",17,19),
         main=paste0("Kgroup",clust," spline"),xlab="last day rosette area (mm2)",ylab="scaled logTPM")
    #draw modLin fitted values for hybrid and inbred
    lines(newdat.sp[newdat.sp$IsHybrid=="Inbred","lastdaymm2"],newdat.sp[newdat.sp$IsHybrid=="Inbred","fitted"],col="purple4",lwd=3)
    lines(newdat.sp[newdat.sp$IsHybrid=="Inbred","lastdaymm2"],newdat.sp[newdat.sp$IsHybrid=="Inbred","upper"],col=alpha("purple4",0.8),lwd=1.2)
    lines(newdat.sp[newdat.sp$IsHybrid=="Inbred","lastdaymm2"],newdat.sp[newdat.sp$IsHybrid=="Inbred","lower"],col=alpha("purple4",0.8),lwd=1.2)
    
    lines(newdat.sp[newdat.sp$IsHybrid=="Hybrid","lastdaymm2"],newdat.sp[newdat.sp$IsHybrid=="Hybrid","fitted"],col="turquoise4",lwd=3)
    lines(newdat.sp[newdat.sp$IsHybrid=="Hybrid","lastdaymm2"],newdat.sp[newdat.sp$IsHybrid=="Hybrid","upper"],col=alpha("turquoise4",0.8),lwd=1.2)
    lines(newdat.sp[newdat.sp$IsHybrid=="Hybrid","lastdaymm2"],newdat.sp[newdat.sp$IsHybrid=="Hybrid","lower"],col=alpha("turquoise4",0.8),lwd=1.2)
    
    #add transparent polygon
    index<-newdat.sp$IsHybrid=="Inbred"
    polygon(c(newdat.sp$lastdaymm2[index],rev(newdat.sp$lastdaymm2[index])),c(newdat.sp$lower[index],rev(newdat.sp$upper[index])),
            border=NA,col=alpha("purple4",0.2))
    index<-newdat.sp$IsHybrid=="Hybrid"
    polygon(c(newdat.sp$lastdaymm2[index],rev(newdat.sp$lastdaymm2[index])),c(newdat.sp$lower[index],rev(newdat.sp$upper[index])),
            border=NA,col=alpha("turquoise4",0.2))
    
    ############use fitmat to calc difference between hybrid and inbred, use Hybrid - Inbred#####
    
    print("diffHvI")
    fitdiffHvI.sp<-fitmat.sp[,seq(1,ncol(fitmat.sp),2)]-fitmat.sp[,seq(2,ncol(fitmat.sp),2)] 
    colnames(fitdiffHvI.sp)<-sapply(colnames(fitdiffHvI.sp),sub,pattern="H",replacement="HvI")
    perclustlist$splDiffHvI<-fitdiffHvI.sp
    
    fitteddiff.sp<-(newdat.sp[seq(1,nrow(newdat.sp),2),]$fitted-newdat.sp[seq(2,nrow(newdat.sp),2),]$fitted)
    fittedlwr.sp<-apply(fitdiffHvI.sp,2,quantile,prob=0.025)
    fittedupr.sp<-apply(fitdiffHvI.sp,2,quantile,prob=0.975)
    
    #plot
    plot(reordermeta$lastdaymm2,perclustmean[clust,],type="n", ylim=c((min(fittedlwr.sp)-0.1),(max(fittedupr.sp)+0.1)),
         main=paste0("Kgroup",clust, " Hybrid - Inbred spline"),xlab="last day rosette area (mm2)",ylab="Scaled logTPM difference")
    #draw modLin fitted values for hybrid minus inbred
    lines(unique(newdat.sp$lastdaymm2),fitteddiff.sp,type="l",col="grey30",lwd=4)
    lines(unique(newdat.sp$lastdaymm2),fittedlwr.sp,col="grey60",lwd=1.2)
    lines(unique(newdat.sp$lastdaymm2),fittedupr.sp,col="grey60",lwd=1.2)
    #add transparent polygon
    polygon(c(unique(newdat.sp$lastdaymm2),rev(unique(newdat.sp$lastdaymm2))),c(fittedlwr.sp,rev(fittedupr.sp)),
            border=NA,col=alpha("grey30",0.2))
    #add zero line
    abline(h=0,col="grey30",lwd=3,lty=2)
    
    perKlist[[clust]]<-perclustlist
    print(c("namesperclustlist", names(perclustlist)))
    print(length(perKlist))
  }
  dev.off()
  sumSplineLMMlist[[op]]<-perKlist  
}
###############
#structure of sumSplineLMMlist
#[[optimal K divisions]] -- length=6, 6 different K numbers tested
#    [[model and Bayesian statistics for all clusters given a K division ]] -- length= K number being tested in a particular iteration
          # [[lmm Summary]]
          # [[Bayesian fitmat based on lmm]]
          # [[Bayesian CrI based on lmm]]
          # [[extapolated difference between hybrid and inbred, based on lmm]]
          # [[splinelmm Summary]]
          # [[Bayesian fitmat based on splinelmm]]
          # [[Bayesian CrI based on splinelmm]]
          # [[extapolated difference between hybrid and inbred, based on splinelmm]]
#############
#correlation between clusters, it is unecessary to plot k=3
myDistpalette<-brewer.pal(7,"BrBG")  
for(op in 2:ncol(reassignclustdf)){
  reassignclust<-reassignclustdf[,op]
  k<-optimalK[op]
  print(paste0("No.K group=",k))
  list2scale<-list1scale[order(reassignclust),order(meta$IsHybrid,meta$lastdaymm2)]
  perclustmean<-apply(list2scale,2,function(x){tapply(x,sort(reassignclust),mean)})
  perclustsize<-table(reassignclust)
  
  #pearson correlation of clusters
  print("calc pearson")
  clustcor<-cor(t(perclustmean),method="pearson")
  clustlabel<-paste0("K",rownames(clustcor)," n=",perclustsize)
  clustcordist<-as.dist(1-clustcor) #measure distance between cluster means as 1-pearson correlation value
  hclust.row<-hclust(clustcordist,method="complete")
  ##manually scale the color scheme for pearson correlation
  colbreaks<-c(-1,-0.8,-0.5,-0.2,0.2,0.5,0.8,1)
  #plotting
  print("plot")
  pdf(paste0(outbasename,"_Kgroup",k,"_clustcor.pdf"),height=15,width=15)
  
  heatmap.2(clustcor,breaks=colbreaks,dendrogram="row",
            Rowv=as.dendrogram(hclust.row),Colv=as.dendrogram(hclust.row),
            labRow=clustlabel,labCol="",trace="none",
            main=paste0("Pearson correlation cluster mean, k=",k),
            cexRow=1.5,cex.main=9,margins=c(3,10),
            density.info="density",key.xlab="Pearson correlation",keysize = 0.8,key.par=list(cex.lab=1,cex.main=1.6),
            col=myDistpalette
  )
  dev.off()
}
