################
#README
#The script takes minimally filtered logTPM file and perform gene-by-gene LMM.
#between expression in F1s (y-axis) and Mid-Parent Values of the corresponding parents.

################
#Config
library(RColorBrewer)
library(scales)
library(lme4)
wd<-"/ebio/abt6_projects7/SHB2/data/1_analysis/6_BTHdiffgenes/"
outbasename<-"run197_RowFilter_addpredict_"
triocol<-colorRampPalette(brewer.pal(8,"Set2"))
metafile<-"Run197_rmCMRTPHBS_rm2MSeqtkDupOL_corr.meta.txt"
logTPMfile<-"run197_BTHdiffgene_rowfiltered_logTPM.txt"
MPHidxfile<-"run197_BTHdiffgene_MetaRowIdx_forpairedMPHcalc_manualcorr.txt"
#################

#data input
setwd(wd)
meta<-read.table(metafile,header=T,sep="\t",row.names=1)
logTPM<-read.table(logTPMfile,header=T,sep="\t",row.names=1)
MPHidx<-read.table(MPHidxfile,header=T,sep="\t")
meta$Treatment<-factor(meta$Treatment,levels=c("Mock","BTH"))
meta$IsHybrid<-factor(meta$IsHybrid,levels=c("Inbred","Hybrid"))
scalelogTPM<-t(scale(t(logTPM)))
MPHidx$Treatment<-factor(MPHidx$Treatment, levels=c("Mock","BTH"))

####################
#calculate MPV(mid-parent value) and do regression with actual F1 values
calc_trioMPVpredict<-function(idx,gene){
  fxp<-gene[as.numeric(idx[7])]
  mxp<-gene[as.numeric(idx[8])]
  hxp<-gene[as.numeric(idx[9])]
  mphxp<-c(hxp,0.5*(fxp+mxp))
  return(mphxp)
}
trioMPVpredict<-t(apply(scalelogTPM,1,function(x){unlist(apply(MPHidx,1,calc_trioMPVpredict,gene=x))}))
#each gene is a row, 
#and on the columns it's the alternate F1 and MPV for each trio*treatment*rep (correspond to the 177 rows of MPHidx)
#odd rows are F1s, and even rows are MPVs

##lmm to identify genes generally showing additive expression in hybrids
lmmdf<-data.frame()
for(i in 1:nrow(trioMPVpredict)){
  gene<-trioMPVpredict[i,]
  f1<-gene[seq(1,length(gene),2)]
  mpv<-gene[seq(2,length(gene),2)]
  fit<-lmer(f1~MPHidx$Treatment+mpv+MPHidx$Treatment*mpv+(1|MPHidx$Traycode)+(1|MPHidx$PlantBatch)+(1|MPHidx$TrioID))
  returnvec<-c(fixef(fit),summary(fit)$sigma)
  lmmdf<-rbind(lmmdf,returnvec)
}
####
names(lmmdf)<-c("Mock","BTH","MPV","BTHxMPV","Sigma")
pdf(paste0(outbasename,"_mpvlmm_additivepredict.pdf"),height=20, width=30)
par(mfrow=c(2,3))
vioplot(lmmdf$Mock,(lmmdf$BTH+lmmdf$Mock),col=c("springgreen4","darkorange2"),
        main="Treatment effect size on F1 expression",names=c("Mock","BTH"),xlab="Treatment",ylab="Effectsize")
abline(h=0,lty=2)

vioplot(lmmdf$MPV,(lmmdf$MPV+lmmdf$BTHxMPV),col=c("springgreen4","darkorange2"),
        main="F1 regression on MPV, treatment:MPV",names=c("Mock","BTH"),xlab="Treatment",ylab="MPV slope")
abline(h=0,lty=2)
#hist(lmmdf$MPV,breaks=100)
#hist(lmmdf$BTHxMPV,breaks=100)
hist(lmmdf$Sigma,breaks=100)
plot(lmmdf$Mock, (lmmdf$Mock+lmmdf$BTH), pch=20,cex=0.5, col="grey30",main="Effectsize Treatment",xlab="Mock effect",ylab="net BTH effect")
plot(lmmdf$MPV, (lmmdf$MPV+lmmdf$BTHxMPV),pch=20,cex=0.5, col="grey30",main="Treatment-specific MPV regression slope",xlab="Mock slope",ylab="BTH slope")
dev.off()
#once controlled for treatment effect, MPV regression is much like in SHB1, that is between -0.5~1.1
#filter then by slope and residual
print(quantile(lmmdf$MPV))
##with original row filter
# 0%        25%        50%        75%       100% 
# -0.6014054  0.1276357  0.2544798  0.3836101  1.0969998

#with stringent row filter
# 0%        25%        50%        75%       100% 
# -0.6014054  0.1170576  0.2434771  0.3751693  1.0142178 

#Run197 row filtered
# 0%        25%        50%        75%       100% 
# -0.5633489  0.1057889  0.2353029  0.3637190  1.0466364 

print(quantile(lmmdf$MPV+lmmdf$BTHxMPV))
##with original row filter
# 0%        25%        50%        75%       100% 
# -0.4928438  0.1149424  0.2501998  0.3831608  1.1222129 

#with stringent row filter
# 0%        25%        50%        75%       100% 
# -0.4928437  0.1061561  0.2420768  0.3770710  1.0491165 

#Run197 row filtered
# 0%        25%        50%        75%       100% 
# -0.6954263  0.1222936  0.2629245  0.4002843  1.0729535 

print(quantile(lmmdf$Sigma))
##with original row filter
# 0%       25%       50%       75%      100% 
# 0.1295077 0.6549362 0.7704809 0.8617275 1.2684903 

#with stringent row filter
# 0%       25%       50%       75%      100% 
# 0.1295077 0.6602045 0.7788758 0.8682892 1.2684903 

#Run197 row filtered
# 0%       25%       50%       75%      100% 
# 0.1325550 0.6335523 0.7569564 0.8549827 1.1528697 

AdditiveGenes<-rownames(trioMPVpredict)[(lmmdf$MPV+lmmdf$BTHxMPV)>0.4 & lmmdf$MPV>0.4 & lmmdf$Sigma<0.63]
AdditiveGenes2<-rownames(trioMPVpredict)[(lmmdf$MPV+lmmdf$BTHxMPV)>0.5 & lmmdf$MPV>0.5 & lmmdf$Sigma<0.6]
cat(AdditiveGenes, file=paste0(outbasename,"slp0.4Sigma0.65_geneID.txt"),sep="\n") #850 genes with run 189, 901 genes in run 197 
cat(AdditiveGenes2, file=paste0(outbasename,"slp0.5Sigma0.6_geneID.txt"),sep="\n") #270 genes with run 189, 309 genes in run 197 

#############PCA diagnostics####

#make PCA of the filtered and scaled genelist
TPMpca<-prcomp(scalelogTPM, center=T,scale.=T)
sizepal<-c("goldenrod3","#c7e9b4","#7fcdbb","#41b6c4","#2c7fb8","#253494")
sizecol<-ifelse(meta$Day13mpp<=quantile(meta$Day13mpp,probs=0.05),sizepal[1],
                ifelse(meta$Day13mpp<=quantile(meta$Day13mpp,probs=0.25),sizepal[2],
                       ifelse(meta$Day13mpp<=quantile(meta$Day13mpp,probs=0.5),sizepal[3],
                              ifelse(meta$Day13mpp<=quantile(meta$Day13mpp,probs=0.75),sizepal[4],
                                     ifelse(meta$Day13mpp<=quantile(meta$Day13mpp,probs=0.95),sizepal[5],sizepal[6])))))
pdf(paste0(outbasename,"_pca.pdf"),height=30,width=15)
par(mfrow=c(6,3))
#%byPC
plot(x=1:ncol(summary(TPMpca)[[6]]),y=summary(TPMpca)[[6]][2,],
     main="Individual PC Variance Explained",
     xlab="PCs",ylab="Proportion Variance Explained",pch=19,xlim=c(0,450),ylim=c(0,1),cex=2)
plot(x=1:ncol(summary(TPMpca)[[6]]),y=summary(TPMpca)[[6]][3,],
     main="cumulative percentage by PC",
     xlab="PCs",ylab=,pch=19,xlim=c(0,450),ylim=c(0,1),cex=2)
h<-hist(meta$Day13mpp,breaks=50, plot=F)
hist(meta$Day13mpp,breaks=50,
     col=ifelse(h$breaks<=quantile(meta$Day13mpp,probs=0.05),sizepal[1],
                ifelse(h$breaks<=quantile(meta$Day13mpp,probs=0.25),sizepal[2],
                      ifelse(h$breaks<=quantile(meta$Day13mpp,probs=0.5),sizepal[3],
                            ifelse(h$breaks<=quantile(meta$Day13mpp,probs=0.75),sizepal[4],
                                  ifelse(h$breaks<=quantile(meta$Day13mpp,probs=0.95),sizepal[5],sizepal[6]))))))

for (j in 1:5){
  #size color, HI pch
  plot(x=summary(TPMpca)[[2]][,j],y=summary(TPMpca)[[2]][,(j+1)],
       xlab=paste0("PC",j,":",summary(TPMpca)[[6]][2,j],"variance explained"),
       ylab=paste0("PC",j+1,":",summary(TPMpca)[[6]][2,j+1],"variance explained"),
       main=paste0("logTPM PC",j,"vs", j+1," Hybrid vs Inbred ~ Lib batch"),cex=1.2,
       col=sizecol,
       pch=ifelse(meta$IsHybrid=="Hybrid",17,19))
  mtext(side=3,line=-1,"Size color, HI pch")
  
  
  #Treatment color, Replicate shade, PlantBatch pch
  plot(x=summary(TPMpca)[[2]][,j],y=summary(TPMpca)[[2]][,(j+1)],
       xlab=paste0("PC",j,":",summary(TPMpca)[[6]][2,j],"variance explained"),
       ylab=paste0("PC",j+1,":",summary(TPMpca)[[6]][2,j+1],"variance explained"),
       main=paste0("logTPM PC",j,"vs", j+1," BTH vs. Mock ~ Plant batch"),cex=1.2,
       col=ifelse(meta$Treatment=="Mock",
                  ifelse(meta$Replicate=="R1",alpha("springgreen4",0.9),ifelse(meta$Replicate=="R2",alpha("springgreen4",0.7),alpha("springgreen4",0.5))),
                  ifelse(meta$Replicate=="R1",alpha("darkorange",0.9),ifelse(meta$Replicate=="R2",alpha("darkorange",0.7),alpha("darkorange",0.5)))),
       pch=ifelse(meta$PlantBatch=="A",17,19))
  mtext(side=3,line=-1,"Treatment color, Replicate shade, PlantBatch pch")
  
  #HI color, Replicate shade, treatment pch
  plot(x=summary(TPMpca)[[2]][,j],y=summary(TPMpca)[[2]][,(j+1)],
       xlab=paste0("PC",j,":",summary(TPMpca)[[6]][2,j],"variance explained"),
       ylab=paste0("PC",j+1,":",summary(TPMpca)[[6]][2,j+1],"variance explained"),
       main=paste0("logTPM PC",j,"vs", j+1," BTH vs. Mock ~ Hybrid Inbred"),cex=1.2,
       col=ifelse(meta$IsHybrid=="Hybrid",
                  ifelse(meta$Replicate=="R1",alpha("turquoise4",0.9),ifelse(meta$Replicate=="R2",alpha("turquoise4",0.7),alpha("turquoise4",0.5))),
                  ifelse(meta$Replicate=="R1",alpha("purple4",0.9),ifelse(meta$Replicate=="R2",alpha("purple4",0.7),alpha("purple4",0.5)))),
       pch=ifelse(meta$Treatment=="Mock",17,19))
  mtext(side=3,line=-1,"HI color, Replicate shade, Treatment pch")
  
}
dev.off()



#make PCA of these additive genes
add_scaledlogTPM<-scalelogTPM[rownames(scalelogTPM)%in%AdditiveGenes,]
TPMpca<-prcomp(add_scaledlogTPM, center=T,scale.=T)
pdf(paste0(outbasename, "_MPVlmmadditive_pca.pdf"),height=30,width=15)
par(mfrow=c(6,3))
#%byPC
plot(x=1:ncol(summary(TPMpca)[[6]]),y=summary(TPMpca)[[6]][2,],
     main="Individual PC Variance Explained",
     xlab="PCs",ylab="Proportion Variance Explained",pch=19,xlim=c(0,450),ylim=c(0,1),cex=2)
plot(x=1:ncol(summary(TPMpca)[[6]]),y=summary(TPMpca)[[6]][3,],
     main="cumulative percentage by PC",
     xlab="PCs",ylab=,pch=19,xlim=c(0,450),ylim=c(0,1),cex=2)
hist(meta$percMapped,breaks=20)

for (j in 1:5){
  #size color, HI pch
  plot(x=summary(TPMpca)[[2]][,j],y=summary(TPMpca)[[2]][,(j+1)],
       xlab=paste0("PC",j,":",summary(TPMpca)[[6]][2,j],"variance explained"),
       ylab=paste0("PC",j+1,":",summary(TPMpca)[[6]][2,j+1],"variance explained"),
       main=paste0("logTPM PC",j,"vs", j+1," Hybrid vs Inbred ~ Lib batch"),cex=1.2,
       col=sizecol,
       pch=ifelse(meta$IsHybrid=="Hybrid",17,19))
  mtext(side=3,line=-1,"Size color, HI pch")
  
  #Treatment color, Replicate shade, PlantBatch pch
  plot(x=summary(TPMpca)[[2]][,j],y=summary(TPMpca)[[2]][,(j+1)],
       xlab=paste0("PC",j,":",summary(TPMpca)[[6]][2,j],"variance explained"),
       ylab=paste0("PC",j+1,":",summary(TPMpca)[[6]][2,j+1],"variance explained"),
       main=paste0("logTPM PC",j,"vs", j+1," BTH vs. Mock ~ Plant batch"),cex=1.2,
       col=ifelse(meta$Treatment=="Mock",
                  ifelse(meta$Replicate=="R1",alpha("springgreen4",0.9),ifelse(meta$Replicate=="R2",alpha("springgreen4",0.7),alpha("springgreen4",0.5))),
                  ifelse(meta$Replicate=="R1",alpha("darkorange",0.9),ifelse(meta$Replicate=="R2",alpha("darkorange",0.7),alpha("darkorange",0.5)))),
       pch=ifelse(meta$PlantBatch=="A",17,19))
  mtext(side=3,line=-1,"Treatment color, Replicate shade, PlantBatch pch")
  
  #HI color, Replicate shade, treatment pch
  plot(x=summary(TPMpca)[[2]][,j],y=summary(TPMpca)[[2]][,(j+1)],
       xlab=paste0("PC",j,":",summary(TPMpca)[[6]][2,j],"variance explained"),
       ylab=paste0("PC",j+1,":",summary(TPMpca)[[6]][2,j+1],"variance explained"),
       main=paste0("logTPM PC",j,"vs", j+1," BTH vs. Mock ~ Hybrid Inbred"),cex=1.2,
       col=ifelse(meta$IsHybrid=="Hybrid",
                  ifelse(meta$Replicate=="R1",alpha("turquoise4",0.9),ifelse(meta$Replicate=="R2",alpha("turquoise4",0.7),alpha("turquoise4",0.5))),
                  ifelse(meta$Replicate=="R1",alpha("purple4",0.9),ifelse(meta$Replicate=="R2",alpha("purple4",0.7),alpha("purple4",0.5)))),
       pch=ifelse(meta$Treatment=="Mock",17,19))
  mtext(side=3,line=-1,"HI color, Replicate shade, Treatment pch")
  
}
dev.off()
#additive PCA is more separated between hybrid and inbred. It is hard not to have size signature in the PCA, though, 
#since SHB2.2 plants are larger, and the transriptome anyways segregate by SHB2.1 and SHB2.2.

#remove the additive genes from clustering analysis. Use the slope 0.4, sigma 0.65 category
addrm_scaledlogTPM<-scalelogTPM[!rownames(scalelogTPM)%in%AdditiveGenes,] #13061 genes
write.table(addrm_scaledlogTPM,file=paste0(outbasename,"_rmadditivegenes_scaledlogTPM.txt"),quote=F,sep="\t")
