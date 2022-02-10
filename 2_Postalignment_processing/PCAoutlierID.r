#take a logTPM file, filter genes with a dropout threshold, performPCA, remove outliers and dup'ed samples
args<-commandArgs(trailingOnly=T)
wd<-args[1]
logTPMin<-args[2]
metain<-args[3]
cutoff<-as.numeric(args[4])
outbasename<-args[5]
TPMin<-args[6]

setwd(wd)
library(scales)
library(RColorBrewer)

print("inputting data")
meta<-read.table(metain, header=T,sep="\t", row.names=1)
logTPM<-read.table(logTPMin, header=T,sep="\t",row.names=1)
TPM<-read.table(TPMin, header=T,sep="\t",row.names=1)

meta$IsHybrid<-factor(meta$IsHybrid,levels=c("Inbred","Hybrid"))
meta$Treatment<-factor(meta$Treatment,levels=c("Mock","BTH"))

print("setting plot parameters")
mappedcol=alpha(brewer.pal(4,"BuPu")[2:4],0.7)
batchlib<-meta$LibBatch
treatcol<-c("springgreen4","darkorange") #green for mock, orange for BTH

#batchRS<-paste0(batchlib,meta$TissueType) #coloring criteria for the batch
#batchRScol=c(alpha(c("red","goldenrod3","forestgreen","turquoise4","purple","grey30"),0.9),alpha(c("red","goldenrod3","forestgreen","turquoise4","purple","grey30"),0.6))
#batchlab<-ifelse(batchRS=="aS",batchRScol[1],
#            ifelse(batchRS=="bS",batchRScol[2],
#              ifelse(batchRS=="cS",batchRScol[3],
#                ifelse(batchRS=="dS",batchRScol[4],
#                  ifelse(batchRS=="eS",batchRScol[5],
#                    ifelse(batchRS=="fS",batchRScol[6],
#		      ifelse(batchRS=="aR",batchRScol[7],
#                        ifelse(batchRS=="bR",batchRScol[8],
#                          ifelse(batchRS=="cR",batchRScol[9],
#                            ifelse(batchRS=="dR",batchRScol[10],
#			      ifelse(batchRS=="eR",batchRScol[11],
#				 batchRScol[12])))))))))))


#filter dropout genes by threshold
print("filtering dropouts")
logTPM<-logTPM[apply(logTPM,1,function(x){sum(x==0)})<=cutoff,]
logTPM<-logTPM[apply(logTPM,1,function(x){!any(x==Inf)}),]

#PCA
print("Generating PCA")
TPMpca<-prcomp(logTPM, center=T,scale.=T)

#outpdf<-paste0(outbasename,".pdf")
#outpdf<-"run189_rmCMRTPHBS_OL_PCA.pdf"
outpdf<-"run197_rmCMRTPHBS_rm2MSeqtkDupOL_PCA.pdf"
pdf(file=outpdf,height=35,width=20)
par(mfrow=c(7,4))

#%byPC
plot(x=1:ncol(summary(TPMpca)[[6]]),y=summary(TPMpca)[[6]][2,],
          main="Individual PC Variance Explained",
	       xlab="PCs",ylab="Proportion Variance Explained",pch=19,xlim=c(0,450),ylim=c(0,1),cex=2)
plot(x=1:ncol(summary(TPMpca)[[6]]),y=summary(TPMpca)[[6]][3,],
          main="cumulative percentage by PC",
	       xlab="PCs",ylab=,pch=19,xlim=c(0,450),ylim=c(0,1),cex=2)
hist(meta$percMapped,breaks=20)
plot(1:10,1:10,axes=F,main="",xlab="",ylab="",type="n")
for (j in 1:5){
  #HI color, LibBatch pch
  plot(x=summary(TPMpca)[[2]][,j],y=summary(TPMpca)[[2]][,(j+1)],
          xlab=paste0("PC",j,":",summary(TPMpca)[[6]][2,j],"variance explained"),
          ylab=paste0("PC",j+1,":",summary(TPMpca)[[6]][2,j+1],"variance explained"),
          main=paste0("logTPM PC",j,"vs", j+1," Hybrid vs Inbred ~ Lib batch"),cex=1.2,
	  col=ifelse(meta$IsHybrid=="Hybrid",alpha("turquoise4",0.5),alpha("purple4",0.5)),
          pch=ifelse(meta$LibBatch=="a",17,19))
  mtext(side=3,line=-1,"HI color, LibBatch pch")

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

 #Treatment color, Replicate shade, HI pch
 plot(x=summary(TPMpca)[[2]][,j],y=summary(TPMpca)[[2]][,(j+1)],
	xlab=paste0("PC",j,":",summary(TPMpca)[[6]][2,j],"variance explained"),
	ylab=paste0("PC",j+1,":",summary(TPMpca)[[6]][2,j+1],"variance explained"),
	main=paste0("logTPM PC",j,"vs", j+1," BTH vs. Mock ~ Hybrid Inbred"),cex=1.2,
	col=ifelse(meta$Treatment=="Mock",
		ifelse(meta$Replicate=="R1",alpha("springgreen4",0.9),ifelse(meta$Replicate=="R2",alpha("springgreen4",0.7),alpha("springgreen4",0.5))),
		ifelse(meta$Replicate=="R1",alpha("darkorange",0.9),ifelse(meta$Replicate=="R2",alpha("darkorange",0.7),alpha("darkorange",0.5)))),
	pch=ifelse(meta$IsHybrid=="Hybrid",17,19))
 mtext(side=3,line=-1,"Treatment color, Replicate shade, HI pch")

 #ttlreads color, LibBatch pch
  plot(x=summary(TPMpca)[[2]][,j],y=summary(TPMpca)[[2]][,(j+1)],
       xlab=paste0("PC",j,":",summary(TPMpca)[[6]][2,j],"variance explained"),
       ylab=paste0("PC",j+1,":",summary(TPMpca)[[6]][2,j+1],"variance explained"),
       main=paste0("logTPM PC",j,"vs", j+1," total reads ~ Lib batch"),cex=1.2,
       col=ifelse(meta$ttlReads<=10^6,mappedcol[1],ifelse(meta$ttlReads>10^6 & meta$ttlReads<=2*10^6 ,mappedcol[2],mappedcol[3])),
       pch=ifelse(meta$LibBatch=="a",17,19))
  mtext(side=3,line=-1,"total reads color, LibBatch pch")
  
}
 
#legend plots
#plot(x=rep(1,4),y=c(4,3,2,1),pch=c(4,4,17,19),col=c("turquoise4","purple4","grey30","grey30"),xlab="",ylab="",axes=F,xlim=c(0.5,2),ylim=c(0,5),cex=4,bty="n")
#text(1.2,2,"Hybrid",cex=3,pos=4)
#text(1.2,1,"Inbred",cex=3,pos=4)

#bar<-barplot(tapply(rep(1,length(batchRS)),batchRS,sum)[c(2,4,6,8,10,12,1,3,5,7,9,11)],pch=15,
#	    col=batchRScol,
#	    names=names(tapply(rep(1,length(batchRS)),batchRS,sum)[c(2,4,6,8,10,12,1,3,5,7,9,11)]),
#	    xlab="",ylab="No.Samples",main="Shoot Root Sample Count per Batch",ylim=c(0,60),cex.lab=1.2)
#text(x = bar, y = tapply(rep(1,length(batchRS)),batchRS,sum)[c(2,4,6,8,10,12,1,3,5,7,9,11)],
#          label = tapply(rep(1,length(batchRS)),batchRS,sum)[c(2,4,6,8,10,12,1,3,5,7,9,11)], pos = 3, cex = 2)

#plot(x=seq(0.5,1,0.1),y=seq(0.5,1,0.1),type="n",xlab="Mapping Stats Intervals",ylab="",main="Color Code for % Reads Mapped",
#          xlim=c(0.4,1),ylim=c(0.4,1))
#rect(0.4,0.4,0.7,0.7,col=mappedcol[1])
#rect(0.7,0.7,0.85,0.85,col=mappedcol[2])
#rect(0.85,0.85,1,1,col=mappedcol[3])
#text(0.4,0.55,"<=70% mapped",pos=4,cex=3)
#text(0.7,0.8,"70-85% mapped",pos=4,cex=3)
#text(0.85,0.95,">85% mapped",pos=4,cex=3)

dev.off()

#######################
#plot to filter outliers on TPMpca PC1
print("Performing outlier analysis")
outlierpdf<-paste0(outbasename,"outlier_boxplot.pdf")
pdf(file=outlierpdf,height=5,width=10)
par(mfrow=c(1,2))

pc1box<-boxplot(TPMpca[[2]][,1],main="All samples on AllsamplePC1",ylim=c(min(TPMpca[[2]][,1])-0.05*diff(range(TPMpca[[2]][,1])),max(TPMpca[[2]][,1])+0.05*diff(range(TPMpca[[2]][,1]))))
pc1out<-pc1box$out

# pc2box<-boxplot(TPMpca[[2]][,2],main="All samples on AllsamplePC2",ylim=c(min(TPMpca[[2]][,2])-0.05*diff(range(TPMpca[[2]][,2])),max(TPMpca[[2]][,2])+0.05*diff(range(TPMpca[[2]][,2]))))
# pc2Sbox<-boxplot(TPMpca[[2]][,2][which(meta$TissueType=="S")],main="Shoot samples on AllsamplePC2",ylim=c(min(TPMpca[[2]][,2])-0.05*diff(range(TPMpca[[2]][,2])),max(TPMpca[[2]][,2])+0.05*diff(range(TPMpca[[2]][,2]))))
# pc2Rbox<-boxplot(TPMpca[[2]][,2][which(meta$TissueType=="R")],main="Root samples on AllsamplePC2",ylim=c(min(TPMpca[[2]][,2])-0.05*diff(range(TPMpca[[2]][,2])),max(TPMpca[[2]][,2])+0.05*diff(range(TPMpca[[2]][,2]))))
# pc2out<-pc2box$out
# pc2Sout<-pc2Sbox$out
# pc2Rout<-pc2Rbox$out

samp<-colnames(logTPM)

for(j in 1:2){
  plot(x=summary(TPMpca)[[2]][,j],y=summary(TPMpca)[[2]][,(j+1)],
      xlab=paste0("PC",j,":",summary(TPMpca)[[6]][2,j],"variance explained"),
      ylab=paste0("PC",j+1,":",summary(TPMpca)[[6]][2,j+1],"variance explained"),
      main=paste0("logTPM PC",j,"vs", j+1," Hybrid vs Inbred"),cex=1.2,
       col=ifelse(samp%in%names(pc1out),"black",
                  alpha("grey60",0.6)))
      #            ifelse(samp%in%names(pc1Sout),"turquoise4",
      #             ifelse(samp%in%names(pc2Sout),alpha("turquoise2",0.7),
      #               ifelse(samp%in%names(pc2Rout),alpha("purple4",0.7),alpha("grey60",0.6))))),
      #pch=ifelse(meta$TissueType=="S",17,19))
}
#legendplot
# plot(rep(1,6),6:1,xlim=c(0.5,2.5),ylim=c(0,7),cex=3,pch=c(17,19,17,17,17,19),col=c("grey60","grey60","black","turquoise4","turquoise2","purple4"),bty="n",axes=F,xlab="",ylab="",main="")
# text(1.2,6,pos=4,cex=2,"Shoot in range")
# text(1.2,5,pos=4,cex=2,"Root in range")
# text(1.2,4,pos=4,cex=2,"PC1 outlier")
# text(1.2,3,pos=4,cex=2,"PC1 Shoot outlier")
# text(1.2,2,pos=4,cex=2,"PC2 Shoot outlier")
# text(1.2,1,pos=4,cex=2,"PC2 Root outlier")
dev.off()

outlierfilename<-paste0(outbasename,"outlierSampleID.txt")
cat(paste("PC1outliers","PC1loading",sep="\t"), file=outlierfilename,sep="\n")
cat(paste(names(pc1out),pc1out,sep="\t"), file=outlierfilename,sep="\n",append=T)
# cat("PC1Shoot_out", file=outlierfilename,sep="\n",append=T)
# cat(names(pc1Sout), file=outlierfilename,sep="\n",append=T)
# cat("PC1Root_out", file=outlierfilename,sep="\n",append=T)
# cat(names(pc1Rout), file=outlierfilename,sep="\n",append=T)
# cat("PC2out", file=outlierfilename,sep="\n",append=T)
# cat(names(pc2out), file=outlierfilename,sep="\n",append=T)
# cat("PC2Shoot_out", file=outlierfilename,sep="\n",append=T)
# cat(names(pc2Sout), file=outlierfilename,sep="\n",append=T)
# cat("PC2Root_out", file=outlierfilename,sep="\n",append=T)
# cat(names(pc2Rout), file=outlierfilename,sep="\n",append=T)
#cat("NTC",file=outlierfilename,sep="\n",append=T)
#cat(c("S383_CMRTPHBS","S384_CMRTPHBS"),file=outlierfilename,sep="\n",append=T)
#filter samples whose PC1 loading<0.035 (aka far outliers)
# TPM_OL<-TPM[,!names(TPM)%in%c(names(pc1out)[pc1out<=0.035],"S383_CMRTPHBS","S384_CMRTPHBS")]
# write.table(TPM_OL,"run189_rmCMRTPHBS_OL.TPM.txt",quote=F,sep="\t")
# filtersampID<-sapply(c(names(pc1out)[pc1out<=0.035],"S383_CMRTPHBS","S384_CMRTPHBS"),function(x){unlist(strsplit(x,"_"))[1]})
# meta_OL<-meta[!rownames(meta)%in%filtersampID,]
# write.table(meta_OL,"Run189_rmCMRTPHBS_OL.meta.txt",quote=F,sep="\t")
# #in order to go back an generate the filtered PCA, modify the logTPM file as well, without writing table
# logTPM<-logTPM[,!names(logTPM)%in%c(names(pc1out)[pc1out<=0.035],"S383_CMRTPHBS","S384_CMRTPHBS")]

#for run197, NTC and low-coverage samples already filtered.
#filter samples whose PC1 loading <0.04 (aka far outliers)
TPM_OL<-TPM[,!names(TPM)%in%names(pc1out)[pc1out<=0.04]]
write.table(TPM_OL,"run197_rmCMRTPHBS_rm2MSeqtkDupOL.TPM.txt",quote=F,sep="\t")
filtersampID<-sapply(c(names(pc1out)[pc1out<=0.04]),function(x){unlist(strsplit(x,"_"))[1]})
meta_OL<-meta[!rownames(meta)%in%filtersampID,]
write.table(meta_OL,"Run197_rmCMRTPHBS_rm2MSeqtkDupOL.meta.txt",quote=F,sep="\t")
logTPM<-logTPM[,!names(logTPM)%in%names(pc1out)[pc1out<=0.04]]
# generate PCA after filter