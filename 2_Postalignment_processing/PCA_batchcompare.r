#evaluate gene batch mean by PCA
args<-commandArgs(trailingOnly=T)
sourcepath<-args[1]
pdfwidth<-as.numeric(args[2])
#print(pdfwidth)
source(sourcepath)
print(length(datalist))
library(scales)

pdfname<-paste0(outbasename,"_batcheval.pdf")
pdf(pdfname,height=70,width=pdfwidth)
par(mfcol=c(7,length(datalist)))
batchScol=alpha(c("red","goldenrod3","forestgreen","turquoise4","purple","grey30"),0.7)
for(i in 1:length(datalist)){
 currlist<-datalist[[i]]
 pc<-prcomp(currlist,center=T,scale.=T)
 batchcol<-batchScol[as.numeric(meta$LibBatch)]

 #%byPC
 plot(x=1:ncol(summary(pc)[[6]]),y=summary(pc)[[6]][2,],
     main=paste0(names(datalist)[i],": Individual PC Variance Explained"),xlab="PCs",ylab="Proportion Variance Explained",
     pch=19,xlim=c(0,240),ylim=c(0,1),cex=2)

 #color by batch, shape by HI
         for(j in 1:5){		 
		 plot(x=summary(pc)[[2]][,j],y=summary(pc)[[2]][,(j+1)],
  	         xlab=paste0("PC",j,":",summary(pc)[[6]][2,j],"variance explained"),
		 ylab=paste0("PC",j+1,":",summary(pc)[[6]][2,j+1],"variance explained"),
		 main=paste0(names(datalist)[i],"PC",j,"vs", j+1),
		 cex=1.2,
		 pch=ifelse(meta$IsHybrid=="Hybrid",17,19),
		 col=batchcol)
	 }
	#boxplot of batch variance
	pergenemean<-as.data.frame(t(apply(currlist,1,function(x){tapply(x,meta$LibBatch, mean)})))
	boxplot(pergenemean$a,pergenemean$b,pergenemean$c,pergenemean$d,pergenemean$e,pergenemean$f,col=batchScol,names=letters[1:6], main=paste0(names(datalist)[i]," Batchmean logTPM"))
  	print(i)
}
dev.off()


