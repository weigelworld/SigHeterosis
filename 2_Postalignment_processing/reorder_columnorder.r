#Rscript /PATH/TO/SCRIPT <PATH/TO/WORKING/DIRECTORY> <input filename>
args=commandArgs(trailingOnly=T)
wd<-args[1]
infilename<-args[2]
outfilename<-args[3]

setwd(wd)
infile<-read.table(file=infilename, header=T,sep="\t")
samplename<-colnames(infile)[2:ncol(infile)]
samplenum<-sapply(
	sapply(samplename,function(x){unlist(strsplit(x,"_"))[1]}), 
	function(x){as.numeric(substr(x,2,nchar(x)))})
infileorder<-infile[,2:ncol(infile)]
infileorder<-infileorder[,order(samplenum)]
infileorder<-apply(infileorder,c(1,2),function(x){ifelse(is.na(x),0,x)})
outfile<-data.frame("GeneID"=infile$gene_id)
outfile<-cbind(outfile,infileorder)

write.table(outfile, file=outfilename,quote=F,sep="\t",row.names=F)
