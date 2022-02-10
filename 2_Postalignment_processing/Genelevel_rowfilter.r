#use this script to filter high drop-out, low-expressed, and in variable genes. 
#the script first check and remove all lines that dropped out on all genes
#for the remaining genes, the script appends dropout rate, expression mean val, and variance coef to the gene list
#the script over writes the logTPM value

args<-commandArgs(trailingOnly=T)
setwd(args[1])
TPM<-read.table(file=args[2], header=T, sep="\t")
TPMrs<-apply(TPM[,2:ncol(TPM)],1,function(x){sum(abs(x))!=0})
TPM<-TPM[TPMrs,] #filter completely dropped out lines
dropout<-apply(TPM[,2:ncol(TPM)],1, function(x){sum(x==0)})
expmean<-apply(TPM[,2:ncol(TPM)],1,function(x){mean(x[x!=0])})
varcoef<-apply(TPM[,2:ncol(TPM)], 1, function(x){sd(x[x!=0])/mean(x[x!=0])})
TPM[,"Dropout_count"]<-dropout
TPM[,"Non0Exp_mean"]<-expmean
TPM[,"VarCoef"]<-varcoef
write.table(TPM, file=args[2],quote=F,sep="\t", row.names=F)
