#the script calculates log2TPM, and re-intercalate zero values. Note that sample rankings for each gene would change.

args<-commandArgs(trailingOnly=T)
infilename<-args[1]
outfilename<-args[2]

print("reading input")
TPM<-read.table(infilename, header=T,sep="\t",row.names=1) #must be TPM files with geneIDs as rownames
print(dim(TPM))
print(head(TPM)[,1:6])

print("performing log transformation")
#logTPM<-t(apply(TPM,1,function(x){unlist(sapply(x, function(y){ifelse(y==0,0,log2(y))}))}))
logTPM<-apply(TPM,c(1,2),function(x){log2(x+1)})
print(dim(logTPM))
print(head(logTPM)[,1:6])

print("creating output data table")
write.table(logTPM, outfilename, quote=F,sep="\t")
