args=commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  stop("input directory not supplied", call.=FALSE)
}

rsemdatadir=args[1]
rsem_ttlsamp=args[2]
print(c(rsemdatadir,rsem_ttlsamp),sep="\n")
datadirname<-unlist(strsplit(rsemdatadir,split="/"))
datadirname<-datadirname[length(datadirname)]
print(datadirname)

setwd(rsemdatadir)
print(getwd())

samplename<-readLines("samplename.txt")
firstline<-samplename[1]
shortname<-unlist(strsplit(firstline,"/"))
shortname<-shortname[length(shortname)-1]
print(shortname)
template<-read.table(file=firstline,header=T,sep="\t")
templatetpm=template[,c(1,3:4,6)]
print(head(templatetpm))
print("################")
names(templatetpm)<-c("gene","length","effective_length", shortname)
for(i in 2:rsem_ttlsamp){
  currsample<-samplename[i]
  currshortname<-unlist(strsplit(currsample,"/"))
  currshortname<-currshortname[length(currshortname)-1]
  currdata<-read.table(file=currsample,header=T,sep="\t")
  templatetpm[,currshortname]<-currdata$TPM
  print(currshortname)
}
outname=paste("rsem",datadirname,"TPM.txt", sep="_")
write.table(templatetpm, file=outname, quote=F,sep="\t", row.names=F)

