#to be ran in the most recent common ancester of all rsem results files

args<-commandArgs(trailingOnly=T)
Routfilename=args[1]
filelistname=args[2]
Rheader=args[3]
SampIDfield=as.numeric(args[4])

print(paste0("Output TPM file: ", Routfilename))
print(paste0("list of individual RSEM results: ",filelistname))
print(paste0("SampleID found in field:",SampIDfield))

###check the Rscript for indexing of sample ID (the SXXX number) from input file name via strsplit(). It varies from run to run
Rout<-read.table(file=Routfilename,header=T,sep="\t")
filelist<-readLines(filelistname)
print("Head of output TPM file:")
print(head(Rout))
print("Examples of individual results file:")
print(filelist[1:5])

for (i in 1:length(filelist)){
	infilename<-filelist[i]
	print(paste(i,infilename,sep=" "))
	res<-read.table(file=infilename,header=T, sep="\t",row.names=1)
	#print(head(res))
	RPK.eff<-res$expected_count*1000/res$effective_length
	sumRPKeff<-sum(RPK.eff[!is.na(RPK.eff)& RPK.eff!="Inf"])
	TPM.eff<-RPK.eff*1000000/sumRPKeff

	sampleID<-unlist(strsplit(infilename, "_"))[SampIDfield] #the SXXX number
	#print(sampleID)
	#filterinfo<-unlist(strsplit(infilename, "_"))[4]
	#sampleapd<-ifelse(grepl("genes.results",filterinfo), "All",
	#	     ifelse(grepl("rmChrCM.results",filterinfo),"CM",
	#		    ifelse(grepl("rDNA.results",filterinfo),"CMR",
	#			   ifelse(grepl("pseudogene.results",filterinfo),"CMRTP","ERROR"))))
	samplename<-paste(sampleID, Rheader, sep="_")
	print(samplename)

	Rout[,samplename]<-TPM.eff
	if(i%%10==1){
		print("Updated output TPM file:")
		print(dim(Rout))
	}
}
write.table(Rout, file=Routfilename, quote=F,sep="\t",row.names=F)
