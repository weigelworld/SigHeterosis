#########
#README
#The script sets up for rare allele analysis (whether rare allele burden at certain genes correlates with phenotype).
#index rare alleles from 1kb-upstream of all genes, exclude singletons
#go to each accession, count its rare alleles upstream
#index the accession-specific rare allele counts to upstream of the gene lists of interest

#The daughter scripts of this script (10A/B) associate rare allele counts upstream of gene lists with expression level /phenotype in inbreeding parents,
#and address how the rare allele status in the accessions affect F1 expression/phenotype
########
#Config and data input
library(vioplot)
setwd("/ebio/abt6_projects7/SHB2/data/1_analysis/10_Subgenome_Genetdist")
freq.bkgd<-read.table("sub1001_BTHBkgdlist_1kbpromoterONLY.filt.frq",header=F,skip=1) #441340 SNPs
BED.bkgd<-read.table("run197_BTHdiffgene_rowfiltered_Bkgdlist_1kbpromoterONLY.BED",header=F)
names(freq.bkgd)<-c("Chr","Pos","No.allele","Cnt.allele","Freq1","Freq2")
#input the genelists of interest
AllBTHneg<-read.table("/ebio/abt6_projects7/SHB2/data/1_analysis/7_MPHsplineGO/run197_MPHsplinepredict_rmflatK61_AllBTHnegative.txt",header=F,sep="\t")[,1]
pospos<-read.table("/ebio/abt6_projects7/SHB2/data/1_analysis/7_MPHsplineGO/run197_MPHsplinepredict_rmflatK61_5_positive-positive.txt",header=F,sep="\t")[,1]
additive<-read.table("/ebio/abt6_projects7/SHB2/data/1_analysis/6_BTHdiffgenes/run125_run197_loosethreshold_additiveintersect.txt",header=F,sep="\t")[,1] #common additive genes
#where accession IDs are stored:
seqmeta<-read.table("/ebio/abt6_projects7/SHB2/data/1_analysis/6_BTHdiffgenes/Run197_rmCMRTPHBS_rm2MSeqtkDupOL_corr.meta.txt",header=T,sep="\t",row.names=1)
#variants table
vcftab<-"sub1001_BTHBkgdlist_1kbpromoterONLY.filt.variants.table"

############
#calculate minor allele frequency (MAF)
freq1<-sapply(freq.bkgd$Freq1,function(x){as.numeric(substring(x,first=3))})
freq2<-sapply(freq.bkgd$Freq2,function(x){as.numeric(substring(x,first=3))})
freq.bkgd[,"MAF"]<-ifelse(freq1<=freq2,freq1,freq2)
# summary(freq.bkgd$MAF)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001000 0.002760 0.008108 0.046146 0.036684 0.500000
pdf("sub1001_BTHBkgdlist_1kbpromoterONLY_filt_MAFhistogram.pdf",height=14,width=14)
par(mfrow=c(2,2))
hist(freq.bkgd$MAF,breaks=1000,main="Histogram of all SNPs",xlab="MAF")
hist(freq.bkgd$MAF[freq.bkgd$MAF<=0.05],breaks=500,main="Histogram of MAF<=0.05",xlab="MAF")
hist(freq.bkgd$MAF[freq.bkgd$MAF<=0.01],breaks=100,main="Histogram of MAF<=0.01",xlab="MAF")
hist(freq.bkgd$MAF[freq.bkgd$MAF<=0.002],breaks=50,main="Histogram of MAF<=0.002",xlab="MAF")
dev.off()

#from the histogram it seems that the singletons are a group of SNPs whose MAF<0.0015
rareAllele<-subset(freq.bkgd,MAF>0.0015 & MAF<=0.05)
#assign rare allele to gene
whichgenelist<-list()
for(chr in 1:5){
  print(chr)
  subrare<-subset(rareAllele,Chr==chr)
  subBED<-subset(BED.bkgd,V1==chr)
  whichgene<-t(sapply(subrare$Pos,function(x){return(c(chr,x,subBED[which(subBED$V2<=x & subBED$V3>=x),]$V4))})) 
  # a small number of SNPs (<5%) hits multiple genes, in which 31 SNPs are assigned to 3 genes (they tend to be consecutive SNPs) 
  whichgenelist[[chr]]<-whichgene
}
#for those SNPs that are assigned to multiple genes, give +1 count for all genes
genevote<-unlist(lapply(whichgenelist,function(x){unlist(lapply(x,function(y){return(y[3:5])}))}))
genevote<-table(genevote[!is.na(genevote)])
quantile(genevote,probs=c(0.05,0.1,0.25,0.5,0.75,0.95))
# 5%  10%  25%  50%  75%  95% 
# 5.0  9.0 16.0 24.0 33.0 56.4 
whichgenelist<-unlist(whichgenelist,recursive = F) #some SNPs will have more than 1 genes assigned

### sub index the SNPs in gene lists of interest
genevote.allBTHneg<-genevote[names(genevote)%in%AllBTHneg]
genevote.pospos<-genevote[names(genevote)%in%pospos]
genevote.additive<-genevote[names(genevote)%in%additive]
#plot distribution of rare alleles in each gene list. 
pdf("run197_subgenelists_pergene_1kbupstreamONLY_MAF0.05_1135rareallelecount.pdf")
par(mar=c(10,4,3,1))
vioplot(genevote,genevote.allBTHneg,genevote.pospos,genevote.additive,
        main="per-gene upstream rare allele count", ylab="Count",
        names=c("All Genes","BTH-Negative","Positive","Additive"),las=2)
dev.off()
write.table(genevote,"run197_BTHdiffgene_1kbupstreamONLY_MAF0.05_1135rareallelecount.txt",quote=F,sep="\t",row.names=F)
#13393 genes had rare allele vote

#############
#count per-accession rare allele
  #get accession IDs of the parents of interest
gt<-unique(seqmeta$Genotype)
accgt<-gt[!(grepl("SH",gt)|grepl("rd",gt))]
rareChrPos<-paste0(rareAllele[,1],"_",rareAllele[,2])
con<-file(vcftab)
open(con)
vcfheader<-readLines(con,n=1)
vcfheader<-unlist(strsplit(vcfheader,"\t"))
vcfacc<-sapply(vcfheader[6:length(vcfheader)],function(x){paste0("p",substr(x,1,(nchar(x)-3)))})
SHB2accidx<-as.numeric(sapply(accgt,function(x){which(vcfacc==x)}))+5 #This gives the column index in the vcf file to retrieve the correct acc genotype

#read in vcf table line-by-line, id if it is a rare allele, id the genotypes of the parents, and whether the parents are carrying the rare allele
SHB2rarealleleGT<-matrix(nrow=nrow(freq.bkgd),ncol=71)
for(i in 2:(nrow(freq.bkgd)+1)){
  if(i%%1000==2){print(i)}
  snp<-readLines(con,n=1)#read in the next SNP
  snp<-unlist(strsplit(snp,"\t"))
  if(!paste0(snp[1],"_",snp[2]) %in% rareChrPos){  #if it is not a rare SNP, skip
    next
  }else{ #get the genotypes of the accessions, together with the REF-ALT info
    rarefreq<-rareAllele[rareChrPos==paste0(snp[1],"_",snp[2]),5:7]
    snpacc<-snp[SHB2accidx]
    snpinfo<-snp[c(1,2,4,5)]
    if(any(snpacc=="./.")&length(unique(snpacc))<=2){#if it's not polymorphic within the SHB2 parents
      next
    }else{
      if(sum(grepl("./.",snpacc))>=16){ #if the SNP has >25% missing info in the SHB2 accessions
        next
      }else{ #if SNP is a rare, polymorphic within SHB2, and has no more than 25% missing info, record it
        returninfo<-as.character(c(snpinfo,rarefreq,snpacc))
        SHB2rarealleleGT[(i-1),]<-returninfo
      }
    }
  }
}
close(con)

SHB2rarealleleGT.NA<-na.omit(SHB2rarealleleGT) #152418 rare SNPs remained that are polymorphic in SHB2 accessions
SHB2rarealleleGT.NA<-as.data.frame(SHB2rarealleleGT.NA)
names(SHB2rarealleleGT.NA)<-c("Chr","Pos","Ref","Alt","Freq1","Freq2","MAF",accgt)

#find which genes these SNPs belong to
whichgene.ChrPos<-unlist(lapply(whichgenelist,function(x){paste0(x[1],"_",x[2])}))
SHB2rareallele.ChrPos<-paste0(SHB2rarealleleGT.NA[,1],"_",SHB2rarealleleGT.NA[,2])
whichgeneidx<-sapply(SHB2rareallele.ChrPos,
                     function(x){which(whichgene.ChrPos==x)})
  #check how many genes are involved
len<-sapply(whichgeneidx,function(x){return(length(whichgenelist[[x]]))})
  # table(len)
  # len
  # 3      4      5 
  # 146444   5961     13

#count how many rare allele each gene in each accession has. 
  #First expand the rare SNP table, allowing duplicated SNPs to make sure each row corresponds to a single gene
gene1<-sapply(whichgeneidx,function(x){return(whichgenelist[[x]][3])}) #first all SNPs get associated with their first hit
gene2<-sapply(whichgeneidx[which(len==4)],function(x){return(whichgenelist[[x]][4])})#then SNPs with 2 hits get associated with their second hit
dupSNPgene2<-SHB2rarealleleGT.NA[sapply(names(gene2),function(x){which(SHB2rareallele.ChrPos==x)}),]
gene3<-sapply(whichgeneidx[which(len==5)],function(x){return(whichgenelist[[x]][5])})#finally SNPs with 3 hits get associated with their last hit
dupSNPgene3<-SHB2rarealleleGT.NA[sapply(names(gene3),function(x){which(SHB2rareallele.ChrPos==x)}),]
SHB2rarealleleGT.NA[,"Gene"]<-gene1
dupSNPgene2[,"Gene"]<-gene2
dupSNPgene3[,"Gene"]<-gene3
SHB2rarealleleGT.NA<-rbind(SHB2rarealleleGT.NA,dupSNPgene2,dupSNPgene3)
SHB2rarealleleGT.NA[,"Freq1"]<-sapply(SHB2rarealleleGT.NA$Freq1,substring,first=3)
SHB2rarealleleGT.NA[,"Freq2"]<-sapply(SHB2rarealleleGT.NA$Freq2,substring,first=3)
write.table(SHB2rarealleleGT.NA,"run197_BTHdiffgene_1kbupstreamONLY_MAF0.05_SHB2accGT.txt",quote=F,sep="\t",row.names=F)

#assign common(C)/rare(R)/missing(M) to each SNP in each accession
assignCRM<-function(x){
  ref<-x[3]
  alt<-x[4]
  fref<-as.numeric(x[5])
  falt<-as.numeric(x[6])
  Maf<-as.numeric(x[7])
  rare<-as.character(ifelse(fref==Maf,ref,alt))
  crm<-rep("C",64)
  crm[grepl(rare,x[8:71])] <-"R"
  crm[grepl("./.",x[8:71])]<-"M"
  return(as.character(c(x[c(1:7,72)],crm)))
}
CRM<-as.data.frame(t(apply(SHB2rarealleleGT.NA,1,assignCRM)))
names(CRM)<-c("Chr","Pos","Ref","Alt","FreqRef","FreqAlt","MAF","Gene",accgt)
write.table(CRM,"run197_BTHdiffgene_1kbupstreamONLY_MAF0.05_SHB2accCommon.Rare.txt",quote=F,sep="\t",row.names=F)

#count nubmer of rare alleles per gene per accession
Acc.rareallele<-apply(CRM[,9:ncol(CRM)],2,function(x){tapply(x,CRM$Gene,function(y){sum(y=="M")})})
write.table(Acc.rareallele,"run197_BTHdiffgene_1kbupstreamONLY_MAF0.005_SHB2acc_pergeneRareallelecount.txt",quote=F,sep="\t")

#index the rare allele counts to gene lists of interest
AllBTHneg.rareallele<-Acc.rareallele[rownames(Acc.rareallele)%in%AllBTHneg,]
pospos.rareallele<-Acc.rareallele[rownames(Acc.rareallele)%in%pospos,]
additive.rareallele<-Acc.rareallele[rownames(Acc.rareallele)%in%additive,]

#plot distribution of rare alleles in each gene list. 
maxRA<-apply(Acc.rareallele,1,max)
meanRA<-apply(Acc.rareallele,1,mean)
pdf("run197_subgenelists_pergene_1kbupstreamONLY_MAF0.05_SHB2rareallelecount.pdf",width=15,height = 8)
par(mfrow=c(1,2),mar=c(10,4,3,1))
vioplot(maxRA,
        maxRA[rownames(Acc.rareallele)%in%AllBTHneg],
        maxRA[rownames(Acc.rareallele)%in%pospos],
        maxRA[rownames(Acc.rareallele)%in%additive],
        main="per-gene max upstream rare allele count", ylab="Count",
        names=c("All Genes","BTH-Negative","Positive","Additive"),las=2)
vioplot(meanRA,
        meanRA[rownames(Acc.rareallele)%in%AllBTHneg],
        meanRA[rownames(Acc.rareallele)%in%pospos],
        meanRA[rownames(Acc.rareallele)%in%additive],
        main="per-gene mean upstream rare allele count", ylab="Count",
        names=c("All Genes","BTH-Negative","Positive","Additive"),las=2)
dev.off()

