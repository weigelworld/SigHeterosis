#provide environmental variables for PCA plotting
setwd("/ebio/abt6_projects9/SigHeterosis_Batch1/data/1_analysis/7_rsem_tair10rmChrCM/rmCMRTPHBS_TPM")
meta<-read.table("run125_RSEM_rmCMRTPHBS_rmColOL3.5M_Shoot_metainfo.txt",header=T,sep="\t",row.names=1)
outbasename<-"run125_rmCMRTPHBSColOL3.5M_ShootGIS_PCA"
GIS<-read.table("run125_rmCMRTPHBSColOL3.5M_ShootGIS0.25_invariantGeneID.txt",header=T,sep="\t")[,1]

logTPM<-read.table("run125_rmCMRTPHBSColOL3.5M_ShootGIS_rowfiltered_logTPM.txt",header=T,sep="\t",row.names=1)[,1:208]
logTPMnew<-read.table("run125_rmCMRTPHBSColOL3.5M_ShootGIS0.25_f0.25_restore0_logTPM.txt",header=T,sep="\t",row.names=1)

datalist<-list()
datalist$logTPM<-logTPM
datalist$logTPMnew<-logTPMnew
datalist$GIS<-logTPM[rownames(logTPM)%in%GIS,]
datalist$GISnew<-logTPMnew[rownames(logTPMnew)%in%GIS,]

#pclist<-list()
#pclist$logTPM<-prcomp(datalist[[1]],center=T,scale.=T)
#pclist$logTPMnew<-prcomp(datalist[[1]],center=T,scale.=T)
#pclist$GIS<-prcomp(datalist[[3]],center=T,scale.=T)
#pclist$GISnew<-prcomp(datalist[[4]],center=T,scale.=T)


