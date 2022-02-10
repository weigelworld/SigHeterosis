#this is to filter col-0, outlier samples, and samples with RSEM (bt2) counted mapped reads<3Mil
#then separate files into shoot and root
meta<-read.table("/ebio/abt6_projects9/SigHeterosis_Batch1/data/1_analysis/7_rsem_tair10rmChrCM/TPMlistcomparison/rsem_run125_metainfo_sorted.txt",header=T,sep="\t")
TPM<-read.table("/ebio/abt6_projects9/SigHeterosis_Batch1/data/1_analysis/7_rsem_tair10rmChrCM/TPMlistcomparison/run125_RSEM_rmCMRTPHBS_TPM_S.txt",header=T,sep="\t",row.names=1)
ol<-read.table("/ebio/abt6_projects9/SigHeterosis_Batch1/data/1_analysis/7_rsem_tair10rmChrCM/TPMlistcomparison//SampleID_2filter.txt",header=F)
ol<-paste0(ol[,1],"CMRTPHBS")
nocol<-which(meta$Genotype!="Col-0")
metancol<-meta[nocol,]
TPMncol<-TPM[,nocol]
olID<-which(colnames(TPMncol)%in%ol)
TPMncol_ol<-TPMncol[,-olID]
metancol_ol<-metancol[-olID,]
TPM3.0<-TPMncol_ol[,metancol_ol$ttlMapped.rsem/1000000>=3]
meta3.0<-metancol_ol[metancol_ol$ttlMapped.rsem/1000000>=3,]

TPM3.0S<-TPM3.0[,meta3.0$TissueType=="S"]
meta3.0S<-meta3.0[meta3.0$TissueType=="S",]

write.table(TPM3.0,"run125_RSEM_rmCMRTPHBS_rmColOL3M_seqtk28_TPM.txt",quote=F,sep="\t")
write.table(meta3.0,"run125_RSEM_rmCMRTPHBS_rmColOL3M_seqtk28_metainfo.txt",quote=F,sep="\t")
write.table(TPM3.0S,"run125_RSEM_rmCMRTPHBS_rmColOL3M_seqtk28_ShootTPM.txt",quote=F,sep="\t")
write.table(meta3.0S,"run125_RSEM_rmCMRTPHBS_rmColOL3M_seqtk28_Shootmetainfo.txt",quote=F,sep="\t")

#write.table(TPMncol,"run125_RSEM_rmCMRTPHL_rmCol_TPM.txt",quote=F,sep="\t")
#write.table(TPMncol_ol,"run125_RSEM_rmCMRTPHL_rmColOL_TPM.txt",quote=F,sep="\t")
#write.table(metancol,"run125_RSEM_rmCol_metainfo.txt",quote=F,sep="\t")
#write.table(metancol_ol,"run125_RSEM_rmCMRTPHL_rmColOL_metainfo.txt",quote=F,sep="\t")
#TPMS<-TPMncol_ol[,metancol_ol$TissueType=="S"]
#metaS<-metancol_ol[metancol_ol$TissueType=="S",]
#write.table(TPMS,"run125_RSEM_rmCMRTPHL_rmColOL_Shoot_TPM.txt",quote=F,sep="\t")
#write.table(metaS,"run125_RSEM_rmCMRTPHL_rmColOL_Shoot_metainfo.txt",quote=F,sep="\t")
