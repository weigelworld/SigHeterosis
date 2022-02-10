#!/bin/bash
####################################
README
#this is a semi-automated pipeline. The adjust the config in the config section, and the commands for running each individual section is at the end of this script.
####################################
#CONFIG SECTION
pipelineoutput='/ebio/abt6_projects7/SHB2/tmp/pipelineoutput/'

#config for Genelist_trimming_rsemgenesresults.sh
sourcedir='/ebio/abt6_projects7/SHB2/data/1_analysis/2_rsem'
operationdir='20_01_2021_run197_seqtk'
rDNAID='/ebio/abt6_projects9/SigHeterosis_Batch1/data/0_reference/TAIR10_rDNAcluster_BLASThits_GeneID.txt'
TE='/ebio/abt6_projects9/SigHeterosis_Batch1/data/0_reference/TAIR10_annotatedTEinGFF3_GeneID.txt'
pseudogene='/ebio/abt6_projects9/SigHeterosis_Batch1/data/0_reference/TAIR10_annotatedPSEUDOGENEinGFF3_GeneID.txt'
othergene='/ebio/abt6_projects9/SigHeterosis_Batch1/data/0_reference/Highexp_ChloroplastLoc_batchPCextremeload_GeneIDrnd2.txt'
smallgene='/ebio/abt6_projects9/SigHeterosis_Batch1/data/0_reference/EffectiveLength_lt150_GeneID.txt'
##moving trimmed genes.results into a separate dir
movedir='/ebio/abt6_projects7/SHB2/data/1_analysis/3_rsem_tair10rmChrCM/210121_RSEMcounts_rmCMRTPHBS'
filesuffix='rmCMRTPHBS.gene.results' #this will also be used in the Calc_rsemTPM_intotable.sh
##

#config for rename_Rsemfiles.sh and Calc_rsemTPM_intotable.sh
basedir='/ebio/abt6_projects7/SHB2/data/1_analysis/3_rsem_tair10rmChrCM/'
operationdir='210121_RSEMcounts_rmCMRTPHBS'

#config for Calc_rsemTPM_intotable.sh
basedir='/ebio/abt6_projects7/SHB2/data/1_analysis/3_rsem_tair10rmChrCM/210121_RSEMcounts_rmCMRTPHBS'
Routfilename='run197_RSEM_rmCMRTPHBS_seqtk_TPM.txt'
filesuffix='rmCMRTPHBS.gene.results' #the genes results file to be processed
filelistname='tmp.filels.tmp'
Rheader='CMRTPHBS'
SampIDfield=2 #which element in the sample name is the sample ID (SXXX)

#config for Calc_log2TPM.r
TPMdir='/ebio/abt6_projects7/SHB2/data/1_analysis/3_rsem_tair10rmChrCM/'
infilename=' run197_RSEM_rmCMRTPHBS_rmSeqtkDup_TPM_S.txt'
outfilename='run197_RSEM_rmCMRTPHBS_Sorted_rmSeqtkDup_log2TPM.txt'

#config for simple PCA (outlier identification)
wd='/ebio/abt6_projects7/SHB2/data/1_analysis/3_rsem_tair10rmChrCM/'
logTPMin='run197_RSEM_rmCMRTPHBS_Sorted_rm2MSeqtkDup_log2TPM.txt'
metain='Run197_SHB2_rm2MSeqtkDup_Meta.txt'
cutoff=350 #max number of dropout per gene
outbasename='run197_rmCMRTPHBS_rm2MSeqtkDup_PCA'
TPMin='run197_RSEM_rmCMRTPHBS_rm2MSeqtkDup_TPM_S.txt'

###################################
#trim gene list from rsem genes results
bash /ebio/abt6_projects7/SHB2/code/Postalignment_processing/Genelist_trimming_rsemgenesresults.sh $sourcedir $operationdir $rDNAID $TE $pseudogene $othergene $smallgene > $pipelineoutput$(date +'%d_%m_%Y')_Genelisttrim.out 2>&1

#copy genes results into a single directory
cd $movedir
find $sourcedir/$operationdir -name '*'$filesuffix -exec cp {} $movedir \;

#convert counts into TPM
bash /ebio/abt6_projects7/SHB2/code/Postalignment_processing/Calc_rsemTPM_intotable.sh $basedir $Routfilename $filesuffix $filelistname $Rheader $SampIDfield > $pipelineoutput$(date +'%d_%m_%Y')_CalcrsemTPM.out 2>&1

#manually replace the first column name ("GeneID") into tab in the TPM file, to conform with rownames=T option in R.

#calculate log2TPM
#the script must be operated on TPM files with gene ID as rownames
#edit Calc_log2TPM.r to choose between log2(TPM) or log2(TPM+1)
cd $TPMdir
Rscript /ebio/abt6_projects7/SHB2/code/Postalignment_processing/Calc_log2TPM.r $infilename $outfilename

#simple row filter and PCA
Rscript /ebio/abt6_projects7/SHB2/code/Postalignment_processing/PCAoutlierID.r $wd $logTPMin $metain $cutoff $outbasename $TPMin

