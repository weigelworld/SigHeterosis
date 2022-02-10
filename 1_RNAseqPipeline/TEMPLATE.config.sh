#!/bin/bash
#setting parameters
scriptpath='/ebio/abt6_projects7/SHB2/code/RNAseqpipeline'
pipelineoutput='/ebio/abt6_projects9/SigHeterosis_Batch1/tmp/pipelineoutput/'
parsefullrun="FALSE" # TRUE/FALSE values indicating whether all data should be ran. a fastq file list must be supplied if FALSE.
rundate=$(date +'%d_%m_%Y')

#1. copy decompress fastq files
SRAsourcedir='/ebio/abt6_sra/years/2018/11_28/SigHeterosis_Batch1'
additionalSRAsourcedir='/ebio/abt6_sra/years/2018/11_27/SigHeterosis_Batch1' #if source fastqs are in different SRA directories
cpfastqtargetdir='/ebio/abt6_projects9/SigHeterosis_Batch1/data/run125_4files' #this is a variable used in subsequent kallisto too
fastqfileprefix="*.fastq.gz"
fastqfilels='run125_4fastqfilels.txt'

#2. merge files
mergefileprefix="L001_R1_001.fastq"
trimmer="1_R1_001.fastq"
lanenum=8

#3. rsem (bowtie2) reads mapping
datadir='/ebio/abt6_projects9/SigHeterosis_Batch1/data/run125_4files/merged_reads' #directory of merged reads, the variable is used in subsequent RSEM and Hisat2
rsem_ttlsamp=4
rsemoutdir='/ebio/abt6_projects9/SigHeterosis_Batch1/data/1_analysis/2_rsem_tair10gff' #this directory is used also for getting mapping stats
rsemrefdir='/ebio/abt6_projects9/SigHeterosis_Batch1/data/1_analysis/2_rsem_tair10gff/0_TAIR10_reference_RSEM'
rsemexeID="${rundate}_run125_testpipe" #this directory is also used for getting mapping stats
rsemrefname=RSEM_TAIR10_reference

#4. sort and index bam
unsortedBAMdir=/ebio/abt6_projects9/SigHeterosis_Batch1/data/1_analysis/2_rsem_tair10gff/${rsemexeID} #this is the directory where all genome .bams are. i.e. genome.bams are within the subdirectories of this directory

#5. acquiring mapping stats from rsem.cnt. 
cntfiledir=$rsemoutdir/$rsemexeID #this directory is dependent on RSEM output in step 4
rsemCntoutputdir=$rsemoutdir/$rsemexeID/sortedgBAM
rsemCntoutfile='mappingstats_frm_RSEMcnt.csv'

