pipeline to perform read mapping and qc filtering from split-lane, SE illumina sequencing
#A. Preprocessing
        #1. copy and decompress fastq files
        #2. merge fastqfiles
#B. Read mapping
        #3. RSEM(bowtie2)
          #flagstat of bam
          #QC using RSEM stats
#C. Next step (not included in this pipeline): post-alignment filtering and processing
####################################################################################################################
BEFORE RUNNING THE PIPELINE SCRIPT,

1. copy the pipeline folder into ~/code

2. create the following data directory structure:
~/tmp
	-d pipelineoutput : $pipelineoutput
	-d qsublog : qsub -o
~/data
	-d runXX_for_raw_reads : $cpfastqtargetdir
	-d 0_reference
		-d reference_fasta&GFF
		-d RSEM_reference
		-d Hisat2_reference
	-d 1_analysis
		-d 1_rsem
		-d 2_filterTPM
		
###########################################################################################################
pipeline assumes that indices have already been made for the reference genome with each mapping program
makinng indices

#rsem
#use a small shell script to generate the RSEM refererence
#e.g. /ebio/abt6_projects9/SigHeterosis_Batch1/code/MiscScript/mk_rsem_ref.sh
#e.g. /ebio/abt6_projects9/SigHeterosis_Batch1/code/MiscScript/mk_rsem_ref_tair10.sh

###########################################################################################################
to execute the pipeline:
1. check if config.sh exists in the subdirectory ~/data/runxx_for_raw_reads
2. if not taking data from a full submission, check if a list of the full path to SRA fastq files exists in ~/data/runxx_for_raw_reads
3. From ~/code/RNAseqpipeline, type in the commandline:
  bash ./pipeline.sh [PATH/TO/CONFIG.sh] > [PATH/TO/PIPELINELOG/pipeline_DATE.out] 2>&1 
the log file can be found in [PATH/TO/PIPELINELOG/pipeline_DATE.out]
