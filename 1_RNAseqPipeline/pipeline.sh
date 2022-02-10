#!/bin/bash
#pipeline to perform read mapping and qc filtering from split-lane, SE illumina sequencing
#MAKESURE TO READ AND COMPLETE INSTURCTIONS IN RNAseqpipelineREADME.txt before running the pipeline
#Double check if config.sh is within the target directory
#if running pipeline on a subset of the samples, check if the fastq file list is within the target directory

############################################################################################
# setting parameters
echo Running RNAseq Read mapping pipeline
echo Pipeline version: 1.0
echo Author: Wei Yuan
echo $'\n'
echo `date`
echo Sourcing configurations...
config=$1
source $config

echo log files stored within $pipelineoutput $'\n'
echo original fastq file locations: 
echo $SRAsourcedir
echo $additionalSRAsourcedir
echo $'\n'

echo local merged reads within $cpfastqtargetdir/merged_reads
echo rsem results in $rsemoutdir/$rsemexeID
echo sorted rsem genome .bam stored in: $rsemoutdir/$rsemexeID/sortedgBAM 

##########################################################################################

#A. Preprocessing
#1. copy and decompress fastq files
echo `date`
echo running SE-RNAseq read mapping pipeline, log files stored in $pipelineoutput
if [ $parsefullrun = TRUE ]; then
	echo "processing full run mode"
else
	echo "processing partial run mode, a fastq file list must be supplied"
fi
echo $'\n'
echo '#######################'
echo 1. copy and decompress fastq files
echo running script "copy_fastq_to_project.sh"
echo SRAsourcedir: $SRAsourcedir
echo AdditionalSRNdir: $additionalSRAsourcedir
echo copying fastqs to: $cpfastqtargetdir
echo fastqfilelist: $fastqfilels
$scriptpath/copy_fastq_to_project.sh "$SRAsourcedir" "$cpfastqtargetdir" "$fastqfileprefix" "$fastqfilels" "$additionalSRAsourcedir" "$parsefullrun" > $pipelineoutput/Copyfastqfiles$(date +'%d_%m_%Y').out 2>&1

#2. merging fastqfiles
echo '#######################'
echo $'\n'
echo `date`
echo merging fastqfiles within $cpfastqtargetdir
echo No.lanes: $lanenum
$scriptpath/merge_fastq_files.sh "$cpfastqtargetdir" "$mergefileprefix" "$trimmer" "$lanenum" > $pipelineoutput/mergefastq$(date +'%d_%m_%Y').out 2>&1
echo fastqmerged
rm ./Lane1filenames.txt
rm -f $cpfastqtargetdir/*.fastq

#B. Read mapping
#3. RSEM(bowtie2)
echo '#######################'
echo $'\n'
echo `date`
echo performing rsem/bowtie2 alignment
echo writing results to $rsemoutdir
echo execution ID: $rsemexeID
echo number of samples: $rsem_ttlsamp
echo using index: $rsemrefdir/$rsemrefname

$scriptpath/multisub_rsem.sh "$rsemoutdir" "$rsemexeID" "$datadir" "$rsem_ttlsamp" "$rsemrefdir" "$rsemrefname" "$scriptpath" > $pipelineoutput/RSEM$(date +'%d_%m_%Y').out 2>&1
#run loop until rsem finishes. The qsub script of rsem is such, that unless scripts are done or killed, qstat would not be zero
until [ `qstat | grep rsem | wc -l` -eq 0 ]; do
	njobs=`qstat | grep rsem | wc -l`
	echo waiting for $njobs qsubs to finish
	echo `date`
	sleep 5m
done
#check if the anticipated files are there
if [ `find $rsemoutdir/$rsemexeID -name *genes.results | wc -l` -eq $rsem_ttlsamp ]; then
	echo RSEM read mapping finished, $rsem_ttlsamp samples
else
	no_finished=`find $rsemoutdir/$rsemexeID -name *genes.results | wc -l`
	echo incorrect sample num: $no_finished
	exit 1
fi

#4. sort and index bam, run flagstat. genome.bam needs to be first sorted and then indexed, in order to be visualized with IGV/IGB. 
#The first half of the script sort and index all genome bam, move sorted bam into a directory (sortedgBAM) on the same level with all sample bam directories.
#The second half of the same script perform flagstat, and output csv file into the sortedgBAM directory
echo '#######################'
echo $'\n'
echo `date`
echo sorting, indexing genome bam, calculating flagstat
echo sorting bam from: $unsortedBAMdir

$scriptpath/batch_sort_index_genomebam.sh "$unsortedBAMdir" > $pipelineoutput/samtools_SI$(date +'%d_%m_%Y').out 2>&1

until [ `qstat | grep samtools | wc -l` -eq 0 ]; do
        echo waiting for samtools qsubs to finish
        echo `date`
        sleep 1m
done

echo finished sorting bams

#5. acquire mapping stats directly from RSEM (bowtie2 mapping), output as .csv file 
echo '######################'
echo $'\n'
echo `date`
echo Acquireing mapping statistics from RSEM .cnt files

$scriptpath/Mapping_stats_from_RSEMcntfile.sh "$cntfiledir" "$rsemCntoutputdir" "$rsemCntoutfile" > $pipelineoutput/MappingStatsRSEM_$(date +'%d_%m_%Y').out 2>&1

