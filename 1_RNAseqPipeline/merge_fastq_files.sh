#!/bin/bash
PROJECT=$1
mergefileprefix=$2
trimmer=$3
lanenum=$4
echo dir=$PROJECT
echo mergefileprefix=$mergefileprefix
echo trimmer=$trimmer
echo total No. lane=$lanenum

mkdir -p ${PROJECT}/merged_reads
cd $PROJECT
#find . -name $mergefileprefix > ./Lane1filenames.txt #make a temporary file that contains all fastq files from one lane
ls | grep $mergefileprefix > ./Lane1filenames.txt

DATA=`cat Lane1filenames.txt` #save this into a variable (hacky way to initiate for loop)
for line in $DATA; do
	#samplename=`echo ${line%$trimmer} | awk -F'/' '{print $2}'` #before the pipe, trim the tail of the character string, after the pipe,trim the head of the character string
	samplename=`echo ${line%$trimmer}` 
	echo $samplename
	if [ `ls ${samplename}*_R1_001.fastq | wc -l` -eq $lanenum ]; then
		cat ${samplename}*_R1_001.fastq > ./merged_reads/${samplename%_L00}_merge_R1.fastq #merge fastq files carrying the same sample header (aka from all lanes) into a single fastq file in the merged_reads folder
	else
		echo ${samplename}_ERROR_filecount
		continue
	fi	
done


