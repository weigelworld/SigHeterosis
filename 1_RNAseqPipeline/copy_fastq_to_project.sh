#!/bin/bash

#######setting variables that pass down from the pipeline####
sourcedir=$1 #where raw data are stored in SRA
targetdir=$2 #to which fastqs are to be copied
fileprefix=$3 #search criteria for the fastq names in SRA
fastfilels=$4 #name of the fastqfilelist
addsourcedir=$5 #if source fastqs are in different SRA directories (different date of demultiplexing)
parsefullrun=$6 #if false, fastfilels must be provided within the directory
######################################
if [ $parsefullrun = TRUE ];then
  find $sourcedir -name $fileprefix > $targetdir/$fastfilels
  find $addsourcedir -name $fileprefix >> $targetdir/$fastfilels
else
  if [ ! -f $targetdir/$fastfilels ];then 
	  exit 1
  else 
	  echo $fastfilels exists in $targetdir
  fi
fi

cd $targetdir

echo copy *.fastq.gz files from SRA

DATA=`cat $fastfilels`
for line in $DATA; do
	cp $line $targetdir
	
done

echo gunzip *.fastq.gz files inside my own data directory
for file in $targetdir/*fastq.gz; do
	gunzip $file 
done

