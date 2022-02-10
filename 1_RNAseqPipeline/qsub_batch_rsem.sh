#!/bin/bash
#$ -l h_vmem=5G
#$ -l h_rt=8:0:0
#$ -N rsem_wei_V1.0
#$ -pe parallel 5
#$ -j y
#$ -S /bin/bash
#$ -o /ebio/abt6_projects7/SHB2/tmp/qsublog
#######set up paths and import variable from wrapper script#####
export PATH=/ebio/abt6_projects9/abt6_software/bin/rsem:$PATH
export PATH=/ebio/abt6_projects9/abt6_software/bin/bowtie2-2.2.3:$PATH
set -xv #this shows all the lines that got executed
idx1=$1
idx2=$2
outputdir=$3
now=$4
# the now variable is the $rsemexeID
i=$5
Referencename=$6
referencedir=$7

################################################################

echo `date`
echo ${idx1}_${idx2}:rsemcalcxpn

DATA=`cat $outputdir/$now/tmp$i.txt`
cd $outputdir/$now
for line in $DATA; do
	samplename=${line##*/}
	samplename=${samplename%.fastq}
	mkdir $samplename
	cd $samplename
	rsem-calculate-expression --bowtie2 \
		--output-genome-bam \
		-p 5 \
		$line \
		$referencedir/$Referencename \
		$samplename
	cd ../
	echo $samplename completed

done

echo `date`
