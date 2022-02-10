#!/bin/bash
#$ -l h_vmem=5G
#$ -l h_rt=6:0:0
#$ -N rsem_wei_V1.0
#$ -pe parallel 5
#$ -S /bin/bash
#$ -o /ebio/abt6_projects7/SHB2/tmp/qsublog
#$ -j y
#######set up paths and import variable from wrapper script#####
export PATH=/ebio/abt6_projects9/abt6_software/bin/rsem:$PATH
export PATH=/ebio/abt6_projects9/abt6_software/bin/bowtie2-2.2.3:$PATH

datadir=$1
outputdir=$2
referencedir=$3
Referencename=$4
now=$5
samplename=$6
line=$7
################################################################

echo `date`
echo $samplename:rsemcalcxpn
mkdir $outputdir/$now/$samplename
cd $outputdir/$now/$samplename
rsem-calculate-expression --bowtie2 \
	--output-genome-bam \
	-p 5 \
	$line \
	$referencedir/$Referencename \
	$samplename

echo `date`
