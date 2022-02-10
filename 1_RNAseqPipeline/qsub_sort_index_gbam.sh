#!/bin/bash
#$ -l h_vmem=3G
#$ -l h_rt=1:0:0
#$ -N samtools_wei_V1.0
#$ -o /ebio/abt6_projects7/SHB2/tmp/qsublog  
#$ -j y

export PATH=/ebio/abt6_projects9/abt6_software/bin/samtools-1.9:$PATH  

line=$1
fn=$2
unsortedBAMdir=$3

echo `date`
echo '$fn'

cd $unsortedBAMdir/sortedgBAM

#sort sample bam
samtools sort $line > $unsortedBAMdir/sortedgBAM/$fn.S.genome.bam
echo $fn sorted
#index sample bam
samtools index $fn.S.genome.bam
echo $fn indexed
echo `date`

