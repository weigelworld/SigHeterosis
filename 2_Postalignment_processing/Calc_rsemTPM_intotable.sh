#!/bin/bash
basedir=$1 #full path, with /, last common ancester of all the .results files
Routfilename=$2 #the proper name of the output calculated tpm file, though it needs to be repetitively read in and 
filesuffix=$3
filelistname=$4 #the tmp file that has the list of all paths to individual results file
Rheader=$5 #the common suffix of r colnames
SampIDfield=$6 #the element in the file name that corresponds to the sample ID (SXXX)

cd $basedir
find . -name '*'$filesuffix | sort > $filelistname
echo $filelistname generated

firstfile=`head -1 $filelistname`
cut -f 1 $firstfile > $Routfilename
echo $Routfilename generated
###check the Rscript for indexing of sample ID (the SXXX number) from input file name via strsplit(). It varies from run to run
Rscript /ebio/abt6_projects7/SHB2/code/Postalignment_processing/Calc_rsemTPM_intotable.r $Routfilename $filelistname $Rheader $SampIDfield

echo 'TPM conversion finished'
sortfilename=${Routfilename%.txt}_S.txt

echo 'sorting table'
Rscript /ebio/abt6_projects7/SHB2/code/Postalignment_processing/reorder_columnorder.r $basedir $Routfilename $sortfilename


