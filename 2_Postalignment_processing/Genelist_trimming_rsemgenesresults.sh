#!/bin/bash

sourcedir=$1 #full path, one above the "most recent common ancester" of all .genes.results files of interest, without /
operationdir=$2 #relative path from $sourcedir, the "most recent common ancester" of all .genes.results
rDNAID=$3 #full path
TE=$4 #full path
pseudogene=$5 #full path
othergene=$6
smallgene=$7

cd $sourcedir
echo "current working directory:"
echo $sourcedir

find $sourcedir/$operationdir -name *genes.results > rsemoutput.ls.tmp
echo `wc -l rsemoutput.ls.tmp`
rawoutput=`cat rsemoutput.ls.tmp`


for line in $rawoutput; do
	dir=${line%/*genes.results}
	echo $dir
	cd $dir
	bname=${line##*/}
      	bname=${bname%.genes.results}
	echo $bname
	#grep -v ATMG $line | grep -v ATCG > ${bname}.rmChrCM.results
	#echo `wc -l  ${bname}.rmChrCM.results`
	#grep -vf $rDNAID ${bname}.rmChrCM.results > ${bname}.rmChrCM.rDNA.results
	#echo `wc -l  ${bname}.rmChrCM.rDNA.results`
	#grep -vf $TE ${bname}.rmChrCM.rDNA.results | grep -vf $pseudogene > ${bname}.rmChrCM.rDNA.TEpseudogene.results
	#echo `wc -l ${bname}.rmChrCM.rDNA.TEpseudogene.results`
grep -v ATMG $line | grep -v ATCG | grep -vf $rDNAID | grep -vf $TE | grep -vf $pseudogene | grep -vf $othergene | grep -vf $smallgene > ${bname}.rmCMRTPHBS.gene.results
	echo `wc -l ${bname}.rmCMRTPHBS.gene.results`
	cd ../..
done
