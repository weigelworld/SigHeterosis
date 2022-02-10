#!/bin/bash
#use this file to run rsem. taking the config of ttl sample number, and determine which rsem script to call (individual or batch). The key is to limit total # of simultaneous RSEM runs to 40
rsemoutdir=$1
rsemexeID=$2
datadir=$3
rsem_ttlsamp=$4
rsemrefdir=$5
rsemrefname=$6
scriptpath=$7
if [ -d $rsemoutdir/$rsemexeID ]; then
	echo "rsem directory exists"
else
	mkdir $rsemoutdir/$rsemexeID
fi	
cd $rsemoutdir/$rsemexeID

if [ -f tmp.txt ] ; then
	echo "tmp file exists"
else
	ls $datadir/*fastq > tmp.txt
fi

if [ $rsem_ttlsamp -lt 40 ]
then
  echo running rsem individual submission
  DATA=`cat tmp.txt`
  for line in $DATA; do
	samplename=${line##*/}
	samplename=${samplename%.fastq}
	mkdir $samplename
	cd $samplename
	qsub $scriptpath/qsub_individual_rsem.sh "$datadir" "$rsemoutdir" "$rsemrefdir" "$rsemrefname" "$rsemexeID" "$samplename" "$line"
	cd ../
	echo $samplename has been submitted for rsem_calc_expression
  done

else
  echo running rsem batch submission
  let loopnum=$rsem_ttlsamp/24 #never let more than 40 samples to run simultaneously--this is subsequently changed to 24 for this particular run
  for((i=1;i<=24;i++));do
	let idx1=$loopnum*$i-$loopnum+1
	let idx2=$loopnum*$i
	sed -ne "$idx1,$idx2 p" $rsemoutdir/$rsemexeID/tmp.txt > $rsemoutdir/$rsemexeID/tmp$i.txt
	qsub $scriptpath/qsub_batch_rsem.sh "$idx1" "$idx2" "$rsemoutdir" "$rsemexeID" "$i" "$rsemrefname" "$rsemrefdir"
        echo ${idx1}_${idx2} has been submitted for rsem_calc_expression
  done
 # submit the last few files
  	let idx1=$loopnum*24+1
  	let idx2=$rsem_ttlsamp
  	let i=25
  	sed -ne "$idx1,$idx2 p" $rsemoutdir/$rsemexeID/tmp.txt > $rsemoutdir/$rsemexeID/tmp$i.txt
   	qsub $scriptpath/qsub_batch_rsem.sh "$idx1" "$idx2" "$rsemoutdir" "$rsemexeID" "$i" "$rsemrefname" "$rsemrefdir"
  	echo ${idx1}_${idx2} has been submitted for rsem_calc_expression

fi


