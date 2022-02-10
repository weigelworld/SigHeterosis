!/bin/bash
datadir=$1 #the working directory. All file path can be releative to this
outputdir=$2
outfilename=$3

cd $datadir
rsemcntfile=`find . -name *.cnt`

echo "sample,total_reads,total_mapped,total_unique,total_multi" > $outputdir/$outfilename
for line in $rsemcntfile; do
  fname=${line##*/}
  fname=`basename $fname .cnt`
  echo Acquiring mapping statistics from $fname
  mapstat=`sed -n '1,2p' $line | paste -sd ' ' | awk 'BEGIN{FS=" ";OFS=","}{print $4,$2,$5,$6}'`
  statsentry=$fname,$mapstat
  echo $statsentry >> $outputdir/$outfilename
done

echo mapping stats written as $outputdir/$outfilename


