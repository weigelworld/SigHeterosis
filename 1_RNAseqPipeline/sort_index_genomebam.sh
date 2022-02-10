#make .bai from .genome.bam, move all the sorted .bam into a directory parallel to all sample bam dirs

unsortedBAMdir=$1 #directory containing all sample bam dirs
export PATH=/ebio/abt6_projects9/abt6_software/bin/samtools-1.9:$PATH 

echo `date`

mkdir $unsortedBAMdir/sortedgBAM
cd $unsortedBAMdir
unsortedBAMls=`find $unsortedBAMdir -name *.genome.bam`
cd $unsortedBAMdir/sortedgBAM
for line in $unsortedBAMls; do
   bamname=${line##*/}
   echo Sorting and indexing:$bamname
   fn=`basename $line .genome.bam` 
   samtools sort $line > $unsortedBAMdir/sortedgBAM/$fn.S.genome.bam
   echo Indexing:$fn.sorted.bam
   samtools index $fn.S.genome.bam
done

echo bam and index written in $unsortedBAMdir/sortedgBAM
echo `date`
####
#run flagstat with the sorted bam, within the sorted bam dir
cd sortedgBAM
mkdir tmpdir #make temporary directory to store intermediate flagstat files
bamlist=`ls *genome.bam`

echo "sample,ttlqcReads,secondaryReads,mappedReads" > tmpdir/flagstatsum.csv
cd tmpdir
for line in $bamlist; do
	 samtools flagstat $unsortedBAMdir/sortedgBAM/$line > $line.tmp.txt
	 flag=`sed -n '1,5 p' $line.tmp.txt | awk -F' +' '{print $1}' | paste -sd, | awk 'BEGIN{FS=",";OFS=",";}{print $1,$2,$5}'`
	 flagentry=$line,$flag
	 echo $flagentry >> flagstatsum.csv
	 #rm $line.tmp.txt
	 echo flagstat for $line calculated
 done
mv flagstatsum.csv ../
cd ..

