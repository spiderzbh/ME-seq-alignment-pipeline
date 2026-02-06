#!/bin/bash

Name=INPUTNAME
rawdir=INPUTDIR/fastqs
dir=INPUTDIR/$Name
Species=INPUTSPECIES
Genomes=$Species
cores=16
Type=TAPS
ReadsPerBarcode=100

PythonCode='Path_to_the_python_scripts' # adjust it based on the location of python scripts
RCode='Path_to_the_R_script' # adjust it based on the location of R scripts
tssFilesPATH='/mnt/Scripts/SHARE-seq-alignment/TSSfiles/' # cleaned TSSs on unknown chr for mm10
genomeBed='/mnt/Scripts/SHARE-seq-alignment/genomeBed/'
picardPATH='/bin/picard.jar'
bowtieGenome='/mnt/Genome/bowtie2/'
starGenome='/mnt/Genome/star/'
bismarkGenome='/mnt/Genome/bismark/'
FAFiles='Path_to_the_fa_files'

# trim fastq files
mkdir $dir/$Name.$Species.$Type
mkdir $dir/trimmed
mkdir $dir/temp
mv $rawdir/$Name*fastq.gz $dir

cd $dir/trimmed
trim_galore --cores 8 --fastqc -o . $dir/$Name*R1*.fastq.gz

if [ $Genomes == "both" ]; then
	Genome1=(both)
	Genome2=(hg38 mm10)
else
	Genome1=$Genomes
	Genome2=$Genomes
fi

# align the fastq files to reference genome
cd $dir/$Name.$Species.$Type
(biscuit align -t $cores $FAFiles/$Species/$Species.fa ../trimmed/$Name.R1_trimmed.fq.gz | samtools view -@ $cores -bS > ./$Name.$Species.bam) 2> $Name.$Species.align.log

#sort the aligned bam file and remove the original bam file
samtools sort -@ 4 -m 8G $Name.$Species.bam -o $Name.$Species.st.bam
samtools index $Name.$Species.st.bam
rm $Name.$Species.bam

if [ $Genomes == "both" ] && [ ! -f $Name.mm10.st.bam ]; then
	echo "Split into hg and mm"
	chrs1=`samtools view -H $Name.$Species.st.bam | grep GRCh38 | cut -f2 | sed 's/SN://g' | awk '{if(length($0)<13)print}'`
	chrs2=`samtools view -H $Name.$Species.st.bam | grep mm10 | cut -f2 | sed 's/SN://g' | awk '{if(length($0)<13)print}'`
	samtools view -@ $cores -h $Name.$Species.st.bam `echo ${chrs1[@]}` | sed 's/GRCh38_//g' | samtools view -@ $cores -b -o $Name.hg38.st.bam &
	samtools view -@ $cores -h $Name.$Species.st.bam `echo ${chrs2[@]}` | sed 's/mm10___//g' | sed 's/mm10_//g' | samtools view -@ $cores -b -o $Name.mm10.st.bam &
	wait
	samtools index -@ $cores $Name.hg38.st.bam &
	samtools index -@ $cores $Name.mm10.st.bam &
	wait
fi

doawk() {
	chr=$1
	Name=$2
	Species2=$3
	zcat "$Name"."$Species2"."$chr".temp.bed.gz | awk '$3-$8 > 9' | awk '$9-$3 > 9' | sed 's/\_R1/\tR1/g' | cut -f2,3,4,5,8 | awk -v OFS='\t' '{if ($5 == "C") {print $2,$3,$3+1,$1,"1",$4} else if ($5 == "R") {print $2,$3,$3+1,$1,"0",$4}}' | gzip > "$Name"."$Species2"."$chr".temp.CpG.bed.gz
}
export -f doawk

echo "Genomes2 is " $Genome2

for Species2 in ${Genome2[@]}; do
	if [ -f $Name.$Species2.nonCpG.bed.gz ]; then
		echo "Found $Name.$Species2.nonCpG.bed.gz"
	else
		echo $Species2
		echo "Extract DNA methylation info, skip extracting DNA methylation"
		chrs=`samtools view -H $Name.$Species2.st.bam | grep chr | cut -f2 | sed 's/SN://g' | grep -v Y | awk '{if(length($0)<6)print}'`
		parallel 'samtools view -bS -q 1 -F 0x104 '$Name'.'$Species2'.st.bam {} > '$Name'.'$Species2'.{}.bam' ::: `echo ${chrs[@]}`

		parallel 'biscuit cinread -t cg  -p QNAME,CHRM,CRPOS,BSSTRAND,CRBASE,CQBASE,CRETENTION,QBEG,QEND $FAFiles/'$Species'/'$Species'.fa '$Name'.'$Species2'.{}.bam | gzip > '$Name'.'$Species2'.{}.temp.bed.gz' ::: `echo ${chrs[@]}` > /dev/null 2>&1
		parallel doawk {} $Name $Species2 ::: `echo ${chrs[@]}`
		cat $Name.$Species2.*.temp.CpG.bed.gz > $Name.$Species2.CpG.bed.gz
		parallel 'biscuit cinread -t ch  -p QNAME,CHRM,CRPOS,BSSTRAND,CRBASE,CQBASE,CRETENTION,QBEG,QEND $FAFiles/'$Species'/'$Species'.fa '$Name'.'$Species2'.{}.bam | gzip > '$Name'.'$Species2'.{}.temp.bed.gz' ::: `echo ${chrs[@]}` > /dev/null 2>&1
		parallel doawk {} $Name $Species2 ::: `echo ${chrs[@]}`
		cat $Name.$Species2.*.temp.CpG.bed.gz > $Name.$Species2.nonCpG.bed.gz
		rm $Name.$Species2.chr*.temp.bed.gz $Name.$Species2.chr*.temp.CpG.bed.gz
		rm $Name.$Species2.chr*.bam
		zcat  $Name.$Species2.CpG.bed.gz | awk -v OFS='\t' '{me+=$5}END{print "CpG", me/NR}' > $Name.$Species2.methy.pct.log
	fi
	if [ -f $Name.$Species2.CpG.bedGraph ]; then
		echo "Found $Name.$Species2.CpG.bedGraph"
	else
		echo "Generate bedgraph"
		echo "track type=bedGraph description=methylation level" >  $Name.$Species2.CpG.bedGraph
		zless $Name.$Species2.CpG.bed.gz | cut -f1,2,3,5 | sort --parallel=$cores -S 40G -k1,1 -k2n,2 -k4,3 | uniq -c | awk -v OFS="\t" 'NR==1 { t1=$2;t2=$3;t3=$4;count=$1; sum=$5*$1} NR > 1 {if(t1==$2 && t2==$3) {count=count+$1; sum=sum+$5*$1} else {print t1,t2,t3,sum/count*100,sum,count-sum; t1=$2;t2=$3;t3=$4;count=$1; sum=$5*$1}} END {print t1,t2,t3,sum/count*100,sum,count-sum}' >> $Name.$Species2.CpG.bedGraph
	fi
	zless $Name.$Species2.CpG.bed.gz | awk -v OFS='\t' '{key=$1 FS $2 FS $3 FS $4; count[key]++; sum[key]+=$5} END {for (k in count) {split(k, arr, FS); print arr[1], arr[2], arr[3], arr[4], sum[k],count[k]-sum[k], sum[k]/count[k]*100}}' | pigz --fast -p $cores > $Name.$Species2.CpG_Me.bed.gz
        zless $Name.$Species2.CpG_Me.bed.gz | awk -v OFS='\t' '{key=$4; count[key]++; sum[key]+=$7}END {for (k in count) {print k, count[k],sum[k]/count[k]}}' > summary_$Name.bed

	zcat $Name.$Species2.CpG.bed.gz | sort -T $dir/temp --parallel=$cores -S 40G -k1,1 -k2,2n -k3,3n -k4,4 | awk '{if ($5 == 1) {print $0 | "bgzip -@ 8 > '$Name'.'$Species2'.CpG.mC.tsv.gz"} else {print $0 | "bgzip -@ 8 > '$Name'.'$Species2'.CpG.C.tsv.gz"}}'
        zcat $Name.$Species2.nonCpG.bed.gz | sort -T $dir/temp --parallel=$cores -S 40G -k1,1 -k2,2n -k3,3n -k4,4 | awk '{if ($5 == 1) {print $0 | "bgzip -@ 8 > '$Name'.'$Species2'.nonCpG.mC.tsv.gz"} else {print $0 | "bgzip -@ 8 > '$Name'.'$Species2'.nonCpG.C.tsv.gz"}}'

done

rm -r $dir/trimmed $dir/temp

echo "TAPS pipeline finished"
