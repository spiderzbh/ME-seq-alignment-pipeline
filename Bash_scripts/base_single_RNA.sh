#!/bin/bash

Name=INPUTNAME
rawdir=INPUTDIR/fastqs
dir=INPUTDIR/$Name
Species=INPUTSPECIES
cores=16
genename=gene_name
refgene=gencode
ReadsPerBarcode=100
keepMultiMapping=T

PythonCode='Path_to_the_python_scripts' # adjust it based on the location of python scripts
RCode='Path_to_the_R_script' # adjust it based on the location of R scripts
tssFilesPATH='Path_to_the_TSS_files' # cleaned TSSs on unknown chr for mm10
genomeBed='Path_to_the_Genomebed_files'
picardPATH='/bin/picard.jar'
bowtieGenome='/mnt/Genome/bowtie2/'
starGenome='/mnt/Genome/star/'
bismarkGenome='/mnt/Genome/bismark/'

# trim fastq files
mkdir $dir/trimmed
cd $dir/trimmed
trim_galore --cores 8 --fastqc -o . $rawdir/$Name*R1*.fastq.gz
mv $rawdir/$Name*fastq.gz $dir

cd ..
mkdir $dir/$Name.$Species
cd $dir/$Name.$Species

STAR --chimOutType WithinBAM --runThreadN $cores --genomeDir $starGenome/$Species/ --readFilesIn $dir/trimmed/$Name.R1_trimmed.fq.gz --outFileNamePrefix $dir/$Name.$Species/$Name.$Species. --outFilterMultimapNmax 20 --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 --outSAMattributes NH HI AS nM MD --limitOutSJcollapsed 5000000 --outSAMtype BAM Unsorted --limitIObufferSize 400000000 --outReadsUnmapped None --readFilesCommand zcat

echo "Add PT tag to label N6"
samtools view -H ./$Name.$Species.Aligned.out.bam > ./$Name.$Species.header.sam
# add PT =1 tag to dT reads; PT = 0 to N6 reads
cat ./$Name.$Species.header.sam <(samtools view -@ $cores ./$Name.$Species.Aligned.out.bam | awk -v OFS='\t' '{if ($1 ~ /NNNNNNNNNN/) print $0, "PT:i:0"; else print $0, "PT:i:1"}') | samtools view -@ $cores -bS > ./$Name.$Species.bam
rm $Name.$Species.header.sam $Name.$Species.Aligned.out.bam
rm *Log.progress.out *SJ.out.tab *.Log.out
mv ./$Name.$Species.Log.final.out ./$Name.$Species.align.log

echo "Sort $Name.$Species.bam"
samtools sort -@ 4 -m 8G $Name.$Species.bam > $Name.$Species.st.bam
samtools index -@ $cores $Name.$Species.st.bam
rm $Name.$Species.bam

samtools view -H $Name.$Species.st.bam | sed 's/chrMT/chrM/g' > $Name.$Species.st.header.sam
if [ $keepMultiMapping == "T" ]; then
       cat $Name.$Species.st.header.sam <(samtools view -@ $cores $Name.$Species.st.bam | sed 's/chrMT/chrM/g') | samtools view -@ $cores -bS -F 256 > $Name.$Species.rigid.reheader.st.bam
else
       cat $Name.$Species.st.header.sam <(samtools view -@ $cores $Name.$Species.st.bam | sed 's/chrMT/chrM/g') | samtools view -@ $cores -bS -q 30 > $Name.$Species.rigid.reheader.st.bam
fi
samtools index -@ $cores $Name.$Species.rigid.reheader.st.bam

rm $Name.$Species.st.header.sam

echo $Species

if [ $Species == 'both' ]; then
	echo "Split into hg and mm"
chrs1=`samtools view -H $Name.$Species.rigid.reheader.st.bam | grep GRCh38 | cut -f2 | sed 's/SN://g' | awk '{if(length($0)<12)print}'`
	chrs2=`samtools view -H $Name.$Species.rigid.reheader.st.bam | grep mm10 | cut -f2 | sed 's/SN://g' | awk '{if(length($0)<12)print}'`
	samtools view -@ $cores -b $Name.$Species.rigid.reheader.st.bam -o $Name.$Species.temp1.bam `echo ${chrs1[@]}`
	samtools view -@ $cores -b $Name.$Species.rigid.reheader.st.bam -o $Name.$Species.temp2.bam `echo ${chrs2[@]}`
	samtools view -@ $cores -h $Name.$Species.temp1.bam | sed 's/GRCh38_//g' | sed 's/chrMT/chrM/g'| samtools view -@ $cores -b -o $Name.hg38.rigid.reheader.st.bam
	samtools view -@ $cores -h $Name.$Species.temp2.bam | sed 's/mm10___//g' | sed 's/chrMT/chrM/g'| samtools view -@ $cores -b -o $Name.mm10.rigid.reheader.st.bam
	samtools index -@ $cores $Name.hg38.rigid.reheader.st.bam &
	samtools index -@ $cores $Name.mm10.rigid.reheader.st.bam &
	wait
	rm $Name.$Species.temp1.bam $Name.$Species.temp2.bam
	Species2=(hg38 mm10)
else
	echo "Single species is aligned"
	Species2=$Species
fi
echo ${Species2[@]}

for Species in ${Species2[@]}; do
	echo "processing $Species part"
	if [ $keepMultiMapping == "T" ]; then
		featureCounts -T $cores -Q 0 -M -a $genomeBed/gtf/$Species.$refgene.gtf -t exon -g $genename -o $Name.$Species.feature.count.txt -R BAM $Name.$Species.rigid.reheader.st.bam
		mv $Name.$Species.rigid.reheader.st.bam.featureCounts.bam $Name.$Species.exon.featureCounts.bam
		featureCounts -T $cores -Q 0 -M -a $genomeBed/gtf/$Species.$refgene.gtf -t gene -g $genename -o $Name.$Species.feature.count.txt -R BAM $Name.$Species.exon.featureCounts.bam
	else 
		featureCounts -T $cores -Q 30 -M -a $genomeBed/gtf/$Species.$refgene.gtf -t exon -g $genename -o $Name.$Species.feature.count.txt -R BAM $Name.$Species.rigid.reheader.st.bam
		mv $Name.$Species.rigid.reheader.st.bam.featureCounts.bam $Name.$Species.exon.featureCounts.bam
		featureCounts -T $cores -Q 30 -M -a $genomeBed/gtf/$Species.$refgene.gtf -t gene -g $genename -o $Name.$Species.feature.count.txt -R BAM $Name.$Species.exon.featureCounts.bam
	fi
	samtools sort -@ $cores -m 2G -o $Name.$Species.wdup.bam  $Name.$Species.exon.featureCounts.bam.featureCounts.bam
	rm $Name.$Species.exon.featureCounts.bam.featureCounts.bam
	samtools index -@ $cores $Name.$Species.wdup.bam
	## group poly T reads by UMI and cutting position
	samtools view -@ $cores $Name.$Species.wdup.bam | grep PT:i:1 | grep XT:Z: | sed 's/Unassigned_Ambiguity/discard/g' | sed 's/Unassigned_MappingQuality/discard/g' | awk 'gsub(/[_]/,"\t", $1)' | awk -v OFS='\t' '{if($NF ~/discard/){$NF=$(NF-1)} print $6, $7, $2, $3, $NF}' | sed 's/XT:Z://g' > $Name.$Species.wdup.dT.bed
	# remove dup reads
	python3 $PythonCode/rm_dup_barcode_UMI_v3.py -i $Name.$Species.wdup.dT.bed -o $Name.$Species.groups.dT.tsv --m 1
	rm $Name.$Species.wdup.dT.bed
	# group N6 reads by gene, barcode and cutting position
	samtools view -@ $cores $Name.$Species.wdup.bam | grep PT:i:0 | grep XT:Z: | sed 's/Unassigned_Ambiguity/discard/g' | sed 's/Unassigned_MappingQuality/discard/g' | awk 'gsub(/[_]/,"\t", $1)' | awk -v OFS='\t' '{if($NF ~/discard/){$NF=$(NF-1)} print $6"_"$7, $2, $NF}' | sed 's/XT:Z://g' > $Name.$Species.wdup.N6.bed
	# remove dup reads
	cat $Name.$Species.wdup.N6.bed | sort --parallel=$cores -S 24G -k1,1 -k2,2 -k3,3 | uniq -c | awk -v OFS='\t' '{print $3, $4, $1, $2}' > $Name.$Species.groups.N6.tsv
	rm $Name.$Species.wdup.N6.bed
	# convert groupped UMI to bed file
	cat ./$Name.$Species.groups.dT.tsv ./$Name.$Species.groups.N6.tsv | sort --parallel=$cores -S 24G -k1,1 -k2,2 | awk -v OFS="\t" 'NR==1 { t1=$1;t2=$2;readsum=0; umisum=0} {if(t1==$1 && t2==$2) {readsum+=$3; umisum+=1} else {print t1, t2, umisum, readsum; t1=$1;t2=$2;umisum=1;readsum=$3}} END {print t1, t2, umisum, readsum}' | pigz --fast -p $cores > $Name.$Species.bed.gz
	# count reads for RNA
	echo "count unfiltered reads"
	zcat $Name.$Species.bed.gz | awk -v OFS='\t' '{a[$1] += $4} END{for (i in a) print a[i], i}' | awk -v OFS='\t' '{if($1 >= '$ReadsPerBarcode') print }'> $Name.$Species.wdup.RG.freq.bed
	Rscript $RCode/sum_reads_v2.R ./ $Name.$Species.wdup.RG.freq.bed --save
	mv $Name.$Species.wdup.RG.freq.bed.csv $Name.$Species.unfiltered.counts.csv
	echo "count filtered reads"
	zcat $Name.$Species.bed.gz | awk -v OFS='\t' '{a[$1] += $3} END{for (i in a) print a[i], i}' | awk -v OFS='\t' '{if($1 >= '$ReadsPerBarcode') print }' > $Name.$Species.rmdup.RG.freq.bed
	Rscript $RCode/sum_reads_v2.R ./ $Name.$Species.rmdup.RG.freq.bed --save
	mv $Name.$Species.rmdup.RG.freq.bed.csv $Name.$Species.filtered.counts.csv
	rm $Name.$Species.wdup.RG.freq.bed $Name.$Species.rmdup.RG.freq.bed
	# remove barcode combination that has less then N reads
	sed -e 's/,/\t/g' $Name.$Species.filtered.counts.csv | awk -v OFS=',' 'NR>=2 {if($6 >= '$ReadsPerBarcode') print $1,$2,$3,$4,$5} '  > $Name.$Species.barcodes.txt
	grep -wFf $Name.$Species.barcodes.txt <(zcat $Name.$Species.bed.gz) | pigz --fast -p $cores > $Name.$Species.cutoff.bed.gz
	# Gene body coverage and reads distribution
	samtools view -@ $cores -s 0.5 -o ./$Name.$Species.temp.bam ./$Name.$Species.wdup.bam
	samtools index ./$Name.$Species.temp.bam
	geneBody_coverage.py -i ./$Name.$Species.temp.bam -r <(head -n 5000 $genomeBed/genomeBed/$Species.UCSC_RefSeq.bed) -o $Name.$Species
	read_distribution.py -i ./$Name.$Species.temp.bam -r $genomeBed/$Species.UCSC_RefSeq.bed > $Name.$Species.read_distribution.txt
	rm $Name.$Species.temp.bam*
	# plot reads disbution
	tail -n +5 $Name.$Species.read_distribution.txt | head -n -1 > temp1.txt
	head -n 3  $Name.$Species.read_distribution.txt | grep Assigned | sed 's/Total Assigned Tags//g' | sed 's/ //g' > temp2.txt
	Rscript $RCode/Read_distribution.R ./ $Name.$Species --save
	rm temp1.txt temp2.txt
	Rscript $RCode/UMI_gene_perCell_plot_v3.R ./ $Name.$Species --save
	Rscript $RCode/lib_size_sc_V5_single_species.R ./ $Name $ReadsPerBarcode $Species $Type --save
done
rm $Name.$Species.exon.featureCounts.bam $Name.$Species.rigid.reheader.st.bam*
rm -r $dir/trimmed

echo "RNA pipeline finished"
