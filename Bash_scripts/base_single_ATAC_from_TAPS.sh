#!/bin/bash

Name=INPUTNAME
dir=INPUTDIR/$Name
Species=INPUTSPECIES
cores=16
ReadsPerBarcode=100
Type=ATAC

PythonCode='Path_to_the_python_scripts' # adjust it based on the location of python scripts
RCode='Path_to_the_R_script' # adjust it based on the location of R scripts
tssFilesPATH='/mnt/Scripts/SHARE-seq-alignment/TSSfiles/' # cleaned TSSs on unknown chr for mm10
genomeBed='/mnt/Scripts/SHARE-seq-alignment/genomeBed/'
picardPATH='/bin/picard.jar'
bowtieGenome='/mnt/Genome/bowtie2/'
starGenome='/mnt/Genome/star/'

mkdir $dir/$Name.$Species.$Type
cd $dir/$Name.$Species.$Type
ln -s $dir/$Name.$Species.TAPS/$Name.$Species.st.bam ./$Name.$Species.st.bam
echo "add barcode information to each line of read"
 
samtools view -H $Name.$Species.st.bam > $Name.$Species.st.header.sam
samtools view -@ $cores $Name.$Species.st.bam | cut -f1 | sed 's/_/\t/g' | cut -f2 | sort --parallel=$cores -S 20G | uniq | awk -v OFS='\t' '{print "@RG", "ID:"$1, "SM:Barcode"NR}' > header.temp.sam

awk '/@SQ/{nr=NR} /@PG/{if(!flag) {system("cat ./header.temp.sam"); flag=1}} 1; NR==nr && !flag' ./$Name.$Species.st.header.sam | uniq > $Name.$Species.rigid.st.header.sam

cat $Name.$Species.rigid.st.header.sam <(samtools view $Name.$Species.st.bam | awk -v OFS='\t' '{temp=substr($1,length($1)-31,32); print $0"\t""RG:Z:"temp}') | samtools view -@ $cores -bS > $Name.$Species.rigid.reheader.st.bam
samtools index -@ $cores $Name.$Species.rigid.reheader.st.bam

rm header.temp.sam $Name.$Species.rigid.st.header.sam $Name.$Species.st.header.sam

chrs=`samtools view -H $Name.$Species.rigid.reheader.st.bam | grep chr | cut -f2 | sed 's/SN://g' | awk '{if(length($0)<6)print}'`
echo $chrs
samtools view -b -F 0x104 $Name.$Species.rigid.reheader.st.bam  `echo $chrs` | samtools sort -@ $cores -m 3G -n -o $Name.$Species.namesort.bam

bedtools bamtobed -i $Name.$Species.namesort.bam | sed 's/_/\t/g' | sed 's/\/1//g' | awk -v OFS="\t" '{if($7=="+"){print $1,$2,$3,$5}else if($7=="-"){print $1,$2,$3,$5}}' | sort --parallel=$cores -S 40G  -k4,4 -k1,1 -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $4, $5, $1}' | pigz --fast -p $cores > $Name.$Species.bed.gz
rm $Name.$Species.namesort.bam

bedToBam -i <(zcat $Name.$Species.bed.gz | awk -v OFS="\t" '{if($2>$3) print $1, $3, $2, $4, $5; else print}') -g $genomeBed/$Species.chrom.sizes | samtools sort -@ $cores -m 2G - > $Name.$Species.rmdup.bam
samtools index -@ $cores $Name.$Species.rmdup.bam

echo "count unfiltered reads"
zcat $Name.$Species.bed.gz | awk -v OFS='\t' '{a[$4] += $5} END{for (i in a) print a[i], i}' | awk -v OFS='\t' -v ReadsPerBarcode="$ReadsPerBarcode" '{if($1 >= ReadsPerBarcode) print }' > $Name.$Species.wdup.RG.freq.bed
Rscript $RCode/sum_reads_v2.R ./ $Name.$Species.wdup.RG.freq.bed --save
mv $Name.$Species.wdup.RG.freq.bed.csv $Name.$Species.unfiltered.counts.csv

echo "count filtered reads"
zcat $Name.$Species.bed.gz | cut -f4 | uniq -c | awk -v OFS='\t' -v ReadsPerBarcode="$ReadsPerBarcode" '{if($1 >= ReadsPerBarcode ) print }' > $Name.$Species.rmdup.RG.freq.bed
Rscript $RCode/sum_reads_v2.R ./ $Name.$Species.rmdup.RG.freq.bed --save
mv $Name.$Species.rmdup.RG.freq.bed.csv $Name.$Species.filtered.counts.csv

rm $Name.$Species.wdup.RG.freq.bed $Name.$Species.rmdup.RG.freq.bed

echo "Remove low counts barcode combination"
sed -e 's/,/\t/g' $Name.$Species.filtered.counts.csv | awk -v OFS=',' -v ReadsPerBarcode="$ReadsPerBarcode" 'NR>=2 {if($6 >= ReadsPerBarcode) print $1,$2,$3,$4,$5} '  > $Name.$Species.barcodes.txt
grep -wFf $Name.$Species.barcodes.txt  <(zcat $Name.$Species.bed.gz) | sort --parallel=$cores -S 40G -k1,1 -k2,2n -k3,3n -k4,4 | bgzip  > $Name.$Species.fragments.tsv.gz
tabix -p bed $Name.$Species.fragments.tsv.gz

# make TSS pileup fig
python2.7 $PythonCode/pyMakeVplot.py -a $Name.$Species.rmdup.bam -b $tssFilesPATH/$Species.TSS.bed -e 2000 -p ends -v -u -o $Name.$Species.RefSeqTSS

# estimate lib size
Rscript $RCode/lib_size_sc_V5_single_species.R ./ $Name $ReadsPerBarcode $Species $Type --save

rm $dir/$Name.$Species.$Type/$Name.$Species.rigid.reheader.st.bam*
echo "pipeline ATAC finished"
