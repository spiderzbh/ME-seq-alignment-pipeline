#!/bin/bash
# input file
# 1) yaml
# 2) BCL or fastq or demultiplexed fastqs

rawdir=path_to_your_raw_data # the path to either bcl files or fastq files
dir=path_to_folder_for_aligned_data
yaml=$dir/config_nova.yaml
Project=(example_RNA example_TAPS) #add the name of your sample here, need to be the same as the one in yaml
Type=(RNA TAPS) # ATAC or TAPS (DNA methylation profile) or RNA
cores=8
Genomes=(mm10 mm10) # hg38 or mm10 or both

PythonCode='Path_to_the_python_scripts' # adjust it based on the location of python scripts
BashCode='Path_to_the_python_scripts' # adjust it based on the location of bash scripts

Start=Fastq # Bcl or Fastq, if start from fastq files, make sure you have I1, I2, R1, and R2 ready in $rawdir/fastqs
Runtype=full # QC or full,  QC only analyze 12M reads
chem=rev # rev or fwd: nova 1.5 & nextseq use rev; nova1.0 uses fwd

# real code
export SHELL=$(type -p bash)
source ~/.bashrc
export LC_COLLATE=C
export LANG=C
## this may speed up sorting by limit output to bytes

echo "the number of projects is" ${#Project[@]}
echo "Running $Runtype pipeline"

if [ ! -d $dir ]; then mkdir $dir; fi
if [ ! -d $dir/fastqs ]; then mkdir $dir/fastqs ; fi
if [ ! -d $dir/temp ]; then mkdir $dir/temp ; fi
if [ -f $dir/Run.log ]; then rm $dir/Run.log; fi

# if start with bcl file
if [ "$Start" = Bcl ]; then
    echo "Bcl2fastq"
    if [ -f $rawdir/fastqs/Undetermined_S0_R1_001.fastq.gz ] || [ -f $rawdir/fastqs/Undetermined_S1_R1_001.fastq.gz ] || [ -f $rawdir/fastqs/Undetermined_S0_L00*_R1_001.fastq.gz ]; then
        echo "Found Undetermined_S0_L001_I1_001.fastq.gz, skip Bcl2Fastq"
    else
        echo "Converting bcl to fastq"
        mkdir $rawdir/fastqs/
        bcl2fastq -p $cores -R $rawdir --mask-short-adapter-reads 0 -o $rawdir/fastqs/ --create-fastq-for-index-reads  2>>$dir/Run.log
        cd $rawdir/fastqs/
    fi
    if [ -f $dir/fastqs/filesi1.xls ]; then rm $dir/fastqs/filler*; fi

    rawdir=$rawdir/fastqs/
    Start=Fastq
fi

if [ "$Start" = Fastq ]; then
    echo "Skip bcltofastq"
    if [ -f $dir/fastqs/filesi1.xls ]; then rm $dir/fastqs/filler*; fi
    cd $rawdir
    if ls *L001_R1_001.fastq.gz 1> /dev/null 2>&1; then
        temp=$(ls *_L001_R1_001.fastq.gz)
        Run=$(echo $temp | sed -e 's/\_S0\_L001\_R1\_001.fastq.gz//')
        singlelane=F
        temp=$(ls *L00*_R1_001.fastq.gz)
        VAR=( $temp )
        nolane=${#VAR[@]}
        echo "Detected $nolane lanes"
    elif ls *S1_R1_001.fastq.gz 1> /dev/null 2>&1; then
        echo "Detected single lane"
        temp=$(ls *S1_R1_001.fastq.gz)
        Run=$(echo $temp | sed -e 's/\_\S1\_\R1\_\001.fastq.gz//')
        singlelane=T
        nolane=1
    else
        echo "No fastq with matched naming format detected; exit..."
        exit
    fi
    echo "Run number is:" $Run

    # split fastqs
    mkdir $dir/smallfastqs/
    if [ -f $dir/smallfastqs/0001.1.$Run.R2.fastq.gz ]; then
        echo "Found 0001.$Run.R2.fastq, skip split fastqs"
    else
        if [ ! -f $dir/R1.1.fastq.gz ]; then
            echo "Link fastqs"
            if [ $singlelane == T ]; then
                ln -s $rawdir/"$Run"_S1_R1_001.fastq.gz $dir/R1.1.fastq.gz
                ln -s $rawdir/"$Run"_S1_R2_001.fastq.gz $dir/R2.1.fastq.gz
                ln -s $rawdir/"$Run"_S1_I1_001.fastq.gz $dir/I1.1.fastq.gz
                ln -s $rawdir/"$Run"_S1_I2_001.fastq.gz $dir/I2.1.fastq.gz
            else
                parallel 'ln -s '$rawdir'/'$Run'_S0_L00{}_R1_001.fastq.gz '$dir'/R1.{}.fastq.gz' ::: $(seq $nolane)
                parallel 'ln -s '$rawdir'/'$Run'_S0_L00{}_R2_001.fastq.gz '$dir'/R2.{}.fastq.gz' ::: $(seq $nolane)
                parallel 'ln -s '$rawdir'/'$Run'_S0_L00{}_I1_001.fastq.gz '$dir'/I1.{}.fastq.gz' ::: $(seq $nolane)
                parallel 'ln -s '$rawdir'/'$Run'_S0_L00{}_I2_001.fastq.gz '$dir'/I2.{}.fastq.gz' ::: $(seq $nolane)
            fi
        fi
	if [ "$Runtype" = full ]; then
            # Runing full pipeline
            echo "Split fastqs to small files"
            dosplitfull(){
                fastp -i $2/R1.$1.fastq.gz -o $2/smallfastqs/$1.$3.R1.fastq.gz -S 80000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
                fastp -i $2/R2.$1.fastq.gz -o $2/smallfastqs/$1.$3.R2.fastq.gz -S 80000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
                fastp -i $2/I1.$1.fastq.gz -o $2/smallfastqs/$1.$3.I1.fastq.gz -S 80000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
                fastp -i $2/I2.$1.fastq.gz -o $2/smallfastqs/$1.$3.I2.fastq.gz -S 80000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
                wait
            }
            export -f dosplitfull
            parallel --delay 1 dosplitfull {} $dir $Run ::: $(seq $nolane)
        elif [ "$Runtype" = QC ]; then
            # Runing QC pipeline
            echo "Split fastqs to small files"
            dosplitQC(){
                fastp -i $2/R1.$1.fastq.gz -o $2/smallfastqs/$1.$3.R1.fastq.gz --reads_to_process $4 -S 4000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
                fastp -i $2/R2.$1.fastq.gz -o $2/smallfastqs/$1.$3.R2.fastq.gz --reads_to_process $4 -S 4000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
                fastp -i $2/I1.$1.fastq.gz -o $2/smallfastqs/$1.$3.I1.fastq.gz --reads_to_process $4 -S 4000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
                fastp -i $2/I2.$1.fastq.gz -o $2/smallfastqs/$1.$3.I2.fastq.gz --reads_to_process $4 -S 4000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
                wait
            }
            export -f dosplitQC
            let reads=12100000/$nolane
            parallel --delay 1 dosplitQC {} $dir $Run $reads ::: $(seq $nolane)
        else
            echo "Unknown sequencer type, exiting" && exit
        fi
    fi
  
  # trim and index fastq
    if [ -f $dir/fastp.json ]; then rm $dir/fastp.json $dir/fastp.html; fi
    ls $dir/smallfastqs | grep R1 > $dir/filesr1.xls
    ls $dir/smallfastqs | grep R2 > $dir/filesr2.xls
    ls $dir/smallfastqs | grep I1 > $dir/filesi1.xls
    ls $dir/smallfastqs | grep I2 > $dir/filesi2.xls
    cd $dir/

    if [ -f $dir/fastqs/Sub.0001.1.$Project.R1.fq.gz ] || [ -f $dir/fastqs/$Project.R1.fastq.gz ]; then
        echo "Found fq.gz, skip updating index"
    else
        echo "Update index and trim fastqs"
        noreadfile=`ls $dir/smallfastqs | grep R1 2>/dev/null | wc -l`
        noindexfile=`ls $dir/smallfastqs | grep I1 2>/dev/null | wc -l`
        if [ $noreadfile == $noindexfile ]; then
            paste filesr1.xls filesr2.xls filesi1.xls filesi2.xls | awk -v OFS='\t' '{print $1, $2, $3, $4, substr($1,1,7)}'> Filelist2.xls
            parallel --jobs $cores --colsep '\t' 'if [ -f '$dir'/fastqs/Sub.{5}.'$Project'.R1.fq.gz ]; then echo "found Sub.{5}.'$Project'.R1.fq.gz"; \
                                          else python3 '$PythonCode'/fastq.process.py3.Novadesign.v031425.py \
                                          -a '$dir'/smallfastqs/{1} -b '$dir'/smallfastqs/{2} \
                                          --c '$dir'/smallfastqs/{3} --d '$dir'/smallfastqs/{4} \
                                          --out '$dir'/fastqs/Sub.{5} \
                                          -t '$chem' -y '$yaml' && pigz --fast -p 4 '$dir'/fastqs/Sub.{5}*fq; fi' :::: Filelist2.xls
        else
            paste filesr1.xls filesr2.xls | awk -v OFS='\t' '{print $1, $2, substr($1,1,7)}'> Filelist2.xls
            parallel --jobs $cores --colsep '\t' 'if [ -f '$dir'/fastqs/Sub.{3}.'$Project'.R1.fq.gz ]; then echo "found Sub.{3}.'$Project'.R1.fq.gz"; \
                                          else python3 '$PythonCode'/fastq.process.py3.Novadesign.v031425.py -a '$dir'/smallfastqs/{1} -b '$dir'/smallfastqs/{2} \
                                          --out '$dir'/fastqs/Sub.{3} \
                                          -t '$chem' -y '$yaml' && pigz --fast -p 4 '$dir'/fastqs/Sub.{3}*fq; fi' :::: Filelist2.xls
        fi
    fi
    rm filesr1.xls filesr2.xls filesi1.xls filesi2.xls
    if [ -f Filelist2.xls ]; then
        rm Filelist2.xls
    fi
fi

#remove small fastq files
find $dir/fastqs -type f -size -10000c -delete

# merge fastq
echo "Merge fastqs"
for i in "${!Type[@]}"; do
	if [[ "${Type[i]}" == "RNA2" ]]; then
		continue  # Skip this iteration
    	fi
	process_name=${Project[i]}
	ls $dir/fastqs/Sub*"$process_name"*R1.fq.gz | xargs cat > $dir/fastqs/"$process_name".R1.fastq.gz && echo "Generated $process_name.R1.fastq.gz"
	ls $dir/fastqs/Sub*"$process_name"*R2.fq.gz | xargs cat > $dir/fastqs/"$process_name".R2.fastq.gz && echo "Generated $process_name.R2.fastq.gz"
done

# start process data
for i in "${!Type[@]}"; do
    if [[ "${Type[i]}" == "RNA2" ]]; then
        continue  # Skip this iteration
    fi
    process_name=${Project[i]}
    process_dir=$dir
    process_species=${Genomes[i]}
    echo "$process_name"
    echo "$process_dir"
    echo "$process_species"
    if [[ "${Type[i]}" == "RNA" ]]; then
            echo "process RNA part"
            mkdir $process_dir/$process_name
	    echo "working dir is: $process_dir/$process_name"
            less '$BashCode'/base_single_RNA.sh | sed "s/INPUTNAME/$process_name/" | sed "s|INPUTDIR|$process_dir|" | sed "s/INPUTSPECIES/$process_species/"  > $process_dir/$process_name/run_script.sh
            cd $process_dir/$process_name
            chmod u+x ./run_script.sh
	    ./run_script.sh
    elif [[ "${Type[i]}" == "TAPS" ]]; then
	    echo "process TAPS part"
	    mkdir $process_dir/$process_name
	    echo "working dir is: $process_dir/$process_name"
	    less '$BashCode'/base_single_TAPS.sh | sed "s/INPUTNAME/$process_name/" | sed "s|INPUTDIR|$process_dir|" | sed "s/INPUTSPECIES/$process_species/"  > $process_dir/$process_name/run_script_TAPS.sh
	    cd $process_dir/$process_name
	    chmod u+x ./run_script_TAPS.sh
	    ./run_script_TAPS.sh
	    echo "process ATAC part"
            less '$BashCode'/base_single_ATAC_from_TAPS.sh | sed "s/INPUTNAME/$process_name/" | sed "s|INPUTDIR|$process_dir|" | sed "s/INPUTSPECIES/$process_species/"  > $process_dir/$process_name/run_script_ATAC.sh
            cd $process_dir/$process_name
            chmod u+x ./run_script_ATAC.sh
            ./run_script_ATAC.sh
    else echo "done with the pipeline" && exit
    fi
done
