#!/bin/bash
#SBATCH -J RNA_Quant                    # job name; 3 '#' means comment, 1 '#SBATCH' can be treated as the setting parameters
#SBATCH -p cpu
###SBATCH -N 1
###SBATCH --ntasks-per-node=1
###SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
##SBATCH -t 100:00:00                # set time limit
###SBATCH -a 1-138%15

### Inputs
inputbam_dir=$1

### Static inputs
counts_dirname="htseq_counts_v37"
htseqscript_fp="/stg2/data2/barton/Code/seq_processing/RNAseq/Quantification/Subscripts/htseq_counts.sh"
gtf_fp="/stg2/data2/barton/Data/Annotation/GTF/gencode.v37.annotation.gtf"
bam_ext=".bam"
counts_ext=".counts"

### Make required folders and files
cd $inputbam_dir
cd ..
countsbase_dir=$(pwd)
counts_dir="$countsbase_dir/$counts_dirname"
if [ ! -d $counts_dir ]; then mkdir $counts_dir; fi

cd $counts_dir

admin_dir="${counts_dir}/admin"
if [ ! -d $admin_dir ]; then mkdir $admin_dir; fi

input_fp="${admin_dir}/baminputs.txt"
if [ ! -a $input_fp ]; then touch $input_fp; fi

### Make file of filepaths to be processed ###
find $inputbam_dir -type f -name "*${bam_ext}" > $input_fp

### Determine array job allocation numbers
filecount=$(wc -l < "$input_fp")
if (( $filecount > 10 )); then
  paralleljobs=10
else
  paralleljobs=$filecount
fi
arraycommand="1-${filecount}%${paralleljobs}"

echo $arraycommand
echo $input_fp

### Run script
sbatch --wait -a $arraycommand $htseqscript_fp $input_fp $gtf_fp
wait
