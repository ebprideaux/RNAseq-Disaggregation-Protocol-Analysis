#!/bin/bash
#SBATCH -J RNAnormPipe                    # '##' means comment, 1 '#SBATCH' can be treated as the setting parameters
#SBATCH -p cpu  # can be compute0, compute2 or gpu, use "sinfo -Nel" to check the crosponding nodes under the spicified partion
###SBATCH -N 1
###SBATCH --ntasks-per-node=1
###SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH -t 100:00:00                # set time limit
###SBATCH -a 1-138%15

### Dynamic inputs

# input the raw quantified counts
input_dir=$1
counts_ext=$2

### Static inputs
gtf_fp="/stg2/data2/barton/Data/Annotation/GTF/gencode.v37.annotation.gtf"
normalizationScript_fp="/stg2/data2/barton/Code/seq_processing/RNAseq/Normalization/subscripts/geTMM_normalization.R"

############ Set up admin and locate raw data ############

### Make admin directory and progress file
cd $input_dir
cd ..

norm_dir="$(pwd)/3.Normalization"
if [ ! -d $norm_dir ]; then mkdir $norm_dir; fi

cd $norm_dir

rpk_dir="$(pwd)/1.RPK"
if [ ! -d $rpk_dir ]; then mkdir $rpk_dir; fi

getmm_dir="$(pwd)/2.geTMM"
if [ ! -d $getmm_dir ]; then mkdir $getmm_dir; fi

admin_dir="$(pwd)/0.admin"
if [ ! -d $admin_dir ]; then mkdir $admin_dir; fi

### Create file with all files to be converted from counts to RPK counts
countsinput_fp="${admin_dir}/counts_inputfile.txt"
if [ ! -a $countsinput_fp ]; then touch $countsinput_fp; fi

find $input_dir -type f -name "*${counts_ext}" >> $countsinput_fp

############ Run script ############
conda run -n R_4.0.2 Rscript $normalizationScript_fp $countsinput_fp $gtf_fp $rpk_dir $getmm_dir
