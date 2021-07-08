#!/bin/bash

#SBATCH -J cellranger_job
#SBATCH --export=ALL
#SBATCH --signal=2
#SBATCH -output="cellranger_job_counts.%J.out"
#SBATCH --mem=40g



echo Beginning Cellranger count Script....

input_dir=data3/GL_SJSU_2021_Interns/Project3_scRNAseq_pipeline
out_dir=data3/GL_SJSU_2021_Interns/Project3_scRNAseq_pipeline


start=$(date +%s)
echo "start time: $start"
echo Hostname is :     $HOSTNAME


string_id=$1
reference_transcriptome=$2
fastq_path=$3
sample=$4
expect_cells=$5
local_cores=$6
local_mem=$7



call="cellranger count --id=$string_id \
--transcriptome=$reference_transcriptome \
--fastqs=$fastq_path \
--sample=$sample \
--expect-cells=$expect_cells \
--localcores=$local_cores \
--localmem=$local_mem \
--jobmode=local"


echo $call
eval $call
echo “”


echo ""
end=$(date +%s)	
echo "end time: $end"

echo ""
runtime=$(expr $end - $start)
echo "total runtime(s): $runtime"

sec_per_min=60
sec_per_hour=3600

runtime_m=$(echo "scale=2; $runtime / $sec_per_min;" | bc)
echo "total runtime(m): $runtime_m"

runtime_h=$(echo "scale=2; $runtime / $sec_per_hour;" | bc)
echo "total runtime(hr): $runtime_h"

