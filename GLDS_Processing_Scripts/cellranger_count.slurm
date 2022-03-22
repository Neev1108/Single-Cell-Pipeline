#!/bin/bash
#
#SBATCH --job-name="cellranger"
#SBATCH --output=cellranger_out_logs/cellranger_%A_%a.out
#SBATCH --partition=priority
#SBATCH --mem=50000
#SBATCH --array=1-16%2
#SBATCH -c 80
#SBATCH --mail-user=lauren.m.sanders@nasa.gov
#SBATCH --mail-type=END


# Example script call: sbatch cellranger_count.slurm GLDS# Organism #ExpectedCellsPerSample
# Example script call: sbatch cellranger_count.slurm 402 Mus_musculus 8000

# Organism must be formatted according to the directories in the scRNAseq reference transcriptomes dir [UPDATE WITH REAL PATH]
# This script runs cellranger-count on each FASTQ sample individually 
# The FASTQ files need to be named correctly: [Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz
# This script will look for the FASTQ files in the 00-RawData/Fastq dir and will output in 01-CellRanger/CellRanger_Output



. ~/.profile
echo "cellranger"
echo ""

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

source activate /global/smf/miniconda38_admin/envs/CellRanger_12-2021

glds=$1
organism=$2
expectedCells=$3

fastqDir=/global/data/Data_Processing/scRNAseq_Datasets/GLDS_Datasets/GLDS-$glds/00-RawData/Fastq
# TODO: replace with permanent path to scRNAseq dir

outDir=/global/data/temp_scratch/lmsande2/scRNAseq/GLDS-$glds/01-CellRanger/CellRanger_Output
# TODO: replace with permanent path to scRNAseq dir

transcriptomePath=/global/data/temp_scratch/lmsande2/scRNAseq/scRNAseq_Reference_Transcriptomes/ensembl_release101/$organism
# TODO: replace with permanent path to scRNAseq transcriptomes eg /global/data/Data_Processing/scRNAseq_Datasets/Reference_Transcriptomes/ensembl_release101

# In SLURM array fashion, iterate through all the FASTQ files in the fastqDir
# grep: Grab only the R1 files since sample names are the same for R1 and R2 files
# sed: get the sample for this array iteration
# awk: cut off everything after _S1 since everything before that is the sample name
sample=$(ls $fastqDir | grep 'R1' | sed -n ${SLURM_ARRAY_TASK_ID}p | awk -F '_S1' '{print $1}')

call2="echo $sample"
eval $call2

# Run Cellranger count
call2="cellranger count --id=$sample \
--transcriptome=$transcriptomePath \
--fastqs=$fastqDir \
--sample=$sample \
--expect-cells=$expectedCells \
--localcores=78 
"

echo $call2
eval $call2
echo “”

# Move each cellranger output to the proper output directory
call3="mv $sample $outDir"
eval $call3

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

