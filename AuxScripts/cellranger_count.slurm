#!/bin/bash
#
#SBATCH --job-name="cellranger-count"
#SBATCH --output=cellranger-count.out
#SBATCH --partition=priority
#SBATCH --mem=50000
#SBATCH -c 64 
#SBATCH --mail-user=lauren.m.sanders@nasa.gov
#SBATCH --mail-type=END


# Example script call: sbatch cellranger_count.slurm OutputDirectory Path/To/Cellranger/Reference/Transcriptome /Path/To/Fastqs SampleID 8000 


. ~/.profile
echo "cellranger-test-1sample"
echo ""

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME

string_id=$1
reference_transcriptome=$2
fastq_path=$3
sample=$4
expect_cells=$5

call="cellranger count --id=$string_id \
--transcriptome=$reference_transcriptome \
--fastqs=$fastq_path \
--sample=$sample \
--expect-cells=$expect_cells \
"


echo $call
eval $call
echo “”

#call1 = "mv $string_id $output_directory"
#eval $call1


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

