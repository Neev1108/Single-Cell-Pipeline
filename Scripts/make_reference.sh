#!/bin/bash

#SBATCH -J cellranger_job_mkref
#SBATCH --export=ALL
#SBATCH --output="cellranger_job_mkref.%J.out"
#SBATCH --mem=40g



echo Beginning Cellranger Make Reference Script....



start=$(date +%s)
echo "start time: $start"
echo Hostname is :     $HOSTNAME

genome=$1
fasta=$2
genes=$3


call="cellranger mkref --genome=$genome --fasta=$fasta --genes=$genes"


echo ""

echo $call
eval $call




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

