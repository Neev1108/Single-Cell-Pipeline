#!/bin/bash
#
#SBATCH --job-name="cellranger-mkref"
#SBATCH --output=cellranger-mkref.out
#SBATCH --partition=priority
#SBATCH --mem=50000
#SBATCH --mail-user=lauren.m.sanders@nasa.gov
#SBATCH --mail-type=END

. ~/.profile

echo Beginning Cellranger Make Reference Script....

start=$(date +%s)
echo "start time: $start"
echo Hostname is :     $HOSTNAME



source activate /global/smf/miniconda38_admin/envs/CellRanger_12-2021

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

