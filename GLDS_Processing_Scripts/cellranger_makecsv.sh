#!/bin/bash


# Example script call: bash cellranger_makecsv.sh GLDS#


# This script creates a CSV file with sample IDs in first column and path to cellranger-count outputs in second column
# This CSV requires column names = "sample_id" and "molecule_h5" (also added in this script)

glds=$1

outDir=/global/data/temp_scratch/lmsande2/scRNAseq/GLDS-$glds/01-CellRanger/CellRanger_Output
# TODO: replace with permanent path to scRNAseq dir

fastqDir=/global/data/Data_Processing/scRNAseq_Datasets/GLDS_Datasets/GLDS-$glds/00-RawData/Fastq
# TODO: replace with permanent path to scRNAseq dir

# Write out column headers and samples and their out dirs  to the aggr CSV files
# Samples and their out dirs
# grep: Grab only the R1 files since sample names are the same for R1 and R2 files
# awk: cut off everything after _S1 since everything before that is the sample name
for sample in $(ls $fastqDir | grep 'R1' | awk -F '_S1' '{print $1}'); do 
echo $sample','$outDir/$sample'/outs/molecule_info.h5' >> $outDir/GLDS-$glds\_aggr_CSV.csv
done
# Column headers
sed -i '1s/^/sample_id,molecule_h5\n/' $outDir/GLDS-$glds\_aggr_CSV.csv
