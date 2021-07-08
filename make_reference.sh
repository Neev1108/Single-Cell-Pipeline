#!/bin/bash

#SBATCH -J cellranger_job_mkref
#SBATCH --export=ALL
#SBATCH --signal=2
#SBATCH -output="cellranger_job_counts.%J.out"
#SBATCH --mem=40g



echo Beginning Cellranger Make Reference Script....

input_dir=data3/GL_SJSU_2021_Interns/Project3_scRNAseq_pipeline
out_dir=data3/GL_SJSU_2021_Interns/Project3_scRNAseq_pipeline


start=$(date +%s)
echo "start time: $start"
echo Hostname is :     $HOSTNAME

genome=$1
fasta=$2
genes=$3


call1="cellranger mkgtf $genes \ 
$genes.filtered \
 --attribute=gene_biotype:protein_coding \
                   --attribute=gene_biotype:lincRNA \
                   --attribute=gene_biotype:antisense \
                   --attribute=gene_biotype:IG_LV_gene \
                   --attribute=gene_biotype:IG_V_gene \
                   --attribute=gene_biotype:IG_V_pseudogene \
                   --attribute=gene_biotype:IG_D_gene \
                   --attribute=gene_biotype:IG_J_gene \
                   --attribute=gene_biotype:IG_J_pseudogene \
                   --attribute=gene_biotype:IG_C_gene \
                   --attribute=gene_biotype:IG_C_pseudogene \
                   --attribute=gene_biotype:TR_V_gene \
                   --attribute=gene_biotype:TR_V_pseudogene \
                   --attribute=gene_biotype:TR_D_gene \
                   --attribute=gene_biotype:TR_J_gene \
                   --attribute=gene_biotype:TR_J_pseudogene \
                   --attribute=gene_biotype:TR_C_gene "

call2="cellranger mkref \
--genome=$genome \
--fasta=$fasta \
--genes=$genes.filtered"


echo $call1

eval $call1
echo “”

echo $call2
eval $call2





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

