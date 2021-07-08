Sample command to run bash script is:

sbatch cellranger_script.sh heart_dataset refdata-gex-mm10-2020-A heart_1k_v3_fastq heart_1k_v3 1000 1 32



Notes: 

In our above example, the slurm script will set a default memory of 40gb and then use 32 gb in the cellranger command. 

For the slurm script, you can change it to sbatch --mem=70g cellranger_script.sh... for the slurm memory, but the cellranger memory must be lower. Therefore using 64 will work if using 70g in slurm memory. 

Do this if you might need more memory for the job.