# Genelab Single-Cell RNA Sequencing Pipeline

> **This repository holds the slurm scripts for each cellranger command, sample Jupyter Notebooks for each Genelab experiment, and a final python script to automate the spaceflight pipeline workflow for NASA's GeneLab. The final python script will optimized for spaceflight experiments, like keeping mitochondrial reads, optimizing QC metrics for spaceflight experiments, and indepth visualization for the user to change these inputs if need be.**

R# Table of contents

- [**Software used**](#software-used)
- [**General Pipeline**]
  - **1. Cellranger - make reference transcriptome**
    - [**1a. ENSEMBL Files**](#1a-ensembl-files)
    - [**1b. Cellranger mkgtf**](#1b-cellranger-mkgtf)
    - [**1c. Cellranger mkref**](#1c-cellranger-mkref)
    - [**1c. Test with Scripts**](#1d-test-with-scripts)
  - **2. Cellranger Counts Setup**
    - [**2a. Datasets**](#2a-datasets)
    - [**2b. Cellranger count**](#2b-cellranger-count)
  - [**3. Cellranger Output**](#3-build-star-reference)
  - [**4. Scanpy Load in**](#4-align-reads-to-reference-genome-with-star)
  - [**5. Quality Control in Scanpy**](#5-quality-control-in-scanpy)
  - [**6. KNN **](#6-count-aligned-reads-with-rsem)
  - [**7. Clustering**](#7-normalize-read-counts-perform-differential-gene-expression-analysis-and-add-gene-annotations-in-r)
    - [**7a. For Datasets with ERCC Spike-In**](#7a-for-datasets-with-ercc-spike-in)
    - [**7b. For Datasets without ERCC Spike-In**](#7b-for-datasets-without-ercc-spike-in)
  - [**6. Visualizations **](#6-count-aligned-reads-with-rsem)
  - [**6. Classification **](#6-count-aligned-reads-with-rsem)

---

## Software Used




| Program | Version | Relevant Links |
| - | - | - |
| Cellranger | 6.0.2 | https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger |
| SRA Toolkit | 2.11.0 | https://github.com/ncbi/sra-tools |
| Scanpy | 1.8 | https://scanpy.readthedocs.io/en/stable/ |

## 1. Cellranger Reference Transcriptome

To count a dataset with cellranger, a reference transcriptome, and the fastq files are required. The first step of this pipeline involves making a reference transcriptome, if not using cellranger's given transcriptomes (which are currently only human and mice transcriptomes).

For a simple example, we will be using a plant dataset, Arabidopsis Thaliana to run analysis on.

---

## 1a ENSEMBL files

> ** Before running mkref, we need 2 important files: a whole genome fasta file, and a gtf file, preferably both from Ensembl. While files from other places might work, there might changes needed in the file, for example contigs needed to be changed in the gtf file. The whole genome fasta file will need to be one file, instead of split of files for each chromosome.**

On EMSEMBL, the whole genome file is normally named as *.primary_assembly.fa. If no primary_assembly file is available, top_level.fa will also work fine.

---

## 1b Cellranger mkgtf

Before running the mkref method, we can also use cellranger's mkgtf method to filter out the gtf file of annotated genes we might not need.The command is below:

```
cellranger mkgtf input.gtf output.gtf  --attribute=gene_biotype:protein_coding

```

This command only has one parameter:

* `--attribute` – Use this command multiple times to filter our genes. An example is below.

> **
> --attribute=gene_biotype:protein_coding 
> --attribute=gene_biotype:lincRNA 
> --attribute=gene_biotype:antisense 
> --attribute=gene_biotype:IG_LV_gene 
> --attribute=gene_biotype:IG_V_gene 
> --attribute=gene_biotype:IG_V_pseudogene 
> --attribute=gene_biotype:IG_D_gene 
> --attribute=gene_biotype:IG_J_gene 
> --attribute=gene_biotype:IG_J_pseudogene 
> --attribute=gene_biotype:IG_C_gene 
> --attribute=gene_biotype:IG_C_pseudogene 
> --attribute=gene_biotype:TR_V_gene 
> --attribute=gene_biotype:TR_V_pseudogene 
> --attribute=gene_biotype:TR_D_gene 
> --attribute=gene_biotype:TR_J_gene 
> --attribute=gene_biotype:TR_J_pseudogene 
> --attribute=gene_biotype:TR_C_gene
> **

Output:

- Filtered gtf file

---

## 1c Cellranger mkref

Cellranger has a command to make transcriptomes called:

```
cellranger mkref --genome --fasta --genes.

```

Explanation of each parameter is below:

* `--genome` – genome name
* `--fasta` – path to fasta file with whole reference genome
* `--genes` – path to gtf file

Output:

- fasta file
- genes file
- reference.json
- STAR index folder (most important)

---

## 1d Test with scripts

Let's test our Arabiodopsis dataset with our Slurm scripts. After getting required files, lets make a reference transcriptome.

- filter_genes_mkgtf.sh (makes filtered gtf file)
- make_reference.sh (makes reference transcriptome directory)

For filtering genes, run:

```
sbatch filter_genes_mkgtf.sh  Arabidopsis_thaliana.TAIR10.51.gtf Arabidopsis_filtered.gtf

```

The slurm filter_genes_mkgtf.sh script already has some filters. Edit for which filters to continue with.

Now that the genes are filtered, let's make the reference transcriptome.

For making a reference transcriptome, run:

```
sbatch make_reference.sh Arabidopsis Araport/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa Arabidopsis_filtered.gtf
```

We now have a folder called Arabidopsis which is our reference transcriptome.


# 2. Cellranger Counts Setup

Our next step will be to make a counts matrix of our dataset.


## 2a Datasets

Cellranger has strict requirements for running count on datasets. The counts function will need fastq files in the format:

```
Cell Ranger requires FASTQ file names to follow the bcl2fastq file naming convention.

[Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz

Where Read Type is one of:

I1: Sample index read (optional)
R1: Read 1
R2: Read 2

incompatible: SRR6334436_1.fastq

compatible: SRR6334436_S1_L001_R1_001.fastq

```

If downloading a dataset from NCBI using SRA, then please first check if the data access has a BAM file. If it does, you can use cellrangers bamtofastq function. Most datasets do not, so we will not discuss this function.

A quick explanation of how SRA accession IDS and numbers compare to Cellrangers definitions::

>
> SRR IDs are run accessions, and the fastq files will normally have 2 fastq files, a read1 and a read 2.
>
> SRX IDS are experiment accessions, and they correspond to the libraries in cellranger.
>
> SRS IDs are sample accessions, and they correspond to cell samples in cellranger.
>
> SRP IDs are study accessions, and they can contain multiple organisms.


.SRA files downloaded will need to be extracted into fastq files using this function from SRA toolkit:

```
fastq-dump --split-files SRR6334436
```

which will output 2 files:

- SRR6334436_1.fastq
- SRR6334436_2.fastq

These will need to be renamed as shown above.

The final file structure for the fastq files will look like below (this is for Arabidopsis):

Arabidopsis_dataset
.
├── SRR13040579_S1_L001_R1_001.fastq.gz
├── SRR13040579_S1_L001_R2_001.fastq.gz
├── SRR13040580_S1_L001_R1_001.fastq.gz
├── SRR13040580_S1_L001_R2_001.fastq.gz
├── SRR13040581_S1_L001_R1_001.fastq.gz
├── SRR13040581_S1_L001_R2_001.fastq.gz
├── SRR13040582_S1_L001_R1_001.fastq.gz
├── SRR13040582_S1_L001_R2_001.fastq.gz
├── SRR13040583_S1_L001_R1_001.fastq.gz
├── SRR13040583_S1_L001_R2_001.fastq.gz
├── SRR13040584_S1_L001_R1_001.fastq.gz
├── SRR13040584_S1_L001_R2_001.fastq.gz
├── SRR13040585_S1_L001_R1_001.fastq.gz
├── SRR13040585_S1_L001_R2_001.fastq.gz
├── SRR13040586_S1_L001_R1_001.fastq.gz
└── SRR13040586_S1_L001_R2_001.fastq.gz

--

## 2b Cellranger Counts

Now we will actually run the cellranger counts method. (Note. This may be computationally intensive, mostly on RAM. Please check your system specs to make sure it works correctly)

Cellranger counts work as the following:

```
cellranger count --id --transcriptome --fastqs --sample --expect-cells --localcores --localmem

```

Explanation of each paramater is below:

* `--id` – A string id of the job. This will be the name of the output folder.
* `--transcriptome` – path to the reference transcriptome directory
* `--fastqs` – path to fastqs directory
* `--sample` – name of samples to run. If multiple samples, can have multiple inputs comma-seperated. ex. SRR13040579, SRR13040580 etc.
* `--expected-cells` – an optimal parameter if number of cells is known beforehand. Default is 3000 cells.
* `--localcores` – number of cores to run the job
* `--localmem` – amount of local memory for the job. (will use 90% of available memory if not specified so be specific)

Output:

- web summary
- metrics summary in csv format
- possorted_genome_bam file (this will be the majority of the file size, might be upwards of 5gb)
- filtered_feature_bc_matrix (this will be the most important for our Scanpy analysis)
- raw feature matrix
- analysis (some cellranger specific analysis, can be used to validate and test some scanpy downstream analysis)
- cloupe file (a file to see visualization on cellranger's cloupe software)

---

cellranger_script.sh will be our slurm script for doing the count function.

Sample command to run bash slurm script is:

```
sbatch cellranger_script.sh plant_dataset Arabidopsis Arabidopsis_dataset SRR13040579 8000 1 32

```

### Notes

> **We only use one sample from the 8 we have on Arabidopsis. To use all 8, put all sample names comma-seperated in the sample parameter.
> In our above example, the slurm script will set a default memory of 40gb and then use 32 gb in the cellranger command.
> For the slurm script, you can change it to sbatch --mem=70g cellranger_script.sh... for the slurm memory, but the cellranger memory must be lower. Therefore using 64 will work if using 70g in slurm memory.
> Do this if you might need more memory for the job. **
