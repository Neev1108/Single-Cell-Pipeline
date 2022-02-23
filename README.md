# Genelab Single-Cell RNA Sequencing Pipeline

> **This repository holds the slurm scripts for each cellranger command, sample Jupyter Notebooks for each Genelab experiment, and a final python script to automate the spaceflight pipeline workflow for NASA's GeneLab. The final python script will optimized for spaceflight experiments, like keeping mitochondrial reads, optimizing QC metrics for spaceflight experiments, and in-depth visualization for the user to change these inputs if need be.**

# Table of contents

- [**Software used**](#software-used)
- [**Cell Ranger**](#Cell Ranger)
  - **0. Directory Structure**
  - **1. Cellranger Reference Transcriptome**
    - [**1a. ENSEMBL Files**](#1a-ensembl-files)
    - [**1b. Cellranger mkgtf**](#1b-cellranger-mkgtf)
    - [**1c. Cellranger mkref**](#1c-cellranger-mkref)
  - **2. Cellranger Count Setup**
    - [**2a. Datasets**](#2a-datasets)
    - [**2b. Cell Ranger count**](#2b-cellranger-count)
    - [**2c. Cell Ranger aggr**](#2c-cellranger-aggr)
- [**scRCT and annotate Python Script**](#scRCT-and-annotate-Python-Script)
  - [**1. Setup and Introduction**](#1-setup-and-introduction)
    - [**1a. Packages**](#1a-packages)
    - [**1b. Prompts**](#1b-prompts)
    - [**1c. Arguments**](#1c-arguments)
    - [**1d. Marker Gene File Format**](#1d-marker-gene-file-format)

  - [**2. Quality control and User analysis tools**](#2-quality-control-and-user-analysis-tools)
    - [**2a. Highest expressed genes**](#2a-highest-expressed-genes)
    - [**2b. Gene and Cell distribution**](#2b-gene-and-cell-distribution)
    - [**2c. Remove Highly Variable Genes**](#2c-remove-highly-variable-genes)
  - [**3. Clustering**](#3-clustering)
    - [**3a. Clustering Revisted: Resolution**](#3a-clustering-revisited-resolution)
  - [**4. Gene Expression**](#4-gene-expression)
    - [**4a. Additional Visualization**](#4a-additional-visualization)
  - [**5. Data Annotation and ```annotate.py```**](#5-data-annotation)
    - [**5a. Visualization**](#5a-visualization)
      - [**5ai Annotation Revisited: Parameter tuning**](#5ai-annotation-revisited-parameter-tuning)
    - [**5b. Export annotation**](#5b-export-annotation)
    - [**5c. Arguments**](#5c-arguments)

---

## Software Used


| Program | Version | Relevant Links |
| - | - | - |
| Cellranger | 6.0.2 | https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger |
| SRA Toolkit | 2.11.0 | https://github.com/ncbi/sra-tools |
| Scanpy | 1.8 | https://scanpy.readthedocs.io/en/stable/ |
| ScoreCT | 1.0 | https://github.com/LucasESBS/scoreCT |

## 0. Directory Structure

Create a directory for your dataset with the GLDS number called GLDS-#. Within this directory, create the directory structure below by running /Scripts/scRNAseq_mkdir.sh

> 00-Raw_Data
>
> ├── [Fastq](#2b-cellranger-counts)
>
> ├── FastQC_Reports
>
> 01-CellRanger
>
> ├── [CellRanger_Output](#2b-cellranger-counts)
>
> ├── [Reference_Annotation](#2b-cellranger-counts)
>
> 02-Scanpy
>
> ├── [Adata](#5-data-annotation)
>
> ├── [Images](#2-quality-control-and-user-analysis-tools)
>
> ├── VV_logs
>
> 03-ScoreCT
>
> ├── [Annotation_Exports](#5b-export-annotation)
>
> ├── [Marker_genes](#1d-marker-gene-file-format)
>
> ├── [tSNE_UMAP](#5a-visualization)

# Cell Ranger

Cell Ranger _count_ command uses FASTQ files to generate a cell by gene count matrix. 

## 1. Cellranger Reference Transcriptome

Cell Ranger _count_ requires a reference transcriptome. Below are instructions for creating a Cell Ranger reference transcriptome from a whole genome fasta (.fa) file and a Gene Transfer Format/GTF (.gtf) file. 

The example used is _Arabidopsis Thaliana_.

---

### 1a ENSEMBL files

> ** Before running _mkref_, we need 2 important files: a whole genome fasta file, and a GTF file, preferably both from Ensembl. While files from other places might work, there might changes needed in the file, for example contigs needed to be changed in the GTF file. The whole genome fasta file will need to be one file, instead of split of files for each chromosome.**

On EMSEMBL, the whole genome file is normally named as *.primary_assembly.fa. If no primary_assembly file is available, top_level.fa will also work fine.

---

### 1b Cellranger mkgtf

(OPTIONAL) Before creating a reference transcriptome with _mkref_, we can use the Cell Ranger _mkgtf_ method to filter out annotated genes we might not need from the .gtf file. 

A SLURM script for running this command is here: /Single-Cell-Pipeline/Scripts/filter_genes_mkgtf.slurm 

Example script run:
```
sbatch filter_genes_mkgtf.slurm  Arabidopsis_thaliana.TAIR10.51.gtf Arabidopsis_filtered.gtf
```

This script runs the following command:

```
cellranger mkgtf input.gtf output.gtf  --attribute=gene_biotype:protein_coding
```

Command parameters:

* `--attribute` – Pass this parameter multiple times to filter out specific types of genes. An example is below.

> --attribute=gene_biotype:protein_coding
>
> --attribute=gene_biotype:lincRNA
>
> --attribute=gene_biotype:antisense
>
> --attribute=gene_biotype:IG_LV_gene
>
> --attribute=gene_biotype:IG_V_gene
>
> --attribute=gene_biotype:IG_V_pseudogene
>
> --attribute=gene_biotype:IG_D_gene
>
> --attribute=gene_biotype:IG_J_gene
>
> --attribute=gene_biotype:IG_J_pseudogene
>
> --attribute=gene_biotype:IG_C_gene
>
> --attribute=gene_biotype:IG_C_pseudogene
>
> --attribute=gene_biotype:TR_V_gene
>
> --attribute=gene_biotype:TR_V_pseudogene
>
> --attribute=gene_biotype:TR_D_gene
>
> --attribute=gene_biotype:TR_J_gene
>
> --attribute=gene_biotype:TR_J_pseudogene
>
> --attribute=gene_biotype:TR_C_gene

Output:

- Filtered gtf file

---

### 1c Cellranger mkref

Use Cell Ranger _mkref_ to create a reference transcriptome from a GTF file and a whole genome fasta file. If GeneLab already has an Ensembl RNAseq reference transcriptome for this organism, use those files.

A SLURM script for running this command is here: /Single-Cell-Pipeline/Scripts/make_reference.slurm 

Example script run: 
```
sbatch make_reference.slurm Arabidopsis Araport/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa Arabidopsis_filtered.gtf
```

This script runs the following command:
```
cellranger mkref --genome --fasta --genes
```

Command parameters:
* `--genome` – genome name
* `--fasta` – path to fasta file with whole reference genome
* `--genes` – path to gtf file

Output:
- fasta file
- genes file
- reference.json
- STAR index folder (most important)


We now have a folder called Arabidopsis which is our reference transcriptome.

Place this folder in the GeneLab scRNAseq Reference Transcriptome folder.

## 2. Cellranger Count Setup

The next step is to use Cell Ranger _count_ to make a cells by genes counts matrix of our dataset.

### 2a Datasets

Place your FASTQ files into the GLDS-#/00-Raw_Data directory. 

Cell Ranger requires FASTQ file names to follow the bcl2fastq file naming convention:

[Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz

Where Read Type is one of:

- I1: Sample index read (optional)
- R1: Read 1
- R2: Read 2

For example: 

- An acceptable FASTQ file name: GLDS-402_scRNA-Seq_RRRM2_Femur_BM_FLT_LAR_OLD_FO1_raw_S1_L001_R1_001.fastq.gz

- What Cell Ranger thinks the "sample name" is: "GLDS-402_scRNA-Seq_RRRM2_Femur_BM_FLT_LAR_OLD_FO1_raw"

<details>
  <summary>Example using SRA data</summary>

	If downloading a dataset from NCBI using SRA, then please first check if the data access has a BAM file. If it does, you can use cellrangers bamtofastq 	function. Most datasets do not, so we will not discuss this function.

	A quick explanation of how SRA accession IDS and numbers compare to Cellrangers definitions:

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

	> Arabidopsis_dataset
	>
	> ├── SRR13040579_S1_L001_R1_001.fastq.gz
	>
	> ├── SRR13040579_S1_L001_R2_001.fastq.gz
	>
	> ├── SRR13040580_S1_L001_R1_001.fastq.gz
	>
	> ├── SRR13040580_S1_L001_R2_001.fastq.gz
	>
	> ├── SRR13040581_S1_L001_R1_001.fastq.gz
	>
	> ├── SRR13040581_S1_L001_R2_001.fastq.gz
	>
	> ├── SRR13040582_S1_L001_R1_001.fastq.gz
	>
	> ├── SRR13040582_S1_L001_R2_001.fastq.gz
	>
	> ├── SRR13040583_S1_L001_R1_001.fastq.gz
	>
	> ├── SRR13040583_S1_L001_R2_001.fastq.gz
	>
	> ├── SRR13040584_S1_L001_R1_001.fastq.gz
	>
	> ├── SRR13040584_S1_L001_R2_001.fastq.gz
	>
	> ├── SRR13040585_S1_L001_R1_001.fastq.gz
	>
	> ├── SRR13040585_S1_L001_R2_001.fastq.gz
	>
	> ├── SRR13040586_S1_L001_R1_001.fastq.gz
	>
	> └── SRR13040586_S1_L001_R2_001.fastq.gz
	
</details>

### 2b Cell Ranger count

Once you have a Cell Ranger reference transcriptome file, and properly named FASTQ files, run Cell Ranger _count_ to generate a cells by genes counts matrix for each individual sample.

A SLURM script for running this command is here: /Single-Cell-Pipeline/Scripts/cellranger_count.slurm 

Example script run:
```
sbatch cellranger_count.slurm GLDS# Organism #ExpectedCellsPerSample
sbatch cellranger_count.slurm 402 Mus_musculus 8000
```
NOTE: the "Organism" attribute must match the formatting of the reference transcriptome folder's name.

This script runs the following command on each FASTQ file in the /00-Raw_Data directory.
```
cellranger count --id --transcriptome --fastqs --sample --expect-cells --localcores --localmem
```

Command parameters:
* `--id` – A string id of the job. This is automatically derived from each sample name.
* `--transcriptome` – path to the reference transcriptome directory. This is automatically derived from the organism name.
* `--fastqs` – path to fastqs directory. This is automatically derived from each sample name.
* `--sample` – name of sample to run. The SLURM script iterates through all FASTQ files in the directory.
* `--expected-cells` – an optimal parameter if number of cells is known beforehand. Default is 3000 cells.
* `--localcores` – number of cores to run the job. Automatically set to 78.
* `--localmem` – amount of local memory for the job. Overridden by settings in the SLURM script.

The output for each sample can be found in ```01-CellRanger/CellRanger_Output``` in a folder named for that sample.
- web summary
- metrics summary in csv format
- possorted_genome_bam file (this will be the majority of the file size, might be upwards of 5gb)
- filtered_feature_bc_matrix (this will be the most important for our Scanpy analysis)
- raw feature matrix
- analysis (some cellranger specific analysis, can be used to validate and test some scanpy downstream analysis)
- cloupe file (a file to see visualization on cellranger's cloupe software)

This script also automatically generates a CSV file in ``01-CellRanger/CellRanger_Output``` called "GLDS-#_aggr_CSV.csv" which is used for running [**2c. Cell Ranger aggr**](#2c-cellranger-aggr).


---

### 2c Cell Ranger aggr

After running Cell Ranger _count_, you can aggregate all samples together into a single cell by gene matrix using Cell Ranger _aggr_.

A SLURM script for running this command is here: /Single-Cell-Pipeline/Scripts/cellranger_aggr.slurm 

Example script run:
```
sbatch cellranger_count.slurm GLDS# 
sbatch cellranger_count.slurm 402 
```

This script runs the following command:
```
cellranger aggr --id --csv --nosecondary
```

Command parameters:
* `--id` – A string id of the job. This is automatically derived from each sample name.
* `--csv` – path to a CSV file with paths to the Cell Ranger _count_ outputs from each FASTQ. This is automatically generated by the cellranger_count.slurm script.
* `--nosecondary` – a flag to skip Cell Ranger secondary analysis.

---

# scRCT and annotate Python Script

---

## 1 Setup and Introduction

---

### 1a Packages

This program requires the following packages: ```NumPy```, ```PANDAS```, ```ScanPy```, ```seaborn```, ```matplotlib```, ```HarmonyPy```, ```scikit-misc```, ```leidenalg```, ```openpyxl```, and ```scorect_api```. However, if these packages are not installed, the program will automatically install them. If the packages still fail to load, try running the program again or install the packages on your own. NOTE: Python3 is required.

---

### 1b Prompts

This program will be an interactive experience and will provide the user opportunities to view the data and pass a parameter at each step of the process. Users can also pass the ```-disable_interrupts``` argument to disable this feature and run the program without any interruptions.

IMPORTANT NOTE: If the program is in an environment where it cannot pop-up a new window to display the plots, remember that plot images are immediately saved in the ```Images``` folder and can be viewed while running the program. (ex. running on a cluster)

### 1c Arguments

The majority of these arguments will be prompted during a program run with interruptions making most command arguments optional. The only required argument is ```--filepath``` for the directory location of the cellranger output. The default marker gene file uses CellMarker which only contains Human and Mouse genes. This would require a file to be provided for other species. Optional arguments are put in brackets ([]).

```
python3 scRCT.py --filepath=CELLRANGER_OUTPUT [--markers MARKER_FILE]  
[-disable_interrupts] [--plot PLOT] [--min_cells THRESHOLD] 
[--min_genes THRESHOLD] [--neighbors K] [--res CLUSTER_RESOLUTION]
[--genes GENES] [--species SPECIES] [--tissue TISSUE]
[--K TOP_K] [--bins M]

--filepath CELLRANGER_OUTPUT (REQUIRED)
	File path of the cellranger output folder. The output must 
	contain ```barcodes.tsv.gz```, ```features.tsv.gz```, and 
	```matrix.mtx.gz```. These must be in ```gz``` format.
-disable_interrupts
	Disable program interruptions for user analysis and 	argument prompt.
--markers MARKER_FILE (default: CellMarker)
	Provide a marker gene file to annotate the data. Format 
	described below.
--plot PLOT (default: umap)
	Visualize the data using ```umap``` or ```t-SNE```. This 
	arugment will NOT be prompted.
--min_cells THRESHOLD (default: 1)
	The minimum number of cells that a gene has to be expressed 
	in to pass the filter.
--min_genes THRESHOLD (default: 1st precentile)
	The minimum number of genes a cell has to express to pass 
	the filter.
--neighbors K (default: 15)
	The number of neighbors to "vote" for a cell in the KNN 
	algorithm. This argument will NOT be prompted.
--res CLUSTER_RESOLUTION (default: 0.5)
	A decimal number determining the amount of clusters. The 
	higher the number, the more clusters in the results.
--genes GENES (default: [])
	A list of comma (,) separated genes to visualize gene 
	expression once the data has been filtered. If not provided, 
	the top 2 genes will be displayed. 
--species SPECIES (default: "")
<<<<<<< HEAD
	Species of the sample (ex. human, mouse). Will be used in 
	naming the output files and directories. Will be used to 
	retrieve marker gene data if file not provided (Human and 
	Mouse only). Users will be prompted for species if marker 
	gene file is not provided even on ```-disable_interrupts```.
=======
	Species of the sample (ex. human, mouse). Will be used in 	naming the output files and directories.
>>>>>>> bd7b98b2111f0a682daa28a9c61dfb42179ad059
--tissue TISSUE (default: "")
	Tissue type of the sample (ex. brain, heart, etc.). Will be 
	used in naming the output files and directories. Will 
	be used to retrieve marker gene data if file not provided 
	(Human and Mouse only). Users will be prompted for tissue 
	if marker gene file is not provided even on 
	```-disable_interrupts```.
--K TOP_K (default: 300)
	Top K genes to include in annotation scoring.
--bins M (default: 15)
	M number of bins to use to divide the gene ranking in 
	annotation scoring
```

---

### 1d Marker Gene File Format

The ScoreCT cell annotation package automatically pulls mouse and human cell type annotations from CellMarkerDB. If your dataset is from another species, you must provide your own marker gene file. Please place your files under ```03-ScoreCT/Marker_genes```. 

```
Marker gene file format example (csv):
CellType, CellType, CellType
gene, gene, gene
gene, gene, gene
gene,     , gene
etc.,      ,
```

![Marker gene file viewed in Microsoft Excel](https://drive.google.com/uc?export=view&id=16mGCDYPs4OOm_fj0oCieUKgqfZRmnfO_)

The file must be in csv format. Cell names should be listed in the first row followed by their corresponding genes listed below. If a cell type has fewer genes than the maximum number of columns, simply leave those entries blank or enter 'NA'.

**Because it is more common to find files that follow this format but tranposed (i.e. cell names in the first *column* with genes list in the following *columns*), the program will automatically tranpose the data and load the data, along with saving it as ```*filename*_corrected.csv```.**

---

For each run of the program, a new ```run#_*species*_*tissue*``` directory will be created in the ```Images``` and ```tSNE_UMAP``` folder to organize each run's image output. For ```annotate.py```, it will be named ```annotate_run#_*species*_*tissue*``` instead. The name is based on the number of runs of the program, the species, and tissue inputed. Other outputs will be named similarly.

---

## 2 Quality control and User analysis tools

This section would let the user analyze the distribution of the data and select a threshold for filtering cells and genes to remove unwanted or contaminated data. Since datasets vary, it is highly encourged users take the time to determine an accurate value than rely on default values.

3 different graphs will be displayed in a new window for the user to analyze. The program will stall until the user closes the window of the graph. Only one graph will be displayed at a time.

```
Input min_genes threshold or leave blank to use default settings: 
Input min_cells threshold or leave blank to use 1:
```

Default values are provided should the user leave the input blank. 

---

### 2a Highest expressed genes

The first graph that will pop up will show boxplot distributions of the highest expressed genes. These plots will be saved as a ```.png``` image under the ```02-Scanpy/Images``` directory. It is recommended that the user note the amount of outliers and the genes listed.

![Highly expressed genes](https://drive.google.com/uc?export=view&id=1kMpqFg5JxWIHqUVD7vokHRU3O5hN07pu)

---

### 2b Gene and Cell distribution

The distribution of the number of genes expressed in cells will be shown. Users are encouraged to use this graph to determine an 'x'--the minimum number of genes--a cell must express to pass the filter. A red line is displayed on the graph that follows the cursor to help with determining this value.

![Gene distribution](https://drive.google.com/uc?export=view&id=1lhVy1kr__g1OjyGsm0fyjFR4LnbZdkvQ)

```
Input min_genes threshold or leave blank to use default settings: 
```

After closing the window for the graph, the user would then be prompted to enter the value for the minimum number of genes threshold. Leaving the prompt blank would use the default value or the argument value if it was passed. 

This would be repeated for the cell distribution, where 'x' is the minimum number of cells a gene needs to be expressed in so that gene can pass the filter.

![Cell distribution](https://drive.google.com/uc?export=view&id=1NQjn4OrRvzJRC5UBmOuIYB49q1VzFrnN)

---

### 2c Remove Highly Variable Genes

The program will then remove genes that are highly variable in the data to reduce variablity in clustering later on. No plots will be displayed and no arguments are required for this step as the values have been predetermined. To reduce the likelihood of removing noteworthy genes, parameters were adjusted to require stricter conditions to be deemed highly variable. 

<details>
  <summary>Explanation of Algorithm</summary>

  ```
  sc.pp.highly_variable_genes(adata, flavor='seurat', min_disp=2)
  ```


  This program uses R's Seurat's algorithm.

  Each gene is put into 20 'bins' based on their mean and variance. Each gene is then normalized based on the other genes in their bin. If a gene's normalized dispersion is greater or equal to a z-score of 2 (~98th percentile) AND the gene has a low mean cell count, it is marked highly variable.

</details>

Output:

- Filtered adata
- Highly variable genes removed in adata

---

## 3 Clustering

Clustering is performed using the Leiden alogirthm, an improved version of the Louvain algorithm. More can be read on the algorithms and their differences [here](https://www.nature.com/articles/s41598-019-41695-z). 

The Leiden algorithm uses the distances calculated by the KNN algorithm to perform its calculation. For those familiar with KNN, 'k' can be adjusted using ```--neighbors``` (this parameter will not be prompt during a normal program run). Note that the program will give a warning while this is performed. This just means it will proceed to calculate the PCA of the data and can be ignored. 

After this is done, a window will pop to display a plot visualizing the clusters. This will be saved under ```02-Scanpy/images```. Close the graph to move on with the program.

![Cluster graph](https://drive.google.com/uc?export=view&id=1nDhSGlVEEOrnaRvnd5zyYjiirhvJuDca)

Output:

- Cluster classification values for each cell in adata

---

### 3a Clustering Revisited: Resolution

At this stage, the user would be able to see the results of the clustering, and it may not be to their liking. Instead of running the program all over again, the user could take this opportunity to repeat the previous steps and adjust the cluster resolution. 

Cluster resolution affects the number of clusters in the output. The higher the value, the more clusters would be in the results. The default setting for cluster resolution is ```0.5```.

The user can repeatedly adjust and view the results of the cluster until they are satisfied. Each visualization will be saved under ```02-Scanpy/Images```. To continue with the program, simply give a blank input.

Output:

- Clustered data

---

## 4 Gene Expression

The program will display another plot to view gene expressions. This plot will be saved under ```02-Scanpy/images```. This color codes the strength of gene expression for each cell in the data on the plot. Purple points indicate a cell with no expression while more green and yellow points indicate higher expression.  Users can input which genes they want to view in the command arguments. If no arguments are passed, the program will sample the first 2 genes in the data.

![Gene Expression](https://drive.google.com/uc?export=view&id=1BB5c9EbUpDHJA2eJ2ri6Voi4W9Kr7X5B)

---

### 4a Additional Visualization

```
Enter the gene(s) you wish to visualize separated by commas (,). Enter "-list page_#" to view a list of genes or "-exit" to proceed with the program: Gm1992, Rp1
```

Following this is an opportunity for the user to view other gene expressions now that the data has been processed. A prompt will appear allowing users to repeatedly view as many genes to their liking. Enter a list of genes separated by commas to be viewed, then a new graph will appear showing those genes' expression. Each visualization will be saved under ```02-Scanpy/Images```.

```
Enter the gene(s) you wish to visualize separated by commas (,). Enter "-list page_#" to view a list of genes or "-exit" to proceed with the program: -list 50
Enter '-list #' to view other pages on the list.
Agrp,	Agt,	Agtpbp1,	Agtr1a
Agtrap,	Ahctf1,	Ahcy,	Ahcyl1
Ahcyl2,	Ahdc1,	Ahi1,	Ahnak
Ahnak2,	Ahrr,	Ahsa1,	Ahsa2
Aicda,	Aida,	Aif1,	Aif1l
Aifm1,	Aifm2,	Aifm3,	Aig1
```

Users could also print out a list of all the genes in alphabetical order a handful at a time. Enter ```-list ``` followed by the desired page number to view the genes.

To exit the program, simply enter ```-exit```

Output:

- Gene expression plots

---

## 5 Data Annotation

The following steps could be repeated using ```annotate.py``` should the user want to experiment with the parameters without having to go through the filtering process again. At this stage, the program will export the processed data as ```run#_*species*_*tissue*_adata_filtered.h5ad``` which can be used for ```annotate.py```. The program will annotate the data using the marker gene file and the cluster results. More about the marker gene file can be read [here](#1d-marker-gene-file-format).

---

### 5a Visualization

Once the program finishes annotating the data, the program will display a final graph of the results to view. The plot will be saved as a ```.png``` image under ```03-ScoreCT/tSNE_UMAP```. Close this window to continue the program.

![Annotated data](https://drive.google.com/uc?export=view&id=19oTGbrAJVibfyxJsisK7PtZ6xY99uR1a)

---

### 5ai Annotation Revisited: Parameter tuning

At this stage, the user would be able to see the results of the annotation, and it may not be to their liking. Instead of running the program all over again, the user could take this opportunity to repeat the previous steps and adjust the the three values that affect the algorithm: cluster resolution, K, and m. 

**Cluster resolution** affects the number of clusters in the output. The higher the value, the more clusters would be in the result. The default setting for cluster resolution is ```0.5```.

**K** represents the top K differentially expressed genes from each cluster to include or "vote" for in the scoring algorithm. A K value too low will result in falsely unidentified cells, while a K value too high would heavily affect performance only to converge towards similar results.

**m** represents the number of bins used to split the data based on their distance (calculated during clustering) so that the data can be annotated in groups with similar properties. A m value too low will result in poor annotation, while a m value too high would heavily affect performance only to converge towards similar results.

The user can repeatedly adjust and view the results of the cluster until they are satisfied. To continue with the program, simply give a blank input.

---

### 5b Export annotation

![Annotation output excel file](https://drive.google.com/uc?export=view&id=19uVckwG7VSTjeI6OY7z817G5TXKhak4r)

<<<<<<< HEAD
The last step the program will do is export an .xlsx file with the annotated data. The file will be exported to ```03-ScoreCT/Annotation_Exports``` The file will have 4 columns: 

- The sequence is the barcodes for each cell
=======
The last step program will do is export an .xlsx file with the annotated data. The file will be exorted to ```03-ScoreCT/Annotation_Exports``` The file will have 4 columns:

- The sequence is the raw cell sequence from the fastq files
>>>>>>> bd7b98b2111f0a682daa28a9c61dfb42179ad059
- n_genes: the number of genes expressed in the cell
- cluster: the cluster group the cell is in
- cell_type: annotation of the cell

Additionally, a ```run#_*species*_*tissue*_adata_filtered.h5ad``` will be exported that can be used in future projects.

Output:

- .xlsx file for the annotation
- .h5ad file of the annotated adata

---

### 5c Arguments

```annotate.py``` will run exactly the same as ```scRCT.py```'s annotation step. This can be used to re-annotate the data without having to go through the filtered process again. The only requirement is having the ```adata_filtered.h5ad``` file that is exported by ```scRCT.py``` once filtering has finished. This is the post-filtered version of the counted data and is saved under ```02-Scanpy/Adata```. The arguments for the program match ```scRCT.py```, with the exception of ```--adata```.

```
python3 annotate.py --adata=FILEPATH [-disable_interrupts]
[--markers MARKER_FILE] [--plot PLOT] [--species SPECIES] 
[--tissue TISSUE] [--K TOP_K] [--bins M]

--adata FILEPATH(REQUIRED)
	File path of the post-filtered adata_filtered.h5ad file.
-disable_interrupts
	Disable program interruptions for user analysis and 	argument prompt.
--markers MARKER_FILE (default: Human and Mouse tissues)
	Provide a marker gene file to annotate the data. Format 
	described below.
--plot PLOT (default: umap)
	Visualize the data using ```umap``` or ```t-SNE```. This 
	arugment will NOT be prompted.
--species SPECIES (default: "")
<<<<<<< HEAD
	Species of the sample (ex. human, mouse). Will be used in 
	naming the output files and directories. Will be used to 
	retrieve marker gene data if file not provided (Human and 
	Mouse only). Users will be prompted for species if marker 
	gene file is not provided even on ```-disable_interrupts```.
=======
	Species of the sample (ex. human, mouse). Will be used in 	naming the output files and directories.
>>>>>>> bd7b98b2111f0a682daa28a9c61dfb42179ad059
--tissue TISSUE (default: "")
	Tissue type of the sample (ex. brain, heart, etc.). Will be 
	used in naming the output files and directories. Will 
	be used to retrieve marker gene data if file not provided 
	(Human and Mouse only). Users will be prompted for tissue 
	if marker gene file is not provided even on 
	```-disable_interrupts```.
--K TOP_K (default: 450)
	Top K genes to include in annotation scoring.
--bins M (default: 20)
	M number of bins to use to divide the gene ranking in 
	annotation scoring
```
