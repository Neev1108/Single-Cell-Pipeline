# Genelab Single-Cell RNA Sequencing Pipeline

> **This repository holds the slurm scripts for each cellranger command, sample Jupyter Notebooks for each Genelab experiment, and a final python script to automate the spaceflight pipeline workflow for NASA's GeneLab. The final python script will optimized for spaceflight experiments, like keeping mitochondrial reads, optimizing QC metrics for spaceflight experiments, and indepth visualization for the user to change these inputs if need be.**

R# Table of contents

- [**Software used**](#software-used)
- [**General Pipeline**]
- [**Cellranger**]
  - **1. Cellranger - make reference transcriptome**
    - [**1a. ENSEMBL Files**](#1a-ensembl-files)
    - [**1b. Cellranger mkgtf**](#1b-cellranger-mkgtf)
    - [**1c. Cellranger mkref**](#1c-cellranger-mkref)
    - [**1c. Test with Scripts**](#1d-test-with-scripts)
  - **2. Cellranger Counts Setup**
    - [**2a. Datasets**](#2a-datasets)
    - [**2b. Cellranger count**](#2b-cellranger-count)
  - [**3. Cellranger Output**](#3-build-star-reference)
- [**scrnaseq.py**] 
  - [**1. Quality Control in Scanpy**](#1-quality-control-in-scanpy)
    - [**1a. User analysis plots**](#1a-user-analysis-tools)
      - [**5ai Highest expressed genes overview**](#5ai-Highest-expressed-genes)
      - [**1aii Gene distribution**](#1aii-gene-distribution)
      - [**1aiii Cell distribution**](#1aiii-cell-distribution)
    - [**1b Filtering**](#1b-Filtering)
    - [**1c Batch Effect Correction**](#1c-batch-effect-correction)
  - [**2 Remove Highly Variable Genes**](#2-remove-highly-variable-genes)
    - [**2a Preparation**](#2a-normalization-and-logarithmization)
    - [**2b Remove Highly Variable Genes**](#2b-remove-highly-variable-genes)
  - [**3 Clustering**](#3-clustering)
    - [**3a K-Nearest Neighbors**](#3a-knn)
    - [**3b Leiden Clustering](#3b-leiden-clustering)
  - [**4 Visualization**](#4-visualization)
    - [**4a UMAP**](#4a-umap)
    - [**4b Cluster Overview](#4b-clusters)
      - [**4bi
    - [**5c Gene Expression](#4c-gene-expression)
      - [**5ci


---

## Software Used


| Program | Version | Relevant Links |
| - | - | - |
| Cellranger | 6.0.2 | https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger |
| SRA Toolkit | 2.11.0 | https://github.com/ncbi/sra-tools |
| Scanpy | 1.8 | https://scanpy.readthedocs.io/en/stable/ |

# Cellranger

## 1. Cellranger Reference Transcriptome

To count a dataset with cellranger, a reference transcriptome, and the fastq files are required. The first step of this pipeline involves making a reference transcriptome, if not using cellranger's given transcriptomes (which are currently only human and mice transcriptomes).

For a simple example, we will be using a plant dataset, Arabidopsis Thaliana to run analysis on.

---

### 1a ENSEMBL files

> ** Before running mkref, we need 2 important files: a whole genome fasta file, and a gtf file, preferably both from Ensembl. While files from other places might work, there might changes needed in the file, for example contigs needed to be changed in the gtf file. The whole genome fasta file will need to be one file, instead of split of files for each chromosome.**

On EMSEMBL, the whole genome file is normally named as *.primary_assembly.fa. If no primary_assembly file is available, top_level.fa will also work fine.

---

### 1b Cellranger mkgtf

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

### 1c Cellranger mkref

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

### 1d Test with scripts

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

## 2. Cellranger Counts Setup

Our next step will be to make a counts matrix of our dataset.

### 2a Datasets

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
> .
> ├── SRR13040579_S1_L001_R1_001.fastq.gz
> ├── SRR13040579_S1_L001_R2_001.fastq.gz
> ├── SRR13040580_S1_L001_R1_001.fastq.gz
> ├── SRR13040580_S1_L001_R2_001.fastq.gz
> ├── SRR13040581_S1_L001_R1_001.fastq.gz
> ├── SRR13040581_S1_L001_R2_001.fastq.gz
> ├── SRR13040582_S1_L001_R1_001.fastq.gz
> ├── SRR13040582_S1_L001_R2_001.fastq.gz
> ├── SRR13040583_S1_L001_R1_001.fastq.gz
> ├── SRR13040583_S1_L001_R2_001.fastq.gz
> ├── SRR13040584_S1_L001_R1_001.fastq.gz
> ├── SRR13040584_S1_L001_R2_001.fastq.gz
> ├── SRR13040585_S1_L001_R1_001.fastq.gz
> ├── SRR13040585_S1_L001_R2_001.fastq.gz
> ├── SRR13040586_S1_L001_R1_001.fastq.gz
> └── SRR13040586_S1_L001_R2_001.fastq.gz


### 2b Cellranger Counts

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

---

# SCRNASEQ Python Script

---

## 1 Setup and Introduction

---

## 1a Packages

This program requires the following packages: ```NumPy```, ```PANDAS```, ```ScanPy```, ```seaborn```, ```matplotlib```, ```HarmonyPy```, ```scikit-misc```, and ```leidenalg```. However, if these packages are not install, the program will automatically install them.

---

## 1b Prompts

This program will be an interactive experience and will provide the user opportunities to view the data and pass an argument at each step of the pipeline process. Users can also pass the ```-disable_interrupts``` argument to disable this feature and run the program without any interruptions. 

## 1c Arguments

The majority of these arguments will be prompted during a program run with interruptions making most command arguments optional. The only required argument is ```--filepath``` for the location and directory name of the cellranger output.

```
scrnaseq.py --filepath CELLRANGER_OUTPUT [-disable_interrupts]
[--min_cells THRESHOLD] [--min_genes THRESHOLD] [--neighbors K]
[--res CLUSTER_RESOLUTION] [--genes GENES] [--species SPECIES]
[--tissue TISSUE]

--filepath CELLRANGER_OUTPUT (REQUIRED)
	File path of the cellranger output folder.
-disable_interrupts
	Disable program interruptions for user analysis and argument 
	prompt. 
--min_cells THRESHOLD (default: 1)
	The minimum number of cells that a gene has to be expressed 
	in to pass the filter.
--min_genes THRESHOLD (default: 1st precentile)
	The minimum number of genes a cell has to express to pass 
	the filter.
--neighbors K (default: sqrt(n))
	The number of neighbors to "vote" for a cell in the KNN 
	algorithm. This argument will NOT be prompted.
--res CLUSTER_RESOLUTION (default: 1.0)
	A decimal number determining the amount of clusters. The 
	higher the number, the more clusters in the results.
--genes GENES
	A list of comma (,) separated genes to visualize gene 
	expression once the data has been filtered.
--markers MARKER_FILE (default: Provided marker gene files)
	Provide a marker gene file to annotate the data. Format 
	described below.
--species SPECIES (default: "")
	Species of the sample (ex. human, mouse, etc.). Not 	required if marker gene file was provided.
--tissue TISSUE (default: "")
	Tissue type of the sample (ex. brain, heart, etc.). Not 	required if marker gene file was provided.
--K TOP_K (default: 450)
	Top K genes to include in annotation scoring.
--bins M (default: 20)
	M number of bins to use to divide the gene ranking in 
	annotation scoring
```

---

## 2 Quality control and User analysis tools

This section would let the user analyze the distribution of the data and select a threshold for filtering cells and genes to remove unwanted or contaminated data. Since datasets vary, it is highly encourged users take the time to determine an accurate value than rely on default values.

3 different graphs will be displayed in a new window for the user to analyze. The program will stall until the user closes the window of the graph. Only one graph will be displayed at a time. Note that none of the images will be automatically saved! 

```
Input min_genes threshold or leave blank to use default settings: 
Input min_cells threshold or leave blank to use 1:
```

Default values are provided should the user not provide an input and leave the input blank. 

---

## 2a Highest expressed genes

The first graph that will pop up will show boxplox distributions of the highest expressed genes. It is recommended that the user note the amount of outliers and the genes listed.

---

## 2b Gene distribution

The distribution of the number of genes expressed in cells will be shown. Users are encouraged to use this graph to determine an 'x' value as the minimum number of genes required to pass the filter. A red line is displayed on the graph that follows the cursor to help with determining this value.

```
Input min_genes threshold or leave blank to use default settings: 
```

After closing the window for the graph, the user would then be prompted to enter the value for the minimum number of genes threshold. Leaving the prompt blank would use the default value or the argument value if it was passed.



---

## 2c Cell distribution

Next, the distribution of the number of cells that has a certain gene expressed. Users are encouraged to use this graph to determine an 'x' value as the minimum number of cells required to pass the filter. A red line is displayed on the graph that follows the cursor to help with determining this value.

```
Input min_cells threshold or leave blank to use 1:
```

After closing the window for the graph, the user would then be prompted to enter the value for the minimum number of cells threshold. Leaving the prompt blank would use the default value or the argument value if it was passed.

Output: 
-Filtered adata
-List of genes that have been filtered

---

## 3 Remove Highly Variable Genes

The program will then remove genes that are highly variable in the data to reduce variablity in clustering later on. No arguments are required for this step and the values have been predetermined. To reduce the liklihood of removing noteworthy genes, parameters were adjusted to require stricter conditions to be deemed highly variable. 

```
sc.pp.highly_variable_genes(adata, flavor='seurat', min_disp=2)
adata = adata[:, adata.var.highly_variable==False]
```

<details>
  <summary>Explanation of the Algorithm</summary>

  This program use's R's Seurat's algorithm.

  Each gene is put into 20 'bins' based and their mean and variance. Each gene is then normalized based on the other genes in their bin. If a gene's normalized dispersion is greater or equal to a z-score of 2 (~98th percentile) AND the gene has a low mean cell count, it is marked highly variable.

</details>

Output
-Highly variable genes removed in adata

---

## 4 Clustering

Once the data has been cleaned through filtering, the data is ready to be clustered to determine similarities and pattersns within the data.

---

## 4a KNN

Before Leiden clustering could be performed, the distance must be calculated for each cell. This is done using the K-Nearest Neighbors algorithm (KNN) which calculates euclidean distances based on the "votes" of its neighbors. The top number of neighbors to vote (also known as "k") could be provided by the user using ```--neighbors```. There will not a prompt for this and providing a new 'k' is entirely optional.

Note: the program will give a warning while this is performed. This just means it will proceed to calculate the PCA of the data and can be ignored. It should also be noted this is the most computationally intense portion of the program.

Output:
-Distances caculated for each cell in adata
-PCA calculated in adata

---

## 4b Leiden Clustering

Clustering is performed using the Leiden alogirthm, an improved version of the Louvain algorithm. More can be read on the algorithms and theirs differences [here](https://www.nature.com/articles/s41598-019-41695-z#:~:text=Unlike%20the%20Louvain%20algorithm%2C%20the,moved%20to%20a%20different%20community.). The Leiden algorithm uses the distances calculated in ```neighbors``` to perform its calculation. 

Output:
-Cluster classification values for each cell in adata

---

## 4c Visualization

A window will pop to display a UMAP graph visualizing the clusters. Close the graphs to move on with the program. Note that the images will be automatically saved!

---

## 4d Cluster Revisited: Resolution

At this stage, the user would be able to see the results of the clustering and it may not be to their liking. Instead of running the program all over again, the user could take this opportunity to repeat the previous steps and adjust the *cluster resolution. 

Cluster resolution affects the number of clusters in the output. The higher the value, the more clusters would be in the result. The default setting for cluster resolution is ```1.0```.

The user can repeatedly adjust and view the results of the cluster until they are satisfied. To continue with the program, simply give a blank input.

## 5 Visualization

## 5a Gene Expression

The program will display another UMAP graph, this time to view gene expression. This color codes the strength of gene expression for each cell in the data on the umap. Purple points indicate a cell with no expression while more green and yellow points indicate higher expression.  Users can input which genes they want to view in the command arguments. If no arguments are passed, the program will display the first 2 genes' plot in adata as a sample.

---

## 5b Additional Visualization

```
Enter the gene(s) you wish to visualize separated by commas (,). Enter "-list page_#" to view a list of genes or "-exit" to proceed with the program: Gm1992, Rp1
```

Following this is an opportunity for the user to view other genes now that the data has been processed. A prompt will appear allowing users to repeatedly view as many genes to their liking. Enter a list of genes separated commas to be view, then a new graph will appear showing those genes' expressing

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
-Gene expression umap plots

---

## 6 Data Annotation

The following steps could be repeated using ```annotate.py``` should the user want to experiment with the parameters without having to go through the filtering process again. At this stage, the program will export the processed data as ``adata.h5ad``` which is used in ```annotate.py```.

---

## 6a Marker Gene File

---

## 6ai Default Marker Gene Files

---

## 6b Visualization

Once the program finishes annotating the data, the program will display a UMAP of the clustered results for the user to view. Close this window to continue the program. Note that this image will not be automatically saved!

---

## 6c Parameter tuning (```annotate.py``` only)

## 6d Export annotation
