# -*- coding: utf-8 -*-
"""scRNAseq Tutorial Neuron Dataset.ipynb
Automatically generated by Colaboratory.
Original file is located at
    https://colab.research.google.com/drive/1Ir3kLJnSJcWnXfZZZ6qX4oIDIbU6IZ79
"""

#Import or Install Libraries
while True:
    try:
        print("Loading packages...")
        import sys
        import os
        import numpy as np
        import pandas as pd
        import scanpy as sc
        import seaborn as sb
        import matplotlib.pyplot as plt
        import subprocess
        break
    except:
        print("***Some packages have not been installed. Installing now...***")
        import urllib.request

        # Retrieve installer if not available
        remove = False
        if not os.path.exists("get-pip.py"):
            urllib.request.urlretrieve("https://bootstrap.pypa.io/get-pip.py", "get-pip.py")
            remove = True
        subprocess.check_call([sys.executable, "get-pip.py"])
        # Download and install packages if not installed
        subprocess.check_call([sys.executable, "-m", "pip", "install", "numpy"])
        subprocess.check_call([sys.executable, "-m", "pip", "install", "pandas"])
        subprocess.check_call([sys.executable, "-m", "pip", "install", "seaborn"])
        subprocess.check_call([sys.executable, "-m", "pip", "install", "scanpy"])
        subprocess.check_call([sys.executable, "-m", "pip", "install", "matplotlib"])

        # Packages used but not included in scanpy package
        subprocess.check_call([sys.executable, "-m", "pip", "install", "harmonypy"])
        subprocess.check_call([sys.executable, "-m", "pip", "install", "scikit-misc"])
        subprocess.check_call([sys.executable, "-m", "pip", "install", "leidenalg"])

        # Remove installer if it wasn't available before for discretion
        if remove: os.remove("get-pip.py")
print("Package import success!")

"""Variables
These will be default numbers if the user does not change these inputs. Much of these are variable throughout experiments so the defaults will be basic at best.
"""
#These arguments will be dependent on user input.
#User will check visualizations and other experiment variables to decide input to optimize the experiment to their needs.
n=-1                    #Initialized later
cluster_res = 1.0
min_cells = 1
min_genes = 1
genes = []

#Cell ranger arguements: These are not optional
cellranger_script =""
string_id=""
reference_transcriptome=""
fastqs=""
sample=""
expected_cells=0
local_cores=0
local_mem=0

#SCSA arguments
species = "Human"
tissue = "All"      #SCSA's default parameter

#Other arguments
skip_cellranger = False
file_path = ""
interrupt=True

"""Command line argument syntax"""
#Scanpy
batch_arg = "--batch"
n_arg = "--neighbors"
resoution_arg = "--res"
cell_arg = "--min_cells"
gene_arg = "--min_genes"
display_arg = "--genes"

#Cellranger: These are not optional
cellranger_script_arg = "--script"
string_id_arg = "--id"
reference_transcriptome_arg = "--transcriptome"
fastqs_arg = "--fastqs"
sample_arg = "--sample"
expected_cells_arg = "--min_genes"
local_cores_arg = "--cores"
local_mem_arg = "--mem"

#SCSA
species_arg = "--species"
tissie_arg = "--tissue"

#Other
disable_interrupt_arg = "-disable_interrupts"
skip_cellranger_arg = "-skip_cellranger"
file_path_arg = "--filepath"        #Required if skipping cell ranger


"""Parse arguments"""
args = sys.argv[1:]
for arg in args:
    if batch_arg+"=" in arg:
        batch_effect = arg[arg.index("=") + 1:]
    elif n_arg+"=" in arg:
        n = int(arg[arg.index("=") + 1:])
    elif resoution_arg+"=" in arg:
        cluster_res = float(arg[arg.index("=") + 1:])
    elif cell_arg+"=" in arg:
        min_cells = float(arg[arg.index("=") + 1:])
    elif gene_arg+"=" in arg:
        min_genes = float(arg[arg.index("=") + 1:])

    #Cell ranger arguments
    elif cellranger_script_arg+"=" in arg:
        cellranger_script = arg[arg.index("=") + 1:]
    elif string_id_arg+"=" in arg:
        string_id = arg[arg.index("=") + 1:]
    elif reference_transcriptome_arg+"=" in arg:
        reference_transcriptome = arg[arg.index("=") + 1:]
    elif fastqs_arg+"=" in arg:
        fastqs = arg[arg.index("=") + 1:]
    elif sample_arg+"=" in arg:
        sample = arg[arg.index("=") + 1:]
    elif expected_cells_arg+"=" in arg:
        expected_cells = float(arg[arg.index("=") + 1:])
    elif local_cores_arg+"=" in arg:
        local_cores = float(arg[arg.index("=") + 1:])
    elif local_mem_arg+"=" in arg:
        local_mem = float(arg[arg.index("=") + 1:])

    #SCSA
    elif species_arg+"=" in arg:
        species = arg[arg.index("=") + 1:]
    elif tissie_arg+"=" in arg:
        tissue = arg[arg.index("=") + 1:]

    #Other
    elif file_path_arg in arg:
        file_path = arg[arg.index("=") + 1:]
        if file_path[-1] != "/" or file_path[-1] != "\\": file_path+="/"
    elif disable_interrupt_arg in arg:
        interrupt = False
    elif skip_cellranger_arg in arg:
        skip_cellranger = True
    elif display_arg+"=" in arg:                 #Genes listed must be comma separated
        genes = arg[arg.index("=") + 1:].split(",")
    else:
        sys.exit(arg+" is not a valid argument!")

"""Run cell ranger
Run the cell ranger using the passed arguments. Let program crash naturally if missing arguments or errors occur. 
Program continues only once process finishes.

Cellranger outputs a features.tsv file, a barcode file, and the matrix file. 
Features are the genes, barcode are the cells, and the matrix is the actual data. 
Anndata is a Scanpy object that can hold all these important variables and data.
"""

if (not skip_cellranger) and os.path.isfile(file_path+"matrix.mtx.gz") and os.path.isfile(file_path+"features.tsv.gz") and os.path.isfile(file_path+"barcodes.tsv.gz"):
    while True:
        run_prompt = input("Cell ranger has already be run on this file path. Would you like to run it again? [y/n]")
        if run_prompt == 'y' or run_prompt == 'Y':
            break
        elif run_prompt =='n' or run_prompt == 'N':
            skip_cellranger = True
            break
if not skip_cellranger:
    if not os.path.exists(string_id):
        os.mkdir(string_id)
    file_path = string_id + "/outs/filtered_features_bc_matrix/"        #Reset file_path if passed
    subprocess.run(["sbatch", cellranger_script, string_id, reference_transcriptome, fastqs, sample, expected_cells, local_cores, local_mem])



"""Load data"""
adata = sc.read_10x_mtx(file_path)

#Initialize n now that data length is defined
if n == -1:
  n = round(np.sqrt(adata.n_obs))

"""Overview Data"""
print("Number of cells: "+str(adata.n_obs))
print("Numebr of genes: "+str(adata.n_vars))

#Display a graph for the highest expressed genes
"""
In our test, MT genes seem to be highly expressed, and thus will be a cause of concern for our experiment. 
In standard sCRNA-seq experiments, the common thought is to remove MT genes as a sign of low quality cells due to cell perforation of cytoplasmic RNA loss. 
However for NASA's Genelab, spaceflight seems to affect mitochondrial gene function, according to our Genelab researcher Dr. Afshin Beheshti.
Therefore we have made the decision to not remove mitochondrial reads because that would filter our goal of analyzing biological processes from spaceflight.

**Inclusion of mitochondrial reads is vital as spaceflight has demonstrated its involvement in the following:**
- innate immunity 
- lipid metabolism
  => can contribute to greater risk of cardiovascular issues
- gene regulation 
- ETC, ATP synthesis (from increased radiation)
  - increased levels of oxidative damage and stress
  - muscle loss from metabolic flux changes 
"""
if interrupt:
    sc.pl.highest_expr_genes(adata)

"""Quality Control"""
#Basic Filtering
"""
Our first important step in our analysis will be to do some quality control. Our 2 most important steps of quality control are:
- Basic filtering 
- Removal of highly expressed genes

In a standard scRNAseq experiment, there will be 3 steps, including removing MT reads. 
But of course, of spaceflight experiments, we at Genelab have decided to not include this step in our analysis.
"""

#Calculate QC metrics and collect dataframe results for cell and gene
stats = sc.pp.calculate_qc_metrics(adata)
cell_qc_dataframe = stats[0]
gene_qc_dataframe = stats[1]

##Graph quality control graphs to find thresholds.

#Helper class that follows cursor and displays a vertical red line when viewing plot to help determine a threshold.
#This is implementation is a modified version of Matplotlib's open source code: https://matplotlib.org/stable/gallery/misc/cursor_demo.html
class VerticalCursor:
    #crosshair cursor.
    def __init__(self, ax, x_name):
        self.ax = ax
        self.vertical_line = ax.axvline(color='red', lw=0.8)
        # text location in axes coordinates
        self.text = ax.text(0.72, 0.9, '', transform=ax.transAxes)
        self.x_name = x_name

    def set_cross_hair_visible(self, visible):
        need_redraw = self.vertical_line.get_visible() != visible
        self.vertical_line.set_visible(visible)
        self.text.set_visible(visible)
        return need_redraw

    def on_mouse_move(self, event):
        if not event.inaxes:
            need_redraw = self.set_cross_hair_visible(False)
            if need_redraw:
                self.ax.figure.canvas.draw()
        else:
            self.set_cross_hair_visible(True)
            x, y = event.xdata, event.ydata
            # update the line position
            self.vertical_line.set_xdata(x)
            self.text.set_text(self.x_name+'=%1.2f' % x)
            self.ax.figure.canvas.draw()

"""Plot cell and  gene distribution
Looking at this graph, users must decide where to do the minimum cutoff for minimum genes. 
It is necessary to filter cells based on minimum genes because with that threshold, we can filter cells that have may have been contaminated.
It is important to do gene filtering after cell filtering because some genes may be detected only in low quality cells.
Users can zoom into graph by clicking onto the magnifying glass and drawing a rectangle where they wish to have a closer look. 
"""

#Cell distribution
if interrupt:
    fig, ax = plt.subplots()
    cursor = VerticalCursor(ax, "min_genes")
    fig.canvas.mpl_connect('motion_notify_event', cursor.on_mouse_move)

    plt.hist(cell_qc_dataframe['n_genes_by_counts'], bins="sqrt")
    plt.xlabel('N genes')
    plt.ylabel('N cells')
    plt.title("Use the graph to determine a threshold for filtering the dataset.")
    plt.show()

    #Collect user input for the minimum amount of genes to filter cells. Can be left blank if arguments already given or for default parameters.
    while True:
      try:
        g_value = input("Input min_genes threshold or leave blank to use \'"+str(min_genes)+"\': ")
        if g_value == "": break
        min_genes = int(g_value)
        break
      except:
        print("\nError: Please enter an integer value or leave blank!")

#Gene distribution
if interrupt:
    fig, ax = plt.subplots()
    cursor = VerticalCursor(ax, "min_cells")
    fig.canvas.mpl_connect('motion_notify_event', cursor.on_mouse_move)

    plt.hist(gene_qc_dataframe['n_cells_by_counts'], bins="sqrt")
    plt.xlabel('N cells expressing > 0')
    plt.ylabel('log(N genes)')
    plt.yscale('log')
    plt.title("Use the graph to determine a threshold for filtering the dataset.")
    plt.show()

    #Collect user input for the minimum amount of cells to filter genes. Can be left blank if arguments already given or for default parameters.
    while True:
      try:
        c_value = input("Input min_cells threshold or leave blank to use \'"+str(min_cells)+"\': ")
        if c_value == "": break
        min_cells = int(c_value)
        break
      except:
        print("\nError: Please enter an integer value or leave blank!")

"""Filtering
For our defaults, the minimum number of genes for filtering cells is a set number of 1.
The minimum number of cells for filtering genes is instead going to be percentile based.
"""
unfiltered_genes = adata.var_names
cell_filter_percentile = 0.01

if min_cells == 1 and min_genes == 1:             #Percentile-based filtering (default)
  print("Filtering using default settings.")
  stats = sc.pp.calculate_qc_metrics(adata)
  gene_counts_mean = len(stats[0]['n_genes_by_counts'])
  min_genes = round(gene_counts_mean*cell_filter_percentile)

print("Filtering using min_genes="+str(int(min_genes))+" and min_cells="+str(int(min_cells))+".")
sc.pp.filter_cells(adata, min_genes = min_genes)
sc.pp.filter_genes(adata, min_cells = min_cells)

#Display filtering results
filtered_genes = np.setdiff1d(unfiltered_genes, adata.var_names)    #This will also be used later
print(str(len(filtered_genes))+" genes filtered")

"""Normalize and Logarithmize Data
Normalization will help to preserve biological heterogeneity without the influence of any technical noise like sequencing depth and gene abundance.
Note: target_sum = 1e6 here refers to counts per million. 
Although different methods of normalization does exist, that are more accurate and have better performance, the CPM method is more flexible and 
scalable to all datasets and pipelines making it the best choice for our Genelab pipeline.

sc.pp.log1p helps to logarithmize the data to improve data "symmetry" on a linear scale for more relevant and accurate data. 
For further information on this topic, feel free to check out the following link: https://blog.qbaseplus.com/seven-tips-for-bio-statistical-analysis-of-gene-expression-data
"""
sc.pp.normalize_total(adata,target_sum=1e6)
sc.pp.log1p(adata)
print("Data normalized and logarithmized.")

"""Remove highly variable genes
Determine and remove highly variable genes based on each genes' mean and variance. Removing these genes eliminates noise due to high variability in the data.

Algorithm Description: Each gene is put into 20 'bins' based and their mean and variance. Each gene is then normalized based on the other genes in their bin. 
If a gene's normalized dispersion is greater or equal to a z-score of 2 (~98th percentile) AND the gene has a low mean cell count, it is marked highly variable.
(Decribed in the 'Identification of highly variable genes.' section of https://www.nature.com/articles/nbt.3192)
"""
sc.pp.highly_variable_genes(adata, flavor='seurat', min_disp=2)
highly_variable = adata[:, adata.var.highly_variable==True].var_names
print(str(len(adata.var[adata.var['highly_variable']==True]))+"/"+str(adata.n_vars)+" genes are highly variable and removed.")
adata = adata[:, adata.var.highly_variable==False]

"""### K-Nearest Neighbors
Calculate the distance between each cell using the KNN algorithm. Distances will be used to cluster the cells in the next step.
Note: warning just means it will proceed to automatically calculate PCA since it was not done beforehand.

Default (n=sqrt(adata.n_obs)): setting n as the square root of the length of the data is the general consensus if n is not provided.
-Sources:
https://towardsdatascience.com/how-to-find-the-optimal-value-of-k-in-knn-35d936e554eb
https://discuss.analyticsvidhya.com/t/how-to-choose-the-value-of-k-in-knn-algorithm/2606/7
https://stackoverflow.com/questions/11568897/value-of-k-in-k-nearest-neighbor-algorithm
"""
print("\n*Do not be concerned about the following warning.")
sc.pp.neighbors(adata, n_neighbors=n)

"""Cluster
Cluster/Group each cell based on the distances calculated in the previous step using the Leiden algorithm.
'resolution' determines the amount of clusters that will be formed (default: 1.0. The higher the resolution, the more clusters in the result)
"""
sc.tl.leiden(adata, resolution=cluster_res)

#Visualize Cluster Results
#Prepares the data to be visualized by simplifying multiple dimensions down to two dimensional coordinates using the UMAP algorithm.
#This algoithm also uses the distances calculated in 'neighbors()'.
sc.tl.umap(adata)

#Plot Clusters
#We will later label these clusters with cell types using marker gene identification.
if interrupt:
    sc.pl.umap(adata, color=['leiden'])

#Repeat the previous steps and allow user to adjust cluster resolution until satisifed.
while interrupt:
    try:
        prompt = input("Enter a decimal to change the cluster resolution (1.0 is the default) or leave blank to keep the results: ")
        if prompt == "": break
        cluster_res = float(prompt)
        sc.tl.leiden(adata, resolution=cluster_res)
        sc.tl.umap(adata)
        sc.pl.umap(adata, color=['leiden'])
    except:
        print("\nYou must enter a decimal value!")
print("Using cluster resolution of: "+str(cluster_res))

"""Visualized data based on selected gene(s)
Color cells based on their expression of a specific gene. 
The first 2 genes will be displayed as a sample if argument not provided, then prompt user to choose which gene(s) to visualize until satisfied.
Users can also view a list of genes by a handful at a time so it doesn't overload the command output.
Users a exit this prompt at anytime to proceed to the next steps.
Note: purple = no expression
"""
#Initial display (argument passed genes or sample first 2 genes)
if interrupt:
    view_genes = []
    if len(genes) != 0:    #Display selected genes' expression
        for gene in genes:
            if gene in adata.var_names:
                view_genes.append(gene)
            elif gene in highly_variable:
                print(gene+" gene was highly variable and was filtered out.")
            elif gene in filtered_genes:
                print(gene+" gene had been filtered out.")
            else:
                print(gene+" gene does not exist in the dataset.")
        sc.pl.umap(adata, color=view_genes)
    else:                   #Display the first 2 gene as a sample
        view_genes = []
        title = []
        i = 2
        for gene in adata.var_names:
            view_genes.append(gene)
            title.append("Sample: "+gene)
            i -= 1
            if i == 0: break
        sc.pl.umap(adata, color=view_genes, title=title)

#Additional gene visualization
PAGE_ROWS = 6
PAGE_COLS = 4
gene_list = adata.var_names.sort_values()
while True and interrupt:
    prompt = input("Enter the gene(s) you wish to visualize separated by commas (,). Enter \"-list page_#\" to view a list of genes or \"-exit\" to proceed with the program: ")
    if "-exit" == prompt[0:5] or prompt == "":
        break
    elif "-list" == prompt[0:5]:
        try: page_num = int(prompt[6:])
        except: page_num = 0

        page = page_num*PAGE_ROWS*PAGE_COLS
        content = ""
        newline = 0
        for i in range(PAGE_ROWS*PAGE_COLS):
            try:
                content += gene_list[page+i]+",\t"
                newline+=1
                if newline == PAGE_COLS:
                  content=content[:-2]+"\n"
                  newline = 0
            except IndexError:
                break
        print("Enter '-list #' to view other pages on the list.\n"+content)
    else:
        prompted_genes = "".join(prompt.split()).split(",")
        view_genes = []
        for gene in prompted_genes:
            if gene in adata.var_names:
                view_genes.append(gene)
            elif gene in highly_variable:
                print("\'"+gene+"\' gene was highly variable and was filtered out.")
            elif gene in filtered_genes:
                print("\'"+gene+"\' gene had been filtered out.")
            else:
                print("\'"+gene+"\' gene does not exist in the dataset.")
        sc.pl.umap(adata, color=view_genes)

"""Annotate data"""
#Prepare data to be annotated on SCSA
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
dat = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'logfoldchanges','scores','pvals']})
dat.to_csv("adata.csv")

print("WARNING: the following can only annotate human or mice cells!")

#Prompt for species argument
if interrupt:
    species_prompt = input("Input species of the dataset or leave blank to use \'"+species+"\': ")
    if not species_prompt == "":
        species = species_prompt

#Prompt for tissue argument
    tissue_prompt = input("Input tissue of the dataset or leave blank to use \'"+tissue+"\': ")
    if not tissue_prompt == "":
        tissue = tissue_prompt

#Run SCSA
print("Annotating with species="+species+" and tissue="+tissue)
subprocess.check_call([sys.executable, "SCSA.py", "-i", "adata.csv", "-o", "output", "-s", "scanpy", "-E", "-g", species, "-p", "1", "-f", "1", "-k", tissue])

"""Revisit cluster resolution"""
while interrupt:
    print("\n*View the output.xlsx file. \n*If you are not satisfied with the results, use this opportunity to adjust the cluster resolution and re-annotate the data.")
    while interrupt:
        try:
            prompt = input("Enter a decimal to change the cluster resolution. Leave blank to keep the results. Enter 'exit' to exit: ")
            if prompt == "" or prompt == "exit": break
            cluster_res = float(prompt)
            sc.tl.leiden(adata, resolution=cluster_res)
            sc.tl.umap(adata)
            sc.pl.umap(adata, color=['leiden'])
        except:
            print("\nYou must enter a decimal value!")
    if prompt == "exit": break

    print("Using cluster resolution of: "+str(cluster_res))
    sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    dat = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'logfoldchanges', 'scores', 'pvals']})
    dat.to_csv("adata.csv")

    subprocess.check_call([sys.executable, "SCSA.py", "-i", "adata.csv", "-o", "output", "-s", "scanpy", "-E", "-g", species, "-p", "1", "-f", "1", "-k", tissue])

"""TO BE REMOVED? Reload output into program"""
df = pd.read_excel("output.xlsx", None)