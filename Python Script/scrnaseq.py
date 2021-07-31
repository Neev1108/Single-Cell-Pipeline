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
        import subprocess

        import numpy as np
        import pandas as pd
        import scanpy as sc
        import seaborn as sb
        import matplotlib.pyplot as plt
        import scorect_api as ct
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

        if not os.path.exists("scorect_api.py"):
            urllib.request.urlretrieve("https://raw.githubusercontent.com/LucasESBS/scoreCT/master/src/scorect_api.py", "scorect_api.py")
            # Required by scorect_api
            subprocess.check_call([sys.executable, "-m", "pip", "install", "requests"])

        # Packages used but not included in scanpy package
        subprocess.check_call([sys.executable, "-m", "pip", "install", "harmonypy"])
        subprocess.check_call([sys.executable, "-m", "pip", "install", "scikit-misc"])
        subprocess.check_call([sys.executable, "-m", "pip", "install", "leidenalg"])

        # Remove installer if it wasn't available before for discretion
        if remove: os.remove("get-pip.py")
print("Package import success!")
input()

"""Variables
These will be default numbers if the user does not change these inputs. Much of these are variable throughout experiments so the defaults will be basic at best.
"""
#These arguments will be dependent on user input.
#User will check visualizations and other experiment variables to decide input to optimize the experiment to their needs.
n=-1                    #Initialized later
cluster_res = 1.0
min_cells = -1
min_genes = -1
genes = []

#Annotation arguments
species = ""
tissue = ""
K = 350
m = 15

#Other arguments
file_path = ""
markers_path = ""
interrupt=True

"""Command line argument syntax"""
#Scanpy
n_arg = "--neighbors"
resoution_arg = "--res"
cell_arg = "--min_cells"
gene_arg = "--min_genes"
display_arg = "--genes"

#Annotation
species_arg = "--species"
tissie_arg = "--tissue"
K_arg = "--K"
m_arg = "--bins"

#Other
file_path_arg = "--filepath"
markers_path_arg = "--markers"
disable_interrupt_arg = "-disable_interrupts"


"""Parse arguments"""
args = sys.argv[1:]
for arg in args:
    #Required
    if file_path_arg in arg:
        file_path = arg[arg.index("=") + 1:]
        if file_path[-1] != "/" or file_path[-1] != "\\": file_path += "/"

    #ScanPy arguments
    elif n_arg+"=" in arg:
        n = int(arg[arg.index("=") + 1:])
    elif resoution_arg+"=" in arg:
        cluster_res = float(arg[arg.index("=") + 1:])
    elif cell_arg+"=" in arg:
        min_cells = float(arg[arg.index("=") + 1:])
    elif gene_arg+"=" in arg:
        min_genes = float(arg[arg.index("=") + 1:])
    elif display_arg+"=" in arg:                 #Genes listed must be comma separated
        genes = arg[arg.index("=") + 1:].split(",")

    #Annotation
    elif species_arg+"=" in arg:
        species = arg[arg.index("=") + 1:]
    elif tissie_arg+"=" in arg:
        tissue = arg[arg.index("=") + 1:]
    elif K_arg+"=" in arg:
        K = float(arg[arg.index("=") + 1:])
    elif m_arg+"=" in arg:
        m = float(arg[arg.index("=") + 1:])
    if markers_path_arg in arg:
        markers_path = arg[arg.index("=") + 1:]
        if markers_path[-1] != "/" or markers_path[-1] != "\\": markers_path += "/"

    #Other
    elif disable_interrupt_arg in arg:
        interrupt = False
    else:
        sys.exit(arg+" is not a valid argument!")

"""Load data"""
if file_path == "":             #File path is required!
    raise FileExistsError("You must input a file path of the cell ranger output!")

#If filepath points directly to data or to cellranger output
if os.path.isfile(file_path+"matrix.mtx.gz") and os.path.isfile(file_path+"features.tsv.gz") and os.path.isfile(file_path+"barcodes.tsv.gz"):
    adata = sc.read_10x_mtx(file_path)
else:
    if "/outs" in file_path:                                                                         #File path points to other paths within cellranger output
        file_path = file_path[:file_path.index("outs")]
    if not os.path.exists(file_path+"outs/filtered_feature_bc_matrix/"):                             #Corrupted cellranger output (that isn't pointed to data)
        raise FileExistsError(file_path+"outs/filtered_feature_bc_matrix/ is missing in cell ranger output!")
    if not (os.path.isfile(file_path+"outs/filtered_feature_bc_matrix/matrix.mtx.gz") and os.path.isfile(file_path+"outs/filtered_feature_bc_matrix/features.tsv.gz") and os.path.isfile(file_path+"outs/filtered_feature_bc_matrix/barcodes.tsv.gz")):
        raise FileNotFoundError("One or more of the following files in "+file_path+"outs/filtered_feature_bc_matrix/ is missing: matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz")   #Missing files in filepath
    adata = sc.read_10x_mtx(file_path+"outs/filtered_feature_bc_matrix/")

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
        min_gene_text = min_genes
        if min_gene_text == -1:
            min_gene_text = "default settings"
        g_value = input("Input min_genes threshold or leave blank to use \'"+str(min_gene_text)+"\': ")
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
        min_cell_text = min_cells
        if min_cell_text == -1:
            min_cell_text = "default settings"
        c_value = input("Input min_cells threshold or leave blank to use \'"+str(min_cell_text)+"\': ")
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

prefiltered_length = adata.n_vars
if min_cells == -1 and min_genes == -1:             #Percentile-based filtering (default)
  print("Filtering using default settings.")
  stats = sc.pp.calculate_qc_metrics(adata)
  gene_counts_mean = len(stats[0]['n_genes_by_counts'])
  min_genes = round(gene_counts_mean*cell_filter_percentile)
  min_cells = 1

print("Filtering using min_genes="+str(int(min_genes))+" and min_cells="+str(int(min_cells))+".")
sc.pp.filter_cells(adata, min_genes = min_genes)
sc.pp.filter_genes(adata, min_cells = min_cells)

#Display filtering results
print(str(adata.n_vars)+"/"+str(prefiltered_length)+" genes filtered")

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
#Initialize n now that filtering has finished
if n == -1:
  n = round(np.sqrt(adata.n_obs))

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
        prompt = input("Enter a decimal to change the cluster resolution (current: "+str(cluster_res)+") or leave blank to keep the results: ")
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
            else:
                print("\'"+gene+"\' gene is not in the dataset.")
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
            else:
                print("\'"+gene+"\' gene is not in the dataset.")
        sc.pl.umap(adata, color=view_genes)

"""Annotate data"""
#Prepare data to be annotated on SCSA
#sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
#result = adata.uns['rank_genes_groups']
#groups = result['names'].dtype.names
#dat = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'logfoldchanges','scores','pvals']})
#dat.to_csv("adata.csv")

#print("WARNING: the following can only annotate human or mice cells!")

#Prompt for species argument
#if interrupt:
#    species_prompt = input("Input species of the dataset or leave blank to use \'"+species+"\': ")
#    if not species_prompt == "":
#        species = species_prompt

#Prompt for tissue argument
#if interrupt:
#    tissue_prompt = input("Input tissue of the dataset or leave blank to use \'"+tissue+"\': ")
#    if not tissue_prompt == "":
#        tissue = tissue_prompt

#Run SCSA
#print("Annotating with species="+species+" and tissue="+tissue)
#subprocess.check_call([sys.executable, "SCSA.py", "-i", "adata.csv", "-o", "output", "-s", "scanpy", "-E", "-g", species, "-p", "1", "-f", "1", "-k", tissue])

"""Revisit cluster resolution"""
#while interrupt:
#    print("\n*View the output.xlsx file. \n*If you are not satisfied with the results, use this opportunity to adjust the cluster resolution and re-annotate the data.")
#    while interrupt:
#        try:
#            prompt = input("Enter a decimal to change the cluster resolution. Leave blank to keep the resolution. Enter 'exit' to exit: ")
#            if prompt == "" or prompt == "exit": break
#            cluster_res = float(prompt)
#            sc.tl.leiden(adata, resolution=cluster_res)
#            sc.tl.umap(adata)
#            sc.pl.umap(adata, color=['leiden'])
#        except:
#            print("\nYou must enter a decimal value!")
#    if prompt == "exit": break

#    print("Using cluster resolution of: "+str(cluster_res))
#    sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
#    result = adata.uns['rank_genes_groups']
#    groups = result['names'].dtype.names
#    dat = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'logfoldchanges', 'scores', 'pvals']})
#   dat.to_csv("adata.csv")
#
#    subprocess.check_call([sys.executable, "SCSA.py", "-i", "adata.csv", "-o", "output", "-s", "scanpy", "-E", "-g", species, "-p", "1", "-f", "1", "-k", tissue])

"""TO BE REMOVED? Reload output into program"""
#df = pd.read_excel("output.xlsx", None)

"""Export adata to use in annotate.py"""
adata.write("adata.h5ad")
print("Exported data to 'adata.h5ad'. The following steps can be repeated using [SCRIPT NAME] and this program can safely exited.")

"""Annotate data
TODO DESCRIPTION
"""
#Load marker file
#Allow users to use own marker gene file
marker_loaded = False
if not markers_path == "":
    try:
        ref_marker = ct.read_markers_from_file(markers_path)
        marker_loaded = True
    except:
        print("Invalid marker file. Try again using annotation.py or continue to use default marker data.")
#Otherwise Uuse default marker data
if not marker_loaded:
    # Prompt for species argument
    while species == "":
        species_prompt = input("Input species of the dataset: ")
        if not species_prompt == "":
            species = species_prompt

    #Prompt for tissue argument
    while tissue == "":
        tissue_prompt = input("Input tissue of the dataset: ")
        if not tissue_prompt == "":
            tissue = tissue_prompt

    #TODO: Use 'species' to collect marker gene file
    if species == "asdlkasjd":
        pass
        remove = 'D_melanogaster_genes_corrected.csv'
    else:
        ref_marker = get_markers_from_db(species, tissue)
    print("Using species="+species+" and tissue="+tissue)

ref_marker = ct.read_markers_from_file('D_melanogaster_genes_corrected.csv')

#Calculate statistics
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
marker_df = ct.wrangle_ranks_from_anndata(adata)

#Calculate p-value and scores needed for 'assign_celltypes'
background = adata.var.index.tolist()
ct_pval, ct_score = ct.celltype_scores(nb_bins=m,
                                        ranked_genes=marker_df,
                                        K_top = K,
                                        marker_ref=ref_marker,
                                        background_genes=background)

#Annotate
adata.obs['cell_type'] = ct.assign_celltypes(cluster_assignment=adata.obs['leiden'], ct_pval_df=ct_pval, ct_score_df=ct_score)

#Visualize results
sc.pl.umap(adata, color=['cell_type'], title=['Cell Type Annotation for '+species+" "+tissue])

#Export results as an excel
adata.obs = adata.obs.rename(columns={"leiden":"cluster"})  #Rename to avoid confusion
adata.obs.to_excel(species+' '+tissue+' annotation.xlsx')

#Message about other python script
print("If you would like to experiment with different parameters regarding only annotation, please use 'annotate.py' and input the post-filtered 'adata.h5ad' file.")