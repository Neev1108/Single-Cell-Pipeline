#Import or Install Libraries
while True:
    try:
        print("Loading packages...")
        import sys
        import os
        import subprocess

        import scanpy as sc
        import scorect_api as ct
        import openpyxl
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
        subprocess.check_call([sys.executable, "-m", "pip", "install", "scanpy"])
        subprocess.check_call([sys.executable, "-m", "pip", "install", "openpyxl"])

        if not os.path.exists("scorect_api.py"):
            urllib.request.urlretrieve("https://raw.githubusercontent.com/LucasESBS/scoreCT/master/src/scorect_api.py", "scorect_api.py")
            # Required by scorect_api
            subprocess.check_call([sys.executable, "-m", "pip", "install", "requests"])

        # Remove installer if it wasn't available before for discretion
        if remove: os.remove("get-pip.py")
print("Package import success!")

"""Variables
These will be default numbers if the user does not change these inputs. Much of these are variable throughout experiments so the defaults will be basic at best.
"""
#These arguments will be dependent on user input.
#User will check visualizations and other experiment variables to decide input to optimize the experiment to their needs.

species = ""
tissue = "\"\""
K = 450
m = 20

#Not required, but highly recommended
adata_path = "adata.h5ad"
markers_path = ""

#Other
interrupt=True

"""Command line argument syntax"""
markers_path_arg = "--markers"
adata_path_arg = "--adata"
species_arg = "--species"
tissie_arg = "--tissue"
K_arg = "--K"
m_arg = "--bins"
disable_interrupt_arg = "-disable_interrupts"

"""Parse arguments"""
args = sys.argv[1:]
for arg in args:
    #Files
    if markers_path_arg in arg:
        markers_path = arg[arg.index("=") + 1:]
    elif adata_path_arg in arg:
        adata_path = arg[arg.index("=") + 1:]

    #Parameters
    elif species_arg+"=" in arg:
        species = arg[arg.index("=") + 1:]
    elif tissie_arg+"=" in arg:
        tissue = arg[arg.index("=") + 1:]
    elif K_arg+"=" in arg:
        K = float(arg[arg.index("=") + 1:])
    elif m_arg+"=" in arg:
        m = float(arg[arg.index("=") + 1:])

    # Other
    elif disable_interrupt_arg in arg:
        interrupt = False
    else:
        sys.exit(arg + " is not a valid argument!")

"""Load adata"""
adata = sc.read(adata_path)
print("Successfully loaded adata!")

"""Annotate data
Use ScoreCT to annotate data. Marker gene file must be provided for the package to work.
If file not provided to this program, pre-loaded marker files will be used (in 'marker genes' folder).
If annotating Human or Mouse, tissue name is required to collect data from Cell Marker.
Using 'adata' and the marker gene file, annotate data and add results onto 'adata.obs', then visualize
results as a UMAP graph. Finally, output adata.obs as a excel to externally view cell data stats and
classification.

Marker gene file format example (csv):
Cell, Cell, Cell
gene, gene, gene
gene, gene, gene
gene,     , gene
etc.,      ,
"""
#Load marker file

#Allow users to use own marker gene file
marker_loaded = False
if not markers_path == "":
    try:
        ref_marker = ct.read_markers_from_file(markers_path)
        print("Using marker file: \""+markers_path+"\"")
        marker_loaded = True
    except:
        print("Invalid marker file. Try again using annotate.py or continue to use provided marker gene data.")

#Otherwise use default marker data
if not marker_loaded:
    # Prompt for species argument
    while species == "":
        species_prompt = input("Input species of the dataset: ")
        if not species_prompt == "":
            species = species_prompt

    #TODO: Use 'species' to collect marker gene file
    species = species[0].upper() + species[1:].lower()
    if species == "Drosophila melanogaster" or species == "Fruitfly" or species == "D_melanogaster":
        markers_path = 'marker genes/D_melanogaster_genes.csv'
    elif species == "Mouse-ear cress" or species == "Mouse ear cress" or species == "Thale cress" or species == "Arabidopsis thaliana" or species == "A_thaliana":
        markers_path = 'marker genes/A_thaliana_genes.csv'
    elif species == "Dario rerio" or species == "Zebrafish":
        markers_path = 'marker genes/Dario_rerio_genes.csv'
    else:
        # Prompt for tissue argument
        while tissue == "\"\"" or tissue == "":
            tissue_prompt = input("Input tissue of the dataset: ")
            if not tissue_prompt == "":
                tissue = tissue_prompt
        tissue = tissue[0].upper() + tissue[1:].lower()

        print("Retrieving data...")
        ref_marker = ct.get_markers_from_db(species, tissue)
        print("Using species="+species+" and tissue="+tissue)
        marker_loaded = True

    if not marker_loaded:
        print("Using species="+species)
        ref_marker = ct.read_markers_from_file(markers_path)

#Repeated prompt for K, m, and cluster_res; recalculate annotation, and display results until user is satisfied.
cluster_res = 0.5
while True:
    #Prompt for cluster resolution
    while interrupt:
        try:
            prompt = input("Enter a decimal to change the cluster resolution (current: " + str(cluster_res) + ") or leave blank to keep the results: ")
            if prompt == "": break
            cluster_res = float(prompt)
            sc.tl.leiden(adata, resolution=cluster_res)
            sc.tl.umap(adata)
            sc.pl.umap(adata, color=['leiden'])
        except:
            print("\nYou must enter a decimal value!")
    print("Using cluster resolution of: " + str(cluster_res))

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
    if interrupt:
        sc.pl.umap(adata, color=['cell_type'], title=['Cell Type Annotation for '+species+" "+tissue])
    else:
        break

    #Repeat the previous steps and allow user to adjust K and m until satisifed.
    while interrupt:
        try:
            prompt = input("Enter a value for the top K genes to include in scoring (current: " + str(K) + ") or leave blank to keep the results." +
                           "\nEnter \"cluster\" to re-adjust cluster resolution (current: " + str(cluster_res) + ") or \"exit\" to exit: ")
            if prompt == "cluster" or prompt == "exit": break
            elif not prompt == "": K = int(prompt)

            prompt = input("Enter a value for the m number of bins to use to divide the gene ranking (current: " + str(m) + ") or leave blank to keep the results: ")
            if not prompt == "": m = int(prompt)

            ct_pval, ct_score = ct.celltype_scores(nb_bins=m,
                                                   ranked_genes=marker_df,
                                                   K_top=K,
                                                   marker_ref=ref_marker,
                                                   background_genes=background)
            adata.obs['cell_type'] = ct.assign_celltypes(cluster_assignment=adata.obs['leiden'], ct_pval_df=ct_pval, ct_score_df=ct_score)
            sc.pl.umap(adata, color=['cell_type'], title=['Cell Type Annotation for ' + species + " " + tissue])
        except:
            print("\nYou must enter integer values!")
    if prompt == "exit": break

#Export results as an excel
adata.obs = adata.obs.rename(columns={"leiden":"cluster"})  #Rename to avoid confusion
output_name = species+" "
if not tissue == "\"\"":
    output_name += tissue +" "
adata.obs.to_excel(output_name+'annotation.xlsx')

#Message about other python script
print("\nAnnotation exported to \'"+output_name+'annotation.xlsx\'')