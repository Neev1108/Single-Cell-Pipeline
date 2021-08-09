"""Helper Functions
Created helper functionss to minimize padding from log related code within the program."""
#Flush file on write in case of forced program exit
def log(msg):
    file.write(msg+"\n")
    file.flush()

#Save figure and V&V
def save_figure(image_name):
    plt.savefig(image_directory + image_name)
    if interrupt: plt.show()
    VV_file_save(image_directory,image_name,"Saved \"" + image_name + "\"")
    plt.close()

#Validate and verify flies were successfully saved
def VV_file_save(dir, filename, success_msg):
    if not os.path.isfile(dir + filename):
        print("ERROR: Could not save \""+filename+"\" file!")
        log("ERROR: Could not save \""+filename+"\" file!")
    else:
        log(success_msg)

#Validate and Verify: Count amount of NA's that were generated in the classification output. Warn if there are too many NA's.
def VnV_classification():
    count = 0
    for classification in adata.obs['cell_type']:
        if classification.lower() == "na" or pd.isnull(classification):
            count+=1
    print(str(count)+"/"+str(adata.n_obs)+" cells were not able to be classified. Try adjusting the parameters.")

    assert not count == adata.n_obs, "None of the cells were able to be classified! Try using a different marker gene file."
    if (count/adata.n_obs) > annotate_limit:
        print("WARNING: More than " + str(annotate_limit * 100) + "% of cells were NOT able to be classified!")

"""Setup Log File
Log major outputs and actions for validation and verification of a successfully run program."""
import os
import traceback

#Directories
log_dir = "logs/"
image_directory = "images/"
marker_dir = "marker genes/"
annotation_dir = "annotations/"

run_num = 1
while os.path.exists(log_dir+"annotate_run"+str(run_num)+".txt") or os.path.exists(image_directory+"annotate_run"+str(run_num)+".png"):
    run_num += 1
file = open(log_dir+"annotate_run"+str(run_num)+".txt", "w")
print("Created log file \"annotate_run"+str(run_num)+".txt\" in "+log_dir)

#Set up Image Directory
image_directory = image_directory + "annotate_run" + str(run_num) + "/"
os.mkdir(image_directory)
print("Images will be automatically saved at \"" + image_directory + "\"")
log("Set up \"images/annotate run" + str(run_num) + "/\" directory")

"""***Program starts here***
Try-catch statement used to log exceptions."""
try:
    """Import or Install Libraries"""
    while True:
        try:
            print("Loading packages...")
            import sys
            import os
            import subprocess

            import pandas as pd
            import matplotlib.pyplot as plt
            import scanpy as sc
            import scorect_api as ct
            import openpyxl
            break
        except:
            print("***Some packages have not been installed. Installing now...***")
            log("Some packages have not been installed. Attempting to install...")

            # Retrieve installer if not available
            import urllib.request
            remove = False
            if not os.path.exists("get-pip.py"):
                urllib.request.urlretrieve("https://bootstrap.pypa.io/get-pip.py", "get-pip.py")
                remove = True
            subprocess.check_call([sys.executable, "get-pip.py"])

            # Download and install packages if not installed
            subprocess.check_call([sys.executable, "-m", "pip", "install", "pandas"])
            subprocess.check_call([sys.executable, "-m", "pip", "install", "matplotlib"])
            subprocess.check_call([sys.executable, "-m", "pip", "install", "scanpy"])
            subprocess.check_call([sys.executable, "-m", "pip", "install", "openpyxl"])

            if not os.path.exists("scorect_api.py"):
                urllib.request.urlretrieve("https://raw.githubusercontent.com/LucasESBS/scoreCT/master/src/scorect_api.py", "scorect_api.py")
                # Required by scorect_api
                subprocess.check_call([sys.executable, "-m", "pip", "install", "requests"])

            # Remove installer if it wasn't available before for discretion
            if remove: os.remove("get-pip.py")
            log("Succesfully installed packages")
    print("Package import success!")
    log("Successfully loaded packages")

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
    annotate_limit = 0.3

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
    log("Successfully parsed command arguments")

    """Load adata"""
    #Check validity of  filepath
    if not os.path.isfile(adata_path):  # Corrupted cellranger output (that isn't pointed to data)
        raise FileExistsError(adata_path + " is not a valid file path!")

    #Load file
    adata = sc.read(adata_path)

    #V&V sequences are in correct format for classification
    for i in range(1000):
        try:
            entry = adata.obs_names[i]
            assert entry[:-2].isalpha() and entry[-2:] == "-1", "Entries in barcodes.tsv are NOT in the correct format: " + entry
        except IndexError:
            break

    print("Successfully loaded adata!")
    log("Successfully loaded adata")

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
            print("Using marker file: \"" + markers_path + "\"")
            log("Loaded marker gene file \"" + markers_path + "\"")
            marker_loaded = True

            # Additional V&V (read_markers_from_file already checks for formating)
            assert ref_marker.empty, ""
        except:
            print("Failed to load marker gene file. Try again using annotate.py or continue to use provided marker gene data.")
            log("Failed to load marker gene file")

    # Otherwise use default marker data
    if not marker_loaded:
        # Prompt for species argument
        while species == "":
            species_prompt = input("Input species of the dataset: ")
            if not species_prompt == "":
                species = species_prompt

        # Use 'species' to collect marker gene file
        species = species[0].upper() + species[1:].lower()
        if species == "Drosophila melanogaster" or species == "Fruitfly" or species == "D_melanogaster":
            markers_path = marker_dir + 'D_melanogaster_genes.csv'
            assert os.path.isfile(markers_path), 'D_melanogaster_genes.csv is missing from ' + markers_path + "!"
        elif species == "Mouse-ear cress" or species == "Mouse ear cress" or species == "Thale cress" or species == "Arabidopsis thaliana" or species == "A_thaliana":
            markers_path = marker_dir + 'A_thaliana_genes.csv'
            assert os.path.isfile(markers_path), 'A_thaliana_genes.csv is missing from ' + markers_path + "!"
        elif species == "Dario rerio" or species == "Zebrafish":
            markers_path = marker_dir + 'Dario_rerio_genes.csv'
            assert os.path.isfile(markers_path), 'Dario_rerio_genes.csv is missing from ' + markers_path + "!"
        # Retrieve data from Cell Marker useing ct.get_markers_from_db(species, tissue)
        else:
            # Change keyword to match format if applies
            if species == "Homo sapian":
                species = "Human"
            elif species == "Mus musculus":
                species = "Mouse"

            # V&V: Assert if species not available
            assert species == "Mouse" or species == "Human", "Could not find information for " + species + ". Please provide a marker gene file for that species in \"annotate.py\"."

            # Prompt for tissue argument
            while tissue == "\"\"" or tissue == "":
                tissue_prompt = input("Input tissue of the dataset: ")
                if not tissue_prompt == "":
                    tissue = tissue_prompt
            tissue = tissue[0].upper() + tissue[1:].lower()

            # Retrieve data from Cell Marker
            print("Retrieving marker gene data...")
            ref_marker = ct.get_markers_from_db(species, tissue)
            print("Using species=" + species + " and tissue=" + tissue)
            log("Using Cell Marker with species=" + species + " tissue=" + tissue + " for marker gene reference")
            marker_loaded = True

            # V&V: Assert if empty (due to tissue)
            # if species == "Human" or species == "Mouse":
            assert not ref_marker.empty, "Could not find information on " + tissue + " for " + species + ". Please try again using a different keyword for \'tissue\' on \"annotate.py\"."

        if not marker_loaded:
            print("Using species=" + species)
            ref_marker = ct.read_markers_from_file(markers_path)
            log("Loaded marker gene file \"" + markers_path + "\"")

    #Calculate statistics
    sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
    log("Ranked gene groups")

    marker_df = ct.wrangle_ranks_from_anndata(adata)
    log("Wraggled ranks form adata")

    #Calculate p-value and scores needed for 'assign_celltypes'
    background = adata.var.index.tolist()
    ct_pval, ct_score = ct.celltype_scores(nb_bins=m,
                                            ranked_genes=marker_df,
                                            K_top = K,
                                            marker_ref=ref_marker,
                                            background_genes=background)
    log("Calculated cell type scores")

    #Annotate
    adata.obs['cell_type'] = ct.assign_celltypes(cluster_assignment=adata.obs['leiden'], ct_pval_df=ct_pval, ct_score_df=ct_score)
    VnV_classification()

    #Visualize results
    if interrupt:
        sc.pl.umap(adata, color=['cell_type'], title=['Cell Type Annotation for '+species+" "+tissue], show=False)
        save_figure("annotation_plot.png")

    # V&V cluster results
    if not '1' in adata.obs['leiden'].to_list():
        print("WARNING: Only one cluster was created! Please increase the cluster resolution!")

    #Repeated prompt for K, m, and cluster_res; recalculate annotation, and display results until user is satisfied.
    cluster_res = 0.5
    res_i = 1
    ann_i = 1
    while True:
        #Prompt for cluster resolution
        while interrupt or not ('1' in adata.obs['leiden'].to_list()):
            try:
                prompt = input("Enter a decimal to change the cluster resolution (current: " + str(cluster_res) + ") or leave blank to keep the results: ")
                if prompt == "":
                    if not '1' in adata.obs['leiden'].to_list():
                        print("WARNING: Only one cluster was created! Please increase the cluster resolution!")
                    else:
                        break
                cluster_res = float(prompt)
                sc.tl.leiden(adata, resolution=cluster_res)
                sc.pl.umap(adata, color=['leiden'], show=False)
                save_figure("cluster_("+str(res_i)+").png")
                res_i += 1
            except:
                if not prompt == "":  # Not a crash due to mono-cluster warning
                    print("You must enter a decimal value!")
        print("Using cluster resolution of: " + str(cluster_res))

        # Calculate statistics
        sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
        log("Ranked gene groups")

        marker_df = ct.wrangle_ranks_from_anndata(adata)
        log("Wraggled ranks form adata")

        background = adata.var.index.tolist()

        #Repeat the previous steps and allow user to adjust K and m until satisifed.
        while interrupt:
            try:
                prompt = input("Enter a value for the top K genes to include in scoring (current: " + str(K) + ") or leave blank to keep the results." +
                               "\nOr enter \"cluster\" to re-adjust cluster resolution or \"exit\" to exit: ")
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
                VnV_classification()
                sc.pl.umap(adata, color=['cell_type'], title=['Cell Type Annotation for ' + species + " " + tissue], show=False)
                save_figure("annotation_plot_("+str(ann_i)+").png")
                ann_i+=1
            except:
                print("\nYou must enter integer values!")
        if (not interrupt) or prompt == "exit": break

    #Export results as an excel
    adata.obs = adata.obs.rename(columns={"leiden":"cluster"})  #Rename to avoid confusion
    output_name = species+"_"
    if not tissue == "\"\"":
        output_name += tissue +"_"
    output_name = "annotate_run" + str(run_num) + "_" + output_name + "annotation.xlsx"
    adata.obs.to_excel(annotation_dir+output_name)
    VV_file_save(annotation_dir, output_name, "Exported \"" + output_name + "\"")

    #Message about other python script
    print("*Annotation exported to \'"+output_name+'\'*')

    log("*Program finished running without any errors*")
except Exception:
    traceback.print_exc()
    traceback.print_exc(file=file)
finally:
    file.close()