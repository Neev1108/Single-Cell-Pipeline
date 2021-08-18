"""Helper Functions
Created helper functionss to minimize padding from log related code within the program."""
#Flush file on write in case of forced program exit
def log(msg):
    file.write(msg+"\n")
    file.flush()

#Save figure and V&V
def save_figure(image_name, dir):
    plt.savefig(dir + image_name)
    if interrupt: plt.show()
    VV_file_save(dir,image_name,"Saved \"" + image_name + "\"")
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
    if count > 0:
        print(str(count)+"/"+str(adata.n_obs)+" cells were not able to be classified. Try adjusting the parameters.")
    else:
        print("All cells were classified!")

    assert not count == adata.n_obs, "None of the cells were able to be classified! Try using a different marker gene file and double check format."
    if (count/adata.n_obs) > annotate_limit:
        print("WARNING: More than " + str(annotate_limit * 100) + "% of cells were NOT able to be classified!")

#Validate and Verify: Marker gene file format
def VV_markers(ref):
    for col in ref.columns:
        if "Unnamed" in col:
            return False
    return True

#Plot UMAP or tSNE
def plot(adata, color, title=""):
    if tsne:
        if title == "":
            sc.pl.tsne(adata, color=color, show=False)
        else:
            sc.pl.tsne(adata, color=color, title=title, show=False)
    else:
        if title == "":
            sc.pl.umap(adata, color=color, show=False)
        else:
            sc.pl.umap(adata, color=color, title=title, show=False)

def plot_annotation(adata, color, title):
    sc.set_figure_params(figsize=(10, 10))
    if tsne:
        sc.pl.tsne(adata, color=color, title=title, legend_loc='on data', show=False)
    else:
        sc.pl.umap(adata, color=color, title=title, legend_loc='on data', show=False)

"""Setup Log File
Log major outputs and actions for validation and verification of a successfully run program."""
import os
import traceback

# Directories
log_dir = "../../02-Scanpy/VV_logs/"
image_directory = "../../02-Scanpy/Images/"
marker_dir = "../../03-ScoreCT/Marker_genes/"
annotation_dir = "../../03-ScoreCT/Annotation_Exports/"
final_plot_dir = "../../03-ScoreCT/tSNE_UMAP/"

run_num = 1
while os.path.exists(log_dir+"annotate_run"+str(run_num)+".txt") or os.path.exists(image_directory+"annotate_run"+str(run_num)+".png"):
    run_num += 1
file = open(log_dir+"annotate_run"+str(run_num)+".txt", "w")
print("Created log file \"annotate_run"+str(run_num)+".txt\" in "+log_dir)

"""***Program starts here***
Try-catch statement used to log exceptions."""
try:
    """Import or Install Libraries"""
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

        # Try importing packages again. If fail, proceed with normal crash protocol
        print("Attempting to loading packages again. If program crashes, please install package on your own or re-run program.")
        import sys
        import os
        import subprocess

        import pandas as pd
        import matplotlib.pyplot as plt
        import scanpy as sc
        import scorect_api as ct
        import openpyxl
    print("Package import success!")
    log("Successfully loaded packages")

    """Variables
    These will be default numbers if the user does not change these inputs. Much of these are variable throughout experiments so the defaults will be basic at best.
    """
    # These arguments will be dependent on user input.
    # User will check visualizations and other experiment variables to decide input to optimize the experiment to their needs.

    species = ""
    tissue = "tissue"
    K = 300
    m = 15

    # Not required, but highly recommended
    adata_path = "adata.h5ad"
    markers_path = ""
    annotate_limit = 0.3

    # Other
    interrupt = True
    tsne = False

    """Command line argument syntax"""
    markers_path_arg = "--markers"
    adata_path_arg = "--adata"
    species_arg = "--species"
    tissie_arg = "--tissue"
    K_arg = "--K"
    m_arg = "--bins"
    plot_arg = "--plot"
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
        elif plot_arg+"=" in arg:
            plot_input = arg[arg.index("=")+1:].lower()
            if plot_input == "tsne" or plot_input == "t-sne": tsne = True
        else:
            sys.exit(arg + " is not a valid argument!")
    log("Successfully parsed command arguments")

    plot_log = "umap"
    if tsne: plot_log = "tsne"
    log("\n***Parameters***\n"
        "adata path: " + adata_path + "\n"
        "marker gene file: " + markers_path + "\n"
        "species: " + species + "\n"
        "tissue: " + tissue + "\n"
        "K: " + str(K) + "\n"
        "m: " + str(m) + "\n"
        "plot: "+plot_log+"\n"
        "interrupt: " + str(interrupt) + "\n"
        "****************\n")

    """Finish setting up directories"""
    output_name = species
    if not tissue == "\"\"":
        output_name += "_"+tissue

    # Set up Image Directory
    image_directory = image_directory + "annotate_run" + str(run_num) + "_"+output_name + "/"
    os.mkdir(image_directory)
    print("Images will be automatically saved at \"" + image_directory + "\"")
    log("Set up \"annotate_run" + str(run_num) + "_"+output_name + "/\" directory in " + image_directory)

    # Set up Final Plot Directory
    final_plot_dir = final_plot_dir + "annotate_run" + str(run_num) + "_"+output_name + "/"
    os.mkdir(final_plot_dir)
    print("Anntation plot will be automatically saved at \"" + final_plot_dir + "\"")
    log("Set up \"annotate_run" + str(run_num) + "_"+output_name + "/\" directory" + final_plot_dir)

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

    # V&V plot parameter
    if (not 'X_umap' in adata.obsm) and 'X_tsne' in adata.obsm: tsne = True

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
    retried = False
    if not markers_path == "":
        try:
            #Check if direct path or in marker directory
            if os.path.isfile(marker_dir+markers_path):
                markers_path = marker_dir+markers_path
            else:
                assert os.path.isfile(markers_path), "Could not find file. Please place marker gene files in "+marker_dir+"."

            ref_marker = ct.read_markers_from_file(markers_path)

            #V&V. Transpose and retry if fail.
            if not VV_markers(ref_marker): #Retry loading file after transposing
                #Transpose data
                df = pd.read_csv(markers_path, header=None)
                df = df.T
                markers_path = markers_path[:-4] + "_corrected.csv"

                #Remove first numerical column
                new_header = df.iloc[0]
                df = df[1:]
                df.columns = new_header
                df.to_csv(markers_path, header=True)

                retried = True
                ref_marker = ct.read_markers_from_file(markers_path)  # Only with this can we properly check if the data is formated correctly
                assert not ref_marker.empty, "Failed to load marker gene file: " + markers_path

                #Upon Success
                print("Marker gene file was in an incorrect format and has been transposed to \""+markers_path+"\".")
                log("Marker gene file was in an incorrect format and has been transposed to \""+markers_path+"\".")
            print("Using marker file: \"" + markers_path + "\"")
            log("Loaded marker gene file \"" + markers_path + "\"")
            marker_loaded = True
        except:
            if retried:
                os.remove(markers_path)  # For discretion if failed, otherwise allow user to keep
                traceback.print_exc()
                print("There is an issue with the file format. Correct the file or continue to use provided marker gene data (Human and Mouse only).")
                log("There is an issue with the file format")
            else:
                print("Failed to load marker gene file. Try again or continue to use provided marker gene data (Human and Mouse only).")
                log("Failed to load marker gene file")

    # Otherwise use default marker data
    if not marker_loaded:
        print("Retrieving information from Cell Marker. Note only Human and Mouse data are available!\n"
              "If your species is neither of those, please provide a marker gene file.")
        # Prompt for species argument
        while species == "":
            species_prompt = input("Input species of the dataset: ")
            if not species_prompt == "":
                species = species_prompt

        species = species[0].upper() + species[1:].lower()
        # Change keyword to match format if applies
        if species == "Homo sapien" or species == "Homo sapiens":
            species = "Human"
        elif species == "Mus musculus":
            species = "Mouse"

        # V&V: Assert if species not available
        assert species == "Mouse" or species == "Human", "Could not find information for " + species + ". Please provide a marker gene file for that species in \"annotate.py\"."

        # Prompt for tissue argument
        while tissue == "tissue" or tissue == "":
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
        plot_annotation(adata, color=['cell_type'], title=['Cell Type Annotation for '+species+" "+tissue])
        save_figure("annotation_plot.png", final_plot_dir)

    # V&V cluster results
    if not '1' in adata.obs['leiden'].to_list():
        print("WARNING: Only one cluster was created! Please increase the cluster resolution!")

    #Repeated prompt for K, m, and cluster_res; recalculate annotation, and display results until user is satisfied.
    if interrupt:
        print("*Take this opportunity to experiment with different parameters regarding just the annotation process.*")
    cluster_res = 0.5
    res_i = 1
    ann_i = 1
    while interrupt or not ('1' in adata.obs['leiden'].to_list()):
        # Prompt for cluster resolution
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
                plot(adata, ['leiden'], "cluster_res=" + str(cluster_res))
                save_figure("cluster_(" + str(res_i) + ").png", image_directory)
                res_i += 1
            except ValueError:
                if not prompt == "":  # Not a crash due to mono-cluster warning
                    print("You must enter a decimal value!")

        print("Using cluster resolution of: " + str(cluster_res))

        # Calculate statistics
        sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
        log("Ranked gene groups")

        marker_df = ct.wrangle_ranks_from_anndata(adata)
        log("Wraggled ranks form adata")

        background = adata.var.index.tolist()

        # Repeat the previous steps and allow user to adjust K and m until satisifed.
        while interrupt:
            try:
                prompt = input("Enter a value for the top K genes to include in scoring (current: " + str(K) + ") or leave blank to keep the value." +
                               "\nOr enter \"cluster\" to re-adjust cluster resolution or \"exit\" to exit: ")
                if prompt == "cluster" or prompt == "exit":
                    break
                elif not prompt == "":
                    K = int(prompt)

                prompt = input("Enter a value for the m number of bins to use to divide the gene ranking (current: " + str(m) + ") or leave blank to keep the value: ")
                if not prompt == "": m = int(prompt)

                ct_pval, ct_score = ct.celltype_scores(nb_bins=m,
                                                       ranked_genes=marker_df,
                                                       K_top=K,
                                                       marker_ref=ref_marker,
                                                       background_genes=background)
                adata.obs['cell_type'] = ct.assign_celltypes(cluster_assignment=adata.obs['leiden'], ct_pval_df=ct_pval,
                                                             ct_score_df=ct_score)
                VnV_classification()
                plot_annotation(adata, ['cell_type'], "K=" + str(K) + ", m=" + str(m) + ", cluster_res=" + str(cluster_res))
                save_figure("annotation_plot_(" + str(ann_i) + ").png", final_plot_dir)
                ann_i += 1
            except:
                print("\nYou must enter integer values!")
        if (not interrupt) or prompt == "exit": break

    #Export results as an excel
    adata.obs = adata.obs.rename(columns={"leiden":"cluster"})  #Rename to avoid confusion
    ann_output_name = "annotate_run" + str(run_num) + "_" + output_name + "_annotation.xlsx"
    adata.obs.to_excel(annotation_dir+ann_output_name)
    VV_file_save(annotation_dir, ann_output_name, "Exported \"" + ann_output_name + "\"")

    #Export adata for future uses
    adata.write(adata_dir + "annotate_run" + str(run_num) + "_" + output_name + "_adata_annotated.h5ad")
    VV_file_save(adata_dir, "annotate_run" + str(run_num) + "_" + output_name + "_adata_annotated.h5ad", "Successfully exported adata_annotated")

    #Message about other python script
    print("*Annotation exported to \'"+output_name+'\'*')
    print("*Exported data to \'annotate_run" + str(run_num) + "_adata_annotated.h5ad\'*")

    # Log final parameters if on interrupt run
    if interrupt:
        log("\n***Final Parameters***\n"
            "adata path: " + adata_path + "\n"
            "marker gene file: " + markers_path + "\n"
            "species: " + species + "\n"
            "tissue: " + tissue + "\n"
            "cluster resolution: " + str(cluster_res) + "\n"
            "K: " + str(K) + "\n"
            "m: " + str(m) + "\n"
            "plot: "+plot_log+"\n"
            "interrupt: " + str(interrupt) + "\n"
            "***********************")

    log("*Program finished running without any errors*")
except Exception:
    traceback.print_exc()
    traceback.print_exc(file=file)
finally:
    file.close()