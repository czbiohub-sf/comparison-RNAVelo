import numpy as np
import pandas as pd
import anndata
import os
import scanpy as sc
import scvelo as scv

#conda environment is vel38zebra

'''Takes anndata objects with velocity methods and creates a dataframe with the following columns:
| cell_id | magnitude of velocity vector | magnitude of transition vector | 
cosine similarity of transition vector with full object | proportion (subset) | iteration | velocity method |

Run for a single method at a time, with inputs as all subset adatas for the velocity method and the full velocity method
'''


def load_anndata(file_path):
    return anndata.read_h5ad(file_path)

def scale_transition_matrices(adata, key1):
    A = adata.uns[key1].copy()
    A_norm = scv.utils.get_transition_matrix(adata=adata, vgraph=A)
    return A_norm

def align_matrices_with_velocity(transition_matrix, fixed_matrix, velocity_matrix, transition_cell_ids, fixed_cell_ids):
    # Create pandas DataFrames for transition and fixed matrices for easy manipulation
    transition_df = pd.DataFrame(transition_matrix.toarray(), index=transition_cell_ids, columns=transition_cell_ids)
    fixed_df = pd.DataFrame(fixed_matrix.toarray(), index=fixed_cell_ids, columns=fixed_cell_ids)
    velocity_df = pd.DataFrame(velocity_matrix, index=transition_cell_ids)
    
    # Identify common cell IDs for transition and fixed matrices
    common_cell_ids = sorted(set(transition_cell_ids).intersection(set(fixed_cell_ids)))

    # Reorder and align both matrices based on common cell IDs
    aligned_transition_df = transition_df.loc[common_cell_ids, common_cell_ids]
    aligned_fixed_df = fixed_df.loc[common_cell_ids, common_cell_ids]
    aligned_velocity_df = velocity_df.loc[common_cell_ids]  # Reorder only the rows
    return aligned_transition_df.values, aligned_fixed_df.values, aligned_velocity_df.values, common_cell_ids

def compute_magnitudes(matrix):
    return np.linalg.norm(matrix, axis=1)


def compute_cosine_similarities(matrix, fixed_matrix):
    dot_products = np.sum(matrix * fixed_matrix, axis=1)
    norms = np.linalg.norm(matrix, axis=1) * np.linalg.norm(fixed_matrix, axis=1)
    return dot_products / norms

def extract_metadata_from_path(path):
    filename = os.path.basename(path)  # Get the base name of the file
    filename_without_extension = filename.replace('.h5ad', '')  # Remove the file extension for easier parsing
    
    # Split by the first hyphen to separate allfish_x from the rest
    first_split = filename_without_extension.split('-', 1)
    dataset_proportion = first_split[0]  # 'allfish_x'
    remainder = first_split[1]  # '1-scvelo-deterministic'

    # Split the dataset_proportion to extract proportion
    proportion = dataset_proportion.replace('allfish_', '')
    
    # Split the remainder to extract iteration and method
    iteration_method_split = remainder.split('-', 1)
    iteration = iteration_method_split[0]
    method = iteration_method_split[1]

    return proportion, iteration, method

def make_dataframe_foradata(adata_path, fixed_matrix, fixed_cell_ids):
    mydata={}
    adata = load_anndata(adata_path)

    velocity_matrix = adata.layers['velocity']
    transition_matrix = scale_transition_matrices(adata, 'velocity_graph')
    transition_cell_ids = adata.obs.index.tolist()

    #align matrices, get correct ids and subset
    a_trans, a_fixed, a_velo, cell_id = align_matrices_with_velocity(transition_matrix, fixed_matrix, velocity_matrix,transition_cell_ids,fixed_cell_ids)
    mydata['cell_id']=cell_id
    mydata['mag_velo']=compute_magnitudes(a_velo)
    mydata['mag_trans']=compute_magnitudes(a_trans)
    mydata['cosine_sim']=compute_cosine_similarities(a_trans, a_fixed)

    df=pd.DataFrame.from_dict(mydata)
    proportion, iteration, method = extract_metadata_from_path(adata_path)
    df['proportion']=proportion
    df['iteration']=iteration
    df['velo_method']=method

    return df


def process_all_adata_files(input_path, save_path, fixed_matrix, fixed_cell_ids):
    # List all files in the input folder that end with '.h5ad'
    file_paths = [os.path.join(input_path, file) for file in os.listdir(input_path) if file.endswith('.h5ad')]
    
        # Filter files to include only those starting with "allfish_0.98" or "allfish_0.95"
    #filtered_paths = [file for file in file_paths if "allfish_0.98" in file or "allfish_0.95" in file]
    
    # List to store each dataframe created for each .h5ad file
    dataframes = []
    
    # Process each .h5ad file
    for file_path in file_paths:
        # Call your custom function to process the file and generate a dataframe
        df = make_dataframe_foradata(file_path, fixed_matrix, fixed_cell_ids)
        dataframes.append(df)
    
    # Concatenate all dataframes into one
    if dataframes:
        combined_df = pd.concat(dataframes, ignore_index=True)
    else:
        combined_df = pd.DataFrame()  # Return an empty dataframe if no files were processed
    
    return combined_df.to_csv(save_path+'cosim_mag_uniTvelo_May30.csv')


#test with a different adata_full (same pipeline)
#input the adata object with 100% of the reads, as the object to compare to
adata_full=sc.read_h5ad('/velocity_subsets_zf24hpf/uniTvelo/allfish_1-1-uniTvelo.h5ad')
fixed_matrix = scale_transition_matrices(adata_full, 'velocity_graph')
fixed_cell_ids = adata_full.obs.index.tolist()
#input path for all the subset adata objects, run with a specific method and save path for the dataframe
input_path='/velocity_subsets_zf24hpf/uniTvelo'
save_path='/velocity_subsets_zf24hpf/cosim_mag_csvs/'

process_all_adata_files(input_path, save_path, fixed_matrix, fixed_cell_ids)





