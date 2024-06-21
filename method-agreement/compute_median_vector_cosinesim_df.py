import pandas as pd
import anndata
import scvelo as scv
import numpy as np

def load_anndata(file_path):
    return anndata.read_h5ad(file_path)

def compute_cosine_similarities(matrix, fixed_matrix):
    dot_products = np.sum(matrix * fixed_matrix, axis=1)
    norms = np.linalg.norm(matrix, axis=1) * np.linalg.norm(fixed_matrix, axis=1)
    return dot_products / norms

def scale_transition_matrices(adata, key1):
    A = adata.uns[key1].copy()
    A_norm = scv.utils.get_transition_matrix(adata=adata, vgraph=A)
    return A_norm

def align_matrices(matrices, cell_ids):
    # Initialize a DataFrame from the first matrix to start the alignment
    
    reference_df = pd.DataFrame(matrices[0], index=cell_ids[0], columns=cell_ids[0])
    aligned_matrices = [reference_df.values]
    
    # Align each subsequent matrix to the reference cell ID order
    for matrix, ids in zip(matrices[1:], cell_ids[1:]):
        df = pd.DataFrame(matrix, index=ids, columns=ids)
        aligned_df = df.reindex(index=reference_df.index, columns=reference_df.index)
        aligned_matrices.append(aligned_df.values)
    
    return np.array(aligned_matrices), reference_df.index.tolist()

def compute_median_matrix(matrices):
    # Stack the aligned matrices and compute the median across the first axis
    stacked_matrices = np.stack(matrices, axis=0)
    median_matrix = np.median(stacked_matrices, axis=0)
    return median_matrix

#example of running for ZF whole embryo 24hpf
adata_path1 = '/zebrahub_24hrs_dataset/deepVelo_Jun22.h5ad'
adata_path2 = '/zebrahub_24hrs_dataset/scVelo_deterministic_Jun22.h5ad'
adata_path3 = '/zebrahub_24hrs_dataset/scVelo_dynamical_Jun22.h5ad'
adata_path4 = '/zebrahub_24hrs_dataset/scVelo_stochastic_Jun22.h5ad'
adata_path5 = '/zebrahub_24hrs_dataset/uniTVelo_Jun22.h5ad'

save_path='/zebrahub_24hrs_dataset/'
file_paths = [adata_path1, adata_path2, adata_path3, adata_path4, adata_path5]
matrices = []
cell_ids = []
methods = ['DeepVelo', 'Velocyto', 'scv-Dyn', 'scv-Sto', 'UniTVelo']

for file_path in file_paths:
    adata = load_anndata(file_path)
    matrix = scale_transition_matrices(adata, 'velocity_graph').toarray()
    ids = adata.obs.index.tolist().copy()
    matrices.append(matrix)
    cell_ids.append(ids)



aligned_matrices, common_cell_ids = align_matrices(matrices, cell_ids)
median_matrix = compute_median_matrix(aligned_matrices)



# Now median_matrix is your median transition matrix and common_cell_ids are the cell IDs in order
print("Median Matrix Shape:", median_matrix.shape)
print("Common Cell IDs:", common_cell_ids[:10])  # Show the first 10

# Compute cosine similarities for each matrix against the median matrix
cosine_similarities = {method: compute_cosine_similarities(matrix, median_matrix) for method, matrix in zip(methods, aligned_matrices)}

# Construct the DataFrame
cosine_similarity_df = pd.DataFrame(cosine_similarities, index=common_cell_ids)
print(cosine_similarity_df.head())

cosine_similarity_df.to_csv(save_path+'zebrahub_24hrs_cosine_sim_medianvec_May14.csv')




