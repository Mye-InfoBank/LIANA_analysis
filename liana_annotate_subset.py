import os
import platform
import scanpy as sc
import liana as li
import pandas as pd

# define variables

# adata object
adata_path = '/nfs/data/COST_IBD/versions/IBD/03_00_00/build/results/finalized/annotated/liana/core_atlas_liana_cleaned_coarse_annotation.h5ad'
# which obs column to use as the label
obs_key = "annotation:coarse:cleaned"
# output directory for plots and results
output_path = '/nfs/data/COST_IBD/versions/IBD/03_00_00/build/results/finalized/annotated/liana/'
# intermediate folder for storing subatlases
folder_for_storing = 'xxx'

# key for subsetting (f.e. condition)
splitting_key = "condition"

# starting creating sub_adatas for subsetting

adata = sc.read_h5ad(adata_path)
adata.obs
print('adata loaded')

unique_key_values = adata.obs[splitting_key].unique()
subadatas = {}
for key_value in unique_key_values:
    subadatas[key_value] = adata[adata.obs[splitting_key] == key_value].copy()
    # write the h5ad of each subadata into folder with sub_key name
    # create path for storing in  folder_for_storing/splitting_key/sub_key.h5ad
    if not os.path.exists(os.path.join(folder_for_storing, splitting_key)):
        os.makedirs(os.path.join(folder_for_storing, splitting_key))
    subadatas[key_value].write_h5ad(os.path.join(folder_for_storing, splitting_key, key_value + '.h5ad'))

print('adata split done')

if (adata.X < 0).nnz == 0:
    sc.pp.log1p(adata)
print('log1p done')
    
# Run LIANA
for key_value in unique_key_values:
    print(f"running LIANA for sub {key_value}")
    #print(os.path.join(folder_for_storing, splitting_key, key_value +'/'+ '.h5ad'))
    adata_path = os.path.join(folder_for_storing, splitting_key+ '/'+ key_value + '.h5ad')
    output_path = os.path.join(folder_for_storing, splitting_key, key_value, key_value)
    
    adata = sc.read_h5ad(adata_path)
    if adata.obs[obs_key].nunique() > 1:
        li.mt.rank_aggregate(adata, obs_key, use_raw=False, verbose=True, n_jobs=4)
        adata.uns["liana_res"] = adata.uns["liana_res"]
        #updated_atlas_path = output_path + "core_atlas_liana.h5ad"
        adata.write_h5ad(adata_path)
        liana_results_df = pd.DataFrame(adata.uns["liana_res"])
        
        # create output_path if it does not exist
        if not os.path.exists(output_path):
            os.makedirs(output_path)
                
        # Save as pickle
        liana_results_df.to_pickle(output_path + "_liana_results_cleaned.pkl")
        # Save as CSV
        liana_results_df.to_csv(output_path + "_liana_results_cleaned.csv", index=False)
    else:
        print(f"Skipping LIANA as {obs_key} has only one unique value.")