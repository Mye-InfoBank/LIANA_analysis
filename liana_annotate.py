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


# starting liana
print('start liana')

adata = sc.read_h5ad(adata_path)
adata.obs
print('adata loaded')

if (adata.X < 0).nnz == 0:
    sc.pp.log1p(adata)
print('log1p done')
    
# Run LIANA
if adata.obs[obs_key].nunique() > 1:
    print('inside liana')
    li.mt.rank_aggregate(adata, obs_key, use_raw=False, verbose=True, n_jobs=4)
    print('liana done')

    adata.uns["liana_res"] = adata.uns["liana_res"]
    
    updated_atlas_path = output_path + "atlas_liana.h5ad"
    # Save the updated AnnData object
    adata.write_h5ad(updated_atlas_path)
    print('write done')
    
    liana_results_df = pd.DataFrame(adata.uns["liana_res"])

    # Save as pickle
    liana_results_df.to_pickle(output_path + "liana_results_cleaned.pkl")

    # Save as CSV
    liana_results_df.to_csv(output_path + "liana_results_cleaned.csv", index=False)
    print('liana_results saving done')
    
else:
    print(f"Skipping LIANA as {obs_key} has only one unique value.")