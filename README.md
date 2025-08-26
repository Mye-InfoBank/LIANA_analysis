# How to use the LIANA analysis

1. Clone the repository and set up all python scripts in one folder

2. Start the script liana_analysis_pipeline.py with parameters:

### Required Parameters:

--adata_path \
--obs_key \
--output_dir

### Optional Parameters:

--splitting_key "condition" (default: None) \
--source_cell_types (default: myeloid) \
--target_cell_types (default: myeloid stromal epithelial T_NK_ILC B_plasma)

--condition_source_cell (default: myeloid) \
--condition_target_cell (default: stromal) \

--logfc_threshold (default: 0.5) \
--lrscore_threshold (default: 0.9) \
--specificity_threshold (default: 0.05) \
--n_jobs (default: 4)

## Command:

```
python liana_analysis_pipeline.py  --adata_path /your/path/to/anndata/object --obs_key annotation_column --output_dir /path/where/to/store/results  --splitting_key obs_column_you_want_to_compare --source_cell_types list_of_cell_types_source  --target_cell_types list_of_cell_types_target  --condition_source_cell cell_type_for_comparing_conditions  --condition_target_cell cell_type_for_comparing_target 
```

### Example usage:
```
python liana_analysis_pipeline.py   --adata_path /nfs/data/COST_IBD/versions/IBD/05_00_00/sub/results/finalized/base.h5ad   --obs_key "annotation:coarse"   --output_dir /nfs/data/COST_IBD/liana/results_full_adata   --splitting_key "condition"   --source_cell_types myeloid   --target_cell_types myeloid stromal epithelial T_NK_ILC B_plasma   --condition_source_cell myeloid   --condition_target_cell stromal 
```
