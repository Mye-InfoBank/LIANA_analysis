# LIANA Cell-Cell Interaction Analysis

This repository contains tools for performing comprehensive cell-cell interaction analysis using LIANA (Ligand-Receptor Analysis).

## Files Overview

- `liana_analysis_pipeline.py` - Main analysis pipeline
- `submit_liana_analysis.sh` - SLURM wrapper script for cluster execution
- `example_liana_job.sh` - Example job submission script with various configurations
- `liana_config_template.conf` - Configuration template for easy job setup

## Usage Options

### Option 1: Direct Python Execution

Start the script liana_analysis_pipeline.py with parameters:

### Required Parameters:

--adata_path \
--obs_key \
--output_dir

### Optional Parameters:

--splitting_key "condition" (default: None) \
--source_cell_types (default: all) \
--target_cell_types (default: all)

--condition_source_cell (default: auto-detected most abundant cell type) \
--condition_target_cell (default: auto-detected second most abundant cell type)

--logfc_threshold (default: 0.5) \
--lrscore_threshold (default: 0.9) \
--specificity_threshold (default: 0.05) \
--n_jobs (default: 4)

## Command:

```
python liana_analysis_pipeline.py  --adata_path /your/path/to/anndata/object --obs_key annotation_column  --output_dir /path/where/to/store/results  --splitting_key obs_column_you_want_to_compare --source_cell_types list_of_cell_types_source  --target_cell_types list_of_cell_types_target  --condition_source_cell cell_type_for_comparing_conditions  --condition_target_cell cell_type_for_comparing_target
```

#### Example usage:

```bash
# Example 1: Using defaults (auto-detected cell types)
python liana_analysis_pipeline.py \
  --adata_path /nfs/data/COST_IBD/versions/IBD/05_00_00/sub/results/finalized/base.h5ad \
  --obs_key "annotation:coarse" \
  --output_dir /nfs/data/COST_IBD/liana/results_full_adata \
  --splitting_key "condition"

# Example 2: Specifying specific cell types
python liana_analysis_pipeline.py \
  --adata_path /nfs/data/COST_IBD/versions/IBD/05_00_00/sub/results/finalized/base.h5ad \
  --obs_key "annotation:coarse" \
  --output_dir /nfs/data/COST_IBD/liana/results_full_adata \
  --splitting_key "condition" \
  --source_cell_types myeloid \
  --target_cell_types myeloid stromal epithelial T_NK_ILC B_plasma \
  --condition_source_cell myeloid \
  --condition_target_cell stromal
```

### Option 2: SLURM Cluster Execution (Recommended)

For large datasets and cluster environments, use the SLURM wrapper script:

#### Basic Usage:

```bash
sbatch submit_liana_analysis.sh \
  /path/to/your/data.h5ad \
  cell_type \
  /path/to/output
```

#### Usage with Condition Splitting:

```bash
sbatch submit_liana_analysis.sh \
  /path/to/your/data.h5ad \
  cell_type \
  /path/to/output \
  condition
```

#### Advanced Usage with Custom Resources:

```bash
sbatch --mem=128G --time=48:00:00 --cpus-per-task=16 submit_liana_analysis.sh \
  /path/to/your/data.h5ad \
  cell_type \
  /path/to/output \
  condition
```

#### Parameter Details:

The SLURM wrapper accepts the following parameters in order:

1. **adata_path** (required): Path to your AnnData (.h5ad) file
2. **obs_key** (required): Column name in adata.obs containing cell type annotations
3. **output_dir** (required): Directory where results will be saved
4. **splitting_key** (optional): Column name for condition comparison

**Examples:**

```bash
# Basic analysis
sbatch submit_liana_analysis.sh data.h5ad cell_type results/

# With condition splitting
sbatch submit_liana_analysis.sh data.h5ad cell_type results/ condition

# With quoted obs_key (for complex column names)
sbatch submit_liana_analysis.sh data.h5ad "annotation:coarse" results/ treatment_status
```

## SLURM Job Management

### Monitor Jobs:

```bash
squeue -u $USER                    # Check job status
scontrol show job JOBID            # Check job details
sacct -j JOBID --format=JobID,JobName,MaxRSS,Elapsed,State  # Resource usage
```

### Cancel Jobs:

```bash
scancel JOBID                      # Cancel specific job
scancel -u $USER                   # Cancel all your jobs
```

## Resource Recommendations

| Dataset Size             | Memory | CPUs  | Time Limit | Partition |
| ------------------------ | ------ | ----- | ---------- | --------- |
| Small (<10k cells)       | 32G    | 4-8   | 12:00:00   | cpu       |
| Medium (10k-50k cells)   | 64G    | 8-16  | 24:00:00   | cpu       |
| Large (50k-100k cells)   | 128G   | 16-32 | 48:00:00   | highmem   |
| Very Large (>100k cells) | 256G+  | 32+   | 72:00:00   | highmem   |

## Output Structure

```
output_dir/
├── liana_results.csv                    # Full dataset results
├── liana_results.pkl                    # Full dataset results (pickle)
├── full_interaction_network.png         # Network visualization
├── full_interaction_heatmap_*.png       # Heatmaps
├── full_detailed_interactions.png       # Detailed interaction plot
├── condition_scatter_*.png              # Condition comparison (if splitting_key provided)
├── dotplot_top10_magnitude_rank.png     # Top interactions dot plot
├── logs/                                # SLURM job logs
│   ├── liana_JOBID.out
│   └── liana_JOBID.err
└── [CONDITION]/                         # Per-condition subdirectories (if splitting)
    ├── liana_results.csv
    ├── liana_results.pkl
    └── [CONDITION]_*.png                # Condition-specific plots
```

## Troubleshooting

### Common Issues:

1. **Out of Memory**: Increase `--mem` parameter or use highmem partition
2. **Time Limit Exceeded**: Increase `--time` parameter
3. **Job Fails to Start**: Check resource availability with `sinfo`
4. **Python Package Missing**: Ensure all dependencies are installed in your environment

### Log Files:

- SLURM output: `logs/liana_JOBID.out`
- SLURM errors: `logs/liana_JOBID.err`
- Pipeline logs: Integrated in SLURM output with timestamps

## Key Features

### Auto-Detection of Cell Types

- **Condition Source/Target Cells**: The pipeline automatically detects the most abundant and second most abundant cell types in your dataset for condition comparisons
- **Cell Type Resolution**: Use `"all"` for source/target cell types to include all cell types from your dataset
- **Smart Defaults**: No need to specify cell types manually unless you want specific comparisons

### Comprehensive Analysis

- Full dataset analysis with network visualizations and heatmaps
- Condition-specific analysis when splitting key is provided
- Multiple visualization types: networks, heatmaps, scatter plots, and dot plots

## Dependencies

- Python 3.8+
- scanpy
- liana
- pandas
- numpy
- matplotlib
- seaborn
- networkx
- anndata

## Example Configurations

See `example_liana_job.sh` for various job submission examples including:

- Basic analysis without condition splitting
- Analysis with condition splitting
- High-memory analysis with custom parameters
- Analysis focusing on specific cell types
