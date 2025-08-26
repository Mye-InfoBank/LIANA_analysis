#!/usr/bin/env python3
"""
LIANA Cell-Cell Interaction Analysis Pipeline
This script performs comprehensive cell-cell interaction analysis using LIANA.
It includes data preprocessing, interaction analysis, and visualization.
"""

import os
import argparse
import scanpy as sc
import liana as li
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from glob import glob
import anndata as ad
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='LIANA Cell-Cell Interaction Analysis Pipeline')
    
    # Required arguments
    parser.add_argument('--adata_path', required=True, help='Path to the AnnData object')
    parser.add_argument('--obs_key', required=True, help='Observation key for cell type annotation')
    parser.add_argument('--output_dir', required=True, help='Output directory for results')
    
    # Optional arguments
    parser.add_argument('--splitting_key', default=None, 
                       help='Key for subsetting data (e.g., condition)')
    parser.add_argument('--logfc_threshold', type=float, default=0.5, 
                       help='LogFC threshold for significant interactions')
    parser.add_argument('--lrscore_threshold', type=float, default=0.9, 
                       help='LRscore threshold for significant interactions')
    parser.add_argument('--specificity_threshold', type=float, default=0.05, 
                       help='Specificity threshold for significant interactions')
    parser.add_argument('--n_jobs', type=int, default=4, 
                       help='Number of jobs for parallel processing')
    
    # Visualization parameters
    parser.add_argument('--source_cell_types', nargs='+', default=['myeloid'],
                       help='Source cell types for detailed visualization')
    parser.add_argument('--target_cell_types', nargs='+', 
                       default=['myeloid', 'stromal', 'epithelial', 'T_NK_ILC', 'B_plasma'],
                       help='Target cell types for detailed visualization')
    parser.add_argument('--condition_source_cell', default='myeloid',
                       help='Source cell type for condition comparison')
    parser.add_argument('--condition_target_cell', default='stromal',
                       help='Target cell type for condition comparison')
    
    return parser.parse_args()

def setup_directories(output_dir, splitting_key=None):
    """Create necessary directories for output."""
    os.makedirs(output_dir, exist_ok=True)
    
    if splitting_key:
        split_dir = os.path.join(output_dir, splitting_key)
        os.makedirs(split_dir, exist_ok=True)
        return split_dir
    return output_dir

def load_and_preprocess_data(adata_path, obs_key):
    """Load and preprocess the AnnData data."""
    logger.info(f"Loading AnnData from {adata_path}")
    adata = sc.read_h5ad(adata_path)
    
    # Check if data needs log transformation
    if (adata.X.data if hasattr(adata.X, 'data') else adata.X < 0).any():
        logger.info("Applying log1p transformation")
        sc.pp.log1p(adata)
    
    logger.info(f"Data shape: {adata.shape}")
    logger.info(f"Cell types: {adata.obs[obs_key].unique()}")
    
    return adata

def split_data_by_condition(adata, splitting_key, output_dir):
    """Split data by condition and save sub-adatas."""
    logger.info(f"Splitting data by {splitting_key}")
    
    unique_key_values = adata.obs[splitting_key].unique()
    subadatas = {}
    
    for key_value in unique_key_values:
        subadatas[key_value] = adata[adata.obs[splitting_key] == key_value].copy()
        
        # Save sub-adatas
        subadata_path = os.path.join(output_dir, f"{key_value}.h5ad")
        subadatas[key_value].write_h5ad(subadata_path)
        logger.info(f"Saved sub-adata for {key_value} to {subadata_path}")
    
    return subadatas, unique_key_values

def run_liana_analysis(adata, obs_key, output_path, n_jobs=4):
    """Run LIANA analysis on the provided data."""
    logger.info(f"Running LIANA analysis with {obs_key} as cell type key")
    
    if adata.obs[obs_key].nunique() <= 1:
        logger.warning(f"Skipping LIANA as {obs_key} has only one unique value.")
        return None, adata
    
    # Run LIANA rank aggregate
    li.mt.rank_aggregate(adata, obs_key, use_raw=False, verbose=True, n_jobs=n_jobs)
    
    # Save results
    liana_results_df = pd.DataFrame(adata.uns["liana_res"])
    
    # Save as CSV and pickle
    liana_results_df.to_pickle(os.path.join(output_path, "liana_results.pkl"))
    liana_results_df.to_csv(os.path.join(output_path, "liana_results.csv"), index=False)
    
    logger.info(f"LIANA results saved to {output_path}")
    return liana_results_df, adata

def filter_significant_interactions(df, logfc_threshold, lrscore_threshold, specificity_threshold):
    """Filter significant interactions based on thresholds."""
    logger.info("Filtering significant interactions")
    
    # Use .copy() to avoid SettingWithCopyWarning
    significant_interactions = df[
        (df["lr_logfc"] > logfc_threshold) &
        (df["lrscore"] > lrscore_threshold) &
        (df["specificity_rank"] < specificity_threshold)
    ].copy()
    
    # Use .loc to avoid SettingWithCopyWarning
    significant_interactions.loc[:, 'source_target'] = significant_interactions['source'] + "_" + significant_interactions['target']
    significant_interactions.loc[:, 'count'] = 1.0
    
    logger.info(f"Number of significant interactions: {len(significant_interactions)}")
    return significant_interactions

def visualize_interaction_network(df_filtered, output_path, plot_name="interaction_network"):
    """Create and save interaction network visualization."""
    logger.info("Creating interaction network visualization")
    
    # Count interactions between source-target pairs
    interaction_counts = df_filtered.groupby(["source", "target"]).size().reset_index(name="count")
    
    if interaction_counts.empty:
        logger.warning("No interactions to visualize in network")
        return
    
    # Create a directed graph
    G = nx.DiGraph()
    
    # Add edges with weights
    for source, target, count in interaction_counts.itertuples(index=False):
        G.add_edge(source, target, weight=count)
    
    # Draw the network
    plt.figure(figsize=(12, 9))
    pos = nx.spring_layout(G, seed=42)
    
    # Draw nodes
    nx.draw_networkx_nodes(G, pos, node_size=2000, node_color="lightblue", edgecolors="black")
    
    # Draw edges with width proportional to interaction count
    max_count = max(interaction_counts["count"]) if not interaction_counts.empty else 1
    edge_widths = [G[u][v]["weight"] / max_count * 5 for u, v in G.edges()]
    
    nx.draw_networkx_edges(
        G, pos, arrowstyle="->", arrowsize=20, width=edge_widths
    )
    
    # Draw labels
    nx.draw_networkx_labels(G, pos, font_size=10, font_weight="bold")
    
    plt.title("Cell-Cell Interaction Network (Significant Interactions)")
    plt.tight_layout()
    
    # Save the plot
    network_plot_path = os.path.join(output_path, f"{plot_name}.png")
    plt.savefig(network_plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Network visualization saved to {network_plot_path}")

def visualize_interaction_heatmap(significant_interactions, value_for_heatmap, output_path, plot_name_prefix=""):
    """Create and save interaction heatmap."""
    logger.info(f"Creating interaction heatmap for {value_for_heatmap}")
    
    # Create a pivot table for the heatmap
    interaction_matrix = significant_interactions.pivot_table(
        index="source", columns="target", values=value_for_heatmap, aggfunc="median"
    ).fillna(0)
    
    if interaction_matrix.empty:
        logger.warning(f"No data to create heatmap for {value_for_heatmap}")
        return
    
    # Plot heatmap
    plt.figure(figsize=(12, 10))
    sns.heatmap(interaction_matrix, cmap="Blues", annot=True, fmt=".2f", linewidths=0.5)
    plt.title(f"Cell-Cell Interaction Strength ({value_for_heatmap}) Heatmap")
    plt.xlabel("Target Cell Type")
    plt.ylabel("Source Cell Type")
    plt.tight_layout()
    
    # Save the plot
    heatmap_name = f"{plot_name_prefix}interaction_heatmap_{value_for_heatmap}.png"
    heatmap_path = os.path.join(output_path, heatmap_name)
    plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Heatmap saved to {heatmap_path}")

def create_detailed_interactions_plot(liana_results_df, source_cell_types, target_cell_types, output_path, plot_name_prefix=""):
    """Create a custom detailed interactions plot using the LIANA results."""
    logger.info("Creating custom detailed interactions plot")
    
    # Filter for the specified source and target cell types
    filtered_df = liana_results_df[
        (liana_results_df['source'].isin(source_cell_types)) &
        (liana_results_df['target'].isin(target_cell_types))
    ]
    
    if filtered_df.empty:
        logger.warning("No interactions found for the specified cell types")
        return
    
    # Get top interactions by magnitude_rank
    top_interactions = filtered_df.sort_values('magnitude_rank').head(20)
    
    if top_interactions.empty:
        logger.warning("No top interactions found after filtering")
        return
    
    # Create a custom dot plot
    plt.figure(figsize=(14, 10))
    
    # Create a scatter plot with size and color encoding
    scatter = plt.scatter(
        x=range(len(top_interactions)),
        y=top_interactions['lrscore'],
        s=100 - top_interactions['specificity_rank'] * 80,  # Inverse size for specificity
        c=top_interactions['magnitude_rank'],
        cmap='plasma',
        alpha=0.7
    )
    
    # Add labels
    interaction_labels = [
        f"{row['ligand_complex']} → {row['receptor_complex']}\n({row['source']}→{row['target']})" 
        for _, row in top_interactions.iterrows()
    ]
    
    plt.xticks(range(len(top_interactions)), interaction_labels, rotation=45, ha='right')
    plt.ylabel('LRscore')
    plt.title('Top Ligand-Receptor Interactions')
    
    # Add colorbar
    plt.colorbar(scatter, label='Magnitude Rank')
    
    # Add legend for size
    sizes = [50, 30, 10]  # Example sizes
    labels = ['High specificity', 'Medium', 'Low']
    for i, size in enumerate(sizes):
        plt.scatter([], [], s=size, c='k', alpha=0.7, label=labels[i])
    plt.legend(scatterpoints=1, frameon=False, labelspacing=1, title='Specificity')
    
    plt.tight_layout()
    
    # Save the plot
    plot_name = f"{plot_name_prefix}detailed_interactions.png"
    plot_path = os.path.join(output_path, plot_name)
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Custom detailed interactions plot saved to {plot_path}")

def create_condition_scatterplots(csv_dir, output_path, source_cell, target_cell, adata, source_cell_type_list, target_cell_type_list):
    """Create scatterplots comparing interactions across conditions."""
    logger.info(f"Creating condition scatterplots for {source_cell} → {target_cell}")
    
    # Find all CSV files in subdirectories
    csv_files = []
    for root, dirs, files in os.walk(csv_dir):
        for file in files:
            if file == "liana_results.csv":
                csv_files.append(os.path.join(root, file))
    
    if not csv_files:
        logger.warning(f"No CSV files found in {csv_dir} or its subdirectories")
        return
    
    # Load and combine all CSVs into a single DataFrame
    dfs = []
    for file_path in csv_files:
        # Extract condition name from directory structure
        condition_dir = os.path.dirname(file_path)
        condition_name = os.path.basename(condition_dir)
        
        df = pd.read_csv(file_path)
        df['condition'] = condition_name
        dfs.append(df)
    
    if not dfs:
        logger.warning("No data loaded from CSV files")
        return
        
    df_combined = pd.concat(dfs, ignore_index=True)
    
    # Create ligand → receptor column safely
    df_combined = df_combined.copy()
    df_combined.loc[:, 'ligand → receptor'] = df_combined['ligand_complex'] + ' → ' + df_combined['receptor_complex']
    
    # Filter for specific cell types
    subset = df_combined[
        (df_combined['source'] == source_cell) &
        (df_combined['target'] == target_cell)
    ]
    
    if subset.empty:
        logger.info(f"No interactions found for {source_cell} → {target_cell}")
        return
    
    # Sort by magnitude and get top pairs
    top = subset.sort_values(by='lr_means', ascending=False).head(25)
    
    if top.empty:
        logger.info(f"No significant interactions after filtering for {source_cell} → {target_cell}")
        return
    
    # Create title
    title_text = f"{source_cell} → {target_cell} Ligand-Receptor Interactions"
    
    # Dynamic figure sizing
    pair_height = 0.25
    condition_width = 1.1
    num_pairs = len(top['ligand → receptor'].unique())
    num_conditions = len(top['condition'].unique())
    
    fig_height = max(6, num_pairs * pair_height)
    fig_width = max(8, num_conditions * condition_width)
    
    # Create the scatterplot
    plt.figure(figsize=(fig_width, fig_height))
    scatter = sns.scatterplot(
        data=top,
        x='condition',
        y='ligand → receptor',
        hue='lr_means',
        size='specificity_rank',
        palette='viridis',
        sizes=(50, 200),
        edgecolor='gray',
        linewidth=0.5
    )
    
    plt.title(title_text)
    plt.xlabel('Condition')
    plt.ylabel('Ligand → Receptor')
    plt.tight_layout()
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Save the plot
    scatterplot_path = os.path.join(output_path, f"condition_scatter_{source_cell}_{target_cell}.png")
    plt.savefig(scatterplot_path, dpi=300, bbox_inches='tight')
    plt.close()

    # Create a custom dot plot instead of using sc.pl.dotplot
    # Filter for top interactions by magnitude_rank
    top_interactions = df_combined.sort_values('magnitude_rank').head(10)
    
    if not top_interactions.empty:
        # Create a custom dot plot
        plt.figure(figsize=(12, 6))
        
        # Create a scatter plot with size and color encoding
        scatter = plt.scatter(
            x=range(len(top_interactions)),
            y=top_interactions['lrscore'],
            s=100 - top_interactions['specificity_rank'] * 80,  # Inverse size for specificity
            c=top_interactions['magnitude_rank'],
            cmap='plasma',
            alpha=0.7
        )
        
        # Add labels
        interaction_labels = [
            f"{row['ligand_complex']} → {row['receptor_complex']}\n({row['source']}→{row['target']})" 
            for _, row in top_interactions.iterrows()
        ]
        
        plt.xticks(range(len(top_interactions)), interaction_labels, rotation=45, ha='right')
        plt.ylabel('LRscore')
        plt.title('Top 10 Ligand-Receptor Interactions by Magnitude Rank')
        
        # Add colorbar
        plt.colorbar(scatter, label='Magnitude Rank')
        
        # Add legend for size
        sizes = [80, 50, 20]  # Example sizes
        labels = ['High specificity', 'Medium', 'Low']
        for i, size in enumerate(sizes):
            plt.scatter([], [], s=size, c='k', alpha=0.7, label=labels[i])
        plt.legend(scatterpoints=1, frameon=False, labelspacing=1, title='Specificity')
        
        plt.tight_layout()
        
        # Save the figure
        dotplot_path = os.path.join(output_path, f"dotplot_top10_magnitude_rank.png")
        plt.savefig(dotplot_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Custom dot plot saved to {dotplot_path}")
    
    logger.info(f"Condition scatterplot saved to {scatterplot_path}")
    logger.info(f"Dot plot saved to {dotplot_path}")

def main():
    """Main function to run the complete LIANA analysis pipeline."""
    args = parse_arguments()
    
    # Set up directories
    output_dir = setup_directories(args.output_dir, args.splitting_key)
    
    # Load and preprocess data
    adata = load_and_preprocess_data(args.adata_path, args.obs_key)
    
    # Run LIANA on the full dataset first (always)
    logger.info("Running LIANA on full dataset")
    liana_results_full, adata_with_results = run_liana_analysis(adata, args.obs_key, output_dir, args.n_jobs)
    
    if liana_results_full is not None:
        # Filter significant interactions for full dataset
        significant_interactions_full = filter_significant_interactions(
            liana_results_full, 
            args.logfc_threshold, 
            args.lrscore_threshold, 
            args.specificity_threshold
        )
        
        # Visualizations for full dataset
        visualize_interaction_network(significant_interactions_full, output_dir, "full_interaction_network")
        visualize_interaction_heatmap(significant_interactions_full, "lr_logfc", output_dir, "full_")
        visualize_interaction_heatmap(significant_interactions_full, "count", output_dir, "full_")
        
        # Create custom detailed interactions plot
        create_detailed_interactions_plot(
            liana_results_full, 
            args.source_cell_types, 
            args.target_cell_types, 
            output_dir, 
            "full_"
        )
    
    # Split data if splitting_key is provided
    if args.splitting_key:
        subadatas, unique_key_values = split_data_by_condition(
            adata, args.splitting_key, output_dir
        )
        
        # Run LIANA on each sub-adata
        for key_value in unique_key_values:
            logger.info(f"Processing {key_value}")
            subadata_path = os.path.join(output_dir, f"{key_value}.h5ad")
            subadata = sc.read_h5ad(subadata_path)
            
            # Create output directory for this condition
            condition_output_dir = os.path.join(output_dir, key_value)
            os.makedirs(condition_output_dir, exist_ok=True)
            
            # Run LIANA
            liana_results, subadata_with_results = run_liana_analysis(
                subadata, args.obs_key, condition_output_dir, args.n_jobs
            )
            
            if liana_results is not None:
                # Filter significant interactions
                significant_interactions = filter_significant_interactions(
                    liana_results, 
                    args.logfc_threshold, 
                    args.lrscore_threshold, 
                    args.specificity_threshold
                )
                
                # Visualizations for this condition
                visualize_interaction_network(significant_interactions, condition_output_dir, f"{key_value}_interaction_network")
                visualize_interaction_heatmap(significant_interactions, "lr_logfc", condition_output_dir, f"{key_value}_")
                visualize_interaction_heatmap(significant_interactions, "count", condition_output_dir, f"{key_value}_")
                
                # Create custom detailed interactions plot for this condition
                create_detailed_interactions_plot(
                    liana_results, 
                    args.source_cell_types, 
                    args.target_cell_types, 
                    condition_output_dir, 
                    f"{key_value}_"
                )
        
        # Create condition comparison scatterplots
        create_condition_scatterplots(
            output_dir, 
            output_dir,
            args.condition_source_cell,
            args.condition_target_cell,
            adata,
            args.source_cell_types,
            args.target_cell_types
        )
    
    logger.info("LIANA analysis pipeline completed successfully")

if __name__ == "__main__":
    main()