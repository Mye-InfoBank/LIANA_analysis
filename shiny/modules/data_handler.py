#!/usr/bin/env python3
"""
Data Handler Module for LIANA Results Explorer
Handles data loading, processing, and filtering operations.
"""

import os
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import logging

logger = logging.getLogger(__name__)

class DataHandler:
    """Handles all data operations for the LIANA Results Explorer."""
    
    def __init__(self):
        """Initialize the data handler."""
        self.results = {}
        self.conditions = []
        self.current_data_dir = None
        self.cell_types = []
        self.interactions_summary = None
        
        # New attributes for hierarchical browsing
        self.splitting_keys = []
        self.current_splitting_key = None
        self.top_level_dir = None
        self.directory_structure = {}
        
    def load_liana_results(self, data_dir: str) -> Tuple[Dict, List[str]]:
        """
        Load LIANA results from a directory structure.
        
        Args:
            data_dir: Path to the directory containing LIANA results
            
        Returns:
            Tuple of (results_dict, conditions_list)
        """
        results = {}
        conditions = []
        
        if not os.path.exists(data_dir):
            raise FileNotFoundError(f"Directory not found: {data_dir}")
        
        # Load main results
        main_csv = os.path.join(data_dir, "liana_results.csv")
        if os.path.exists(main_csv):
            try:
                results['full'] = pd.read_csv(main_csv)
                logger.info(f"Loaded full dataset results: {len(results['full'])} interactions")
            except Exception as e:
                logger.error(f"Error loading main results: {e}")
                raise
        
        # Look for condition-specific results
        for item in os.listdir(data_dir):
            item_path = os.path.join(data_dir, item)
            if os.path.isdir(item_path) and item != "logs":
                csv_path = os.path.join(item_path, "liana_results.csv")
                if os.path.exists(csv_path):
                    try:
                        results[item] = pd.read_csv(csv_path)
                        conditions.append(item)
                        logger.info(f"Loaded {item} results: {len(results[item])} interactions")
                    except Exception as e:
                        logger.warning(f"Error loading {item} results: {e}")
                else:
                    logger.debug(f"No CSV file found at {csv_path}")
            else:
                logger.debug(f"Skipping {item} (not a directory or is logs)")
        
        # Update instance variables
        self.results = results
        self.conditions = ['full'] + conditions
        self.current_data_dir = data_dir
        self.cell_types = self._extract_cell_types(results)
        
        return results, conditions
    
    def discover_hierarchical_structure(self, top_level_dir: str) -> Dict:
        """
        Discover the hierarchical structure of LIANA results.
        
        Args:
            top_level_dir: Top-level directory containing splitting key subdirectories
            
        Returns:
            Dictionary with discovered structure
        """
        if not os.path.exists(top_level_dir):
            raise FileNotFoundError(f"Directory not found: {top_level_dir}")
        
        structure = {
            'top_level_dir': top_level_dir,
            'splitting_keys': [],
            'structure': {}
        }
        
        # Look for splitting key directories
        for item in os.listdir(top_level_dir):
            item_path = os.path.join(top_level_dir, item)
            if os.path.isdir(item_path) and item != "logs":
                # Check if this directory contains LIANA results
                nested_dir = os.path.join(item_path, item)  # e.g., HPV/HPV/
                if os.path.exists(nested_dir):
                    main_results = os.path.join(nested_dir, "liana_results.csv")
                    if os.path.exists(main_results):
                        structure['splitting_keys'].append(item)
                        structure['structure'][item] = self._analyze_splitting_key_structure(nested_dir)
                        logger.info(f"Discovered splitting key: {item}")
        
        # Update instance variables
        self.top_level_dir = top_level_dir
        self.splitting_keys = structure['splitting_keys']
        self.directory_structure = structure
        
        logger.info(f"Discovered {len(structure['splitting_keys'])} splitting keys: {structure['splitting_keys']}")
        return structure
    
    def _analyze_splitting_key_structure(self, splitting_dir: str) -> Dict:
        """
        Analyze the structure of a specific splitting key directory.
        
        Args:
            splitting_dir: Path to splitting key directory (e.g., HPV/HPV/)
            
        Returns:
            Dictionary with structure information
        """
        structure = {
            'main_results': False,
            'conditions': [],
            'total_interactions': 0
        }
        
        # Check for main results
        main_csv = os.path.join(splitting_dir, "liana_results.csv")
        if os.path.exists(main_csv):
            structure['main_results'] = True
            try:
                df = pd.read_csv(main_csv)
                structure['total_interactions'] = len(df)
            except Exception as e:
                logger.warning(f"Error reading main results: {e}")
        
        # Look for condition subdirectories
        for item in os.listdir(splitting_dir):
            item_path = os.path.join(splitting_dir, item)
            if os.path.isdir(item_path) and item != "logs":
                csv_path = os.path.join(item_path, "liana_results.csv")
                if os.path.exists(csv_path):
                    structure['conditions'].append(item)
                    logger.debug(f"Found condition: {item}")
        
        return structure
    
    def load_splitting_key_data(self, splitting_key: str) -> Tuple[Dict, List[str]]:
        """
        Load data for a specific splitting key.
        
        Args:
            splitting_key: Name of the splitting key (e.g., 'HPV', 'classification')
            
        Returns:
            Tuple of (results_dict, conditions_list)
        """
        if not self.top_level_dir or splitting_key not in self.splitting_keys:
            raise ValueError(f"Splitting key '{splitting_key}' not found or not initialized")
        
        splitting_dir = os.path.join(self.top_level_dir, splitting_key, splitting_key)
        results, conditions = self.load_liana_results(splitting_dir)
        
        # Update current state
        self.current_splitting_key = splitting_key
        self.current_data_dir = splitting_dir
        
        return results, conditions
    
    def get_available_contrasts(self, splitting_key: str) -> List[str]:
        """
        Get available contrasts for a specific splitting key.
        
        Args:
            splitting_key: Name of the splitting key
            
        Returns:
            List of available contrast names
        """
        if not self.top_level_dir or splitting_key not in self.splitting_keys:
            return []
        
        splitting_dir = os.path.join(self.top_level_dir, splitting_key, splitting_key)
        if not os.path.exists(splitting_dir):
            return []
        
        contrasts = ['full']  # Always include full results if available
        
        # Check for condition subdirectories
        for item in os.listdir(splitting_dir):
            item_path = os.path.join(splitting_dir, item)
            if os.path.isdir(item_path) and item != "logs":
                csv_path = os.path.join(item_path, "liana_results.csv")
                if os.path.exists(csv_path):
                    contrasts.append(item)
        
        return contrasts
    
    def get_splitting_key_summary(self) -> Dict:
        """
        Get summary information for all discovered splitting keys.
        
        Returns:
            Dictionary with summary information
        """
        if not self.directory_structure:
            return {}
        
        summary = {}
        for key in self.splitting_keys:
            key_info = self.directory_structure['structure'].get(key, {})
            summary[key] = {
                'has_main_results': key_info.get('main_results', False),
                'num_conditions': len(key_info.get('conditions', [])),
                'conditions': key_info.get('conditions', []),
                'total_interactions': key_info.get('total_interactions', 0)
            }
        
        return summary
    
    def _extract_cell_types(self, results: Dict) -> List[str]:
        """Extract unique cell types from results."""
        all_cell_types = set()
        for df in results.values():
            if not df.empty:
                if 'source' in df.columns:
                    all_cell_types.update(df['source'].unique())
                if 'target' in df.columns:
                    all_cell_types.update(df['target'].unique())
        return sorted(list(all_cell_types))
    
    def filter_interactions(self, 
                          df: pd.DataFrame,
                          logfc_threshold: float = 0.5,
                          lrscore_threshold: float = 0.9,
                          specificity_threshold: float = 0.05,
                          min_source_cells: int = 0,
                          min_target_cells: int = 0,
                          source_types: Optional[List[str]] = None,
                          target_types: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Filter interactions based on thresholds and cell types.
        
        Args:
            df: DataFrame to filter
            logfc_threshold: Minimum LogFC threshold
            lrscore_threshold: Minimum LRscore threshold
            specificity_threshold: Maximum specificity rank threshold
            source_types: List of source cell types to include
            target_types: List of target cell types to include
            
        Returns:
            Filtered DataFrame
        """
        if df.empty:
            return df
        
        logger.debug(f"Filtering {len(df)} rows with thresholds: logfc={logfc_threshold}, lrscore={lrscore_threshold}, specificity={specificity_threshold}")
        
        # Apply thresholds
        filtered_df = df.copy()
        
        # LogFC filter - use NOTTTT absolute value since logfc can be negative, but is NOT enriched then, but depleted
        if 'lr_logfc' in df.columns:
            filtered_df = filtered_df[filtered_df['lr_logfc'] > logfc_threshold]
        
        # LRscore filter
        if 'lrscore' in df.columns:
            filtered_df = filtered_df[filtered_df['lrscore'] > lrscore_threshold]
        
        # Specificity filter
        if 'specificity_rank' in df.columns:
            filtered_df = filtered_df[filtered_df['specificity_rank'] < specificity_threshold]
            
        # Source cell-count filter
        if 'source_n_cells' in filtered_df.columns:
            filtered_df = filtered_df[filtered_df['source_n_cells'] >= min_source_cells]
        else:
            logger.debug("Column 'source_n_cells' not found; skipping source cell-count filter")

        # Target cell-count filter
        if 'target_n_cells' in filtered_df.columns:
            filtered_df = filtered_df[filtered_df['target_n_cells'] >= min_target_cells]
        else:
            logger.debug("Column 'target_n_cells' not found; skipping target cell-count filter")
        
        # Filter by cell types if specified
        if source_types and 'source' in filtered_df.columns:
            filtered_df = filtered_df[filtered_df['source'].isin(source_types)]
        
        if target_types and 'target' in filtered_df.columns:
            filtered_df = filtered_df[filtered_df['target'].isin(target_types)]
        
        logger.debug(f"After filtering: {len(filtered_df)} rows remaining")
        return filtered_df
    
    def get_summary_stats(self, contrast: str = 'full') -> Dict:
        """
        Get summary statistics for a given contrast.
        
        Args:
            contrast: Name of the contrast/condition
            
        Returns:
            Dictionary with summary statistics
        """
        if contrast not in self.results:
            return {}
        
        df = self.results[contrast]
        
        if df.empty:
            return {'total_interactions': 0}
        
        stats = {
            'total_interactions': len(df),
            'unique_pairs': len(df[['source', 'target']].drop_duplicates()) if all(col in df.columns for col in ['source', 'target']) else 0,
            'unique_ligands': len(df['ligand_complex'].unique()) if 'ligand_complex' in df.columns else 0,
            'unique_receptors': len(df['receptor_complex'].unique()) if 'receptor_complex' in df.columns else 0,
            'unique_sources': len(df['source'].unique()) if 'source' in df.columns else 0,
            'unique_targets': len(df['target'].unique()) if 'target' in df.columns else 0
        }
        
        # Add score statistics if available
        if 'lrscore' in df.columns:
            stats['mean_lrscore'] = df['lrscore'].mean()
            stats['median_lrscore'] = df['lrscore'].median()
        
        if 'lr_logfc' in df.columns:
            stats['mean_logfc'] = df['lr_logfc'].mean()
            stats['median_logfc'] = df['lr_logfc'].median()
        
        return stats
    
    def get_top_interactions(self, 
                           df: pd.DataFrame, 
                           n: int = 20, 
                           sort_by: str = 'magnitude_rank') -> pd.DataFrame:
        """
        Get top N interactions based on a specified metric.
        
        Args:
            df: DataFrame to sort
            n: Number of top interactions to return
            sort_by: Column to sort by
            
        Returns:
            DataFrame with top interactions
        """
        if df.empty or sort_by not in df.columns:
            return df.head(n)
        
        if sort_by in ['magnitude_rank', 'specificity_rank']:
            # Lower values are better for ranks
            return df.nsmallest(n, sort_by)
        else:
            # Higher values are better for scores
            return df.nlargest(n, sort_by)
    
    def prepare_comparison_data(self, 
                              source_cell: str, 
                              target_cell: str) -> pd.DataFrame:
        """
        Prepare data for cross-condition comparison.
        
        Args:
            source_cell: Source cell type
            target_cell: Target cell type
            
        Returns:
            Combined DataFrame with condition labels
        """
        if len(self.results) < 2:
            return pd.DataFrame()
        
        combined_data = []
        for condition, df in self.results.items():
            if condition != 'full' and not df.empty:
                condition_data = df[
                    (df['source'] == source_cell) & 
                    (df['target'] == target_cell)
                ].copy()
                condition_data['condition'] = condition
                combined_data.append(condition_data)
        
        if not combined_data:
            return pd.DataFrame()
        
        return pd.concat(combined_data, ignore_index=True)
    
    def validate_data_format(self, df: pd.DataFrame) -> Tuple[bool, List[str]]:
        """
        Validate that a DataFrame has the expected LIANA format.
        
        Args:
            df: DataFrame to validate
            
        Returns:
            Tuple of (is_valid, list_of_issues)
        """
        required_columns = ['source', 'target', 'ligand_complex', 'receptor_complex']
        optional_columns = ['lrscore', 'lr_logfc', 'lr_means', 'specificity_rank', 'magnitude_rank']
        
        issues = []
        
        # Check required columns
        missing_required = [col for col in required_columns if col not in df.columns]
        if missing_required:
            issues.append(f"Missing required columns: {missing_required}")
        
        # Check data types and ranges
        if 'lrscore' in df.columns:
            if not pd.api.types.is_numeric_dtype(df['lrscore']):
                issues.append("lrscore column should be numeric")
            elif df['lrscore'].min() < 0 or df['lrscore'].max() > 1:
                issues.append("lrscore values should be between 0 and 1")
        
        if 'specificity_rank' in df.columns:
            if not pd.api.types.is_numeric_dtype(df['specificity_rank']):
                issues.append("specificity_rank column should be numeric")
            elif df['specificity_rank'].min() < 0 or df['specificity_rank'].max() > 1:
                issues.append("specificity_rank values should be between 0 and 1")
        
        # Check for empty values in required columns
        for col in required_columns:
            if col in df.columns and df[col].isna().any():
                issues.append(f"Column {col} contains empty values")
        
        is_valid = len(issues) == 0
        return is_valid, issues
    
    def get_interaction_matrix(self, 
                             df: pd.DataFrame, 
                             value_col: str = 'lrscore',
                             aggfunc: str = 'median') -> pd.DataFrame:
        """
        Create an interaction matrix for heatmap visualization.
        
        Args:
            df: DataFrame with interaction data
            value_col: Column to use for values
            aggfunc: Aggregation function ('median', 'mean', 'sum', 'count')
            
        Returns:
            Pivot table with source as rows and target as columns
        """
        if df.empty or value_col not in df.columns:
            return pd.DataFrame()
        
        try:
            matrix = df.pivot_table(
                index='source',
                columns='target',
                values=value_col,
                aggfunc=aggfunc
            ).fillna(0)
            return matrix
        except Exception as e:
            logger.error(f"Error creating interaction matrix: {e}")
            return pd.DataFrame()
    
    def export_filtered_data(self, 
                           df: pd.DataFrame, 
                           output_path: str,
                           format: str = 'csv') -> bool:
        """
        Export filtered data to file.
        
        Args:
            df: DataFrame to export
            output_path: Path for output file
            format: Export format ('csv', 'excel', 'json')
            
        Returns:
            Success status
        """
        try:
            if format.lower() == 'csv':
                df.to_csv(output_path, index=False)
            elif format.lower() == 'excel':
                df.to_excel(output_path, index=False)
            elif format.lower() == 'json':
                df.to_json(output_path, orient='records')
            else:
                logger.error(f"Unsupported export format: {format}")
                return False
            
            logger.info(f"Data exported to {output_path}")
            return True
        except Exception as e:
            logger.error(f"Error exporting data: {e}")
            return False
