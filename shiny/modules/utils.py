#!/usr/bin/env python3
"""
Utilities Module for LIANA Results Explorer
Contains helper functions, logging setup, and common utilities.
"""

import os
import logging
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import tempfile
import json

def setup_logging(level: str = "INFO", 
                 log_file: Optional[str] = None) -> logging.Logger:
    """
    Set up logging configuration for the application.
    
    Args:
        level: Logging level ('DEBUG', 'INFO', 'WARNING', 'ERROR')
        log_file: Optional log file path
        
    Returns:
        Configured logger instance
    """
    log_level = getattr(logging, level.upper(), logging.INFO)
    
    # Configure logging format
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Set up handlers
    handlers = []
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    handlers.append(console_handler)
    
    # File handler if specified
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        handlers.append(file_handler)
    
    # Configure root logger
    logging.basicConfig(
        level=log_level,
        handlers=handlers,
        force=True
    )
    
    logger = logging.getLogger(__name__)
    logger.info(f"Logging initialized at {level} level")
    
    return logger

def validate_data_format(df: pd.DataFrame) -> Tuple[bool, List[str]]:
    """
    Validate that a DataFrame has the expected LIANA format.
    
    Args:
        df: DataFrame to validate
        
    Returns:
        Tuple of (is_valid, list_of_issues)
    """
    required_columns = ['source', 'target', 'ligand_complex', 'receptor_complex']
    recommended_columns = ['lrscore', 'lr_logfc', 'lr_means', 'specificity_rank', 'magnitude_rank']
    
    issues = []
    
    if df.empty:
        issues.append("DataFrame is empty")
        return False, issues
    
    # Check required columns
    missing_required = [col for col in required_columns if col not in df.columns]
    if missing_required:
        issues.append(f"Missing required columns: {missing_required}")
    
    # Check recommended columns
    missing_recommended = [col for col in recommended_columns if col not in df.columns]
    if missing_recommended:
        issues.append(f"Missing recommended columns (may affect functionality): {missing_recommended}")
    
    # Check data types and ranges for specific columns
    if 'lrscore' in df.columns:
        if not pd.api.types.is_numeric_dtype(df['lrscore']):
            issues.append("lrscore column should be numeric")
        else:
            score_range = (df['lrscore'].min(), df['lrscore'].max())
            if score_range[0] < 0 or score_range[1] > 1:
                issues.append(f"lrscore values should be between 0 and 1 (found range: {score_range})")
    
    if 'specificity_rank' in df.columns:
        if not pd.api.types.is_numeric_dtype(df['specificity_rank']):
            issues.append("specificity_rank column should be numeric")
        else:
            rank_range = (df['specificity_rank'].min(), df['specificity_rank'].max())
            if rank_range[0] < 0 or rank_range[1] > 1:
                issues.append(f"specificity_rank values should be between 0 and 1 (found range: {rank_range})")
    
    # Check for empty values in required columns
    for col in required_columns:
        if col in df.columns and df[col].isna().any():
            na_count = df[col].isna().sum()
            issues.append(f"Column '{col}' contains {na_count} empty values")
    
    # Check for duplicate rows
    if df.duplicated().any():
        dup_count = df.duplicated().sum()
        issues.append(f"Found {dup_count} duplicate rows")
    
    is_valid = len([issue for issue in issues if not issue.startswith("Missing recommended")]) == 0
    return is_valid, issues

def sanitize_filename(filename: str) -> str:
    """
    Sanitize a filename to be safe for file system use.
    
    Args:
        filename: Original filename
        
    Returns:
        Sanitized filename
    """
    # Remove or replace problematic characters
    invalid_chars = '<>:"/\\|?*'
    for char in invalid_chars:
        filename = filename.replace(char, '_')
    
    # Remove leading/trailing spaces and dots
    filename = filename.strip(' .')
    
    # Limit length
    if len(filename) > 200:
        filename = filename[:200]
    
    return filename

def format_number(value: float, 
                 precision: int = 3,
                 use_scientific: bool = False) -> str:
    """
    Format a number for display.
    
    Args:
        value: Number to format
        precision: Number of decimal places
        use_scientific: Whether to use scientific notation for very small/large numbers
        
    Returns:
        Formatted string
    """
    if pd.isna(value):
        return "N/A"
    
    if use_scientific and (abs(value) < 0.001 or abs(value) > 1000):
        return f"{value:.{precision}e}"
    else:
        return f"{value:.{precision}f}"

def calculate_interaction_stats(df: pd.DataFrame) -> Dict[str, Any]:
    """
    Calculate comprehensive statistics for interaction data.
    
    Args:
        df: DataFrame with interaction data
        
    Returns:
        Dictionary with various statistics
    """
    if df.empty:
        return {}
    
    stats = {
        'basic': {
            'total_interactions': len(df),
            'unique_sources': df['source'].nunique() if 'source' in df.columns else 0,
            'unique_targets': df['target'].nunique() if 'target' in df.columns else 0,
            'unique_ligands': df['ligand_complex'].nunique() if 'ligand_complex' in df.columns else 0,
            'unique_receptors': df['receptor_complex'].nunique() if 'receptor_complex' in df.columns else 0,
        },
        'score_stats': {},
        'top_interactions': {},
        'cell_type_stats': {}
    }
    
    # Score statistics
    score_columns = ['lrscore', 'lr_logfc', 'lr_means', 'specificity_rank']
    for col in score_columns:
        if col in df.columns:
            stats['score_stats'][col] = {
                'mean': df[col].mean(),
                'median': df[col].median(),
                'std': df[col].std(),
                'min': df[col].min(),
                'max': df[col].max()
            }
    
    # Top interactions by different metrics
    if 'magnitude_rank' in df.columns:
        top_magnitude = df.nsmallest(10, 'magnitude_rank')[['source', 'target', 'ligand_complex', 'receptor_complex', 'magnitude_rank']]
        stats['top_interactions']['by_magnitude'] = top_magnitude.to_dict('records')
    
    if 'lrscore' in df.columns:
        top_score = df.nlargest(10, 'lrscore')[['source', 'target', 'ligand_complex', 'receptor_complex', 'lrscore']]
        stats['top_interactions']['by_lrscore'] = top_score.to_dict('records')
    
    # Cell type interaction statistics
    if all(col in df.columns for col in ['source', 'target']):
        source_counts = df['source'].value_counts().head(10)
        target_counts = df['target'].value_counts().head(10)
        
        stats['cell_type_stats']['top_sources'] = source_counts.to_dict()
        stats['cell_type_stats']['top_targets'] = target_counts.to_dict()
        
        # Most common cell type pairs
        pair_counts = df.groupby(['source', 'target']).size().nlargest(10)
        stats['cell_type_stats']['top_pairs'] = {f"{pair[0]} → {pair[1]}": count for pair, count in pair_counts.items()}
    
    return stats

def create_export_summary(df: pd.DataFrame, 
                         filters_applied: Dict[str, Any]) -> str:
    """
    Create a summary text for data export.
    
    Args:
        df: DataFrame being exported
        filters_applied: Dictionary of filters that were applied
        
    Returns:
        Summary text
    """
    summary_lines = [
        "LIANA Results Export Summary",
        "=" * 30,
        f"Export Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"Total Interactions: {len(df):,}",
        ""
    ]
    
    if filters_applied:
        summary_lines.append("Filters Applied:")
        for filter_name, filter_value in filters_applied.items():
            if filter_value is not None:
                summary_lines.append(f"  - {filter_name}: {filter_value}")
        summary_lines.append("")
    
    # Basic statistics
    if not df.empty:
        summary_lines.extend([
            "Dataset Statistics:",
            f"  - Unique Sources: {df['source'].nunique() if 'source' in df.columns else 'N/A'}",
            f"  - Unique Targets: {df['target'].nunique() if 'target' in df.columns else 'N/A'}",
            f"  - Unique L-R Pairs: {len(df[['ligand_complex', 'receptor_complex']].drop_duplicates()) if all(col in df.columns for col in ['ligand_complex', 'receptor_complex']) else 'N/A'}",
            ""
        ])
        
        # Score ranges
        score_cols = ['lrscore', 'lr_logfc', 'specificity_rank']
        for col in score_cols:
            if col in df.columns:
                summary_lines.append(f"  - {col}: {df[col].min():.3f} - {df[col].max():.3f} (median: {df[col].median():.3f})")
    
    return "\n".join(summary_lines)

def detect_data_directory_structure(data_dir: str) -> Dict[str, Any]:
    """
    Analyze the structure of a LIANA results directory.
    
    Args:
        data_dir: Path to the data directory
        
    Returns:
        Dictionary with directory structure information
    """
    structure = {
        'main_results': False,
        'conditions': [],
        'total_files': 0,
        'structure_valid': False,
        'issues': []
    }
    
    if not os.path.exists(data_dir):
        structure['issues'].append(f"Directory does not exist: {data_dir}")
        return structure
    
    try:
        # Check for main results file
        main_csv = os.path.join(data_dir, "liana_results.csv")
        if os.path.exists(main_csv):
            structure['main_results'] = True
            structure['total_files'] += 1
        
        # Look for condition subdirectories
        for item in os.listdir(data_dir):
            item_path = os.path.join(data_dir, item)
            if os.path.isdir(item_path) and item != "logs":
                csv_path = os.path.join(item_path, "liana_results.csv")
                if os.path.exists(csv_path):
                    structure['conditions'].append(item)
                    structure['total_files'] += 1
        
        # Determine if structure is valid
        if structure['main_results'] or structure['conditions']:
            structure['structure_valid'] = True
        else:
            structure['issues'].append("No LIANA results files found")
        
    except Exception as e:
        structure['issues'].append(f"Error analyzing directory: {str(e)}")
    
    return structure

def create_temp_export_file(df: pd.DataFrame, 
                          format: str = 'csv',
                          summary: Optional[str] = None) -> str:
    """
    Create a temporary file for data export.
    
    Args:
        df: DataFrame to export
        format: Export format ('csv', 'excel', 'json')
        summary: Optional summary text to include
        
    Returns:
        Path to temporary file
    """
    temp_dir = tempfile.gettempdir()
    timestamp = pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')
    
    if format.lower() == 'csv':
        temp_file = os.path.join(temp_dir, f"liana_export_{timestamp}.csv")
        
        # Write summary as comments if provided
        with open(temp_file, 'w') as f:
            if summary:
                for line in summary.split('\n'):
                    f.write(f"# {line}\n")
                f.write("#\n")
            
            # Write the data
            df.to_csv(f, index=False)
    
    elif format.lower() == 'excel':
        temp_file = os.path.join(temp_dir, f"liana_export_{timestamp}.xlsx")
        
        with pd.ExcelWriter(temp_file, engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name='Interactions', index=False)
            
            if summary:
                # Create a summary sheet
                summary_df = pd.DataFrame({'Summary': summary.split('\n')})
                summary_df.to_excel(writer, sheet_name='Export_Summary', index=False)
    
    elif format.lower() == 'json':
        temp_file = os.path.join(temp_dir, f"liana_export_{timestamp}.json")
        
        export_data = {
            'export_info': {
                'timestamp': pd.Timestamp.now().isoformat(),
                'total_interactions': len(df)
            },
            'interactions': df.to_dict('records')
        }
        
        if summary:
            export_data['summary'] = summary
        
        with open(temp_file, 'w') as f:
            json.dump(export_data, f, indent=2, default=str)
    
    else:
        raise ValueError(f"Unsupported export format: {format}")
    
    return temp_file

def validate_input_parameters(params: Dict[str, Any]) -> Tuple[bool, List[str]]:
    """
    Validate input parameters for the analysis.
    
    Args:
        params: Dictionary of parameters to validate
        
    Returns:
        Tuple of (is_valid, list_of_issues)
    """
    issues = []
    
    # Validate threshold parameters
    if 'logfc_threshold' in params:
        if not isinstance(params['logfc_threshold'], (int, float)) or params['logfc_threshold'] < 0:
            issues.append("logfc_threshold must be a non-negative number")
    
    if 'lrscore_threshold' in params:
        if not isinstance(params['lrscore_threshold'], (int, float)) or not (0 <= params['lrscore_threshold'] <= 1):
            issues.append("lrscore_threshold must be between 0 and 1")
    
    if 'specificity_threshold' in params:
        if not isinstance(params['specificity_threshold'], (int, float)) or not (0 <= params['specificity_threshold'] <= 1):
            issues.append("specificity_threshold must be between 0 and 1")
    
    # Validate data directory
    if 'data_dir' in params:
        if not isinstance(params['data_dir'], str) or not params['data_dir'].strip():
            issues.append("data_dir must be a non-empty string")
        elif not os.path.exists(params['data_dir']):
            issues.append(f"data_dir does not exist: {params['data_dir']}")
    
    # Validate numeric parameters
    numeric_params = ['top_n_interactions', 'table_rows', 'network_top_n']
    for param in numeric_params:
        if param in params:
            if not isinstance(params[param], int) or params[param] <= 0:
                issues.append(f"{param} must be a positive integer")
    
    is_valid = len(issues) == 0
    return is_valid, issues

