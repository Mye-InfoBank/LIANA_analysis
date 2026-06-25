#!/usr/bin/env python3
"""
Server Logic Module for LIANA Results Explorer
Contains all server-side reactive logic and event handlers.
"""

import pandas as pd
import asyncio
from shiny import Inputs, Outputs, Session, render, reactive, ui
from shinywidgets import render_plotly
from typing import Dict, List, Optional
import logging

from .data_handler import DataHandler
from .visualizations import (
    create_network_plot, create_cell_count_network_plot,
    create_heatmap_plot, create_dotplot, 
    create_scatter_comparison, create_volcano_plot, create_summary_barplot,
    create_structure_overview_plot, create_cross_analysis_comparison
)
from .utils import (
    validate_input_parameters, calculate_interaction_stats,
    create_export_summary, create_temp_export_file
)

logger = logging.getLogger(__name__)

def create_server_function(data_handler: DataHandler):
    """
    Create the server function for the Shiny app.
    
    Args:
        data_handler: Instance of DataHandler for data operations
        
    Returns:
        Server function
    """
    
    def server(input: Inputs, output: Outputs, session: Session):
        
        # Reactive values for storing app state
        app_state = reactive.Value({
            'data_discovered': False,
            'data_loaded': False,
            'current_splitting_key': None,
            'current_contrast': None,
            'filtered_data': pd.DataFrame(),
            'summary_stats': {},
            'splitting_keys_summary': {}
        })
        
        # Trigger for auto-discovery
        startup_trigger = reactive.Value(0)
        
        # Track user's desired contrast selection across loads
        desired_contrast = reactive.Value(None)
        

        

        
        @reactive.effect
        @reactive.event(input.discover_data)
        def discover_data():
            """Discover hierarchical data structure when button is clicked."""
            # Skip if data directory is set via CLI (button is hidden)
            if data_handler.top_level_dir:
                return
                
            top_level_dir = input.top_level_dir()
            
            # Validate input
            params = {'data_dir': top_level_dir}
            is_valid, issues = validate_input_parameters(params)
            
            if not is_valid:
                for issue in issues:
                    ui.notification_show(issue, type="error")
                return
            
            try:
                # Discover data structure
                structure = data_handler.discover_hierarchical_structure(top_level_dir)
                
                if not structure['splitting_keys']:
                    ui.notification_show("No LIANA analysis directories found", type="error")
                    return
                
                # Update UI choices
                ui.update_select("splitting_key", choices=structure['splitting_keys'], selected=None)
                
                # Get summary information
                summary = data_handler.get_splitting_key_summary()
                
                # Update app state
                current_state = app_state.get()
                current_state.update({
                    'data_discovered': True,
                    'splitting_keys_summary': summary
                })
                app_state.set(current_state)
                
                ui.notification_show(f"Discovered {len(structure['splitting_keys'])} analysis types", type="success")
                logger.info(f"Data structure discovered from {top_level_dir}")
                
            except Exception as e:
                error_msg = f"Error discovering data structure: {str(e)}"
                ui.notification_show(error_msg, type="error")
                logger.error(error_msg)
        
        @reactive.extended_task
        async def load_splitting_key_task(splitting_key: str, data_discovered: bool):
            """Load data for selected splitting key as an extended task."""
            if not splitting_key or not data_discovered:
                raise ValueError("Please discover data structure first")
            
            # Simulate some async work
            await asyncio.sleep(0.1)
            
            # Load data for the selected splitting key
            results, conditions = data_handler.load_splitting_key_data(splitting_key)
            
            if not results:
                raise ValueError(f"No LIANA results found for {splitting_key}")
            
            return {
                'splitting_key': splitting_key,
                'results': results,
                'conditions': conditions
            }
        
        @reactive.effect
        @reactive.event(input.splitting_key)
        def auto_load_on_splitting_key_change():
            """Automatically load data whenever the analysis type changes."""
            splitting_key = input.splitting_key()
            if not splitting_key:
                return
            data_discovered = app_state.get()['data_discovered']
            logger.info(f"Auto-loading data for splitting key: {splitting_key}")
            load_splitting_key_task(splitting_key, data_discovered)
        
        @reactive.effect
        def handle_load_splitting_key_result():
            """Handle the result of the load splitting key task."""
            if load_splitting_key_task.result() is not None:
                try:
                    result = load_splitting_key_task.result()
                    splitting_key = result['splitting_key']
                    results = result['results']
                    conditions = result['conditions']
                    
                    logger.info(f"Loading results for {splitting_key}: {list(results.keys())}")
                    logger.info(f"Available conditions: {conditions}")
                    
                    # Update UI choices using filesystem-based discovery
                    available_contrasts = data_handler.get_available_contrasts(splitting_key)
                    # Try to preserve user's desired contrast
                    preferred = desired_contrast.get()
                    if preferred not in (available_contrasts or []):
                        preferred = 'full' if 'full' in available_contrasts else (available_contrasts[0] if available_contrasts else None)
                    logger.info(f"Available contrasts for {splitting_key}: {available_contrasts}")
                    logger.info(f"Setting contrast choices: {available_contrasts}, selected: {preferred}")
                    ui.update_select("contrast", choices=available_contrasts, selected=preferred)
                    
                    ui.update_selectize("source_types", choices=data_handler.cell_types)
                    ui.update_selectize("target_types", choices=data_handler.cell_types)
                    ui.update_select("comparison_source", choices=data_handler.cell_types, 
                                   selected=data_handler.cell_types[0] if data_handler.cell_types else None)
                    ui.update_select("comparison_target", choices=data_handler.cell_types,
                                   selected=data_handler.cell_types[1] if len(data_handler.cell_types) > 1 else None)
                    
                    # Update cross-analysis controls
                    ui.update_select("cross_analysis_source", choices=data_handler.cell_types)
                    ui.update_select("cross_analysis_target", choices=data_handler.cell_types)
                    ui.update_selectize("compare_splitting_keys", choices=data_handler.splitting_keys, 
                                      selected=data_handler.splitting_keys[:3])  # Select first 3 by default
                    
                    # Update app state
                    current_state = app_state.get()
                    current_state.update({
                        'data_loaded': True,
                        'current_splitting_key': splitting_key,
                        'current_contrast': preferred
                    })
                    app_state.set(current_state)
                    
                    ui.notification_show(f"✅ Successfully loaded {splitting_key} analysis with {len(results)} datasets", type="success")
                    logger.info(f"Splitting key {splitting_key} loaded successfully")
                    
                except Exception as e:
                    error_msg = f"❌ Error loading splitting key data: {str(e)}"
                    ui.notification_show(error_msg, type="error")
                    logger.error(error_msg)
        
        @reactive.calc
        def get_filtered_data():
            """Get filtered data based on current inputs."""
            if not app_state.get()['data_loaded']:
                logger.debug("Data not loaded yet")
                return pd.DataFrame()
            
            # Guard against stale state: ensure handler's current key matches app state
            state = app_state.get()
            if data_handler.current_splitting_key != state.get('current_splitting_key'):
                logger.debug("Stale state detected: data_handler/current key mismatch")
                return pd.DataFrame()
            
            current_contrast = input.contrast()
            logger.debug(f"Current contrast: {current_contrast}")
            logger.debug(f"Available results keys: {list(data_handler.results.keys())}")
            logger.debug(f"Data handler results: {[(k, len(v) if hasattr(v, '__len__') else 'N/A') for k, v in data_handler.results.items()]}")
            
            if not current_contrast or current_contrast not in data_handler.results:
                logger.warning(f"Contrast '{current_contrast}' not found in results")
                logger.warning(f"Available contrasts: {list(data_handler.results.keys())}")
                return pd.DataFrame()
            
            df = data_handler.results[current_contrast]
            logger.debug(f"Loaded {len(df)} rows for contrast '{current_contrast}'")
            logger.debug(f"DataFrame columns: {list(df.columns)}")
            
            # Apply filters
            filtered_df = data_handler.filter_interactions(
                df,
                logfc_threshold=input.logfc_threshold() if input.logfc_threshold() is not None else 0.0,
                lrscore_threshold=input.lrscore_threshold() if input.lrscore_threshold() is not None else 0.0,
                specificity_threshold=input.specificity_threshold() if input.specificity_threshold() is not None else 0.1,
                min_source_cells=input.min_source_cells() if input.min_source_cells() is not None else 0,
                min_target_cells=input.min_target_cells() if input.min_target_cells() is not None else 0,
                source_types=input.source_types() or None,
                target_types=input.target_types() or None
            )
            logger.debug(f"After filtering: {len(filtered_df)} rows")
            
            # Update app state
            current_state = app_state.get()
            current_state['filtered_data'] = filtered_df
            current_state['current_contrast'] = current_contrast
            app_state.set(current_state)
            
            return filtered_df
        
        def get_filtered_data_no_state_update():
            """Get filtered data without modifying app_state.

            Use this for summary/helper outputs to avoid reactive side effects.
            """
            if not app_state.get()['data_loaded']:
                return pd.DataFrame()

            state = app_state.get()
            if data_handler.current_splitting_key != state.get('current_splitting_key'):
                return pd.DataFrame()

            current_contrast = input.contrast()

            if not current_contrast or current_contrast not in data_handler.results:
                return pd.DataFrame()

            df = data_handler.results[current_contrast]

            filtered_df = data_handler.filter_interactions(
                df,
                logfc_threshold=input.logfc_threshold() if input.logfc_threshold() is not None else 0.0,
                lrscore_threshold=input.lrscore_threshold() if input.lrscore_threshold() is not None else 0.0,
                specificity_threshold=input.specificity_threshold() if input.specificity_threshold() is not None else 0.1,
                min_source_cells=input.min_source_cells() if input.min_source_cells() is not None else 0,
                min_target_cells=input.min_target_cells() if input.min_target_cells() is not None else 0,
                source_types=input.source_types() or None,
                target_types=input.target_types() or None
            )

            return filtered_df
        
        @output
        @render_plotly
        def network_plot():
            """Render network plot."""
            df = get_filtered_data()
            
            # Apply top N filter if selected
            if input.network_options() == "top" and not df.empty:
                df = data_handler.get_top_interactions(df, n=input.network_top_n())
            
            layout = input.network_layout()
            return create_network_plot(df, layout_algorithm=layout)
        
        @output
        @render_plotly
        def network_cell_count_plot():
            """Render cell-count-aware network plot."""
            df = get_filtered_data()

            layout = input.network_layout()
            top_n_nodes = input.cell_count_network_top_n()

            return create_cell_count_network_plot(
                df,
                layout_algorithm=layout,
                top_n_nodes=top_n_nodes
            )
            
        @output
        @render_plotly
        def heatmap_plot():
            """Render heatmap plot."""
            df = get_filtered_data()
            metric = input.heatmap_metric()
            colorscale = input.colorscale()
            show_values = input.show_heatmap_values()
            return create_heatmap_plot(df, metric, colorscale, show_values)
        
        @output
        @render_plotly
        def dotplot():
            """Render dot plot."""
            df = get_filtered_data()
            top_n = input.top_n_interactions()
            color_col = input.dotplot_color()
            return create_dotplot(df, top_n, color_col=color_col)
        
        @output
        @render_plotly
        def comparison_plot():
            """Render comparison plot."""
            source_cell = input.comparison_source()
            target_cell = input.comparison_target()
            metric = input.comparison_metric()
            max_interactions = input.comparison_max_interactions()
            
            if not source_cell or not target_cell:
                return create_scatter_comparison({}, "", "")
            
            # Pass the raw results and let the viz enforce the cap consistently per condition
            return create_scatter_comparison(
                data_handler.results, source_cell, target_cell, value_col=metric,
                max_points_per_condition=max_interactions
            )
        
        @output
        @render_plotly
        def volcano_plot():
            """Render volcano plot."""
            df = get_filtered_data()
            x_col = input.volcano_x()
            y_col = input.volcano_y()
            threshold = input.volcano_threshold()
            return create_volcano_plot(df, x_col, y_col, threshold)
        
        @output
        @render_plotly
        def summary_plot():
            """Render summary statistics plot based on filtered data."""
            df = get_filtered_data_no_state_update()

            if df.empty:
                return create_summary_barplot({})

            stats = {}

            stats["filtered_interactions"] = len(df)

            if all(col in df.columns for col in ["source", "target"]):
                stats["unique_cell_type_pairs"] = len(
                    df[["source", "target"]].drop_duplicates()
                )

            if "ligand_complex" in df.columns:
                stats["unique_ligands"] = df["ligand_complex"].nunique()

            if "receptor_complex" in df.columns:
                stats["unique_receptors"] = df["receptor_complex"].nunique()

            if "source" in df.columns:
                stats["unique_sources"] = df["source"].nunique()

            if "target" in df.columns:
                stats["unique_targets"] = df["target"].nunique()

            if "lrscore" in df.columns:
                stats["median_lrscore_x100"] = df["lrscore"].median() * 100

            if "lr_logfc" in df.columns:
                stats["median_logfc_x100"] = df["lr_logfc"].median() * 100

            return create_summary_barplot(stats)
        
        @output
        @render.data_frame
        def interactions_table():
            """Render interactions data table."""
            df = get_filtered_data()
            max_rows = input.table_rows()
            search_term = input.table_search()
            
            if df.empty:
                return pd.DataFrame({"Message": ["No data available with current filters"]})
            
            # Apply search filter if provided
            if search_term and search_term.strip():
                search_term = search_term.lower()
                search_mask = (
                    df['ligand_complex'].str.lower().str.contains(search_term, na=False) |
                    df['receptor_complex'].str.lower().str.contains(search_term, na=False) |
                    df['source'].str.lower().str.contains(search_term, na=False) |
                    df['target'].str.lower().str.contains(search_term, na=False)
                )
                df = df[search_mask]
            
            # Select relevant columns for display
            display_cols = [
                'source', 'target',
                'source_n_cells', 'target_n_cells',
                'ligand_complex', 'receptor_complex',
                'lrscore', 'lr_logfc', 'lr_means',
                'specificity_rank', 'magnitude_rank'
            ]
            available_cols = [col for col in display_cols if col in df.columns]
            
            return df[available_cols].head(max_rows)
        
        @output
        @render.data_frame
        def top_pairs_table():
            """Render top cell type pairs table based on filtered data."""
            df = get_filtered_data_no_state_update()

            if df.empty:
                return pd.DataFrame({"Message": ["No data available with current filters"]})

            if not all(col in df.columns for col in ["source", "target"]):
                return pd.DataFrame({"Message": ["source/target columns missing"]})

            pair_counts = (
                df.groupby(["source", "target"])
                .size()
                .reset_index(name="filtered_interactions")
                .sort_values("filtered_interactions", ascending=False)
                .head(10)
            )

            pair_counts["cell_pair"] = (
                pair_counts["source"] + " → " + pair_counts["target"]
            )

            return pair_counts[["cell_pair", "filtered_interactions"]]
        
        @output
        @render.data_frame
        def top_lr_pairs_table():
            """Render top ligand-receptor pairs table based on filtered data."""
            df = get_filtered_data_no_state_update()

            if df.empty:
                return pd.DataFrame({"Message": ["No data available with current filters"]})

            if not all(col in df.columns for col in ["ligand_complex", "receptor_complex"]):
                return pd.DataFrame({"Message": ["ligand/receptor columns missing"]})

            lr_counts = (
                df.groupby(["ligand_complex", "receptor_complex"])
                .size()
                .reset_index(name="filtered_interactions")
                .sort_values("filtered_interactions", ascending=False)
                .head(10)
            )

            lr_counts["lr_pair"] = (
                lr_counts["ligand_complex"] + " → " + lr_counts["receptor_complex"]
            )

            return lr_counts[["lr_pair", "filtered_interactions"]]
        
        # Auto-discover data on startup when data directory is set via CLI
        @reactive.effect
        def auto_discover_on_startup():
            """Auto-discover data structure on app startup."""
            # Trigger the effect
            startup_trigger.get()
            
            if data_handler.top_level_dir:
                try:
                    # Discover data structure
                    structure = data_handler.discover_hierarchical_structure(data_handler.top_level_dir)
                    
                    if structure['splitting_keys']:
                        # Update UI choices
                        ui.update_select("splitting_key", choices=structure['splitting_keys'], selected=None)
                        
                        # Update app state
                        current_state = app_state.get()
                        current_state.update({
                            'data_discovered': True,
                            'splitting_keys_summary': data_handler.get_splitting_key_summary()
                        })
                        app_state.set(current_state)
                        
                        logger.info(f"Auto-discovered {len(structure['splitting_keys'])} analysis types")
                        
                except Exception as e:
                    logger.warning(f"Auto-discover failed: {e}")
        
        # Trigger auto-discovery on startup and auto-select first analysis type
        startup_trigger.set(1)
        
        @reactive.effect
        def auto_select_first_splitting_key():
            """Select the first available analysis type on startup to trigger loading."""
            if app_state.get()['data_discovered'] and data_handler.splitting_keys:
                current = input.splitting_key()
                if not current:
                    first = data_handler.splitting_keys[0]
                    logger.info(f"Auto-selecting first analysis type: {first}")
                    ui.update_select("splitting_key", choices=data_handler.splitting_keys, selected=first)
        
        @reactive.effect
        @reactive.event(input.export_data)
        def export_filtered_data():
            """Export filtered data when button is clicked."""
            df = get_filtered_data()
            
            if df.empty:
                ui.notification_show("No data to export", type="warning")
                return
            
            try:
                # Create export summary
                filters_applied = {
                    'contrast': input.contrast(),
                    'logfc_threshold': input.logfc_threshold(),
                    'lrscore_threshold': input.lrscore_threshold(),
                    'specificity_threshold': input.specificity_threshold(),
                    'min_source_cells': input.min_source_cells(),
                    'min_target_cells': input.min_target_cells(),
                    'source_types': input.source_types(),
                    'target_types': input.target_types()
                }
                
                summary = create_export_summary(df, filters_applied)
                
                # Create temporary file
                temp_file = create_temp_export_file(df, format='csv', summary=summary)
                
                ui.notification_show(
                    f"Data exported to: {temp_file}", 
                    type="success", 
                    duration=10
                )
                logger.info(f"Data exported to {temp_file}")
                
            except Exception as e:
                error_msg = f"Error exporting data: {str(e)}"
                ui.notification_show(error_msg, type="error")
                logger.error(error_msg)
        
        # Reactive effects for updating UI elements based on data changes
        @reactive.effect
        def update_ui_on_data_load():
            """Update UI elements when data is loaded."""
            if app_state.get()['data_loaded']:
                # Update sidebar with available cell types and conditions
                ui.update_selectize("source_types", choices=data_handler.cell_types)
                ui.update_selectize("target_types", choices=data_handler.cell_types)
                ui.update_select("comparison_source", choices=data_handler.cell_types,
                                 selected=(data_handler.cell_types[0] if data_handler.cell_types else None))
                ui.update_select("comparison_target", choices=data_handler.cell_types,
                                 selected=(data_handler.cell_types[1] if len(data_handler.cell_types) > 1 else None))
                
                # Update contrast choices based on current splitting key
                current_splitting_key = app_state.get().get('current_splitting_key')
                if current_splitting_key:
                    available_contrasts = data_handler.get_available_contrasts(current_splitting_key)
                    current_contrast = app_state.get().get('current_contrast', 'full')
                    if current_contrast not in available_contrasts:
                        current_contrast = 'full' if 'full' in available_contrasts else (available_contrasts[0] if available_contrasts else None)
                    ui.update_select("contrast", choices=available_contrasts, selected=current_contrast)
        
        @reactive.effect
        @reactive.event(input.splitting_key)
        def update_contrast_on_splitting_key_change():
            """Update contrast choices when splitting key changes."""
            splitting_key = input.splitting_key()
            if splitting_key and app_state.get()['data_discovered']:
                logger.info(f"Splitting key changed to: {splitting_key}")
                available_contrasts = data_handler.get_available_contrasts(splitting_key)
                logger.info(f"Available contrasts for {splitting_key}: {available_contrasts}")
                
                # Reset contrast to 'full' or first available
                selected_contrast = 'full' if 'full' in available_contrasts else (available_contrasts[0] if available_contrasts else None)
                ui.update_select("contrast", choices=available_contrasts, selected=selected_contrast)
                
                # Remember user's intent
                desired_contrast.set(selected_contrast)
                
                # Reset app state since we haven't loaded data for this splitting key yet
                current_state = app_state.get()
                current_state.update({
                    'data_loaded': False,
                    'current_splitting_key': splitting_key,
                    'current_contrast': selected_contrast
                })
                app_state.set(current_state)

        @reactive.effect
        @reactive.event(input.contrast)
        def validate_contrast_change():
            """Ensure selected contrast exists for the current splitting key; remember the desired contrast."""
            current_state = app_state.get()
            splitting_key = current_state.get('current_splitting_key')
            selected = input.contrast()
            if not splitting_key or not selected:
                return
            available = data_handler.get_available_contrasts(splitting_key)
            if selected not in available:
                logger.warning(f"Selected contrast '{selected}' not in available for {splitting_key}. Resetting.")
                fallback = 'full' if 'full' in available else (available[0] if available else None)
                ui.update_select("contrast", choices=available, selected=fallback)
                desired_contrast.set(fallback)
            else:
                desired_contrast.set(selected)
        
        @reactive.effect
        def validate_inputs():
            """Validate inputs and show warnings if needed."""
            if app_state.get()['data_loaded']:
                # Validate threshold combinations
                if (input.lrscore_threshold() > 0.95 and 
                    input.logfc_threshold() > 1.5 and 
                    input.specificity_threshold() < 0.01):
                    ui.notification_show(
                        "Very strict filtering may result in no data", 
                        type="warning", 
                        duration=5
                    )
        
        # Data Explorer Tab Outputs
        @output
        @render.data_frame
        def splitting_keys_table():
            """Render table of discovered splitting keys."""
            summary = app_state.get().get('splitting_keys_summary', {})
            
            if not summary:
                return pd.DataFrame({"Status": ["🔍 Discovering data structure..."]})
            
            # Create summary table
            rows = []
            for key, info in summary.items():
                rows.append({
                    'Analysis Type': key,
                    'Has Main Results': '✅' if info['has_main_results'] else '❌',
                    'Number of Conditions': info['num_conditions'],
                    'Total Interactions': f"{info['total_interactions']:,}",
                    'Conditions': ', '.join(info['conditions'][:3]) + ('...' if len(info['conditions']) > 3 else '')
                })
            
            return pd.DataFrame(rows)
        
        @output
        @render_plotly
        def structure_overview_plot():
            """Render structure overview plot."""
            summary = app_state.get().get('splitting_keys_summary', {})
            return create_structure_overview_plot(summary)
        
        @output
        @render_plotly
        def cross_analysis_plot():
            """Render cross-analysis comparison plot."""
            source_cell = input.cross_analysis_source()
            target_cell = input.cross_analysis_target()
            compare_keys = input.compare_splitting_keys()
            
            if not source_cell or not target_cell or not compare_keys:
                return create_cross_analysis_comparison(pd.DataFrame(), "", "")
            
            try:
                combined_df = data_handler.compare_across_splitting_keys(
                    source_cell, target_cell, compare_keys
                )
                return create_cross_analysis_comparison(combined_df, source_cell, target_cell)
            except Exception as e:
                logger.error(f"Error creating cross-analysis comparison: {e}")
                return create_cross_analysis_comparison(pd.DataFrame(), source_cell, target_cell)
        
    return server
