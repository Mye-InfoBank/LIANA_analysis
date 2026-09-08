#!/usr/bin/env python3
"""
Server Logic Module for LIANA Results Explorer
Contains all server-side reactive logic and event handlers.
"""

import json
import pandas as pd
import asyncio
import plotly.graph_objects as go
from pathlib import Path
from shiny import Inputs, Outputs, Session, render, reactive, ui
from shinywidgets import render_plotly, output_widget
from typing import Dict, List, Optional
import logging

from .data_handler import DataHandler
from .visualizations import (
    create_network_plot,
    create_cell_count_network_plot,
    create_heatmap_plot,
    create_dotplot,
    create_lr_boxplot,
    create_volcano_plot,
    create_summary_barplot,
    create_structure_overview_plot,

    # Overview plots
    create_global_condition_dotplot,
    create_cell_pair_difference_heatmap,
    create_lr_difference_barplot
)
from .utils import (
    validate_input_parameters, calculate_interaction_stats,
    create_export_summary, create_temp_export_file
)

from .ai_snapshots import (
    save_ai_snapshot,
    figure_to_dict,
    cleanup_old_snapshots
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
    # Remove expired snapshots whenever the application starts.
    cleanup_old_snapshots()
    
    def server(input: Inputs, output: Outputs, session: Session):
        
        # Reactive values for storing app state
        app_state = reactive.Value({
            'data_discovered': False,
            'data_loaded': False,
            'current_splitting_key': None,
            'current_contrast': None,
            'filtered_data': pd.DataFrame(),
            'summary_stats': {},
            'splitting_keys_summary': {},
        })
        
        # Trigger for auto-discovery
        startup_trigger = reactive.Value(0)
        
        # Track user's desired contrast selection across loads
        desired_contrast = reactive.Value(None)
        
        # Dedicated reactive trigger/data for the Data Explorer count plot
        interaction_counts_state = reactive.Value({})
        
        # Last generated AI snapshot URL
        ai_snapshot_relative_url = reactive.Value(None)
        
        # Temporarily holds uploaded state while a new splitting key loads
        pending_restore_state = reactive.Value(None)
        
        
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
        @reactive.event(load_splitting_key_task.result)
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
                    
                    # =====================================================
                    # Overview: condition comparison choices
                    # =====================================================

                    available_conditions = [
                        c
                        for c in conditions
                        if c != "full"
                    ]

                    comparison_choices = {}

                    # Preferred IBD comparisons first
                    preferred_pairs = [
                        ("CD", "HC"),
                        ("UC", "HC"),
                        ("CD", "UC")
                    ]

                    used_pairs = set()

                    for condition_a, condition_b in preferred_pairs:

                        if (
                            condition_a in available_conditions
                            and condition_b in available_conditions
                        ):

                            key = (
                                f"{condition_a}"
                                f"|||"
                                f"{condition_b}"
                            )

                            comparison_choices[key] = (
                                f"{condition_a} vs "
                                f"{condition_b}"
                            )

                            used_pairs.add(
                                frozenset(
                                    [
                                        condition_a,
                                        condition_b
                                    ]
                                )
                            )

                    # Also support any other condition combinations
                    for i in range(
                        len(available_conditions)
                    ):

                        for j in range(
                            i + 1,
                            len(available_conditions)
                        ):

                            condition_a = (
                                available_conditions[i]
                            )

                            condition_b = (
                                available_conditions[j]
                            )

                            pair_key = frozenset(
                                [
                                    condition_a,
                                    condition_b
                                ]
                            )

                            if pair_key in used_pairs:
                                continue

                            key = (
                                f"{condition_a}"
                                f"|||"
                                f"{condition_b}"
                            )

                            comparison_choices[key] = (
                                f"{condition_a} vs "
                                f"{condition_b}"
                            )

                    if comparison_choices:

                        first_choice = next(
                            iter(
                                comparison_choices
                            )
                        )

                        ui.update_select(
                            "overview_condition_comparison",
                            choices=comparison_choices,
                            selected=first_choice
                        )
                        
                    # Count rows per loaded dataset/contrast
                    counts_for_plot = {
                        name: int(len(df)) if df is not None else 0
                        for name, df in results.items()
                    }

                    interaction_counts_state.set(counts_for_plot)
                    logger.info(f"Interaction counts for Data Explorer plot: {counts_for_plot}")
                    # Update app state
                    current_state = app_state.get()
                    current_state.update({
                        'data_loaded': True,
                        'current_splitting_key': splitting_key,
                        'current_contrast': preferred,
                    })
                    app_state.set(current_state)
                    
                    ui.notification_show(f"✅ Successfully loaded {splitting_key} analysis with {len(results)} datasets", type="success")
                    logger.info(f"Splitting key {splitting_key} loaded successfully")
                                        
                    # ============================================================
                    # FINISH PENDING SETTINGS RESTORE
                    # ============================================================

                    pending_state = (
                        pending_restore_state.get()
                    )

                    if pending_state:

                        wanted_key = (
                            pending_state
                            .get("data_navigation", {})
                            .get("splitting_key")
                        )

                        if wanted_key == splitting_key:

                            apply_restored_state(
                                pending_state
                            )

                            pending_restore_state.set(
                                None
                            )

                            ui.notification_show(
                                "✓ Saved dashboard state restored.",
                                type="success"
                            )
                    
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
        
        def get_overview_condition_results(
            use_sidebar_filters=False
        ):
            """
            Return condition -> DataFrame for the Overview plots.

            Excludes the pooled 'full' dataset.

            If use_sidebar_filters=True, the current sidebar filters
            are applied independently to each condition.
            """

            overview_results = {}

            for condition, df in data_handler.results.items():

                if condition == "full":
                    continue

                if df is None or df.empty:
                    continue

                if use_sidebar_filters:
                    current_df = df.copy()
                else:
                    current_df = df

                if use_sidebar_filters:

                    current_df = data_handler.filter_interactions(
                        current_df,

                        logfc_threshold=(
                            input.logfc_threshold()
                            if input.logfc_threshold() is not None
                            else 0.0
                        ),

                        lrscore_threshold=(
                            input.lrscore_threshold()
                            if input.lrscore_threshold() is not None
                            else 0.0
                        ),

                        specificity_threshold=(
                            input.specificity_threshold()
                            if input.specificity_threshold() is not None
                            else 1.0
                        ),

                        min_source_cells=(
                            input.min_source_cells()
                            if input.min_source_cells() is not None
                            else 0
                        ),

                        min_target_cells=(
                            input.min_target_cells()
                            if input.min_target_cells() is not None
                            else 0
                        ),

                        source_types=(
                            input.source_types()
                            or None
                        ),

                        target_types=(
                            input.target_types()
                            or None
                        )
                    )

                overview_results[
                    condition
                ] = current_df

            return overview_results
        
        def get_selected_overview_conditions():
            """
            Read selected condition comparison.

            Stored format:
                conditionA|||conditionB
            """

            comparison = (
                input.overview_condition_comparison()
            )

            if (
                not comparison
                or "|||" not in comparison
            ):
                return None, None

            condition_a, condition_b = (
                comparison.split(
                    "|||",
                    1
                )
            )

            return condition_a, condition_b
        
        def build_restorable_state():
            """Build a machine-readable representation of the
            Exact dashboard settings required to recreate the current state."""

            return {
                "kind": "liana_results_explorer_state",
                "schema_version": "1.0",

                "data_navigation": {
                    "splitting_key": input.splitting_key(),
                    "contrast": input.contrast()
                },

                "filters": {
                    "logfc_threshold": input.logfc_threshold(),
                    "lrscore_threshold": input.lrscore_threshold(),
                    "specificity_threshold": input.specificity_threshold(),
                    "min_source_cells": input.min_source_cells(),
                    "min_target_cells": input.min_target_cells(),
                    "source_types": list(input.source_types() or []),
                    "target_types": list(input.target_types() or [])
                },

                "overview": {
                    "interaction_top_n": input.overview_interaction_top_n(),
                    "apply_filters": input.overview_apply_filters(),
                    "condition_comparison": input.overview_condition_comparison(),
                    "cell_metric": input.overview_cell_metric(),
                    "cell_mode": input.overview_cell_mode(),
                    "cell_top_n": input.overview_cell_top_n(),
                    "lr_metric": input.overview_lr_diff_metric(),
                    "lr_top_n": input.overview_lr_diff_top_n(),
                    "lr_apply_filters": input.overview_lr_diff_apply_filters()
                },

                "network": {
                    "mode": input.network_options(),
                    "top_n": input.network_top_n(),
                    "layout": input.network_layout(),
                    "cell_count_top_n": input.cell_count_network_top_n()
                },

                "heatmap": {
                    "metric": input.heatmap_metric(),
                    "colorscale": input.colorscale(),
                    "show_values": input.show_heatmap_values()
                },

                "ligand_receptor": {
                    "top_n": input.top_n_interactions(),
                    "dotplot_color": input.dotplot_color(),
                    "boxplot_metric": input.lr_boxplot_metric(),
                    "boxplot_top_n": input.lr_boxplot_top_n()
                },

                "data_table": {
                    "rows": input.table_rows(),
                    "search": input.table_search()
                }
            }
        
        def build_ai_snapshot_payload():
            """
            Build a machine-readable representation of the
            current LIANA dashboard state.

            This includes:
            - dataset metadata
            - sidebar filters
            - plot settings
            - filtered-data summary
            - top filtered interactions
            - Data Table preview
            - current Plotly outputs
            """

            # ========================================================
            # CURRENT DATASET
            # ========================================================

            data_dir = (
                data_handler.top_level_dir
                or ""
            )

            dataset_name = (
                Path(data_dir).name
                if data_dir
                else None
            )

            (
                condition_a,
                condition_b
            ) = get_selected_overview_conditions()

            # ========================================================
            # BASE PAYLOAD
            # ========================================================

            payload = {

                "schema_version": "1.0",

                "description": (
                    "Machine-readable snapshot of the current "
                    "IBD LIANA Results Explorer state."
                ),

                "metadata": {

                    "analysis":
                        "LIANA cell-cell communication",

                    "dataset":
                        dataset_name,

                    "splitting_key":
                        input.splitting_key(),

                    "selected_contrast":
                        input.contrast(),

                    "condition_comparison": {

                        "condition_a":
                            condition_a,

                        "condition_b":
                            condition_b
                    }
                },

                # ====================================================
                # SIDEBAR FILTERS
                # ====================================================

                "filters": {

                    "logfc_threshold":
                        input.logfc_threshold(),

                    "lrscore_threshold":
                        input.lrscore_threshold(),

                    "specificity_rank_threshold":
                        input.specificity_threshold(),

                    "min_source_cells":
                        input.min_source_cells(),

                    "min_target_cells":
                        input.min_target_cells(),

                    "source_cell_types":
                        list(
                            input.source_types()
                            or []
                        ),

                    "target_cell_types":
                        list(
                            input.target_types()
                            or []
                        )
                },

                # ====================================================
                # CURRENT PLOT SETTINGS
                # ====================================================

                "plot_settings": {

                    "interaction_landscape": {

                        "top_n":
                            input.overview_interaction_top_n(),

                        "apply_sidebar_filters":
                            input.overview_apply_filters()
                    },

                    "cell_cell_changes": {

                        "metric":
                            input.overview_cell_metric(),

                        "mode":
                            input.overview_cell_mode(),

                        "top_n":
                            input.overview_cell_top_n()
                    },

                    "ligand_receptor_changes": {

                        "metric":
                            input.overview_lr_diff_metric(),

                        "top_n":
                            input.overview_lr_diff_top_n(),

                        "apply_sidebar_filters":
                            input.overview_lr_diff_apply_filters()
                    },

                    "network": {

                        "mode":
                            input.network_options(),

                        "top_n":
                            input.network_top_n(),

                        "layout":
                            input.network_layout(),

                        "cell_count_network_top_n":
                            input.cell_count_network_top_n()
                    },

                    "heatmap": {

                        "metric":
                            input.heatmap_metric(),

                        "colorscale":
                            input.colorscale(),

                        "show_values":
                            input.show_heatmap_values()
                    },

                    "ligand_receptor_dotplot": {

                        "top_n":
                            input.top_n_interactions(),

                        "color_by":
                            input.dotplot_color()
                    },

                    "ligand_receptor_boxplot": {

                        "metric":
                            input.lr_boxplot_metric(),

                        "top_n":
                            input.lr_boxplot_top_n()
                    }
                },

                # ====================================================
                # IMPORTANT INTERPRETATION NOTES FOR LLMS
                # ====================================================

                "interpretation_notes": [

                    (
                        "LRscore is used as a LIANA interaction "
                        "magnitude score; higher values indicate "
                        "stronger inferred communication."
                    ),

                    (
                        "specificity_rank is a consensus specificity "
                        "rank; lower values indicate greater specificity."
                    ),

                    (
                        "magnitude_rank is a consensus magnitude rank; "
                        "lower values indicate stronger consensus ranking."
                    ),

                    (
                        "lr_logfc is a LIANA specificity-related score "
                        "within a condition and must not be interpreted "
                        "as a formal CD-vs-HC differential log fold change."
                    ),

                    (
                        "Condition-difference plots are descriptive "
                        "comparisons of separately inferred LIANA "
                        "communication networks."
                    ),

                    (
                        "LRscore differences use only interaction "
                        "identities present in both compared conditions."
                    ),

                    (
                        "A LIANA interaction score must not automatically "
                        "be interpreted as statistical significance or causality."
                    )
                ]
            }
            
            payload["restorable_state"] = build_restorable_state()

            # ========================================================
            # FILTERED DATA
            # ========================================================

            filtered_df = (
                get_filtered_data_no_state_update()
            )

            payload["filtered_dataset"] = {
                "n_rows": 0
            }

            if (
                filtered_df is not None
                and not filtered_df.empty
            ):

                # ----------------------------------------------------
                # Summary
                # ----------------------------------------------------

                filtered_summary = {

                    "n_rows":
                        int(
                            len(filtered_df)
                        )
                }

                if all(
                    col in filtered_df.columns
                    for col in [
                        "source",
                        "target"
                    ]
                ):

                    filtered_summary[
                        "n_unique_cell_pairs"
                    ] = int(
                        filtered_df[
                            ["source", "target"]
                        ]
                        .drop_duplicates()
                        .shape[0]
                    )

                if all(
                    col in filtered_df.columns
                    for col in [
                        "ligand_complex",
                        "receptor_complex"
                    ]
                ):

                    filtered_summary[
                        "n_unique_ligand_receptor_pairs"
                    ] = int(
                        filtered_df[
                            [
                                "ligand_complex",
                                "receptor_complex"
                            ]
                        ]
                        .drop_duplicates()
                        .shape[0]
                    )

                if "lrscore" in filtered_df.columns:

                    filtered_summary[
                        "median_lrscore"
                    ] = float(
                        filtered_df[
                            "lrscore"
                        ].median()
                    )

                if "lr_logfc" in filtered_df.columns:

                    filtered_summary[
                        "median_lr_logfc"
                    ] = float(
                        filtered_df[
                            "lr_logfc"
                        ].median()
                    )

                payload[
                    "filtered_dataset"
                ] = filtered_summary

                # ----------------------------------------------------
                # Top 500 filtered interactions
                # ----------------------------------------------------

                top_df = (
                    filtered_df.copy()
                )

                if (
                    "magnitude_rank"
                    in top_df.columns
                ):

                    top_df = (
                        top_df
                        .sort_values(
                            "magnitude_rank",
                            ascending=True
                        )
                        .head(500)
                    )

                elif (
                    "lrscore"
                    in top_df.columns
                ):

                    top_df = (
                        top_df
                        .sort_values(
                            "lrscore",
                            ascending=False
                        )
                        .head(500)
                    )

                else:

                    top_df = (
                        top_df.head(500)
                    )

                snapshot_columns = [

                    "source",
                    "target",

                    "source_n_cells",
                    "target_n_cells",

                    "ligand_complex",
                    "receptor_complex",

                    "lrscore",
                    "lr_logfc",
                    "lr_means",

                    "specificity_rank",
                    "magnitude_rank"
                ]

                snapshot_columns = [
                    col
                    for col in snapshot_columns
                    if col in top_df.columns
                ]

                top_df = (
                    top_df[
                        snapshot_columns
                    ]
                    .replace(
                        [float("inf"), float("-inf")],
                        None
                    )
                )

                payload[
                    "top_filtered_interactions"
                ] = (
                    top_df
                    .where(
                        pd.notna(top_df),
                        None
                    )
                    .to_dict(
                        orient="records"
                    )
                )

                # ----------------------------------------------------
                # Same preview used conceptually by Data Table
                # ----------------------------------------------------

                table_df = (
                    filtered_df.copy()
                )

                search_term = (
                    input.table_search()
                )

                if (
                    search_term
                    and search_term.strip()
                ):

                    search_term = (
                        search_term
                        .lower()
                    )

                    search_mask = (

                        table_df[
                            "ligand_complex"
                        ]
                        .astype(str)
                        .str.lower()
                        .str.contains(
                            search_term,
                            na=False
                        )

                        |

                        table_df[
                            "receptor_complex"
                        ]
                        .astype(str)
                        .str.lower()
                        .str.contains(
                            search_term,
                            na=False
                        )

                        |

                        table_df[
                            "source"
                        ]
                        .astype(str)
                        .str.lower()
                        .str.contains(
                            search_term,
                            na=False
                        )

                        |

                        table_df[
                            "target"
                        ]
                        .astype(str)
                        .str.lower()
                        .str.contains(
                            search_term,
                            na=False
                        )
                    )

                    table_df = (
                        table_df[
                            search_mask
                        ]
                    )

                table_columns = [
                    col
                    for col in snapshot_columns
                    if col in table_df.columns
                ]

                table_limit = min(
                    input.table_rows()
                    or 100,
                    500
                )

                table_df = (
                    table_df[
                        table_columns
                    ]
                    .head(table_limit)
                )

                table_df = (
                    table_df
                    .where(
                        pd.notna(table_df),
                        None
                    )
                )

                payload[
                    "data_table_preview"
                ] = (
                    table_df.to_dict(
                        orient="records"
                    )
                )

            # ========================================================
            # CREATE EXACT CURRENT PLOT OUTPUTS
            #
            # These are machine-readable Plotly specifications.
            # ========================================================

            figures = {}

            # --------------------------------------------------------
            # Dataset interaction-count overview
            # --------------------------------------------------------

            counts = (
                interaction_counts_state.get()
            )

            if counts:

                count_df = pd.DataFrame([
                    {
                        "contrast": contrast,
                        "interactions": n
                    }
                    for contrast, n
                    in counts.items()
                ])

                figures[
                    "dataset_interaction_counts"
                ] = figure_to_dict(
                    create_structure_overview_plot(
                        count_df
                    )
                )

            # --------------------------------------------------------
            # Overview Plot 1
            # --------------------------------------------------------

            overview_results = (
                get_overview_condition_results(
                    use_sidebar_filters=(
                        input.overview_apply_filters()
                    )
                )
            )

            figures[
                "interaction_landscape"
            ] = figure_to_dict(
                create_global_condition_dotplot(
                    overview_results,
                    top_n=(
                        input.overview_interaction_top_n()
                        or 50
                    )
                )
            )

            # --------------------------------------------------------
            # Overview Plot 2
            # --------------------------------------------------------

            if (
                condition_a
                and condition_b
            ):

                cell_mode = (
                    input.overview_cell_mode()
                )

                cell_results = (
                    get_overview_condition_results(
                        use_sidebar_filters=(
                            cell_mode
                            == "filtered"
                        )
                    )
                )

                if (
                    condition_a in cell_results
                    and condition_b in cell_results
                ):

                    cell_top_n = None

                    if (
                        cell_mode
                        == "top"
                    ):

                        cell_top_n = (
                            input.overview_cell_top_n()
                        )

                    figures[
                        "cell_cell_changes"
                    ] = figure_to_dict(

                        create_cell_pair_difference_heatmap(

                            cell_results[
                                condition_a
                            ],

                            cell_results[
                                condition_b
                            ],

                            condition_a=
                                condition_a,

                            condition_b=
                                condition_b,

                            metric=(
                                input.overview_cell_metric()
                                or "lrscore"
                            ),

                            top_n=
                                cell_top_n
                        )
                    )

                # ----------------------------------------------------
                # Overview Plot 3
                # ----------------------------------------------------

                lr_results = (
                    get_overview_condition_results(
                        use_sidebar_filters=(
                            input.overview_lr_diff_apply_filters()
                        )
                    )
                )

                if (
                    condition_a in lr_results
                    and condition_b in lr_results
                ):

                    figures[
                        "ligand_receptor_changes"
                    ] = figure_to_dict(

                        create_lr_difference_barplot(

                            lr_results[
                                condition_a
                            ],

                            lr_results[
                                condition_b
                            ],

                            condition_a=
                                condition_a,

                            condition_b=
                                condition_b,

                            metric=(
                                input.overview_lr_diff_metric()
                                or "lrscore"
                            ),

                            top_n=(
                                input.overview_lr_diff_top_n()
                                or 30
                            )
                        )
                    )

            # ========================================================
            # PLOTS BASED ON CURRENT FILTERED DATA
            # ========================================================

            if (
                filtered_df is not None
                and not filtered_df.empty
            ):

                # ----------------------------------------------------
                # Network
                # ----------------------------------------------------

                network_df = (
                    filtered_df
                )

                if (
                    input.network_options()
                    == "top"
                ):

                    network_df = (
                        data_handler
                        .get_top_interactions(
                            network_df,
                            n=(
                                input.network_top_n()
                                or 50
                            )
                        )
                    )

                figures[
                    "network"
                ] = figure_to_dict(

                    create_network_plot(

                        network_df,

                        layout_algorithm=(
                            input.network_layout()
                            or "spring"
                        )
                    )
                )

                # ----------------------------------------------------
                # Cell-count network
                # ----------------------------------------------------

                figures[
                    "cell_count_network"
                ] = figure_to_dict(

                    create_cell_count_network_plot(

                        filtered_df,

                        layout_algorithm=(
                            input.network_layout()
                            or "spring"
                        ),

                        top_n_nodes=(
                            input.cell_count_network_top_n()
                            or 30
                        )
                    )
                )

                # ----------------------------------------------------
                # Heatmap
                # ----------------------------------------------------

                figures[
                    "heatmap"
                ] = figure_to_dict(

                    create_heatmap_plot(

                        filtered_df,

                        value_col=(
                            input.heatmap_metric()
                            or "interaction_count"
                        ),

                        colorscale=(
                            input.colorscale()
                            or "Blues"
                        ),

                        show_values=(
                            input.show_heatmap_values()
                        )
                    )
                )

                # ----------------------------------------------------
                # LR dotplot
                # ----------------------------------------------------

                figures[
                    "ligand_receptor_dotplot"
                ] = figure_to_dict(

                    create_dotplot(

                        filtered_df,

                        top_n=(
                            input.top_n_interactions()
                            or 20
                        ),

                        color_col=(
                            input.dotplot_color()
                            or "specificity_rank"
                        )
                    )
                )

                # ----------------------------------------------------
                # LR boxplot
                # ----------------------------------------------------

                figures[
                    "ligand_receptor_boxplot"
                ] = figure_to_dict(

                    create_lr_boxplot(

                        filtered_df,

                        metric=(
                            input.lr_boxplot_metric()
                            or "lrscore"
                        ),

                        top_n=(
                            input.lr_boxplot_top_n()
                            or 20
                        )
                    )
                )

            payload[
                "plot_outputs"
            ] = figures

            return payload
        
        @reactive.effect
        @reactive.event(
            input.create_ai_snapshot
        )
        def create_ai_snapshot():

            if not app_state.get()[
                "data_loaded"
            ]:

                ui.notification_show(
                    "Please wait until the LIANA data are loaded.",
                    type="warning"
                )

                return

            try:

                payload = (
                    build_ai_snapshot_payload()
                )

                snapshot_id = (
                    save_ai_snapshot(
                        payload
                    )
                )

                relative_url = (
                    f"ai_snapshots/"
                    f"{snapshot_id}.json"
                )

                ai_snapshot_relative_url.set(
                    relative_url
                )
                

                ui.update_action_button(
                    "create_ai_snapshot",
                    label="✓ Snapshot Created!"
                )

                ui.notification_show(
                    "AI snapshot created successfully.",
                    type="success"
                )

            except Exception as e:

                logger.exception(
                    "Failed to create AI snapshot"
                )

                ui.notification_show(
                    f"Could not create AI snapshot: {e}",
                    type="error"
                )
        
        @output
        @render.download(
            filename=lambda: (
                f"IBD_LIANA_settings_"
                f"{input.splitting_key() or 'analysis'}_"
                f"{input.contrast() or 'contrast'}.json"
            )
        )
        def download_app_settings():

            state = build_restorable_state()

            yield json.dumps(
                state,
                indent=2,
                ensure_ascii=False
            )
        
        def apply_restored_state(state):
            """Apply a previously saved dashboard state."""

            filters = state.get("filters", {})
            overview = state.get("overview", {})
            network = state.get("network", {})
            heatmap = state.get("heatmap", {})
            lr = state.get("ligand_receptor", {})
            table = state.get("data_table", {})
            navigation = state.get("data_navigation", {})

            # --------------------------------------------------------
            # Contrast
            # --------------------------------------------------------

            contrast = navigation.get("contrast")

            if contrast in data_handler.results:
                desired_contrast.set(contrast)
                ui.update_select(
                    "contrast",
                    selected=contrast
                )

            # --------------------------------------------------------
            # Sidebar filters
            # --------------------------------------------------------

            ui.update_slider(
                "logfc_threshold",
                value=filters.get("logfc_threshold", 0.5)
            )

            ui.update_slider(
                "lrscore_threshold",
                value=filters.get("lrscore_threshold", 0.9)
            )

            ui.update_slider(
                "specificity_threshold",
                value=filters.get("specificity_threshold", 0.05)
            )

            ui.update_slider(
                "min_source_cells",
                value=filters.get("min_source_cells", 30)
            )

            ui.update_slider(
                "min_target_cells",
                value=filters.get("min_target_cells", 30)
            )

            ui.update_selectize(
                "source_types",
                choices=data_handler.cell_types,
                selected=filters.get("source_types", [])
            )

            ui.update_selectize(
                "target_types",
                choices=data_handler.cell_types,
                selected=filters.get("target_types", [])
            )

            # --------------------------------------------------------
            # Overview
            # --------------------------------------------------------

            ui.update_numeric(
                "overview_interaction_top_n",
                value=overview.get("interaction_top_n", 50)
            )

            ui.update_checkbox(
                "overview_apply_filters",
                value=overview.get("apply_filters", False)
            )

            comparison = overview.get("condition_comparison")

            if comparison:
                ui.update_select(
                    "overview_condition_comparison",
                    selected=comparison
                )

            ui.update_radio_buttons(
                "overview_cell_metric",
                selected=overview.get("cell_metric", "lrscore")
            )

            ui.update_radio_buttons(
                "overview_cell_mode",
                selected=overview.get("cell_mode", "top")
            )

            ui.update_numeric(
                "overview_cell_top_n",
                value=overview.get("cell_top_n", 30)
            )

            ui.update_radio_buttons(
                "overview_lr_diff_metric",
                selected=overview.get("lr_metric", "lrscore")
            )

            ui.update_numeric(
                "overview_lr_diff_top_n",
                value=overview.get("lr_top_n", 30)
            )

            ui.update_checkbox(
                "overview_lr_diff_apply_filters",
                value=overview.get("lr_apply_filters", False)
            )

            # --------------------------------------------------------
            # Network
            # --------------------------------------------------------

            ui.update_radio_buttons(
                "network_options",
                selected=network.get("mode", "all")
            )

            ui.update_numeric(
                "network_top_n",
                value=network.get("top_n", 50)
            )

            ui.update_select(
                "network_layout",
                selected=network.get("layout", "spring")
            )

            ui.update_slider(
                "cell_count_network_top_n",
                value=network.get("cell_count_top_n", 30)
            )

            # --------------------------------------------------------
            # Heatmap
            # --------------------------------------------------------

            ui.update_radio_buttons(
                "heatmap_metric",
                selected=heatmap.get("metric", "interaction_count")
            )

            ui.update_select(
                "colorscale",
                selected=heatmap.get("colorscale", "Blues")
            )

            ui.update_checkbox(
                "show_heatmap_values",
                value=heatmap.get("show_values", False)
            )

            # --------------------------------------------------------
            # Ligand–receptor tab
            # --------------------------------------------------------

            ui.update_slider(
                "top_n_interactions",
                value=lr.get("top_n", 20)
            )

            ui.update_radio_buttons(
                "dotplot_color",
                selected=lr.get("dotplot_color", "specificity_rank")
            )

            ui.update_select(
                "lr_boxplot_metric",
                selected=lr.get("boxplot_metric", "lrscore")
            )

            ui.update_numeric(
                "lr_boxplot_top_n",
                value=lr.get("boxplot_top_n", 20)
            )

            # --------------------------------------------------------
            # Data table
            # --------------------------------------------------------

            ui.update_numeric(
                "table_rows",
                value=table.get("rows", 100)
            )

            ui.update_text(
                "table_search",
                value=table.get("search", "")
            )
        
        @reactive.effect
        @reactive.event(input.restore_app_settings)
        def restore_app_settings():

            uploaded = input.restore_app_settings()

            if not uploaded:
                return

            try:

                file_path = uploaded[0]["datapath"]

                with open(
                    file_path,
                    "r",
                    encoding="utf-8"
                ) as file:
                    uploaded_json = json.load(file)

                # ----------------------------------------------------
                # Support BOTH file types:
                #
                # 1. small downloaded settings file
                # 2. complete AI snapshot containing restorable_state
                # ----------------------------------------------------

                if "restorable_state" in uploaded_json:
                    state = uploaded_json["restorable_state"]
                else:
                    state = uploaded_json

                if (
                    state.get("kind")
                    != "liana_results_explorer_state"
                ):
                    raise ValueError(
                        "This is not a valid LIANA settings file."
                    )

                target_splitting_key = (
                    state
                    .get("data_navigation", {})
                    .get("splitting_key")
                )

                # Save until correct dataset has finished loading
                pending_restore_state.set(state)

                # ----------------------------------------------------
                # Different analysis type:
                # first load it, THEN restore the remaining settings.
                # ----------------------------------------------------

                if (
                    target_splitting_key
                    and target_splitting_key
                    != input.splitting_key()
                ):

                    if (
                        target_splitting_key
                        not in data_handler.splitting_keys
                    ):
                        raise ValueError(
                            f"Analysis type '{target_splitting_key}' "
                            "is not available in this app."
                        )

                    ui.update_select(
                        "splitting_key",
                        selected=target_splitting_key
                    )

                    ui.notification_show(
                        "Loading saved analysis state...",
                        type="message"
                    )

                    return

                # Same splitting key is already loaded
                apply_restored_state(state)
                pending_restore_state.set(None)

                ui.notification_show(
                    "✓ Settings restored successfully.",
                    type="success"
                )

            except Exception as e:

                pending_restore_state.set(None)

                logger.exception(
                    "Failed to restore app settings"
                )

                ui.notification_show(
                    f"Could not restore settings: {e}",
                    type="error"
                )
           
        @output
        @render.ui
        def ai_snapshot_link():

            relative_url = (
                ai_snapshot_relative_url.get()
            )

            if not relative_url:
                return None

            return ui.div(

                ui.p(
                    "Snapshot created:",
                    style="font-weight: 600;"
                ),

                ui.tags.a(
                    "Open AI snapshot",
                    href=relative_url,
                    target="_blank"
                ),

                ui.br(),
                ui.br(),

                # ----------------------------------------------------
                # URL field
                #
                # JavaScript fills this with the complete absolute URL
                # such as:
                # http://127.0.0.1:8080/ai_snapshots/xxx.json
                # ----------------------------------------------------

                ui.tags.input(
                    id="ai_snapshot_url_field",
                    type="text",
                    readonly="readonly",
                    value=relative_url,
                    class_="form-control ai-snapshot-url",
                    onclick="this.select();"
                ),

                ui.br(),

                # ----------------------------------------------------
                # Copy button
                # ----------------------------------------------------

                ui.tags.button(
                    "📋 Copy URL",

                    id="copy_ai_snapshot_button",

                    type="button",

                    class_="btn btn-secondary",

                    onclick="""
                        const field =
                            document.getElementById(
                                'ai_snapshot_url_field'
                            );

                        if (!field) {
                            alert('Snapshot URL field not found.');
                            return;
                        }

                        field.focus();
                        field.select();
                        field.setSelectionRange(
                            0,
                            field.value.length
                        );

                        const successful =
                            document.execCommand('copy');

                        if (successful) {
                            this.innerText = '✓ Copied!';
                        } else {
                            this.innerText = 'Select URL and press Ctrl+C';
                        }
                    """
                ),

                ui.br(),
                ui.br(),

                ui.p(
                    "Paste this URL into ChatGPT, Claude, Gemini "
                    "or another LLM. The snapshot expires after "
                    "14 days or when the app container is rebuilt.",
                    class_="text-muted"
                ),

                # ----------------------------------------------------
                # Convert relative URL into complete URL
                # ----------------------------------------------------

                ui.tags.script(
                    f"""
                    (function() {{

                        const field =
                            document.getElementById(
                                'ai_snapshot_url_field'
                            );

                        if (!field) return;

                        field.value =
                            new URL(
                                '{relative_url}',
                                window.location.href
                            ).href;

                    }})();
                    """
                )
            )
            
        @reactive.effect
        def reset_ai_snapshot_when_state_changes():
            """
            Remove the displayed snapshot whenever an input that
            affects the represented analysis state changes.

            The already-created JSON file is NOT deleted.
            """

            # ========================================================
            # Register all relevant reactive dependencies
            # ========================================================

            input.splitting_key()
            input.contrast()

            # Sidebar filters
            input.logfc_threshold()
            input.lrscore_threshold()
            input.specificity_threshold()
            input.min_source_cells()
            input.min_target_cells()
            input.source_types()
            input.target_types()

            # Overview
            input.overview_interaction_top_n()
            input.overview_apply_filters()

            input.overview_condition_comparison()

            input.overview_cell_metric()
            input.overview_cell_mode()
            input.overview_cell_top_n()

            input.overview_lr_diff_metric()
            input.overview_lr_diff_top_n()
            input.overview_lr_diff_apply_filters()

            # Network
            input.network_options()
            input.network_top_n()
            input.network_layout()
            input.cell_count_network_top_n()

            # Heatmap
            input.heatmap_metric()
            input.show_heatmap_values()
            input.colorscale()

            # Ligand–receptor tab
            input.top_n_interactions()
            input.dotplot_color()

            input.lr_boxplot_metric()
            input.lr_boxplot_top_n()

            # Data table
            input.table_rows()
            input.table_search()

            # ========================================================
            # Reset the CURRENTLY DISPLAYED snapshot
            # ========================================================

            ai_snapshot_relative_url.set(
                None
            )

            ui.update_action_button(
                "create_ai_snapshot",
                label="🔗 Create AI Snapshot"
            )
        
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
        def lr_boxplot():
            """Render ligand-receptor metric distributions."""

            df = get_filtered_data()

            if df.empty:
                return go.Figure()

            metric = (
                input.lr_boxplot_metric()
                or "lrscore"
            )

            top_n = (
                input.lr_boxplot_top_n()
                or 20
            )

            return create_lr_boxplot(
                df,
                metric=metric,
                top_n=top_n
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
            """Render simple interaction counts per loaded dataset/contrast."""

            counts = interaction_counts_state.get()

            logger.info(f"Rendering structure_overview_plot with counts: {counts}")

            if not counts:
                return create_structure_overview_plot(pd.DataFrame())

            count_df = pd.DataFrame([
                {"contrast": contrast, "interactions": n}
                for contrast, n in counts.items()
            ])

            preferred_order = ["full", "HC", "UC", "CD"]
            count_df["order"] = count_df["contrast"].apply(
                lambda x: preferred_order.index(x) if x in preferred_order else 999
            )
            count_df = count_df.sort_values("order").drop(columns="order")

            return create_structure_overview_plot(count_df)
        
        
        @output
        @render_plotly
        def overview_condition_dotplot():
            """
            Plot 1:
            global Top N condition interaction landscape.

            Top N = LIANA consensus magnitude_rank
            Dot size = LRscore
            Dot color = specificity_rank
            """

            # Reactive dependency that changes after datasets finish loading
            counts = interaction_counts_state.get()

            if not counts:
                fig = go.Figure()
                fig.add_annotation(
                    text="Loading condition data...",
                    x=0.5,
                    y=0.5,
                    xref="paper",
                    yref="paper",
                    showarrow=False
                )
                return fig

            # Checkbox controls whether sidebar filters are applied
            apply_filters = input.overview_apply_filters()

            condition_results = get_overview_condition_results(
                use_sidebar_filters=apply_filters
            )

            top_n = input.overview_interaction_top_n()

            logger.info(
                f"Rendering overview condition plot: "
                f"top_n={top_n}, "
                f"conditions={[(k, len(v)) for k, v in condition_results.items()]}"
            )

            return create_global_condition_dotplot(
                condition_results,
                top_n=top_n
            )
            
        @output
        @render_plotly
        def overview_cell_difference_plot():
            """
            Plot 2:
            cell-cell communication differences.

            Metric options:
                lrscore
                interaction_count
            """

            (
                condition_a,
                condition_b
            ) = get_selected_overview_conditions()

            if (
                not condition_a
                or not condition_b
            ):
                return go.Figure()

            mode = input.overview_cell_mode()

            use_filters = (
                mode == "filtered"
            )

            condition_results = (
                get_overview_condition_results(
                    use_sidebar_filters=use_filters
                )
            )

            if (
                condition_a not in condition_results
                or condition_b not in condition_results
            ):
                return go.Figure()

            metric = (
                input.overview_cell_metric()
                or "lrscore"
            )

            top_n = None

            if mode == "top":
                top_n = (
                    input.overview_cell_top_n()
                )

            return create_cell_pair_difference_heatmap(

                condition_results[
                    condition_a
                ],

                condition_results[
                    condition_b
                ],

                condition_a=condition_a,
                condition_b=condition_b,

                metric=metric,

                top_n=top_n
            )
        @output
        @render_plotly
        def overview_lr_difference_plot():
            """
            Plot 3:
            Top N ligand-receptor differences.

            Metric options:
                lrscore
                interaction_count

            Sidebar filtering is independent from Top N.
            """

            (
                condition_a,
                condition_b
            ) = get_selected_overview_conditions()

            if (
                not condition_a
                or not condition_b
            ):
                return go.Figure()

            # Independent checkbox:
            # Top N is ALWAYS applied.
            use_filters = (
                input.overview_lr_diff_apply_filters()
            )

            condition_results = (
                get_overview_condition_results(
                    use_sidebar_filters=use_filters
                )
            )

            if (
                condition_a not in condition_results
                or condition_b not in condition_results
            ):
                return go.Figure()

            metric = (
                input.overview_lr_diff_metric()
                or "lrscore"
            )

            top_n = (
                input.overview_lr_diff_top_n()
                or 30
            )

            return create_lr_difference_barplot(

                condition_results[
                    condition_a
                ],

                condition_results[
                    condition_b
                ],

                condition_a=condition_a,
                condition_b=condition_b,

                metric=metric,

                top_n=top_n
            )
            
        @output
        @render.ui
        def overview_condition_dotplot_container():

            top_n = input.overview_interaction_top_n() or 50

            plot_height = max(
                600,
                top_n * 24
            )

            return output_widget(
                "overview_condition_dotplot",
                height=f"{plot_height}px"
            )


        @output
        @render.ui
        def overview_cell_difference_plot_container():

            return output_widget(
                "overview_cell_difference_plot",
                height="900px"
            )


        @output
        @render.ui
        def overview_lr_difference_plot_container():

            top_n = (
                input.overview_lr_diff_top_n()
                or 30
            )

            plot_height = max(
                600,
                top_n * 24
            )

            return output_widget(
                "overview_lr_difference_plot",
                height=f"{plot_height}px"
            )
                
    return server
