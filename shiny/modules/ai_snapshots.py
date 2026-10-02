#!/usr/bin/env python3

"""
AI Snapshot utilities for the LIANA Results Explorer.

A snapshot is a machine-readable JSON export of the current dashboard
state. The user downloads it and uploads it to an LLM of their choice.
"""

import json


def figure_to_dict(fig):
    """
    Convert a Plotly figure into a JSON-safe Python dictionary.

    This stores the actual information behind the visible plot:
    axes, labels, traces, values, hover data, etc.
    """

    if fig is None:
        return None

    try:
        return json.loads(fig.to_json())

    except Exception:
        return None


def llm_context():
    """
    Static context so an LLM can interpret a snapshot without further
    input from the user.
    """

    return {
        "purpose": (
            "This file is a snapshot of the COST IBD MyeInfoBank LIANA "
            "Results Explorer (Shiny app). It contains what the user was "
            "looking at when it was exported."
        ),

        "instructions_for_assistant": [
            "Base your answers on this file. If something is not in it, "
            "say so instead of guessing.",
            "Start by briefly stating the dataset, analysis type, contrast "
            "and active filters the snapshot shows.",
            "Treat top-N lists and plots as truncated views: absence from "
            "them does not mean an interaction is absent.",
            "Rank-based scores are relative within one LIANA run; do not "
            "compare raw values across runs without saying so.",
            "Present biological interpretation as hypotheses, not "
            "conclusions; LIANA predicts potential communication from "
            "expression, not confirmed signalling.",
            "When citing interactions, quote ligand, receptor, source and "
            "target cell types with their lrscore and ranks from this file.",
            "If the user's question is vague, summarise the main findings "
            "and ask what they want to explore."
        ],

        "study_context": {
            "project": (
                "COST IBD MyeInfoBank (inflammatory bowel disease "
                "single-cell atlas)"
            ),
            "dataset": (
                "Tissue-specific sub-atlas (e.g. colon or ileum); see "
                "metadata.dataset. Cell types are the tier_2 annotation."
            ),
            "analysis_types": (
                "metadata.splitting_key is the column the data were split "
                "by (e.g. 'condition'). LIANA was run once on all cells and "
                "once per group; metadata.selected_contrast is the run shown."
            ),
            "condition_abbreviations": (
                "Group abbreviations are dataset-specific; do not assume "
                "their meaning, ask the user if unclear."
            )
        },

        "methods": {
            "liana": (
                "LIANA rank_aggregate (liana-py) on normalised expression "
                "(use_raw=False), run separately for each group of the "
                "splitting key. It combines several ligand-receptor methods "
                "into consensus ranks."
            ),
            "per_group_runs": (
                "Each condition is inferred independently, so differences "
                "between conditions are descriptive comparisons of separate "
                "networks, not formal differential tests."
            ),
            "app_filters": (
                "The sidebar filters (lr_logfc, lrscore, specificity_rank, "
                "minimum cells, source/target cell types) are applied after "
                "inference; 'filtered_dataset' and the top interaction list "
                "reflect them."
            )
        },

        "field_glossary": {
            "source": "Cell type expressing the ligand (sender).",
            "target": "Cell type expressing the receptor (receiver).",
            "source_n_cells": "Number of source cells behind the estimate.",
            "target_n_cells": "Number of target cells behind the estimate.",
            "ligand_complex": "Ligand (may be a protein complex).",
            "receptor_complex": "Receptor (may be a protein complex).",
            "lrscore": (
                "Interaction magnitude score (0-1); higher means stronger "
                "inferred communication."
            ),
            "lr_logfc": (
                "Specificity-related log fold change of ligand/receptor "
                "expression for the cell types within one run. Not a "
                "between-condition differential log fold change."
            ),
            "lr_means": "Mean ligand and receptor expression.",
            "specificity_rank": (
                "Consensus specificity rank (0-1); lower means more "
                "cell-type-specific. It is not a p-value."
            ),
            "magnitude_rank": (
                "Consensus magnitude rank (0-1); lower means stronger."
            )
        },

        "truncation": {
            "top_filtered_interactions": (
                "up to 500 interactions passing the filters, ordered by "
                "magnitude_rank"
            ),
            "data_table_preview": (
                "the first rows of the Data Table as displayed, at most 500"
            ),
            "plot_outputs": (
                "Plotly figure specifications of the currently selected "
                "plot settings; they have their own top-N limits. Network "
                "graphs are not included (see top_filtered_interactions)"
            )
        },

        "suggested_questions": [
            "Summarise the strongest cell-cell communication in this view.",
            "Which ligand-receptor pairs stand out between the source and "
            "target cell types, and what might they mean in IBD?",
            "How do the two compared conditions differ in this view?",
            "What are the limitations of this result?"
        ]
    }
