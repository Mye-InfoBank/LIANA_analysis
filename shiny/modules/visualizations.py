#!/usr/bin/env python3
"""
Visualization Module for LIANA Results Explorer
Contains all plotting functions for different visualization types.
"""

import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import networkx as nx
from typing import Dict, List, Optional
import logging

logger = logging.getLogger(__name__)

def create_network_plot(df: pd.DataFrame, 
                       layout_algorithm: str = 'spring',
                       node_size_factor: float = 20,
                       edge_width_factor: float = 5) -> go.Figure:
    """
    Create an interactive network plot of cell-cell interactions.
    
    Args:
        df: DataFrame with interaction data
        layout_algorithm: NetworkX layout algorithm ('spring', 'circular', 'kamada_kawai')
        node_size_factor: Factor for node size scaling
        edge_width_factor: Factor for edge width scaling
        
    Returns:
        Plotly Figure object
    """
    if df.empty:
        fig = go.Figure()
        fig.add_annotation(
            text="No data to display", 
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=16)
        )
        return fig
    
    try:
        # Count interactions between source-target pairs
        interaction_counts = df.groupby(["source", "target"]).size().reset_index(name="count")
        
        if interaction_counts.empty:
            logger.warning("No interactions to visualize in network")
            return go.Figure()
        
        # Create NetworkX graph
        G = nx.DiGraph()
        
        # Add edges with weights
        for _, row in interaction_counts.iterrows():
            G.add_edge(row['source'], row['target'], weight=row['count'])
        
        # Choose layout algorithm
        layout_functions = {
            'spring': nx.spring_layout,
            'circular': nx.circular_layout,
            'kamada_kawai': nx.kamada_kawai_layout,
            'shell': nx.shell_layout
        }
        
        # Compute positions; only spring layout supports 'seed' reliably across versions
        if layout_algorithm == 'spring':
            pos = nx.spring_layout(G, seed=42)
        elif layout_algorithm == 'kamada_kawai':
            pos = nx.kamada_kawai_layout(G)
        elif layout_algorithm == 'circular':
            pos = nx.circular_layout(G)
        elif layout_algorithm == 'shell':
            pos = nx.shell_layout(G)
        else:
            pos = nx.spring_layout(G, seed=42)
        
        # Create node trace
        node_x = [pos[node][0] for node in G.nodes()]
        node_y = [pos[node][1] for node in G.nodes()]
        node_text = list(G.nodes())
        
        # Calculate node sizes based on degree
        node_degrees = [G.degree(node) for node in G.nodes()]
        max_degree = max(node_degrees) if node_degrees else 1
        node_sizes = [degree / max_degree * node_size_factor + 10 for degree in node_degrees]
        
        node_trace = go.Scatter(
            x=node_x, y=node_y,
            mode='markers+text',
            text=node_text,
            textposition='top center',
            textfont=dict(size=12, color='black'),
            marker=dict(
                size=node_sizes,
                color='lightblue',
                line=dict(width=1, color='black'),
                colorscale='Blues'
            ),
            hovertemplate='<b>%{text}</b><br>Degree: %{marker.size}<extra></extra>',
            name='Cell Types'
        )
        
        # Create edge traces
        edge_traces = []
        max_weight = max([G[u][v]['weight'] for u, v in G.edges()]) if G.edges() else 1
        
        for edge in G.edges():
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            weight = G[edge[0]][edge[1]]['weight']
            
            # Calculate edge width
            edge_width = (weight / max_weight) * edge_width_factor + 1
            
            edge_trace = go.Scatter(
                x=[x0, x1, None],
                y=[y0, y1, None],
                mode='lines',
                line=dict(width=edge_width, color='gray'),
                opacity=0.25,
                hovertemplate=f'{edge[0]} → {edge[1]}<br>Interactions: {weight}<extra></extra>',
                showlegend=False
            )
            edge_traces.append(edge_trace)
        
        # Create figure
        fig = go.Figure(data=[node_trace] + edge_traces)
        
        fig.update_layout(
            title=dict(text="Cell-Cell Interaction Network", font=dict(size=16)),
            showlegend=False,
            hovermode='closest',
            margin=dict(b=20, l=5, r=5, t=40),
            annotations=[
                dict(
                    text="Network shows interactions between cell types.<br>Node size = connectivity, Edge thickness = interaction strength.",
                    showarrow=False, xref="paper", yref="paper",
                    x=0.005, y=-0.002, xanchor='left', yanchor='bottom',
                    font=dict(size=10)
                )
            ],
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            plot_bgcolor='white',
            uniformtext_minsize=10,
            uniformtext_mode='hide'
        )
        
        return fig
        
    except Exception as e:
        logger.error(f"Error creating network plot: {e}")
        fig = go.Figure()
        fig.add_annotation(
            text=f"Error creating network plot: {str(e)}", 
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig
    
def create_cell_count_network_plot(
    df: pd.DataFrame,
    layout_algorithm: str = "spring",
    top_n_nodes: Optional[int] = None,
) -> go.Figure:
    """
    Create a cell-count-aware network plot.

    Node size = number of cells in source/target cell type.
    Edge width = number of filtered LR interaction rows between source and target.
    Optional top_n_nodes = keep only the N most connected cell types,
    where connectivity = total number of filtered LR rows touching the node.
    """

    if df.empty:
        fig = go.Figure()
        fig.add_annotation(
            text="No data to display",
            xref="paper", yref="paper",
            x=0.5, y=0.5,
            showarrow=False,
            font=dict(size=16)
        )
        return fig

    required = {"source", "target"}
    if not required.issubset(df.columns):
        fig = go.Figure()
        fig.add_annotation(
            text="Missing source/target columns",
            xref="paper", yref="paper",
            x=0.5, y=0.5,
            showarrow=False,
            font=dict(size=16)
        )
        return fig

    try:
        # --------------------------------------------------
        # Edge weight = number of LR rows after filtering
        # --------------------------------------------------
        edge_df = (
            df.groupby(["source", "target"])
            .size()
            .reset_index(name="interaction_count")
        )

        if edge_df.empty:
            fig = go.Figure()
            fig.add_annotation(
                text="No interactions after filtering",
                xref="paper", yref="paper",
                x=0.5, y=0.5,
                showarrow=False,
                font=dict(size=16)
            )
            return fig

        # --------------------------------------------------
        # Cell counts per node from source_n_cells / target_n_cells
        # --------------------------------------------------
        cell_count_map = {}

        if "source_n_cells" in df.columns:
            source_counts = (
                df[["source", "source_n_cells"]]
                .dropna()
                .drop_duplicates()
            )
            for _, row in source_counts.iterrows():
                cell_count_map[row["source"]] = int(row["source_n_cells"])

        if "target_n_cells" in df.columns:
            target_counts = (
                df[["target", "target_n_cells"]]
                .dropna()
                .drop_duplicates()
            )
            for _, row in target_counts.iterrows():
                celltype = row["target"]
                n_cells = int(row["target_n_cells"])
                cell_count_map[celltype] = max(cell_count_map.get(celltype, 0), n_cells)

        all_nodes = sorted(set(edge_df["source"]).union(edge_df["target"]))

        # fallback if old CSVs do not yet have source_n_cells/target_n_cells
        for node in all_nodes:
            cell_count_map.setdefault(node, 1)

        # --------------------------------------------------
        # Connectivity = total number of LR rows touching node
        # --------------------------------------------------
        connectivity_map = {node: 0 for node in all_nodes}

        for _, row in edge_df.iterrows():
            source = row["source"]
            target = row["target"]
            count = int(row["interaction_count"])

            connectivity_map[source] += count
            connectivity_map[target] += count

        # --------------------------------------------------
        # Optional: keep only top N connected cell types
        # --------------------------------------------------
        if top_n_nodes is not None and top_n_nodes > 0:
            top_nodes = sorted(
                connectivity_map,
                key=connectivity_map.get,
                reverse=True
            )[:top_n_nodes]

            edge_df = edge_df[
                edge_df["source"].isin(top_nodes) &
                edge_df["target"].isin(top_nodes)
            ].copy()

            all_nodes = sorted(set(edge_df["source"]).union(edge_df["target"]))

            if edge_df.empty or not all_nodes:
                fig = go.Figure()
                fig.add_annotation(
                    text="No edges left after Top N connected-cell filtering",
                    xref="paper", yref="paper",
                    x=0.5, y=0.5,
                    showarrow=False,
                    font=dict(size=16)
                )
                return fig

            # recalculate connectivity inside the shown subgraph
            connectivity_map = {node: 0 for node in all_nodes}
            for _, row in edge_df.iterrows():
                source = row["source"]
                target = row["target"]
                count = int(row["interaction_count"])

                connectivity_map[source] += count
                connectivity_map[target] += count

        # --------------------------------------------------
        # Build graph
        # --------------------------------------------------
        G = nx.DiGraph()

        for node in all_nodes:
            G.add_node(
                node,
                n_cells=cell_count_map.get(node, 1),
                connected_interactions=connectivity_map.get(node, 0)
            )

        for _, row in edge_df.iterrows():
            G.add_edge(
                row["source"],
                row["target"],
                weight=int(row["interaction_count"])
            )

        # --------------------------------------------------
        # Layout
        # --------------------------------------------------
        if layout_algorithm == "spring":
            pos = nx.spring_layout(G, seed=42, k=0.9)
        elif layout_algorithm == "kamada_kawai":
            pos = nx.kamada_kawai_layout(G)
        elif layout_algorithm == "circular":
            pos = nx.circular_layout(G)
        else:
            pos = nx.spring_layout(G, seed=42, k=0.9)

        # --------------------------------------------------
        # Helpers for scaling
        # --------------------------------------------------
        def scale_value(x, min_x, max_x, min_out, max_out):
            if max_x == min_x:
                return (min_out + max_out) / 2
            return min_out + ((x - min_x) / (max_x - min_x)) * (max_out - min_out)

        # --------------------------------------------------
        # Nodes: size = cell count
        # --------------------------------------------------
        node_x = []
        node_y = []
        node_labels = []
        node_sizes = []
        node_customdata = []

        n_cells_values = [G.nodes[n]["n_cells"] for n in G.nodes()]
        min_cells = min(n_cells_values)
        max_cells = max(n_cells_values)

        for node in G.nodes():
            n_cells = G.nodes[node]["n_cells"]
            connected_interactions = G.nodes[node]["connected_interactions"]
            n_neighbors = G.degree(node)

            node_x.append(pos[node][0])
            node_y.append(pos[node][1])
            node_labels.append(node)
            node_sizes.append(scale_value(n_cells, min_cells, max_cells, 12, 55))
            node_customdata.append([n_cells, connected_interactions, n_neighbors])

        node_trace = go.Scatter(
            x=node_x,
            y=node_y,
            mode="markers+text",
            text=node_labels,
            textposition="top center",
            textfont=dict(size=11, color="black"),
            marker=dict(
                size=node_sizes,
                color="lightblue",
                line=dict(width=1, color="black"),
                opacity=0.85
            ),
            customdata=node_customdata,
            hovertemplate=(
                "<b>%{text}</b><br>"
                "Cells: %{customdata[0]}<br>"
                "Filtered LR interactions connected to this node: %{customdata[1]}<br>"
                "Network neighbors: %{customdata[2]}"
                "<extra></extra>"
            ),
            name="Cell types"
        )

        # --------------------------------------------------
        # Edges: width = interaction count
        # --------------------------------------------------
        edge_traces = []

        weights = [G[u][v]["weight"] for u, v in G.edges()]
        min_weight = min(weights)
        max_weight = max(weights)

        for source, target in G.edges():
            x0, y0 = pos[source]
            x1, y1 = pos[target]
            weight = G[source][target]["weight"]

            edge_width = scale_value(weight, min_weight, max_weight, 1, 9)

            edge_trace = go.Scatter(
                x=[x0, x1, None],
                y=[y0, y1, None],
                mode="lines",
                line=dict(width=edge_width, color="gray"),
                opacity=0.30,
                hovertemplate=(
                    f"{source} → {target}<br>"
                    f"Filtered LR interactions: {weight}"
                    "<extra></extra>"
                ),
                showlegend=False
            )
            edge_traces.append(edge_trace)

        fig = go.Figure(data=edge_traces + [node_trace])

        title = "Cell-count-aware interaction network"
        if top_n_nodes is not None and top_n_nodes > 0:
            title += f" — top {top_n_nodes} connected cell types"

        fig.update_layout(
            title=dict(text=title, font=dict(size=16)),
            showlegend=False,
            hovermode="closest",
            margin=dict(b=40, l=5, r=5, t=50),
            annotations=[
                dict(
                    text=(
                        "Node size = number of cells. "
                        "Edge thickness = number of filtered ligand–receptor rows."
                    ),
                    showarrow=False,
                    xref="paper",
                    yref="paper",
                    x=0.005,
                    y=-0.04,
                    xanchor="left",
                    yanchor="bottom",
                    font=dict(size=10)
                )
            ],
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            plot_bgcolor="white"
        )

        return fig

    except Exception as e:
        logger.error(f"Error creating cell-count network plot: {e}")
        fig = go.Figure()
        fig.add_annotation(
            text=f"Error creating cell-count network plot: {str(e)}",
            xref="paper", yref="paper",
            x=0.5, y=0.5,
            showarrow=False
        )
        return fig
    
def create_heatmap_plot(df: pd.DataFrame, 
                       value_col: str = "interaction_count",
                       colorscale: str = 'Blues',
                       show_values: bool = True) -> go.Figure:
    """
    Create an interactive heatmap of source-target cell-type interactions.

    If value_col == "interaction_count":
        heatmap value = number of filtered LR rows per source-target pair.

    Otherwise:
        heatmap value = median of selected LIANA metric per source-target pair.
    """

    if df.empty:
        fig = go.Figure()
        fig.add_annotation(
            text="No data to display", 
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig

    if not all(col in df.columns for col in ["source", "target"]):
        fig = go.Figure()
        fig.add_annotation(
            text="Missing source/target columns for heatmap", 
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig

    try:
        # --------------------------------------------------
        # Case 1: count filtered LR rows per source-target pair
        # --------------------------------------------------
        if value_col == "interaction_count":
            heatmap_data = (
                df.groupby(["source", "target"])
                .size()
                .reset_index(name="interaction_count")
                .pivot(index="source", columns="target", values="interaction_count")
                .fillna(0)
            )

            title = "Cell-Cell Interaction Heatmap (filtered LR interaction count)"
            colorbar_title = "Filtered LR count"
            hovertemplate = (
                "Source: %{y}<br>"
                "Target: %{x}<br>"
                "Filtered LR interactions: %{z:.0f}"
                "<extra></extra>"
            )

        # --------------------------------------------------
        # Case 2: median LIANA metric per source-target pair
        # --------------------------------------------------
        else:
            if value_col not in df.columns:
                fig = go.Figure()
                fig.add_annotation(
                    text=f"Column '{value_col}' not found in data", 
                    xref="paper", yref="paper",
                    x=0.5, y=0.5, showarrow=False
                )
                return fig

            heatmap_data = df.pivot_table(
                index="source", 
                columns="target", 
                values=value_col, 
                aggfunc="median"
            ).fillna(0)

            title = f"Cell-Cell Interaction Heatmap (median {value_col})"
            colorbar_title = value_col
            hovertemplate = (
                "Source: %{y}<br>"
                "Target: %{x}<br>"
                f"Median {value_col}: " + "%{z:.3f}"
                "<extra></extra>"
            )

        if heatmap_data.empty:
            fig = go.Figure()
            fig.add_annotation(
                text="No data to create heatmap", 
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False
            )
            return fig

        # Optional: sort rows/columns by total signal so strongest pairs are easier to see
        row_order = heatmap_data.sum(axis=1).sort_values(ascending=False).index
        col_order = heatmap_data.sum(axis=0).sort_values(ascending=False).index
        heatmap_data = heatmap_data.loc[row_order, col_order]

        fig = go.Figure(data=go.Heatmap(
            z=heatmap_data.values,
            x=heatmap_data.columns,
            y=heatmap_data.index,
            colorscale=colorscale,
            hoverongaps=False,
            hovertemplate=hovertemplate,
            showscale=True,
            colorbar=dict(title=colorbar_title)
        ))

        # Add text annotations only if heatmap is small enough
        if show_values and heatmap_data.shape[0] * heatmap_data.shape[1] <= 100:
            annotations = []
            max_value = heatmap_data.values.max()

            for i, row in enumerate(heatmap_data.index):
                for j, col in enumerate(heatmap_data.columns):
                    value = heatmap_data.iloc[i, j]

                    if value != 0:
                        if value_col == "interaction_count":
                            text_value = f"{int(value)}"
                        else:
                            text_value = f"{value:.2f}"

                        annotations.append(
                            dict(
                                x=col,
                                y=row,
                                text=text_value,
                                showarrow=False,
                                font=dict(
                                    color="white" if value > max_value / 2 else "black",
                                    size=9
                                )
                            )
                        )

            fig.update_layout(annotations=annotations)

        fig.update_layout(
            title=title,
            xaxis_title="Target Cell Type",
            yaxis_title="Source Cell Type",
            font=dict(size=10),
            height=max(500, len(heatmap_data.index) * 14),
            margin=dict(l=180, r=80, t=70, b=180)
        )

        fig.update_xaxes(tickangle=90)
        fig.update_yaxes(automargin=True)

        return fig

    except Exception as e:
        logger.error(f"Error creating heatmap: {e}")
        fig = go.Figure()
        fig.add_annotation(
            text=f"Error creating heatmap: {str(e)}", 
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig

def create_dotplot(df: pd.DataFrame, 
                  top_n: int = 20,
                  size_col: str = 'lrscore',
                  color_col: str = 'specificity_rank',
                  x_col: str = 'lr_logfc') -> go.Figure:
    """
    Create a dot plot of top interactions.
    
    Args:
        df: DataFrame with interaction data
        top_n: Number of top interactions to show
        size_col: Column to use for dot size
        color_col: Column to use for dot color
        x_col: Column to use for x-axis
        
    Returns:
        Plotly Figure object
    """
    if df.empty:
        fig = go.Figure()
        fig.add_annotation(
            text="No data to display", 
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig
    
    try:
        # Get top interactions
        if 'magnitude_rank' in df.columns:
            top_interactions = df.nsmallest(top_n, 'magnitude_rank').copy()
        else:
            top_interactions = df.head(top_n).copy()
        
        if top_interactions.empty:
            logger.warning("No top interactions found after filtering")
            return go.Figure()
        
        # Create interaction labels
        top_interactions['interaction_label'] = (
            top_interactions['ligand_complex'] + ' → ' + 
            top_interactions['receptor_complex']
        )
        
        top_interactions['cell_pair'] = (
            top_interactions['source'] + ' → ' + 
            top_interactions['target']
        )
        
        # Calculate sizes with conservative scaling (approx 6–18 px)
        def _normalize(series: pd.Series) -> pd.Series:
            try:
                s = pd.to_numeric(series, errors='coerce')
            except Exception:
                s = pd.Series([np.nan] * len(series), index=series.index)
            s_min = s.min(skipna=True)
            s_max = s.max(skipna=True)
            if pd.isna(s_min) or pd.isna(s_max) or s_max == s_min:
                return pd.Series([0.5] * len(series), index=series.index)
            return (s - s_min) / (s_max - s_min)

        if size_col in top_interactions.columns:
            if size_col == 'specificity_rank':
                # Smaller rank -> larger dot. Normalize then invert.
                norm = _normalize(top_interactions[size_col])
                sizes = 6 + (1 - norm) * 12
            else:
                norm = _normalize(top_interactions[size_col])
                sizes = 6 + norm * 12
        else:
            sizes = [8] * len(top_interactions)
        
        # Create the scatter plot
        fig = go.Figure(data=go.Scatter(
            x=top_interactions[x_col].tolist() if x_col in top_interactions.columns else list(range(len(top_interactions))),
            y=list(range(len(top_interactions))),
            mode='markers',
            marker=dict(
                size=sizes,
                color=top_interactions[color_col] if color_col in top_interactions.columns else 'blue',
                colorscale='Plasma',
                showscale=True,
                colorbar=dict(title=color_col.replace('_', ' ').title()),
                line=dict(width=1, color='black'),
                opacity=0.8
            ),
            text=top_interactions['interaction_label'],
            customdata=top_interactions['cell_pair'],
            hovertemplate='<b>%{text}</b><br>%{customdata}<br>' + 
                         f'{x_col}: %{{x:.3f}}<br>' +
                         f'{color_col}: %{{marker.color}}<br>' +
                         '<extra></extra>',
            showlegend=False
        ))
        
        # Update layout
        fig.update_layout(
            title=f"Top {len(top_interactions)} Ligand-Receptor Interactions",
            xaxis_title=x_col.replace('_', ' ').title(),
            yaxis=dict(
                tickmode='array',
                tickvals=list(range(len(top_interactions))),
                ticktext=[f"{row['ligand_complex']} → {row['receptor_complex']}" 
                         for _, row in top_interactions.iterrows()]
            ),
            height=max(400, len(top_interactions) * 25),
            margin=dict(l=200, r=50, t=50, b=50)
        )
        
        return fig
        
    except Exception as e:
        logger.error(f"Error creating dot plot: {e}")
        fig = go.Figure()
        fig.add_annotation(
            text=f"Error creating dot plot: {str(e)}", 
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig

def create_scatter_comparison(results_dict: Dict, 
                            source_cell: str, 
                            target_cell: str,
                            value_col: str = 'lr_means',
                            size_col: str = 'lrscore',
                            max_points_per_condition: Optional[int] = None) -> go.Figure:
    """
    Create scatter plot comparing interactions across conditions.
    
    Args:
        results_dict: Dictionary of condition -> DataFrame
        source_cell: Source cell type to focus on
        target_cell: Target cell type to focus on
        value_col: Column to use for color
        size_col: Column to use for size
        
    Returns:
        Plotly Figure object
    """
    if len(results_dict) < 2:
        fig = go.Figure()
        fig.add_annotation(
            text="Need at least 2 conditions for comparison", 
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig
    
    try:
        # Combine all condition data
        combined_data = []
        for condition, df in results_dict.items():
            if condition != 'full' and not df.empty:
                condition_data = df[
                    (df['source'] == source_cell) & 
                    (df['target'] == target_cell)
                ].copy()
                condition_data['condition'] = condition
                # Enforce per-condition cap here if requested
                if max_points_per_condition is not None and max_points_per_condition > 0:
                    if 'magnitude_rank' in condition_data.columns:
                        condition_data = condition_data.nsmallest(max_points_per_condition, 'magnitude_rank')
                    else:
                        condition_data = condition_data.head(max_points_per_condition)
                combined_data.append(condition_data)
        
        if not combined_data:
            fig = go.Figure()
            fig.add_annotation(
                text="No data for specified cell types", 
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False
            )
            return fig
        
        combined_df = pd.concat(combined_data, ignore_index=True)
        
        # Create interaction labels
        combined_df['interaction_label'] = (
            combined_df['ligand_complex'] + ' → ' + combined_df['receptor_complex']
        )
        
        # Create the scatter plot
        fig = px.scatter(
            combined_df,
            x='condition',
            y='interaction_label',
            size=size_col if size_col in combined_df.columns else None,
            color=value_col if value_col in combined_df.columns else None,
            color_continuous_scale='Viridis',
            hover_data=['specificity_rank', 'magnitude_rank'] if all(col in combined_df.columns for col in ['specificity_rank', 'magnitude_rank']) else None,
            title=f"{source_cell} → {target_cell} Interactions Across Conditions"
        )
        
        fig.update_layout(
            xaxis_title="Condition",
            yaxis_title="Ligand → Receptor",
            height=max(400, len(combined_df['interaction_label'].unique()) * 25)
        )
        
        return fig
        
    except Exception as e:
        logger.error(f"Error creating comparison plot: {e}")
        fig = go.Figure()
        fig.add_annotation(
            text=f"Error creating comparison plot: {str(e)}", 
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig

def create_volcano_plot(df: pd.DataFrame,
                       x_col: str = 'lr_logfc',
                       y_col: str = 'lrscore',
                       significance_threshold: float = 0.05) -> go.Figure:
    """
    Create a volcano plot showing effect size vs significance.
    
    Args:
        df: DataFrame with interaction data
        x_col: Column for x-axis (effect size)
        y_col: Column for y-axis (significance)
        significance_threshold: Threshold for significance
        
    Returns:
        Plotly Figure object
    """
    if df.empty or x_col not in df.columns or y_col not in df.columns:
        fig = go.Figure()
        fig.add_annotation(
            text="Insufficient data for volcano plot", 
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig
    
    try:
        # Create significance categories
        df_plot = df.copy()
        df_plot['significant'] = df_plot['specificity_rank'] < significance_threshold if 'specificity_rank' in df.columns else False
        df_plot['interaction_label'] = df_plot['ligand_complex'] + ' → ' + df_plot['receptor_complex']
        
        # Create the plot
        fig = px.scatter(
            df_plot,
            x=x_col,
            y=y_col,
            color='significant',
            hover_data=['interaction_label', 'source', 'target'],
            title="Volcano Plot: Effect Size vs Significance",
            color_discrete_map={True: 'red', False: 'gray'}
        )
        
        # Add threshold lines
        if 'specificity_rank' in df.columns:
            fig.add_hline(y=significance_threshold, line_dash="dash", line_color="red")
        
        fig.update_layout(
            xaxis_title=x_col.replace('_', ' ').title(),
            yaxis_title=y_col.replace('_', ' ').title()
        )
        
        return fig
        
    except Exception as e:
        logger.error(f"Error creating volcano plot: {e}")
        fig = go.Figure()
        fig.add_annotation(
            text=f"Error creating volcano plot: {str(e)}", 
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig

def create_summary_barplot(summary_stats: Dict) -> go.Figure:
    """
    Create a bar plot of summary statistics.
    
    Args:
        summary_stats: Dictionary with summary statistics
        
    Returns:
        Plotly Figure object
    """
    if not summary_stats:
        fig = go.Figure()
        fig.add_annotation(
            text="No summary data available", 
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig
    
    try:
        # Prepare data for plotting
        metrics = []
        values = []
        
        for key, value in summary_stats.items():
            if isinstance(value, (int, float)) and key != 'total_interactions':
                metrics.append(key.replace('_', ' ').title())
                values.append(value)
        
        if not metrics:
            return go.Figure()
        
        fig = go.Figure(data=[
            go.Bar(
                x=metrics,
                y=values,
                marker_color='steelblue',
                text=values,
                textposition='auto'
            )
        ])
        
        fig.update_layout(
            title="Dataset Summary Statistics",
            xaxis_title="Metrics",
            yaxis_title="Count",
            showlegend=False
        )
        
        return fig
        
    except Exception as e:
        logger.error(f"Error creating summary plot: {e}")
        return go.Figure()

def create_structure_overview_plot(splitting_keys_summary: Dict) -> go.Figure:
    """
    Create an overview plot of the data structure across splitting keys.
    
    Args:
        splitting_keys_summary: Dictionary with summary for each splitting key
        
    Returns:
        Plotly Figure object
    """
    if not splitting_keys_summary:
        fig = go.Figure()
        fig.add_annotation(
            text="No structure data available", 
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig
    
    try:
        # Prepare data
        splitting_keys = list(splitting_keys_summary.keys())
        num_conditions = [info['num_conditions'] for info in splitting_keys_summary.values()]
        total_interactions = [info['total_interactions'] for info in splitting_keys_summary.values()]
        
        # Create subplot with secondary y-axis
        fig = make_subplots(specs=[[{"secondary_y": True}]])
        
        # Add bar chart for number of conditions
        fig.add_trace(
            go.Bar(
                x=splitting_keys,
                y=num_conditions,
                name="Number of Conditions",
                marker_color='lightblue',
                yaxis='y1'
            ),
            secondary_y=False,
        )
        
        # Add line chart for total interactions
        fig.add_trace(
            go.Scatter(
                x=splitting_keys,
                y=total_interactions,
                mode='lines+markers',
                name="Total Interactions",
                line=dict(color='red', width=3),
                marker=dict(size=8),
                yaxis='y2'
            ),
            secondary_y=True,
        )
        
        # Update layout
        fig.update_xaxes(title_text="Analysis Type")
        fig.update_yaxes(title_text="Number of Conditions", secondary_y=False)
        fig.update_yaxes(title_text="Total Interactions", secondary_y=True)
        
        fig.update_layout(
            title="Data Structure Overview Across Analysis Types",
            legend=dict(x=0.01, y=0.99),
            height=400
        )
        
        return fig
        
    except Exception as e:
        logger.error(f"Error creating structure overview plot: {e}")
        return go.Figure()

def create_cross_analysis_comparison(combined_df: pd.DataFrame, 
                                   source_cell: str, 
                                   target_cell: str) -> go.Figure:
    """
    Create a comparison plot across different analysis types (splitting keys).
    
    Args:
        combined_df: DataFrame with data from multiple splitting keys
        source_cell: Source cell type
        target_cell: Target cell type
        
    Returns:
        Plotly Figure object
    """
    if combined_df.empty:
        fig = go.Figure()
        fig.add_annotation(
            text="No data available for cross-analysis comparison", 
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig
    
    try:
        # Create interaction labels
        combined_df = combined_df.copy()
        combined_df['interaction_label'] = (
            combined_df['ligand_complex'] + ' → ' + combined_df['receptor_complex']
        )
        
        # Create the comparison plot
        fig = px.scatter(
            combined_df,
            x="splitting_key",
            y="interaction_label",
            size="lrscore" if "lrscore" in combined_df.columns else None,
            color="lr_means" if "lr_means" in combined_df.columns else None,
            color_continuous_scale="Viridis",
            hover_data=[
                col for col in [
                    "source",
                    "target",
                    "specificity_rank",
                    "magnitude_rank",
                    "lr_logfc",
                    "lrscore",
                    "lr_means"
                ]
                if col in combined_df.columns
            ],
            title=f"Cross-Analysis Comparison: {source_cell} → {target_cell}",
        )

        n_pairs = combined_df["interaction_label"].nunique()

        fig.update_layout(
            xaxis_title="Analysis Type",
            yaxis_title="Ligand → Receptor",
            height=max(520, n_pairs * 24),
            margin=dict(l=230, r=130, t=60, b=50),
        )

        fig.update_yaxes(
            automargin=False,
            tickfont=dict(size=10),
        )

        fig.update_xaxes(
            automargin=False,
            tickfont=dict(size=11),
        )

        return fig
        
    except Exception as e:
        logger.error(f"Error creating cross-analysis comparison: {e}")
        fig = go.Figure()
        fig.add_annotation(
            text=f"Error creating comparison plot: {str(e)}", 
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig
