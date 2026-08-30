#!/usr/bin/env python3
"""
UI Components Module for LIANA Results Explorer
Contains all UI component definitions and layout functions.
"""

from shiny import ui
from shinywidgets import output_widget
from typing import List, Optional

def create_app_header() -> ui.TagChild:
    """Create the application header section."""
    return ui.div(
        ui.h1("🔬 LIANA Results Explorer", class_="text-center"),
        ui.p("Interactive dashboard for exploring cell-cell interaction analysis results", 
             class_="text-center text-muted"),
        class_="content-header"
    )

def create_sidebar(data_dir: str = "/nfs/data/COST_IBD/downstream_tasks/interactions/output/HNC",
                  conditions: Optional[List[str]] = None,
                  splitting_keys: Optional[List[str]] = None,
                  cli_mode: bool = False) -> ui.TagChild:
    """
    Create the sidebar with core navigation controls.
    
    Args:
        conditions: List of available conditions
        splitting_keys: List of available splitting keys
        
    Returns:
        Sidebar UI component
    """
    conditions = conditions or []
    splitting_keys = splitting_keys or []
    
    return ui.sidebar(
        ui.div(
            # Core Navigation Section
            ui.h4("🗂️ Data Navigation"),
            ui.input_select("splitting_key", "Analysis Type:", 
                           choices=splitting_keys, selected=None),
            ui.br(),
            ui.input_select("contrast", "Select Contrast:", 
                           choices=conditions, selected=None),
            ui.br(),
            
            # Data Filters Section
            ui.h4("🔍 Data Filters"),
            ui.input_slider("logfc_threshold", "LogFC Threshold:", 
                           min=0, max=2.5, value=0.5, step=0.1),
            ui.input_slider("lrscore_threshold", "LRscore Threshold:", 
                           min=0, max=1, value=0.9, step=0.01),
            ui.input_slider("specificity_threshold", "Specificity Rank Threshold:", 
                           min=0, max=0.2, value=0.05, step=0.01),
            ui.input_slider("min_source_cells", "Min Source Cells:", 
                min=0, max=500, value=30, step=10),
            ui.input_slider("min_target_cells", "Min Target Cells:", 
                min=0, max=500, value=30, step=10),
            ui.br(),
            
            # Cell Type Selection Section
            ui.h4("🎯 Cell Type Selection"),
            ui.input_selectize("source_types", "Source Cell Types:", 
                              choices=[], selected=[], multiple=True),
            ui.input_selectize("target_types", "Target Cell Types:", 
                              choices=[], selected=[], multiple=True),
            ui.br(),
            
            class_="sidebar-content"
        ),
        width="320px"
    )


def create_network_tab() -> ui.TagChild:
    """Create the network visualization tab."""
    return ui.nav_panel(
        "🌐 Network", 
        ui.div(
            ui.div(
                ui.h5("Network Visualization Options"),
                ui.div(
                    ui.div(
                        ui.input_radio_buttons("network_options", "Display Options:",
                                             choices={
                                                 "all": "All Interactions",
                                                 "top": "Top Interactions Only"
                                             }, selected="all", inline=True),
                        class_="col-md-6"
                    ),
                    ui.div(
                        ui.input_numeric("network_top_n", "Top N (if selected):", 
                                       value=50, min=10, max=200, step=10),
                        class_="col-md-3"
                    ),
                    ui.div(
                        ui.input_select("network_layout", "Layout:",
                                       choices={
                                           "spring": "Spring Layout",
                                           "circular": "Circular Layout", 
                                           "kamada_kawai": "Kamada-Kawai Layout"
                                       }, selected="spring"),
                        class_="col-md-3"
                    ),
                    class_="row mb-3"
                ),
                class_="mb-3 p-2 border rounded"
            ),
            ui.h5("Original network"),
            ui.div(
                output_widget("network_plot"),
                style="height: 650px;"
            ),

            ui.hr(),

            ui.h5("Cell-count-aware network"),
            ui.p(
                "Node size = number of cells. Edge thickness = number of filtered LR interactions. "
                "This network shows the top N most connected cell types after filtering.",
                class_="text-muted"
            ),
            ui.input_slider(
                "cell_count_network_top_n",
                "Top N connected cell types to show:",
                min=5,
                max=100,
                value=30,
                step=5
            ),
            ui.div(
                output_widget("network_cell_count_plot"),
                style="height: 650px;"
            ),   
        )
    )

def create_heatmap_tab() -> ui.TagChild:
    """Create the heatmap visualization tab."""
    return ui.nav_panel(
        "🔥 Heatmap", 
        ui.div(
            ui.div(
                ui.h5("Heatmap Options"),
                ui.div(
                    ui.div(
                        ui.input_radio_buttons("heatmap_metric", "Metric:", 
                                             choices={
                                                 "interaction_count": "Interaction Count",
                                                 "lrscore": "LRscore", 
                                                 "lr_logfc": "LogFC",
                                                 "lr_means": "LR Means"
                                             }, selected="interaction_count", inline=True),
                        class_="col-md-6"
                    ),
                    ui.div(
                        ui.input_checkbox("show_heatmap_values", "Show Values", value=False),
                        class_="col-md-3"
                    ),
                    ui.div(
                        ui.input_select("colorscale", "Color Scheme:",
                                       choices={
                                           "Blues": "Blues",
                                           "Viridis": "Viridis",
                                           "Plasma": "Plasma",
                                           "Cividis": "Cividis"
                                       }, selected="Blues"),
                        class_="col-md-3"
                    ),
                    class_="row mb-3"
                ),
                class_="mb-3 p-2 border rounded"
            ),
            output_widget("heatmap_plot"),
            style="height: 650px;"
        )
    )

def create_dotplot_tab() -> ui.TagChild:
    """Create the dot plot visualization tab."""
    return ui.nav_panel(
        "🎯 Dot Plot", 
        ui.div(
            ui.div(
                ui.input_slider(
                    "top_n_interactions",
                    "Top N Interactions:",
                    min=5,
                    max=200,
                    value=20,
                    step=5
                ),
                ui.input_radio_buttons("dotplot_color", "Color By:",
                                     choices={
                                         "magnitude_rank": "Magnitude Rank",
                                         "specificity_rank": "Specificity Rank",
                                         "lrscore": "LRscore"
                                     }, selected="specificity_rank", inline=True),
                class_="mb-3"
            ),
            output_widget("dotplot"),
            style="height: 700px;"
        )
    )

def create_comparison_tab() -> ui.TagChild:
    """Create the comparison visualization tab."""
    return ui.nav_panel(
        "📈 Comparison", 
        ui.div(
            ui.div(
                ui.h5("Comparison Settings"),
                ui.div(
                    ui.div(
                        ui.input_select("comparison_source", "Source Cell Type:", 
                                       choices=[]),
                        class_="col-md-6"
                    ),
                    ui.div(
                        ui.input_select("comparison_target", "Target Cell Type:", 
                                       choices=[]),
                        class_="col-md-6"
                    ),
                    class_="row mb-3"
                ),
                ui.div(
                    ui.input_radio_buttons("comparison_metric", "Compare By:",
                                         choices={
                                             "lr_means": "LR Means",
                                             "lrscore": "LRscore",
                                             "lr_logfc": "LogFC"
                                         }, selected="lr_means", inline=True),
                    ui.input_numeric("comparison_max_interactions", "Max interactions per condition:", 
                                   value=100, min=10, max=1000, step=10),
                    class_="mb-3"
                ),
                class_="mb-3 p-2 border rounded"
            ),
            output_widget("comparison_plot"),
            style="height: 600px;"
        )
    )

def create_volcano_tab() -> ui.TagChild:
    """Create the volcano plot visualization tab."""
    return ui.nav_panel(
        "🌋 Volcano Plot",
        ui.div(
            ui.div(
                ui.input_select("volcano_x", "X-axis (Effect Size):",
                               choices={
                                   "lr_logfc": "LogFC",
                                   "lr_means": "LR Means"
                               }, selected="lr_logfc"),
                ui.input_select("volcano_y", "Y-axis (Significance):",
                               choices={
                                   "lrscore": "LRscore",
                                   "specificity_rank": "Specificity Rank"
                               }, selected="lrscore"),
                ui.input_numeric("volcano_threshold", "Significance Threshold:",
                               value=0.05, min=0.001, max=0.1, step=0.001),
                class_="mb-3"
            ),
            output_widget("volcano_plot"),
            style="height: 600px;"
        )
    )

def create_data_table_tab() -> ui.TagChild:
    """Create the data table tab."""
    return ui.nav_panel(
        "📋 Data Table", 
        ui.div(
            ui.div(
                ui.input_numeric("table_rows", "Rows to show:", 
                               value=100, min=10, max=1000, step=10),
                ui.input_text("table_search", "Search interactions:",
                            placeholder="Enter ligand, receptor, or cell type"),
                ui.input_action_button("export_data", "Export Filtered Data", 
                                     class_="btn-secondary"),
                class_="mb-3"
            ),
            ui.output_data_frame("interactions_table"),
            style="height: 600px; overflow-y: auto;"
        )
    )

def create_summary_tab() -> ui.TagChild:
    """Create the summary statistics tab."""
    return ui.nav_panel(
        "📊 Summary",
        ui.div(
            ui.div(
                ui.h4("Dataset Overview"),
                output_widget("summary_plot"),
                ui.br(),
                ui.h4("Top Cell Type Pairs"),
                ui.output_data_frame("top_pairs_table"),
                ui.br(),
                ui.h4("Top Ligand-Receptor Pairs"),
                ui.output_data_frame("top_lr_pairs_table"),
                class_="mb-3"
            ),
            style="height: 800px; overflow-y: auto;"
        )
    )

def create_data_explorer_tab() -> ui.TagChild:
    """Create the main overview / data explorer tab."""

    return ui.nav_panel(
        "📊 Overview",
        ui.div(
            ui.div(

                # ============================================================
                # DATA STRUCTURE
                # ============================================================

                ui.h4("Discovered Analysis Types"),
                ui.output_data_frame("splitting_keys_table"),

                ui.hr(),

                # ============================================================
                # EXISTING PLOT: RAW INTERACTION COUNTS
                # ============================================================

                ui.h4("Number of Interactions per Dataset"),

                ui.p(
                    "Raw number of LIANA result rows loaded for each dataset / condition.",
                    class_="text-muted"
                ),

                output_widget("structure_overview_plot"),

                ui.hr(),

                # ============================================================
                # PLOT 1
                # ============================================================

                ui.h4("Interaction Landscape Across Conditions"),

                ui.p(
                    "Top ligand–receptor interactions across all conditions, "
                    "ranked globally by magnitude.",
                    class_="text-muted"
                ),

                ui.input_numeric(
                    "overview_interaction_top_n",
                    "Top N interactions:",
                    value=50,
                    min=10,
                    max=200,
                    step=10
                ),
                
                ui.input_checkbox(
                    "overview_apply_filters",
                    "Apply sidebar filters",
                    value=False
                ),

                ui.output_ui("overview_condition_dotplot_container"),

                ui.hr(),

                # ============================================================
                # SHARED CONDITION COMPARISON FOR PLOTS 2 + 3
                # ============================================================

                ui.h4("Condition Differences"),

                ui.p(
                    "Select the two conditions used for the differential "
                    "cell-cell and ligand–receptor summaries below.",
                    class_="text-muted"
                ),

                ui.input_select(
                    "overview_condition_comparison",
                    "Compare conditions:",
                    choices={}
                ),

                ui.hr(),

                # ============================================================
                # PLOT 2
                # ============================================================

                ui.h4("1. Cell–Cell Communication Changes"),

                ui.p(
                    "Difference in LRscore between conditions, aggregated across "
                    "matched ligand–receptor interactions for each source → target pair.",
                    class_="text-muted"
                ),

                ui.input_radio_buttons(
                    "overview_cell_mode",
                    "Interactions to include:",
                    choices={
                        "all": "All matched interactions",
                        "filtered": "Apply sidebar filters",
                        "top": "Top N most changed cell pairs"
                    },
                    selected="top",
                    inline=True
                ),

                ui.panel_conditional(
                    "input.overview_cell_mode === 'top'",

                    ui.input_numeric(
                        "overview_cell_top_n",
                        "Top N cell pairs:",
                        value=30,
                        min=10,
                        max=1000,
                        step=10
                    )
                ),

                ui.output_ui("overview_cell_difference_plot_container"),

                ui.hr(),

                # ============================================================
                # PLOT 3
                # ============================================================

                ui.h4("2. Ligand–Receptor Changes"),

                ui.p(
                    "Difference in LRscore between conditions, aggregated across "
                    "all matched source → target cell-type contexts.",
                    class_="text-muted"
                ),

                ui.input_radio_buttons(
                    "overview_lr_diff_mode",
                    "Interactions to include:",
                    choices={
                        "all": "All matched interactions",
                        "filtered": "Apply sidebar filters",
                        "top": "Top N most changed LR pairs"
                    },
                    selected="top",
                    inline=True
                ),

                ui.panel_conditional(
                    "input.overview_lr_diff_mode === 'top'",

                    ui.input_numeric(
                        "overview_lr_diff_top_n",
                        "Top N LR pairs:",
                        value=30,
                        min=10,
                        max=1000,
                        step=10
                    )
                ),

                ui.output_ui("overview_lr_difference_plot_container"),

                ui.br(),
                ui.br(),

                class_="mb-3"
            ),

            style="padding-right: 10px;"
        )
    )

def create_main_content() -> ui.TagChild:
    """Create the main content area with all tabs."""
    return ui.div(
        ui.div(
            # Main tabset with all visualizations
            ui.navset_tab(
                create_data_explorer_tab(),
                create_network_tab(),
                create_heatmap_tab(),
                create_dotplot_tab(),
                create_comparison_tab(),
                create_data_table_tab(),
                create_summary_tab()
            ),
            class_="main-content"
        )
    )

def create_app_styles() -> ui.TagChild:
    """Create custom CSS styles for the app."""
    return ui.tags.head(
        ui.tags.style("""
            .content-header { 
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
                padding: 20px; 
                border-radius: 10px; 
                margin-bottom: 20px;
                box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
            }
            .content-header h1 {
                margin-bottom: 10px;
                font-weight: 300;
            }
            .metric-box { 
                background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
                padding: 20px; 
                border-radius: 10px; 
                text-align: center;
                margin: 5px;
                box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
                transition: transform 0.2s ease;
            }
            .metric-box:hover {
                transform: translateY(-2px);
            }
            .metric-box h5 {
                color: #2c3e50;
                font-weight: bold;
                font-size: 1.5em;
                margin-bottom: 5px;
            }
            .metric-box p {
                color: #7f8c8d;
                margin: 0;
                font-size: 0.9em;
            }
            .sidebar {
                background-color: #f8f9fa;
                border-right: 1px solid #dee2e6;
                height: 100vh;
                overflow: hidden;
            }
            .sidebar-content {
                height: 100%;
                overflow-y: auto;
                padding: 20px;
                padding-bottom: 40px;
            }
            .main-content {
                height: calc(100vh - 120px);
                overflow: hidden;
                padding: 20px;
                padding-bottom: 40px;
                display: flex;
                flex-direction: column;
            }

            .main-content .tab-content {
                overflow-y: auto;
                overflow-x: hidden;
                flex: 1 1 auto;
                min-height: 0;
            }
            .sidebar h4 {
                color: #495057;
                border-bottom: 2px solid #007bff;
                padding-bottom: 5px;
                margin-top: 20px;
                margin-bottom: 15px;
            }
            .btn-primary {
                background: linear-gradient(135deg, #007bff 0%, #0056b3 100%);
                border: none;
                border-radius: 20px;
                padding: 8px 20px;
                transition: all 0.3s ease;
            }
            .btn-primary:hover {
                transform: translateY(-1px);
                box-shadow: 0 4px 8px rgba(0, 123, 255, 0.3);
            }
            .btn-secondary {
                background: linear-gradient(135deg, #6c757d 0%, #495057 100%);
                border: none;
                border-radius: 20px;
                padding: 6px 15px;
                color: white;
            }
            .nav-tabs .nav-link {
                border-radius: 20px 20px 0 0;
                margin-right: 5px;
            }
            .nav-tabs .nav-link.active {
                background: linear-gradient(135deg, #007bff 0%, #0056b3 100%);
                color: white;
                border-color: #007bff;
            }
            .form-control, .form-select {
                border-radius: 10px;
                border: 1px solid #ced4da;
                transition: border-color 0.3s ease;
            }
            .form-control:focus, .form-select:focus {
                border-color: #007bff;
                box-shadow: 0 0 0 0.2rem rgba(0, 123, 255, 0.25);
            }
            .alert {
                border-radius: 10px;
                border: none;
            }
            .alert-warning {
                background: linear-gradient(135deg, #fff3cd 0%, #ffeaa7 100%);
                color: #856404;
            }
            .alert-info {
                background: linear-gradient(135deg, #d1ecf1 0%, #a8dadc 100%);
                color: #0c5460;
            }
            .alert-success {
                background: linear-gradient(135deg, #d4edda 0%, #95e1a3 100%);
                color: #155724;
            }
            .alert-danger {
                background: linear-gradient(135deg, #f8d7da 0%, #f5b7b1 100%);
                color: #721c24;
            }
            .plotly-graph-div {
                border-radius: 10px;
                box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
            }
            /* Custom scrollbar styling */
            .sidebar-content::-webkit-scrollbar,
            .main-content::-webkit-scrollbar {
                width: 8px;
            }
            .sidebar-content::-webkit-scrollbar-track,
            .main-content::-webkit-scrollbar-track {
                background: #f1f1f1;
                border-radius: 4px;
            }
            .sidebar-content::-webkit-scrollbar-thumb,
            .main-content::-webkit-scrollbar-thumb {
                background: #c1c1c1;
                border-radius: 4px;
            }
            .sidebar-content::-webkit-scrollbar-thumb:hover,
            .main-content::-webkit-scrollbar-thumb:hover {
                background: #a8a8a8;
            }
            /* Ensure proper layout */
            .shiny-layout-sidebar {
                height: 100vh;
                overflow: hidden;
            }
        """)
    )

def create_full_ui(data_dir: str = "/nfs/data/COST_IBD/downstream_tasks/interactions/output/HNC", cli_mode: bool = False) -> ui.TagChild:
    """Create the complete UI layout."""
    return ui.page_fluid(
        create_app_styles(),
        
        # Main layout with sidebar and content
        ui.layout_sidebar(
            create_sidebar(data_dir=data_dir, cli_mode=cli_mode),
            create_main_content()
        )
    )
