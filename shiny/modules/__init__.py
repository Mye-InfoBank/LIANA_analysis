"""
LIANA Results Explorer - Modular Components
This package contains the modular components for the LIANA Results Explorer Shiny app.
"""

__version__ = "1.0.0"
__author__ = "COST IBD Project"

# Import main components for easy access
from .data_handler import DataHandler
from .data_handler import DataHandler
from .visualizations import (
    create_network_plot,
    create_cell_count_network_plot,
    create_heatmap_plot,
    create_dotplot,
    create_scatter_comparison,
    create_volcano_plot,
    create_summary_barplot,
    create_structure_overview_plot,
    create_cross_analysis_comparison,
)
from .ui_components import create_sidebar, create_main_content, create_full_ui
from .server_logic import create_server_function
from .utils import setup_logging, validate_data_format

__all__ = [
    "DataHandler",
    "create_network_plot",
    "create_cell_count_network_plot",
    "create_heatmap_plot",
    "create_dotplot",
    "create_scatter_comparison",
    "create_volcano_plot",
    "create_summary_barplot",
    "create_structure_overview_plot",
    "create_cross_analysis_comparison",
    "create_sidebar",
    "create_main_content",
    "create_full_ui",
    "create_server_function",
    "setup_logging",
    "validate_data_format",
]
