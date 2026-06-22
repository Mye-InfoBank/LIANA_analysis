# LIANA Results Explorer 🔬

An interactive Python Shiny dashboard for exploring cell-cell interaction analysis results from the LIANA pipeline.

## Features

### 📊 Interactive Visualizations

- **Network Plot**: Interactive network visualization of cell-cell interactions
- **Heatmap**: Customizable heatmaps showing interaction strength across cell types
- **Dot Plot**: Top interactions ranked by magnitude with size and color encoding
- **Comparison Plot**: Cross-condition comparison of specific cell type interactions
- **Data Table**: Searchable and filterable table of all interactions

### 🎛️ Dynamic Controls

- **Data Loading**: Load results from any LIANA pipeline output directory
- **Contrast Selection**: Switch between full dataset and condition-specific results
- **Filtering**: Adjust thresholds for LogFC, LRscore, and specificity
- **Cell Type Selection**: Focus on specific source and target cell types
- **Interactive Parameters**: Customize plot parameters in real-time

### 📈 Analysis Features

- Real-time filtering and visualization updates
- Summary statistics for loaded datasets
- Cross-condition interaction comparisons
- Export-ready high-quality plots
- Responsive design that works on different screen sizes

## Installation

### Prerequisites

- Python 3.8 or higher
- LIANA pipeline results (CSV files and directory structure)

### Setup

1. **Clone or navigate to the shiny directory:**

   ```bash
   cd /path/to/interactions/shiny
   ```

2. **Create a virtual environment (recommended):**

   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

## Usage

### Starting the App

```bash
python app.py
```

The app will start and display the URL (typically `http://127.0.0.1:8000`).

### Loading Data

1. **Enter the path** to your LIANA results directory in the "Results Directory Path" field
2. **Click "Load Data"** to load the results
3. The app will automatically detect:
   - Main dataset results (`liana_results.csv`)
   - Condition-specific results (subdirectories with `liana_results.csv`)
   - Available cell types
   - Summary statistics

### Expected Directory Structure

```
your_liana_results/
├── liana_results.csv                    # Full dataset results
├── liana_results.pkl                    # Full dataset (pickle format)
├── full_interaction_network.png         # Generated plots
├── full_interaction_heatmap_*.png
├── full_detailed_interactions.png
├── condition_scatter_*.png
├── dotplot_top10_magnitude_rank.png
└── [CONDITION_NAME]/                    # Condition-specific results
    ├── liana_results.csv
    ├── liana_results.pkl
    └── [CONDITION_NAME]_*.png
```

### Navigation Guide

#### 🌐 Network Tab

- Visualizes cell-cell interaction networks
- Node size represents cell type importance
- Edge thickness shows interaction strength
- Interactive: hover for details, zoom, pan

#### 🔥 Heatmap Tab

- Choose metrics: LRscore, LogFC, or LR Means
- Rows = source cell types, Columns = target cell types
- Color intensity shows interaction strength
- Interactive hover shows exact values

#### 🎯 Dot Plot Tab

- Shows top N interactions (adjustable)
- Dot size = specificity (larger = more specific)
- Color = magnitude rank (darker = higher rank)
- Y-axis shows ligand-receptor pairs

#### 📈 Comparison Tab

- Compare interactions across conditions
- Select specific source and target cell types
- Scatter plot shows condition differences
- Useful for identifying condition-specific interactions

#### 📋 Data Table Tab

- Searchable and sortable table
- Adjust number of rows displayed
- Shows all interaction details
- Filtered by current settings

### Filtering and Controls

#### Data Loading Controls

- **Results Directory Path**: Path to LIANA output directory
- **Load Data**: Button to load/reload data

#### Analysis Controls

- **Select Contrast**: Choose between full dataset or specific conditions
- **LogFC Threshold**: Minimum log fold change (default: 0.5)
- **LRscore Threshold**: Minimum LR score (default: 0.9)
- **Specificity Threshold**: Maximum specificity rank (default: 0.05)

#### Cell Type Selection

- **Source Cell Types**: Filter interactions by source cell types
- **Target Cell Types**: Filter interactions by target cell types
- Leave empty to include all cell types

#### Comparison Settings

- **Comparison Source/Target**: Select cell types for cross-condition comparison

## Data Format

The app expects LIANA results with the following columns:

- `source`: Source cell type
- `target`: Target cell type
- `ligand_complex`: Ligand name/complex
- `receptor_complex`: Receptor name/complex
- `lrscore`: LIANA score
- `lr_logfc`: Log fold change
- `lr_means`: Mean expression
- `specificity_rank`: Specificity ranking
- `magnitude_rank`: Magnitude ranking

## Troubleshooting

### Common Issues

1. **"Directory not found" error**

   - Check the path is correct and accessible
   - Ensure the directory contains `liana_results.csv`

2. **"No LIANA results found" error**

   - Verify the directory structure matches expected format
   - Check that CSV files are not corrupted

3. **Empty plots**

   - Adjust filtering thresholds (they might be too stringent)
   - Check that the selected contrast has data
   - Verify cell type selections include available types

4. **App won't start**
   - Check Python version (3.8+ required)
   - Verify all dependencies are installed: `pip install -r requirements.txt`
   - Check for port conflicts (try different port)

### Performance Tips

- Large datasets (>10k interactions) may take time to load
- Consider filtering data before visualization for better performance
- Close unused browser tabs to free memory

## Development

### Extending the App

The app is modular and can be extended with:

- Additional plot types
- New filtering options
- Export functionality
- Statistical tests
- Custom analysis modules

### File Structure

```
shiny/
├── app.py                          # Main application entry point
├── app_original.py                 # Backup of original monolithic version
├── requirements.txt                # Dependencies
├── README.md                      # This file
├── run_app.sh                     # Launch script
├── validate_setup.py              # Setup validation
├── example_config.py              # Sample data generator
├── OVERVIEW.md                    # Project overview
├── MODULAR_STRUCTURE.md           # Modular architecture documentation
└── modules/                       # Modular components
    ├── __init__.py               # Package initialization
    ├── data_handler.py           # Data loading and processing
    ├── visualizations.py         # Plot generation functions
    ├── ui_components.py          # UI layout components
    ├── server_logic.py           # Server-side reactive logic
    └── utils.py                  # Utility functions and helpers
```

## Support

For issues related to:

- **LIANA pipeline**: Check the main pipeline documentation
- **Shiny app**: Review this README and check the console for error messages
- **Data format**: Ensure your LIANA results match the expected format

## Version Information

- **Python Shiny**: 0.6.0+
- **Plotly**: 5.15.0+
- **Pandas**: 2.0.0+
- **NetworkX**: 3.0+

---

**Happy exploring! 🚀**
