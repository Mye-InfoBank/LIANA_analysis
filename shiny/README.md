# LIANA Results Explorer 🔬

An interactive Python Shiny dashboard for exploring cell-cell interaction results from the LIANA pipeline. The app is designed to visualize ligand-receptor interactions across source and target cell types, compare condition-specific LIANA results, and interactively filter interactions using LIANA-derived ranking and magnitude metrics.

## Features

### 📊 Interactive Visualizations

* **Data Explorer**: Overview of loaded LIANA result files and number of result rows per dataset/condition run.
* **Network Plot**: Interactive cell-cell interaction networks based on filtered ligand-receptor interactions.
* **Cell-count-aware Network**: Network where node size reflects the number of cells in each cell type and edge thickness reflects the number of filtered ligand-receptor interaction rows.
* **Heatmap**: Source-target interaction heatmaps using either interaction counts or selected LIANA metrics.
* **Dot Plot**: Top ligand-receptor interactions with enrichment, score, and specificity/ranking information.
* **Comparison Plot**: Cross-condition comparison for a selected source-target cell-type pair.
* **Data Table**: Searchable table of filtered ligand-receptor interactions.
* **Summary Tab**: Summary statistics and top cell-type / ligand-receptor pairs after filtering.

### 🎛️ Dynamic Controls

* **Analysis Type**: Select the available LIANA analysis type, for example `condition`.
* **Contrast Selection**: Switch between the full dataset run and condition-specific runs such as `HC`, `UC`, and `CD`.
* **Metric Filtering**: Adjust LogFC, LRscore, and specificity rank thresholds.
* **Cell Count Filtering**: Exclude interactions involving very small source or target cell-type groups.
* **Cell Type Selection**: Focus on specific source and target cell types.
* **Interactive Parameters**: Customize plot-specific options such as top N interactions, network layout, heatmap metric, and comparison metric.

## Navigation Guide

### 🗂️ Data Explorer

The Data Explorer shows:

* Discovered analysis type(s)
* Available condition-specific result files
* Number of LIANA result rows per loaded dataset/condition run

The bar plot shows the number of rows in each loaded `liana_results.csv` file. These counts describe result-file size, not biological interaction strength.

### 🌐 Network Tab

The Network tab contains two network views.

#### Original network

* Nodes represent cell types.
* Edges represent source-target cell-type pairs.
* Edge thickness reflects the number of filtered ligand-receptor rows between the source and target cell type.
* The plot updates after changing the global filters.

#### Cell-count-aware network

* Node size represents the number of cells in the corresponding cell type.
* Edge thickness represents the number of filtered ligand-receptor interaction rows between source and target.
* The Top N slider keeps the most connected cell types, where connectivity is the number of filtered interaction rows attached to the node.
* Hovering over nodes shows the cell type, cell count, and number of connected filtered ligand-receptor interactions.

This plot is useful for checking whether highly connected interaction patterns are driven by large cell populations or by very small cell-type groups.

### 🔥 Heatmap Tab

The heatmap summarizes filtered source-target interactions.

Available metrics:

* **Interaction Count**: Number of filtered ligand-receptor rows per source-target cell-type pair.
* **LRscore**: Median LRscore per source-target pair.
* **LogFC**: Median ligand-receptor log fold change per source-target pair.
* **LR Means**: Median ligand-receptor mean expression score per source-target pair.

Recommended default: **Interaction Count**, because it gives the clearest overview of how many filtered ligand-receptor pairs connect each source-target cell-type pair.

### 🎯 Dot Plot Tab

The dot plot shows the top N ligand-receptor interactions after applying the global filters.

Current interpretation:

* **X-axis**: `lr_logfc`

  * Higher values indicate stronger enrichment of the ligand-receptor pair in the selected source-target context.
* **Dot size**: `lrscore`

  * Larger dots indicate stronger LIANA interaction score.
* **Dot color / ranking option**: selected ranking or score metric.

  * `specificity_rank`: lower values are more specific.
  * `magnitude_rank`: lower values indicate stronger magnitude rank.
  * `lrscore`: higher values indicate stronger LIANA score.
* **Y-axis**: source cell type → target cell type | ligand → receptor

Recommended default:

```text
Rank / Color By: Specificity Rank
Top N Interactions: 20
```

This shows interactions that are specific, strong, and enriched.

### 📈 Comparison Tab

The Comparison tab compares ligand-receptor interactions across condition-specific LIANA runs for one selected source-target cell-type pair.

Controls:

* **Source Cell Type**: ligand-producing cell type.
* **Target Cell Type**: receptor-expressing cell type.
* **Compare By**:

  * `lr_means`
  * `lrscore`
  * `lr_logfc`
* **Max interactions per condition**: limits the number of displayed ligand-receptor pairs per condition.

This plot is useful for comparing whether a selected cell-cell interaction context differs between conditions such as `HC`, `UC`, and `CD`.

### 📋 Data Table Tab

The Data Table shows the currently filtered ligand-receptor interaction rows.

The table can be searched by ligand, receptor, source cell type, or target cell type.

### 📊 Summary Tab

The Summary tab updates after applying filters and shows:

* Number of filtered interactions
* Number of unique source-target cell-type pairs
* Number of unique ligands and receptors
* Number of unique source and target cell types
* Top source-target cell-type pairs
* Top ligand-receptor pairs

## Filtering and Metric Interpretation

The app applies global filters before generating most plots. These filters are intended to keep interactions that are enriched, strong, specific, and supported by enough cells.

### Recommended Default Filter Settings

| Filter                     | Default | Direction                           | Interpretation                             |
| -------------------------- | ------: | ----------------------------------- | ------------------------------------------ |
| LogFC Threshold            |   `0.5` | keep `lr_logfc > threshold`         | Keeps enriched ligand-receptor pairs       |
| LRscore Threshold          |   `0.9` | keep `lrscore > threshold`          | Keeps high-scoring ligand-receptor pairs   |
| Specificity Rank Threshold |  `0.05` | keep `specificity_rank < threshold` | Keeps specific ligand-receptor pairs       |
| Min Source Cells           |    `30` | keep `source_n_cells >= threshold`  | Removes poorly supported source cell types |
| Min Target Cells           |    `30` | keep `target_n_cells >= threshold`  | Removes poorly supported target cell types |

### Suggested Threshold Presets

| Preset       | LogFC | LRscore | Specificity Rank | Min Source Cells | Min Target Cells | Use case                                 |
| ------------ | ----: | ------: | ---------------: | ---------------: | ---------------: | ---------------------------------------- |
| Exploratory  | `0.5` |   `0.8` |           `0.10` |      `0` or `30` |      `0` or `30` | Broad exploration                        |
| Default      | `0.5` |   `0.9` |           `0.05` |             `30` |             `30` | General app use                          |
| Strict       | `1.0` |   `0.9` |           `0.02` |     `30` or `50` |     `30` or `50` | More selective candidate interactions    |
| Figure-level | `1.0` |  `0.95` |           `0.01` |             `50` |             `50` | Highly filtered interactions for figures |

These are practical defaults, not universal biological cutoffs. The appropriate threshold depends on dataset size, cell-type granularity, and the biological question.

### Metric Details

#### `lr_logfc`

`lr_logfc` describes ligand-receptor enrichment in the selected source-target context.

* Higher is better for enriched interactions.
* Positive values indicate enrichment.
* Negative values indicate depletion or lower relative enrichment.
* The app filters with:

```text
lr_logfc > threshold
```

Recommended slider range:

```text
0 to 2.5, default 0.5, step 0.1
```

#### `lrscore`

`lrscore` is a LIANA interaction score. It is commonly interpreted as an interaction strength or magnitude score.

* Higher is better.
* Values are usually between 0 and 1.
* The app filters with:

```text
lrscore > threshold
```

Recommended slider range:

```text
0 to 1, default 0.9, step 0.01
```

#### `specificity_rank`

`specificity_rank` is a rank-based specificity metric.

* Lower is better.
* Values closer to 0 indicate more specific interactions.
* The app filters with:

```text
specificity_rank < threshold
```

Recommended slider range:

```text
0 to 0.2, default 0.05, step 0.01
```

#### `magnitude_rank`

`magnitude_rank` is a rank-based magnitude metric.

* Lower is better.
* Values closer to 0 indicate stronger magnitude ranking.
* It is used for ranking and visualization, but not currently used as a global sidebar filter.

#### `lr_means`

`lr_means` summarizes ligand and receptor expression magnitude.

* Higher values indicate higher average ligand-receptor expression.
* It is useful for comparison plots and heatmaps.
* It is dataset-dependent and does not have a fixed universal maximum.

#### `interaction_count`

`interaction_count` is not a LIANA score. It is the number of ligand-receptor result rows after filtering.

* Higher means more filtered ligand-receptor pairs connect a source-target cell-type pair.
* It is useful for overview plots, networks, and heatmaps.
* It should not be interpreted as molecular interaction strength by itself.

#### `source_n_cells` and `target_n_cells`

These columns give the number of cells in the source or target cell type for the selected dataset/condition.

* Higher values mean better cell-count support.
* Filtering small cell groups helps avoid interactions driven by very rare clusters.
* The app filters with:

```text
source_n_cells >= threshold
target_n_cells >= threshold
```

Recommended slider range:

```text
0 to 500, default 30, step 10
```

A value of `0` effectively turns off the cell-count filter.

## Data Format

The app expects LIANA result CSV files with the following core columns:

| Column             | Description                                                         |
| ------------------ | ------------------------------------------------------------------- |
| `source`           | Source cell type, interpreted as ligand-producing cell type         |
| `target`           | Target cell type, interpreted as receptor-expressing cell type      |
| `ligand_complex`   | Ligand or ligand complex                                            |
| `receptor_complex` | Receptor or receptor complex                                        |
| `lrscore`          | LIANA interaction score; higher is stronger                         |
| `lr_logfc`         | Ligand-receptor log fold change/enrichment; higher is more enriched |
| `lr_means`         | Mean expression-based ligand-receptor magnitude metric              |
| `specificity_rank` | Specificity rank; lower is more specific                            |
| `magnitude_rank`   | Magnitude rank; lower is stronger                                   |
| `source_n_cells`   | Number of cells in the source cell type                             |
| `target_n_cells`   | Number of cells in the target cell type                             |

The `source_n_cells` and `target_n_cells` columns are recommended for robust filtering, especially when working with fine-grained cell-type annotations.


Happy exploring! 🚀


## Installation and Developement

### Prerequisites

* Python 3.8 or higher
* LIANA pipeline result files in CSV format
* A directory structure containing one main `liana_results.csv` and optional condition-specific subdirectories

### Setup

1. Clone the repository or navigate to the Shiny app directory:

   ```bash
   cd /path/to/interactions/shiny
   ```

2. Create and activate a virtual environment:

   ```bash
   python -m venv venv
   source venv/bin/activate
   ```

   On Windows:

   ```bash
   venv\Scripts\activate
   ```

3. Install dependencies:

   ```bash
   pip install -r requirements.txt
   ```

## Usage

### Starting the App

From the `shiny/` directory:

```bash
python app.py /path/to/liana/output
```

or using the provided shell script:

```bash
./run_app.sh /path/to/liana/output
```

The app will start and display a local URL, usually:

```text
http://127.0.0.1:8080
```

### File Structure

```text
shiny/
├── app.py
├── requirements.txt
├── README.md
├── run_app.sh
├── modules/
│   ├── __init__.py
│   ├── data_handler.py
│   ├── visualizations.py
│   ├── ui_components.py
│   ├── server_logic.py
│   └── utils.py
└── data/
    └── ...
```

### Main Components

* `data_handler.py`: data discovery, loading, filtering, and cell-type extraction
* `visualizations.py`: Plotly visualization functions
* `ui_components.py`: Shiny UI layout
* `server_logic.py`: Shiny reactive server logic
* `utils.py`: helper functions, logging, validation, and export utilities

When running on a remote server, this address is local to the server session and is usually accessed through SSH port forwarding or the VS Code forwarded ports interface.

## Expected Directory Structure (for LIANA runs and developers)

The app expects a LIANA output directory with one analysis type folder, for example:

```text
IBD_v11/
└── colon/
    └── condition/
        └── condition/
            ├── liana_results.csv
            ├── HC/
            │   └── liana_results.csv
            ├── UC/
            │   └── liana_results.csv
            └── CD/
                └── liana_results.csv
```

The top-level `liana_results.csv` is treated as the `full` run. Condition folders such as `HC`, `UC`, and `CD` are treated as separate condition-specific LIANA runs.

Important: the `full` run is not the union of the condition-specific result rows. It is an independently computed LIANA result on the full dataset. Therefore, the number of result rows in a condition-specific run can be higher or lower than in the full run.