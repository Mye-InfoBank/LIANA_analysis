# 🚀 Command-Line Interface for LIANA Results Explorer

## 🎯 Overview

The LIANA Results Explorer now supports a powerful command-line interface that allows you to specify the top-level data directory directly when launching the app. **When a data directory is provided via command line, the UI removes the data directory input field entirely, showing only the path as a display element. This ensures data integrity and provides a cleaner interface.**

## 📋 Usage

### **Basic Usage**

```bash
# Launch with data directory
python app.py /nfs/data/COST_IBD/downstream_tasks/interactions/output/HNC

# Using the run script
./run_app.sh /nfs/data/COST_IBD/downstream_tasks/interactions/output/HNC
```

### **Advanced Options**

```bash
# Custom port and host
python app.py --port 8080 --host 0.0.0.0 /path/to/liana/results

# Debug mode
python app.py --debug /path/to/liana/results
```

## 🔧 Command-Line Arguments

### **Required Arguments**

| Argument   | Description                                           | Example                    |
| ---------- | ----------------------------------------------------- | -------------------------- |
| `DATA_DIR` | Top-level directory containing LIANA analysis results | `/nfs/data/.../output/HNC` |

### **Optional Arguments**

| Option    | Short | Description            | Default     |
| --------- | ----- | ---------------------- | ----------- |
| `--port`  | `-p`  | Port to run the app on | `8000`      |
| `--host`  | `-H`  | Host to bind to        | `127.0.0.1` |
| `--debug` |       | Enable debug logging   | `False`     |
| `--help`  |       | Show help message      |             |

## 🎯 Examples

### **Example 1: Basic HNC Analysis**

```bash
./run_app.sh /nfs/data/COST_IBD/downstream_tasks/interactions/output/HNC
```

- Launches app on `http://127.0.0.1:8000`
- Pre-fills the data directory path
- User clicks "Discover Data Structure" to find HPV, EBV, classification, sex, site analyses

### **Example 2: Custom Network Configuration**

```bash
./run_app.sh --host 0.0.0.0 --port 8080 /nfs/data/COST_IBD/downstream_tasks/interactions/output/HNC
```

- Makes app accessible from other machines on the network
- Runs on port 8080
- Useful for remote access or sharing with colleagues

### **Example 3: Development Mode**

```bash
python app.py --debug /path/to/liana/results
```

- Enables detailed debug logging
- Useful for troubleshooting and development

## 🖥️ UI Behavior

### **When Data Directory is Set via Command Line:**

- ✅ **Display Only**: The data directory is shown as a display element (not editable)
- ✅ **No Input Field**: The data directory input field is completely removed
- ✅ **No Discovery Button**: The "Discover Data Structure" button is hidden
- ✅ **Clean Interface**: Simplified UI with no unnecessary input elements
- ✅ **Automatic Discovery**: Data structure is always discovered on startup

### **When No Data Directory is Provided:**

- 📝 **Editable Input**: Users can manually enter the data directory path
- 🔍 **Discovery Button**: Users can click "Discover Data Structure" to find analyses
- 📋 **Manual Process**: Standard workflow for interactive data exploration

## 🔄 Workflow Comparison

### **Before (UI-Only)**

1. Launch app: `python app.py`
2. Open browser: `http://127.0.0.1:8000`
3. Manually enter data directory path
4. Click "Discover Data Structure"
5. Select analysis type
6. Start exploring

### **After (Command-Line)**

1. Launch with path: `./run_app.sh /path/to/data`
2. Open browser: `http://127.0.0.1:8000`
3. Path is displayed (not editable), data automatically discovered
4. Select analysis type
5. Start exploring

## 🚀 Benefits

### **For Users**

- **Faster Setup**: No need to manually enter paths
- **Clean Interface**: No unnecessary input fields when path is pre-set
- **Data Integrity**: No possibility of accidental path changes
- **Automatic Discovery**: Data structure is always discovered on startup
- **Batch Processing**: Easy to script for multiple datasets
- **Remote Access**: Can bind to different hosts/ports

### **For Automation**

- **Scriptable**: Can be integrated into analysis pipelines
- **Configurable**: Flexible host/port configuration
- **Error Handling**: Validates data directory before launch
- **Logging**: Debug mode for troubleshooting

### **For Development**

- **Quick Testing**: Easy to test with different datasets
- **Debug Mode**: Detailed logging for development
- **Flexible Configuration**: All options available via command line

## 🔧 Integration Examples

### **Bash Script for Multiple Datasets**

```bash
#!/bin/bash
# Process multiple LIANA result directories

datasets=(
    "/nfs/data/COST_IBD/downstream_tasks/interactions/output/HNC"
    "/nfs/data/COST_IBD/downstream_tasks/interactions/output/CRC"
    "/nfs/data/COST_IBD/downstream_tasks/interactions/output/LUAD"
)

for dataset in "${datasets[@]}"; do
    echo "Processing: $dataset"
    ./run_app.sh --auto-discover "$dataset" &
    sleep 5  # Wait for app to start
    # Add your automation logic here
done
```

### **Python Script Integration**

```python
import subprocess
import time

def launch_liana_explorer(data_dir, port=8000, auto_discover=True):
    """Launch LIANA Results Explorer for a dataset."""
    cmd = [
        "python", "app.py",
        "--port", str(port),
        "--auto-discover" if auto_discover else "",
        data_dir
    ]

    # Remove empty strings
    cmd = [arg for arg in cmd if arg]

    process = subprocess.Popen(cmd)
    time.sleep(3)  # Wait for startup
    return process

# Usage
launch_liana_explorer("/path/to/liana/results", port=8080)
```

## 🛠️ Troubleshooting

### **Common Issues**

1. **Data directory not found**

   ```bash
   ❌ Error: Data directory does not exist: /path/to/data
   ```

   **Solution**: Check the path and ensure it exists

2. **Port already in use**

   ```bash
   ❌ Error: Port 8000 is already in use
   ```

   **Solution**: Use a different port with `--port 8080`

3. **Permission denied**
   ```bash
   ❌ Error: Permission denied
   ```
   **Solution**: Ensure the script is executable: `chmod +x run_app.sh`

### **Debug Mode**

```bash
python app.py --debug /path/to/data
```

- Enables detailed logging
- Shows data discovery process
- Useful for troubleshooting issues

## 📚 Related Documentation

- [Hierarchical Browsing Update](HIERARCHICAL_BROWSING_UPDATE.md) - Details on the hierarchical browsing features
- [Modular Structure](MODULAR_STRUCTURE.md) - Information about the modular architecture
- [Plotly Fix](PLOTLY_FIX.md) - Details on Plotly rendering fixes

---

## 🎉 **Ready to Use!**

The command-line interface makes the LIANA Results Explorer much more convenient and powerful. You can now:

- **Launch quickly** with pre-filled data paths
- **Automate workflows** with scripting
- **Share easily** with remote access options
- **Debug effectively** with detailed logging

**Perfect for exploring your HNC dataset with HPV, EBV, classification, sex, and site analyses! 🚀**
