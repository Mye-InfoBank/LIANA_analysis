#!/bin/bash
# Launch script for LIANA Results Explorer
# This script sets up the environment and launches the Shiny app.

# Default values
PORT=8080
HOST="127.0.0.1"
DATA_DIR=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -p|--port)
            PORT="$2"
            shift 2
            ;;
        -h|--host)
            HOST="$2"
            shift 2
            ;;
        -d|--data-dir)
            DATA_DIR="$2"
            shift 2
            ;;
        --help)
            echo "Usage: $0 [OPTIONS] DATA_DIR"
            echo ""
            echo "Arguments:"
            echo "  DATA_DIR              Top-level directory containing LIANA analysis results"
            echo ""
            echo "Options:"
            echo "  -p, --port PORT       Port to run the app on (default: 8000)"
            echo "  -h, --host HOST       Host to bind to (default: 127.0.0.1)"
            echo "  -d, --data-dir DIR    Data directory (alternative to positional argument)"
            echo "  --help                Show this help message"
            echo ""
            echo "Examples:"
            echo "  $0 /nfs/data/COST_IBD/downstream_tasks/interactions/output/HNC"
            echo "  $0 --port 8080 /path/to/liana/results"
            exit 0
            ;;
        *)
            # If DATA_DIR is not set, use this as the data directory
            if [ -z "$DATA_DIR" ]; then
                DATA_DIR="$1"
            else
                echo "Unknown option: $1"
                echo "Use --help for usage information"
                exit 1
            fi
            shift
            ;;
    esac
done

# Check if data directory is provided
if [ -z "$DATA_DIR" ]; then
    echo "❌ Error: Data directory is required"
    echo "Usage: $0 [OPTIONS] DATA_DIR"
    echo "Use --help for usage information"
    exit 1
fi

# Check if data directory exists
if [ ! -d "$DATA_DIR" ]; then
    echo "❌ Error: Data directory does not exist: $DATA_DIR"
    exit 1
fi

echo "🔬 Starting LIANA Results Explorer..."
echo "📊 Data directory: $DATA_DIR"
echo "🌐 Dashboard will be available at: http://${HOST}:${PORT}"
echo ""

# Check if virtual environment exists
if [ -d "venv" ]; then
    echo "📦 Activating virtual environment..."
    source venv/bin/activate
fi

# Check if requirements are installed
echo "🔍 Checking dependencies..."
python -c "import shiny, pandas, plotly, networkx, shinywidgets" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "❌ Missing dependencies. Installing..."
    pip install -r requirements.txt
fi

# Build command
CMD="python app.py --host $HOST --port $PORT \"$DATA_DIR\""

# Launch the app
echo "🚀 Launching app..."
echo "Command: $CMD"
eval $CMD
