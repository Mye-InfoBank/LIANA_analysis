#!/usr/bin/env python3
"""
LIANA Results Explorer - Modular Interactive Shiny Dashboard
This app provides interactive exploration of LIANA cell-cell interaction analysis results.

This is the main entry point that orchestrates all modular components.
"""

import sys
import argparse
from pathlib import Path
from shiny import App
import logging

# Import modular components
from modules.data_handler import DataHandler
from modules.ui_components import create_full_ui
from modules.server_logic import create_server_function
from modules.utils import setup_logging
from modules.ai_snapshots import SNAPSHOT_DIR

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="LIANA Results Explorer - Interactive Dashboard for cell-cell interaction analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python app.py /nfs/data/COST_IBD/downstream_tasks/interactions/output/HNC
  python app.py --port 8080 /path/to/liana/results
  python app.py --host 0.0.0.0 --port 8000 /path/to/results
  python app.py --auto-discover /path/to/results
        """
    )
    
    parser.add_argument(
        "data_dir",
        help="Top-level directory containing LIANA analysis results (e.g., /path/to/output/HNC)"
    )
    
    parser.add_argument(
        "--port", "-p",
        type=int,
        default=8000,
        help="Port to run the app on (default: 8000)"
    )
    
    parser.add_argument(
        "--host", "-H",
        default="127.0.0.1",
        help="Host to bind to (default: 127.0.0.1)"
    )
    

    
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug logging"
    )
    
    return parser.parse_args()

def create_app(data_dir: str) -> App:
    """
    Create and configure the Shiny application.
    
    Args:
        data_dir: Path to the top-level data directory
        
    Returns:
        Configured Shiny App instance
    """
    logger = setup_logging(level="INFO")
    logger.info("Initializing LIANA Results Explorer")
    
    # Initialize data handler
    data_handler = DataHandler()
    
    # Store the data directory for the server to use
    data_handler.top_level_dir = data_dir
    
    # Create UI with CLI mode when data directory is provided via command line
    app_ui = create_full_ui(data_dir=data_dir, cli_mode=True)
    
    # Create server function with automatic discovery
    server_function = create_server_function(data_handler)
    
    # Create and return the app
    app = App(app_ui, server_function, static_assets={"/ai_snapshots": SNAPSHOT_DIR})
    
    logger.info("LIANA Results Explorer initialized successfully")
    return app

def main():
    """Main function to run the application."""
    args = parse_arguments()
    
    # Validate data directory
    data_dir = Path(args.data_dir)
    if not data_dir.exists():
        print(f"❌ Error: Data directory does not exist: {data_dir}")
        sys.exit(1)
    
    if not data_dir.is_dir():
        print(f"❌ Error: Path is not a directory: {data_dir}")
        sys.exit(1)
    
    try:
        # Print startup information
        print(f"🔬 LIANA Results Explorer")
        print(f"📊 Data directory: {data_dir}")
        print(f"🌐 Dashboard will be available at: http://{args.host}:{args.port}")
        print(f"🚀 Starting server...")
        
        # Create and run the app
        app = create_app(str(data_dir))
        app.run(host=args.host, port=args.port)
        
    except KeyboardInterrupt:
        print("\n🛑 Server stopped by user")
    except Exception as e:
        print(f"❌ Error running server: {e}")
        raise

if __name__ == "__main__":
    main()