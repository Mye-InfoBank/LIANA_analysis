#!/usr/bin/env python3
"""
Validation script for LIANA Results Explorer setup
This script checks if all dependencies are installed and the app can start.
"""

import sys
import os
import importlib
from pathlib import Path

def check_python_version():
    """Check if Python version is compatible."""
    version = sys.version_info
    if version.major < 3 or (version.major == 3 and version.minor < 8):
        return False, f"Python {version.major}.{version.minor}.{version.micro}"
    return True, f"Python {version.major}.{version.minor}.{version.micro}"

def check_dependencies():
    """Check if all required dependencies are installed."""
    required_packages = [
        'shiny',
        'pandas', 
        'numpy',
        'plotly',
        'networkx'
    ]
    
    results = {}
    for package in required_packages:
        try:
            module = importlib.import_module(package)
            version = getattr(module, '__version__', 'unknown')
            results[package] = {'installed': True, 'version': version}
        except ImportError:
            results[package] = {'installed': False, 'version': None}
    
    return results

def check_file_structure():
    """Check if all required files exist."""
    required_files = [
        'app.py',
        'requirements.txt',
        'README.md',
        'run_app.sh'
    ]
    
    current_dir = Path(__file__).parent
    results = {}
    
    for file in required_files:
        file_path = current_dir / file
        results[file] = file_path.exists()
    
    return results

def test_app_import():
    """Test if the app can be imported without errors."""
    try:
        # Add current directory to path
        current_dir = Path(__file__).parent
        sys.path.insert(0, str(current_dir))
        
        # Try to import the app
        import app
        return True, "App imported successfully"
    except Exception as e:
        return False, str(e)

def main():
    """Run all validation checks."""
    print("🔬 LIANA Results Explorer - Setup Validation")
    print("=" * 50)
    
    # Check Python version
    print("\n🐍 Python Version Check:")
    python_ok, python_version = check_python_version()
    status = "✅" if python_ok else "❌"
    print(f"   {status} {python_version}")
    if not python_ok:
        print("   ⚠️  Python 3.8 or higher required")
    
    # Check dependencies
    print("\n📦 Dependency Check:")
    deps = check_dependencies()
    all_deps_ok = True
    
    for package, info in deps.items():
        if info['installed']:
            print(f"   ✅ {package} v{info['version']}")
        else:
            print(f"   ❌ {package} - NOT INSTALLED")
            all_deps_ok = False
    
    if not all_deps_ok:
        print("   💡 Install missing packages: pip install -r requirements.txt")
    
    # Check file structure
    print("\n📁 File Structure Check:")
    files = check_file_structure()
    all_files_ok = True
    
    for file, exists in files.items():
        status = "✅" if exists else "❌"
        print(f"   {status} {file}")
        if not exists:
            all_files_ok = False
    
    # Test app import
    print("\n🚀 App Import Test:")
    app_ok, app_message = test_app_import()
    status = "✅" if app_ok else "❌"
    print(f"   {status} {app_message}")
    
    # Overall status
    print("\n" + "=" * 50)
    overall_ok = python_ok and all_deps_ok and all_files_ok and app_ok
    
    if overall_ok:
        print("🎉 All checks passed! The app is ready to run.")
        print("\n🚀 To start the app:")
        print("   Option 1: ./run_app.sh")
        print("   Option 2: python app.py")
        print("   Option 3: python -c 'from shiny import run_app; run_app(\"app:app\")'")
        print("\n📊 To create test data:")
        print("   python example_config.py")
    else:
        print("❌ Some checks failed. Please fix the issues above.")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())
