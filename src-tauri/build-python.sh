#!/bin/bash
set -e

# Build Python bundle for macOS and Linux
# This script creates a self-contained Python environment with baseline correction dependencies

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Determine platform
if [[ "$OSTYPE" == "darwin"* ]]; then
    PLATFORM="macos"
    PYTHON_BIN="bin/python"
elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
    PLATFORM="linux"
    PYTHON_BIN="bin/python"
else
    echo "❌ Error: Unsupported platform: $OSTYPE"
    echo "This script only supports macOS and Linux"
    exit 1
fi

BUNDLE_DIR="$SCRIPT_DIR/resources/python/$PLATFORM"

echo "Building Python bundle for $PLATFORM..."

# Check if uv is installed
if ! command -v uv &> /dev/null; then
    echo "❌ Error: uv is not installed"
    echo ""
    echo "Please install uv first:"
    echo "  curl -LsSf https://astral.sh/uv/install.sh | sh"
    echo ""
    if [[ "$PLATFORM" == "macos" ]]; then
        echo "Or via Homebrew:"
        echo "  brew install uv"
        echo ""
    fi
    exit 1
fi

# Clean up any existing bundle
if [ -d "$BUNDLE_DIR" ]; then
    echo "Removing existing bundle..."
    rm -rf "$BUNDLE_DIR"
fi

# Create virtual environment with Python 3.13
# Using --relocatable to make the venv work when moved/bundled
echo "Creating virtual environment with Python 3.13..."
uv venv --python 3.13 --relocatable "$BUNDLE_DIR"

# Install required packages from requirements.txt
echo "Installing packages from requirements.txt..."
uv pip install --python "$BUNDLE_DIR" -r "$SCRIPT_DIR/python/requirements.txt"

# Copy the baseline correction script
echo "Copying baseline correction script..."
cp "$SCRIPT_DIR/python/baseline_correction.py" "$BUNDLE_DIR/baseline_correction.py"

# Verify the installation
echo "Verifying installation..."
"$BUNDLE_DIR/$PYTHON_BIN" -c "import numpy, scipy; print('✅ numpy version:', numpy.__version__); print('✅ scipy version:', scipy.__version__)"

# Test the baseline correction script
echo "Testing baseline correction..."
echo '{"spectrum": [1,2,3,4,5,4,3,2,1], "lambda_param": 1e7, "p": 0.01, "d": 2}' \
    | "$BUNDLE_DIR/$PYTHON_BIN" "$BUNDLE_DIR/baseline_correction.py" > /dev/null

if [ $? -eq 0 ]; then
    echo "✅ Baseline correction test passed"
else
    echo "❌ Baseline correction test failed"
    exit 1
fi

echo ""
echo "✅ Python bundle built successfully!"
echo "   Location: $BUNDLE_DIR"
echo "   Python: $("$BUNDLE_DIR/$PYTHON_BIN" --version)"
echo ""
echo "You can now run: bun run tauri dev"
