# 2D-3D-PSD-tool

This repository contains a tool that transforms a 2D pore radius distribution (Cross-section radius distribution) to a 3D pore radius distribution (Pore Size Distribution). The method is detailed in the article "A New Method to Measure Pore Radius Distribution of Powders".

## Installation

It is recommended to use a Conda environment to manage dependencies.

1. Create a new Conda environment:
   ```bash
   conda create -n psd_tool python=3.10
   ```

2. Activate the environment:
   ```bash
   conda activate psd_tool
   ```

3. Clone the repository and install the package using pip:
   ```bash
   pip install .
   ```

For development (editable install):
```bash
pip install -e .
```

## Usage

The core functionality is provided by the `psd_tool` package.

### Notebooks

The `notebooks/` directory contains a Jupyter notebook demonstrating the usage:

- **`tutorial.ipynb`**: A step-by-step tutorial on how to use the `psd_tool` to estimate PSD from synthetic data.

### Python API

You can also use the package directly in your Python scripts:

```python
from psd_tool import estimate_c
import numpy as np

# Load your chord length distribution data
# chords = ...

# Estimate Pore Size Distribution
c, bins, residual = estimate_c(chords, r_min=0.05, cutoff=1, num_bins=20, linear=False)

print("Estimated PSD coefficients:", c)
```

#### `estimate_c` Function Reference

```python
def estimate_c(csd, r_min=0.5, cutoff=0.3, num_bins=20, linear=True)
```

Estimate the 3D pore radius distribution from chord size distribution.

**Parameters:**

- `csd` (array_like): Chord size distribution data.
- `r_min` (float, optional): Minimum radius to consider. Default is 0.5.
- `cutoff` (float, optional): Cutoff parameter for the phi function. Default is 0.3.
- `num_bins` (int, optional): Number of bins for the histogram. Default is 20.
- `linear` (bool, optional): If True, use linear binning. If False, use logarithmic binning. Default is True.

**Returns:**

- `c` (array): Estimated coefficients (probability density).
- `bins` (array): Bin edges.
- `residual` (float): Residual of the NNLS solution.

## Directory Structure

- `src/psd_tool/`: Source code for the Python package.
  - `core.py`: Contains the main estimation logic (`estimate_c`).
  - `saltykov.py`: Implementation of the Saltykov method for comparison.
- `notebooks/`: Jupyter notebook for tutorial.
- `data/`: Data files used in the notebooks.

## Dependencies

- numpy
- scipy
- matplotlib
- seaborn
- pandas
