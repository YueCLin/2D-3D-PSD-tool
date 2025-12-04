# 2D-3D-PSD-tool

This repository contains a tool that transforms a 2D pore radius distribution (Chord Length Distribution) to a 3D pore radius distribution (Pore Size Distribution). The method is detailed in the article "A New Method to Measure Pore Radius Distribution of Powders".

## Installation

To install the package, clone the repository and install it using pip:

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

The `notebooks/` directory contains Jupyter notebooks demonstrating the usage and reproduction of results:

- **`tutorial.ipynb`**: A step-by-step tutorial on how to use the `psd_tool` to estimate PSD from synthetic data.
- **`reproduction.ipynb`**: Reproduces the figures from the research paper using the packaged code.
- **`comparison.ipynb`**: Compares the new method with the Saltykov method.

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

## Directory Structure

- `src/psd_tool/`: Source code for the Python package.
  - `core.py`: Contains the main estimation logic (`estimate_c`).
  - `saltykov.py`: Implementation of the Saltykov method for comparison.
- `notebooks/`: Jupyter notebooks for tutorial, reproduction, and comparison.
- `data/`: Data files used in the notebooks.
- `input/`: Original input files (archived).

## Dependencies

- numpy
- scipy
- matplotlib
- seaborn
- pandas
