# cif-site-analyzer

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/emiljaffal/Site-Analysis/blob/main/LICENSE)
![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)

This script parses .cif files based on structure type and creates a CSV table with site compositions, and heat maps for each site and a heat map for compositions. Run the main.py file with (i) a folder containing the same-structure-type CIFs in the same directory as the main.py file (cif-site-analyzer).

> The current README.md serves as a tutorial and documentation - last updated June 24, 2025

## Demo

The code is designed for interactive use without the need to write any code.

![SA-demo-gif](https://github.com/EmilJaffal/Site-Analysis/blob/main/assets/siteanalysis_DEMO.gif)

## Scope

Any `.cif` files.

## How to install `cif-site-analyzer` locally

Download or clone the repository.

`cd` into the project directory:

```bash
cd cif-site-analyzer
```

Create and activate a new conda environment:

```bash
conda create -n cif_site_analyzer_env python=3.13
conda activate cif_site_analyzer_env
```

### Method 1: Install your package with dependencies sourced from pip

It's simple. The only command required is the following:

```bash
pip install .
```

> The above command will automatically install the dependencies listed in `requirements/pip.txt`.

Verify the installation:

```bash
pip list
```

## How to run?

While `cif_site_analyzer_env` is activated in terminal window, type `site_analyzer` and answer the prompts!

## Citations

If you use `CIF Site Analyzer` in your ressearch please cite the following :)
https://link.springer.com/article/10.1007/s40192-025-00400-x#citeas
