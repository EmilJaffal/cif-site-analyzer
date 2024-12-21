# Site-Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/emiljaffal/Site-Analysis/blob/main/LICENSE)
![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)

This script parses .cif files based on structure type and creates a CSV table with site compositions, and heat maps for each site (up to 5 sites) and a heat map for compositions. Run the create_occu_table_and_heatmaps.py file with the path of the folder containing CIFs belonging to the same structure type as argument.

> The current README.md serves as a tutorial and documentation - last update December 21, 2024

## Demo

The code is designed for interactive use without the need to write any code.

![SA-demo-gif](https://github.com/EmilJaffal/Site-Analysis/blob/main/assets/siteanalysis_DEMO.gif)

## Scope

Any `.cif` files.

## Getting started

Copy each line into your command-line applications. The code needs (i) a file containing list of filename/PCD numbers and (ii) path to folder containing all cifs. This will make it easy when you have CIFs with several structure types in one folder - no need to separate them.

```bash
$ git clone https://github.com/EmilJaffal/Site-Analysis
$ cd Site-Analysis
$ pip install -r requirements.txt
$ python main.py [name].csv [cif_folder_name]
```

Once the code is executed using `python main.py`, the following prompt will
appear, asking you to choose one of the available structure types:

```text

File cif_files/26147078.cif not found!
Error reading cif_files/26147078.cif file.
More than one structure types are fould in the input list.
Please select one structure type form the list below.

No  Count Structure type
(1) 19    CuTi,tP4,129
(2) 5     Ca14AlSb11,tI208,142
Please enter the number corresponding to the selected structure type: 

```

You may then choose to process whichever structure type you would like, and it will process the sites.

When a structure type has more than 5 sites, a prompt will be given to select the sites for further analysis. For any option, SA will ask you to choose a structure type from the folder containing `.cif` files:

#### Output 1 CSV

Data for each folder is saved in `[structuretype].csv`. Below is an example of a .csv of your structure type containing: filename, formula, notes (synthesis conditions), # of elements and site occupations.

```csv
Filename	Formula	Notes	Num Elements	2c1	2c2
        0.250 0.250 0.150(1) 	0.250 0.250 0.600(1) 
527049.cif	Ag1.5Hf0.5		2	Ag0.5Hf0.5	Ag
1815856.cif	Ag1.02Hf0.98		2	Hf0.98Ag0.02	Ag
1826520.cif	AgHf		2	Hf	Ag
528891.cif	AgHf		2	Hf	Ag
1826507.cif	AgZr		2	Zr	Ag
260733.cif	AgZr		2	Zr	Ag
527770.cif	AgZr		2	Zr	Ag
528846.cif	AuHf	rt	2	Hf	Au
261356.cif	AuTi	rt	2	Ti	Au
528849.cif	AuTi	rt	2	Ti	Au
1010265.cif	CuTi	tet1	2	Ti	Cu
```

#### Output 2 Heatmaps

Periodic table heat maps showing element distribution of the overall structure type in your folder

![SA-heatmap](https://github.com/EmilJaffal/Site-Analysis/blob/main/ElemDist_CuTi-tP4.png)

and element distribution of sites 

![SA-site-heatmap](https://github.com/EmilJaffal/Site-Analysis/blob/main/ElemDist_CuTi-tP4_2c(1).png)

## Installation

```bash
$ git clone https://github.com/EmilJaffal/Site-Analysis
$ cd Site-Analysis
$ pip install -r requirements.txt
$ python main.py [name].csv [cif_folder_name]
```

## Contributors

- Anton Oliynyk
- Balaranjan Selvaratnam
- Emil Jaffal

## How to ask for help

`CBA` is also designed for experimental materials scientists and chemists.

- If you have any issues or questions, please feel free to reach out or
  [leave an issue](https://github.com/emiljaffal/Site-Analysis/issues).