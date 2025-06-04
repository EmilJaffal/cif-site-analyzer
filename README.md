# CIF Site Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/emiljaffal/Site-Analysis/blob/main/LICENSE)
![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)

This script parses .cif files based on structure type and creates a CSV table with site compositions, and heat maps for each site (up to 5 sites) and a heat map for compositions. Run the main.py file with (i) a text file containing the cif file names to analyze (belonging to the same structure type) and (ii) the path of the folder containing CIFs as arguments.

> The current README.md serves as a tutorial and documentation - last update December 21, 2024

## Demo

The code is designed for interactive use without the need to write any code.

![SA-demo-gif](https://github.com/EmilJaffal/Site-Analysis/blob/main/assets/siteanalysis_DEMO.gif)

## Scope

Any `.cif` files.

## Getting started

Copy each line into your command-line applications. The code needs (i) a file containing a list of filenames/PCD numbers and (ii) a path to the folder containing all cifs. This will make it easy when you have CIFs with several structure types in one folder - no need to separate them.

```bash
$ git clone https://github.com/EmilJaffal/Site-Analysis
$ cd Site-Analysis
$ pip install -r requirements.txt
$ python main.py [cif_names].csv [cif_folder_name]
```

Once the code is executed using `python main.py`, the following prompt will
appear, asking you to choose one of the available structure types:

```text

More than one structure type is found in the input list.
Please select one structure type from the list below.

No  Count Structure type
(1) 19    CuTi,tP4,129
(2) 5     Ca14AlSb11,tI208,142
Please enter the number corresponding to the selected structure type: 
```

You may then choose to process whichever structure type you would like, and it will process the sites.

When a structure type has more than 5 sites, a prompt will be given to select the sites for further analysis. For any option, SA will ask you to choose a structure type from the folder containing `.cif` files:

```
There are 9 sites present for this structure type.
Please select a maximum of five sites from the list below.
Enter the numbers separated by a space. e.g. 1 3 4 6

No Site
(1) 8a
(2) 8b
(3) 16e
(4) 16f
(5) 32g1
(6) 32g2
(7) 32g3
(8) 32g4
(9) 32g5
Please enter the numbers corresponding to the selected sites: 
```

#### Output 1 - .csv

Data for each folder is saved in `[structuretype].csv`. Below is an example of a .csv of your structure type containing: filename, formula, notes (synthesis conditions), # of elements and site occupations.
The first row shows the sites' average coordinates and the associated standard deviation for all the files processed with the selected structure type.  

You can find an example of our test set here: [SA-csv](https://github.com/EmilJaffal/Site-Analysis/blob/main/CuTi-tP4.csv)

#### Output 2 - Heatmaps

Periodic table heat maps showing element distribution of the compounds in the selected structure type:

![SA-heatmap](https://github.com/EmilJaffal/Site-Analysis/blob/main/ElemDist_CuTi-tP4.png)

and the element distribution of sites:

![SA-site-heatmap](https://github.com/EmilJaffal/Site-Analysis/blob/main/ElemDist_CuTi-tP4_2c(1).png)

## Installation

```bash
$ git clone https://github.com/EmilJaffal/Site-Analysis
$ cd Site-Analysis
$ pip install -r requirements.txt
$ python main.py [name].csv [cif_folder_name]
```

## Contributors

- [Balaranjan Selvaratnam](https://github.com/balaranjan)
- [Emil Jaffal](https://github.com/EmilJaffal)
- Anton Oliynyk

## How to ask for help

`Site-Analysis` is also designed for experimental materials scientists and chemists.

- If you have any issues or questions, please feel free to reach out or
  [leave an issue](https://github.com/emiljaffal/Site-Analysis/issues).
