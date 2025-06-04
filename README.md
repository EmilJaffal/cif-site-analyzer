# CIF Site Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/emiljaffal/Site-Analysis/blob/main/LICENSE)
![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)

This script parses .cif files based on structure type and creates a CSV table with site compositions, and heat maps for each site (up to 5 sites) and a heat map for compositions. Run the main.py file with (i) a folder containing same-structure-type CIFs in the same directory as the main.py file (cif-site-analyzer).

> The current README.md serves as a tutorial and documentation - last update December 21, 2024

## Demo

The code is designed for interactive use without the need to write any code.

![SA-demo-gif](https://github.com/EmilJaffal/Site-Analysis/blob/main/assets/siteanalysis_DEMO.gif)

## Scope

Any `.cif` files.

## Getting started

Copy each line into your command-line applications. The code needs (i) a folder containing same-structure-type CIFs in the same directory as the main.py file (cif-site-analyzer).

```bash
$ git clone https://github.com/EmilJaffal/cif-site-analyzer
$ cd cif-site-analyzer
$ pip install -r requirements.txt
$ python main.py
```

Once the code is executed using `python main.py`, the following prompt will
appear, asking you to choose one of the available structure types:

```text

Folders with .cif files:
1. NaCl,cF8,225 (20 files)

Would you like to process each folder above sequentially? (Y/n): 
```

You may then choose to process whichever structure type folder you would like, and it will process the sites.

When a structure type has more than 5 sites, a prompt will be given to select the sites for further analysis. For any option, SA will ask you to choose a structure type from the folder containing `.cif` files:

```
There are 9 sites present for this structure type.
Please select a maximum of five sites from the list below.
Enter the numbers separated by a space. e.g. 1 3 4 6
Enter 0 to plot all sites or enter the numbers corresponding to the selected sites

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

You can find an example of our test set here: [demo_cifs](https://github.com/EmilJaffal/cif-site-analyzer/tree/main/NaCl%2CcF8%2C225)

#### Output 2 - Heatmaps

Periodic table heat maps showing element distribution of the compounds in the selected structure type:

![SA-heatmap](https://github.com/EmilJaffal/Site-Analysis/blob/main/ElemDist_CuTi-tP4.png)

and the element distribution of sites:

![SA-site-heatmap](https://github.com/EmilJaffal/cif-site-analyzer/blob/main/ElemDist_NaCl-cF8.png)

When done, the following prompt will display:

```text
Saved 3 PNG files:
  ElemDist_NaCl-cF8_4a.png
  ElemDist_NaCl-cF8_4b.png
  ElemDist_NaCl-cF8.png
Saved CSV file: NaCl-cF8_NaCl,cF8,225.csv
Thanks for using cif-site-analyzer! teehee
```

## Contributors

- [Balaranjan Selvaratnam](https://github.com/balaranjan)
- [Emil Jaffal](https://github.com/EmilJaffal)
- Anton Oliynyk

## How to ask for help

`cif-site-analyzer` is also designed for experimental materials scientists and chemists.

- If you have any issues or questions, please feel free to reach out or
  [leave an issue](https://github.com/emiljaffal/Site-Analysis/issues).
