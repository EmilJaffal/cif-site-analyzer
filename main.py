import os
import re
import sys
import matplotlib
import matplotlib.patches
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Any
from collections import defaultdict
from mpl_toolkits.axes_grid1 import make_axes_locatable


############################################################################
#
#  Parse Cif and gather required information
#
############################################################################

def cif_to_dict(path: str) -> dict:
    if not os.path.isfile(path):
        print(f"File {path} not found!")
        return
    
    attributes_to_read = [
        '_chemical_formula_structural',
        '_chemical_formula_sum',
        '_chemical_name_structure_type',
        '_chemical_formula_weight'
        '_cell_length_a ',
        '_cell_length_b',
        '_cell_length_c',
        '_cell_angle_alpha',
        '_cell_angle_beta',
        '_cell_angle_gamma',
        '_cell_volume',
        '_cell_formula_units_Z',
        '_space_group_IT_number',
        '_space_group_name_H-M_alt',
        # '#_database_code_PCD'
    ]
    
    data = defaultdict(Any)
    with open(path, 'r') as f:
        lines = f.readlines()
        
        # read conditions, if any
        conds = lines[2].replace('#', "").lstrip().split()
        if len(conds) > 3:
            conds = conds[2]
        else:
            conds = ""
        data['conditions'] = conds
        ln = 0
        while ln < len(lines):
            line = lines[ln].lstrip()
            if not len(line.split()):
                ln += 1
                continue
            
            if line.split()[0] in attributes_to_read:
                next_line = lines[ln+1].lstrip()
                if next_line.startswith('_') or next_line.startswith('loop_'):
                    data[line.split()[0]] = line.split()[1:]
                else:
                    line_data = ''
                    while not next_line.startswith('_'):
                        line_data += next_line.strip().replace(';', '').replace(' ', '')
                        ln += 1
                        next_line = lines[ln+1].lstrip()

                    data[line.split()[0]] = line_data.strip().replace(';', '').replace(' ', '')
                
            if line.startswith('loop_') and lines[ln+1].lstrip().startswith('_atom_site'):
                site_data = []
                keys = []
                ln += 1
                while lines[ln].lstrip().startswith('_atom_site'):
                    keys.append(lines[ln].lstrip().replace('\n', ''))
                    ln += 1
                
                while not lines[ln].lstrip().startswith('_'):
                    if len(lines[ln].strip()):
                        site_data.append(dict(zip(keys, lines[ln].lstrip().split())))
                    ln += 1
                data['atom_site_data'] = site_data
                
            ln += 1
    
    data = format_cif_data(data)

    return dict(data)
    
    
def format_cif_data(cif_data: dict) -> dict:
    numeric_arribs = [
        '_chemical_formula_weight'
        '_cell_length_a ',
        '_cell_length_b',
        '_cell_length_c',
        '_cell_angle_alpha',
        '_cell_angle_beta',
        '_cell_angle_gamma',
        '_cell_volume',
        '_cell_formula_units_Z',
        '_space_group_IT_number',
        # '#_database_code_PCD',
        ]
    
    for k in numeric_arribs:
        if k in cif_data:
            if k in ['_cell_formula_units_Z', '_space_group_IT_number', 
                    #  '#_database_code_PCD'
                     ]:
                try:
                    cif_data[k] = int(cif_data[k][0])
                except:
                    cif_data[k] = -1
            else:
                cif_data[k] = float(cif_data[k][0])
    
    string_arribs = [
        '_chemical_formula_sum',
        '_chemical_name_structure_type',
        '_space_group_name_H-M_alt',
        'conditions'
        ]
    
    for k in string_arribs:
        if k in cif_data:
            if isinstance(cif_data[k], list):
                cif_data[k] = ''.join(cif_data[k]).replace("'", '').replace("~", "")
                
    site_data = []
    if cif_data.get('atom_site_data', None) is not None:
        for site in cif_data['atom_site_data']:
            sdict = {}
            for k, v in site.items():
                if k not in ['_atom_site_label', '_atom_site_type_symbol', '_atom_site_Wyckoff_symbol']:
                    sdict[k] = float(v.split('(')[0])
                elif k == "_atom_site_type_symbol":
                    symbol = v[:2]
                    if str(symbol[-1]).isnumeric():
                        symbol = symbol[:-1]
                    sdict[k] = symbol
                else:
                    sdict[k] = v
            sdict['coordinates'] = [float(c) for c in [site['_atom_site_fract_x'], site['_atom_site_fract_y'], site['_atom_site_fract_z']]]
            site_data.append(sdict)
        cif_data['atom_site_data'] = site_data
    else:
        cif_data['atom_site_data'] = []
    

    cif_data['formula'] = _parse_formula(cif_data['_chemical_formula_sum'])
    return cif_data


def _parse_formula(formula: str, strict: bool = True) -> dict[str, float]:
    """
    copied from pymatgen
    
    Args:
        formula (str): A string formula, e.g. Fe2O3, Li3Fe2(PO4)3.
        strict (bool): Whether to throw an error if formula string is invalid (e.g. empty).
            Defaults to True.

    Returns:
        Composition with that formula.

    Notes:
        In the case of Metallofullerene formula (e.g. Y3N@C80),
        the @ mark will be dropped and passed to parser.
    """
    # Raise error if formula contains special characters or only spaces and/or numbers
    if "'" in formula:
        formula = formula.replace("'", "")

    if strict and re.match(r"[\s\d.*/]*$", formula):
        print(formula)
        raise ValueError(f"Invalid {formula=}")

    # For Metallofullerene like "Y3N@C80"
    formula = formula.replace("@", "")
    # Square brackets are used in formulas to denote coordination complexes (gh-3583)
    formula = formula.replace("[", "(")
    formula = formula.replace("]", ")")
    
    def get_sym_dict(form: str, factor: float) -> dict[str, float]:
        sym_dict: dict[str, float] = defaultdict(float)
        for match in re.finditer(r"([A-Z][a-z]*)\s*([-*\.e\d]*)", form):
            el = match[1]
            amt = 1.0
            if match[2].strip() != "":
                amt = float(match[2])
            sym_dict[el] += amt * factor
            form = form.replace(match.group(), "", 1)
        if form.strip():
            raise ValueError(f"{form} is an invalid formula!")
        return sym_dict

    match = re.search(r"\(([^\(\)]+)\)\s*([\.e\d]*)", formula)
    while match:
        factor = 1.0
        if match[2] != "":
            factor = float(match[2])
        unit_sym_dict = get_sym_dict(match[1], factor)
        expanded_sym = "".join(f"{el}{amt}" for el, amt in unit_sym_dict.items())
        expanded_formula = formula.replace(match.group(), expanded_sym, 1)
        formula = expanded_formula
        match = re.search(r"\(([^\(\)]+)\)\s*([\.e\d]*)", formula)
    return get_sym_dict(formula, 1)


##################################################################################################
# 
# END OF CIF PARSING
##################################################################################################

##################################################################################################
#
# SITE DATA ANALYSIS
##################################################################################################

def find_sites_with_same_wyckoff_symbols(cifpath):
    shared = {}
    sites = {}
    
    cif = cif_to_dict(cifpath)
    for site in cif['atom_site_data']:
        wyckoff = f"{int(site['_atom_site_symmetry_multiplicity'])}{site['_atom_site_Wyckoff_symbol']}"
        if wyckoff not in sites:
            sites[wyckoff] = [site['coordinates']]
        else:
            current_coords = sites[wyckoff]
            duplicate = False
            for coord in current_coords:
                if np.allclose(coord, site['coordinates']):
                    duplicate = True
                    
            if not duplicate:
                current_coords.append(site['coordinates'])
                sites[wyckoff] = current_coords
                
    for k, v in sites.items():
        if len(v) > 1:
            shared[k] = np.array(v)
    
    return shared


def get_site_label(site, shared):
    wyckoff = f"{int(site['_atom_site_symmetry_multiplicity'])}{site['_atom_site_Wyckoff_symbol']}"

    if wyckoff in shared:
        if len(shared):
            order = shared[wyckoff]
            ind = [[i+1, np.linalg.norm(np.array(site['coordinates']) - order[i])] for i in range(order.shape[0])]
            ind = sorted(ind, key=lambda x: x[1])
            wyckoff = f'{wyckoff}{ind[0][0]}'
    #     return wyckoff
    # else:
    return wyckoff



def get_data_row(cifpath, shared_wyckoffs=None):
    
    cif = cif_to_dict(cifpath)
    site_formula = []
    for site in cif['atom_site_data']:
        comp = site['_atom_site_type_symbol']
        wyckoff = get_site_label(site, shared_wyckoffs)
        if site['_atom_site_occupancy'] != 1.0:
            comp += str(site['_atom_site_occupancy'])
        site_formula.append([comp, wyckoff,  site['coordinates']])
        
    # same site
    site_formula_r = {}
    for i in range(len(site_formula)):
        icoordinates = site_formula[i][2]
        iformula = site_formula[i][0]
        iwyckoff = site_formula[i][1]
        
        if iwyckoff not in site_formula_r:
            site_formula_r[iwyckoff] = [iformula, icoordinates]
        else:
            jformula, jcoordinates = site_formula_r[iwyckoff]
            if np.allclose(icoordinates, jcoordinates):
                site_formula_r[iwyckoff] = [jformula+iformula, icoordinates]
                
    data = {
            # 'PCD Entry No.': cif['#_database_code_PCD'],
            'Filename': cifpath.split(os.sep)[-1],
            'Formula': cif['_chemical_formula_sum'], 'Notes': cif['conditions'],
            'Num Elements': len(cif['formula'])}
    for k, v in site_formula_r.items():
        data[k] = v
    return data


###############################################################################################
#
# END OF SITE DATA
###############################################################################################


###############################################################################################
#
# HEATMAP
###############################################################################################


def get_site_element_dist(dataframe: pd.DataFrame, site: str = None) -> pd.DataFrame:
    element_stats = defaultdict(float)
    
    for _, row in dataframe.iterrows():
        if row['Formula'] == "":
            continue
        
        if site is None:
            formula = _parse_formula(str(row['Formula']))
        else:
            formula = _parse_formula(str(row[site]))
        
        for k, v in formula.items():
            element_stats[k] += v

    return element_stats


def format_wyckoff_site_label(site_label, italisize=False):
    ws = [v for v in re.split(r'\d+', site_label) if v !=''][0]
    nums = [s for s in site_label.split(ws) if s != ""]
    if italisize:
        site_label = nums[0] + f'${ws}$'
    else:
        site_label = nums[0] + f'{ws}'
    if len(nums) > 1:
        site_label += f" ({''.join(nums[1:])})"
    return site_label


def ptable_heatmap_mpl(vals_dict: dict, site: str, stype: str, cmap: str):
    
    # format site
    elements_in_data = list(vals_dict.keys())
    if site is not None:
        site = format_wyckoff_site_label(site, italisize=True)
    
    elements = {'H': [1, 1], 'He': [1, 18], 'Li': [2, 1], 'Be': [2, 2], 'B': [2, 13], 'C': [2, 14], 'N': [2, 15], 'O': [2, 16], 'F': [2, 17], 'Ne': [2, 18], 
                'Na': [3, 1], 'Mg': [3, 2], 'Al': [3, 13], 'Si': [3, 14], 'P': [3, 15], 'S': [3, 16], 'Cl': [3, 17], 'Ar': [3, 18], 'K': [4, 1], 'Ca': [4, 2], 
                'Sc': [4, 3], 'Ti': [4, 4], 'V': [4, 5], 'Cr': [4, 6], 'Mn': [4, 7], 'Fe': [4, 8], 'Co': [4, 9], 'Ni': [4, 10], 'Cu': [4, 11], 'Zn': [4, 12], 
                'Ga': [4, 13], 'Ge': [4, 14], 'As': [4, 15], 'Se': [4, 16], 'Br': [4, 17], 'Kr': [4, 18], 'Rb': [5, 1], 'Sr': [5, 2], 'Y': [5, 3], 'Zr': [5, 4], 
                'Nb': [5, 5], 'Mo': [5, 6], 'Tc': [5, 7], 'Ru': [5, 8], 'Rh': [5, 9], 'Pd': [5, 10], 'Ag': [5, 11], 'Cd': [5, 12], 'In': [5, 13], 'Sn': [5, 14], 
                'Sb': [5, 15], 'Te': [5, 16], 'I': [5, 17], 'Xe': [5, 18], 'Cs': [6, 1], 'Ba': [6, 2], 'La': [8, 4], 'Ce': [8, 5], 'Pr': [8, 6], 'Nd': [8, 7], 
                'Pm': [8, 8], 'Sm': [8, 9], 'Eu': [8, 10], 'Gd': [8, 11], 'Tb': [8, 12], 'Dy': [8, 13], 'Ho': [8, 14], 'Er': [8, 15], 'Tm': [8, 16], 'Yb': [8, 17], 
                'Lu': [8, 18], 'Hf': [6, 4], 'Ta': [6, 5], 'W': [6, 6], 'Re': [6, 7], 'Os': [6, 8], 'Ir': [6, 9], 'Pt': [6, 10], 'Au': [6, 11], 'Hg': [6, 12], 
                'Tl': [6, 13], 'Pb': [6, 14], 'Bi': [6, 15], 'Po': [6, 16], 'At': [6, 17], 'Rn': [6, 18], 'Ac': [9, 4], 'Th': [9, 5], 'Pa': [9, 6], 'U': [9, 7], 
                'Np': [9, 8], 'Pu': [9, 9], 'Am': [9, 10], 'Cm': [9, 11], 'Bk': [9, 12], 'Cf': [9, 13], 'Es': [9, 14], 'Fm': [9, 15], 'Md': [9, 16], 'No': [9, 17],
                'Lr': [9, 18]}
    
    pmask = np.zeros(shape=(9, 18))
    for k, v in elements.items():
        p, g = v
        pmask[p-1, g-1] = 1
    
    min_value = min(vals_dict.values())
    norm_const = max(vals_dict.values()) - min_value
    heat_map = np.full(shape=(9, 18), fill_value=-(min_value))
    for k, v in elements.items():
        if k in vals_dict:
            heat_map[v[0]-1, v[1]-1] = vals_dict[k]
    heat_map *= pmask
    heat_map[heat_map==0.] = np.nan

    fig, ax = plt.subplots(figsize=(18, 9))
    
    cmap = plt.get_cmap(cmap)
    cmap.set_under('w', alpha=0.1)
    
    im = plt.imshow(heat_map, cmap=cmap, vmin=min_value)
    # Minor ticks
    ax.set_xticks(np.arange(-.5, 18, 1), minor=True)
    ax.set_yticks(np.arange(-.5, 9, 1), minor=True)

    # Gridlines based on minor ticks
    # ax.grid(which='minor', color='w', linestyle='-', linewidth=3)

    ticks = sorted(set([int(i) for i in np.linspace(min_value, max(vals_dict.values()), 4)]))
    im.set_clim(vmin=ticks[0], vmax=ticks[-1])
    divider = make_axes_locatable(ax)
    
    cbar = plt.colorbar(im, ticks=ticks, location='top', 
                        shrink=0.25, cax=ax.inset_axes((0.19, 0.72, 0.4, 0.02)), orientation='horizontal')
    cbar.ax.tick_params(labelsize=25) 
    
    # element symbols
    kw = dict(horizontalalignment="center",
              verticalalignment="center",
              size=25)
    
    im.axes.text(2, 5.1, '*', alpha=0.6, **kw)
    im.axes.text(2, 7.1, '*', alpha=0.6, **kw)
    for k, v in elements.items():
        c = 'w' if ((vals_dict[k] - min_value) / norm_const) > 0.6 else 'k'
        if k in elements_in_data:
            im.axes.text(v[1]-1, v[0]-1, k, c=c, **kw)
        else:
            im.axes.text(v[1]-1, v[0]-1, k, c=c, alpha=0.6, **kw)
        #  path_effects=[patheffects.withStroke(linewidth=3, foreground='white')], 
                     
        rect = matplotlib.patches.Rectangle(xy=(v[1]-1.5, v[0]-1.5), width=1, height=1, 
                                            facecolor=(0, 0, 0, 0), edgecolor=(0,0,0, 1), lw=1)
        im.axes.add_patch(rect)
    
    # placeholder for title
    rectangle = matplotlib.patches.Rectangle(xy=(2., -0.59), width=9., height=1.7, facecolor='w', edgecolor=None)
    ax.add_patch(rectangle)
    rx, ry = rectangle.get_xy()
    cx = rx + rectangle.get_width() / 2.
    cy = ry + rectangle.get_height() /2.
    
    # compose title and wrap
    if site is not None:
        text_str = f"element distribution of {site} site in {stype.split(',')[0]}-type compounds"
    else:
        text_str = f"element distribution in {stype.split(',')[0]}-type compounds"
    words = text_str.split()
    wrapped_text = ""
    line = ""
    max_chars_per_line = 40  # Adjust based on font and rectangle width

    for word in words:
        if len(line + word) + 1 <= max_chars_per_line:
            line += word + " "
        else:
            wrapped_text += line.rstrip() + "\n"
            line = word + " "
    wrapped_text += line.rstrip()

    ax.annotate(wrapped_text,
                xy=(cx, cy), color='black', fontsize=25, ha='center', va='center', wrap=False)
    
    # remove ticks
    ax.axes.get_xaxis().set_ticks([])
    ax.axes.get_yaxis().set_ticks([])
    # Remove minor ticks
    ax.tick_params(which='minor', bottom=False, left=False)
    
    # Remove all spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    stype = '-'.join(stype.split(',')[:2])
    if site is not None:
        plt.savefig(f"ElemDist_{stype}_{site.replace('$', '').replace(' ', '')}.png", dpi=300)
    else:
        plt.savefig(f"ElemDist_{stype}.png", dpi=300)
        

def get_ws(s):
    s = re.split(r'\d', s)
    return [_s for _s in s if _s != ''][0]


if __name__ == "__main__":
    
    input_cifs = []
    if len(sys.argv) != 3:
        print("This script requires \n\t(i) file contatining names of cif files to process, and \n\t(ii) path to the folder contating .cif files.")
        exit(0)
    
    cif_filenames = []
    root = sys.argv[2]
    with open(sys.argv[1], 'r') as f:
        for line in f.readlines():
            line = line.replace('\n', "")
            if not line.endswith('.cif'):
                line += '.cif'
            cif_filenames.append(line)

    # read cif and get formula and structure type
    # TODO: read cif only once
    cif_file_data = []
    for cif_name in cif_filenames:
        cif_path = os.path.join(root, cif_name)
        try:
            cif = cif_to_dict(cif_path)
            formula = cif['_chemical_formula_sum']
            stype = cif['_chemical_name_structure_type']
            cif_file_data.append({'Filename': cif_name, 'Formula': formula, 'Structure type': stype})
        except Exception as e:
            print(f"Error reading {cif_path} file.")
            # print(e)
    
    cif_file_data = pd.DataFrame(cif_file_data)
    stypes = cif_file_data['Structure type'].value_counts()
    
    # Get user input for structure type
    if len(stypes) > 1: 
        print("More than one structure types are fould in the input list.")
        print("Please select one structure type form the list below.")

        print(f"\nNo  {'Count':<5} Structure type")
        for i, (index, r) in enumerate(stypes.items(), 1):
            print(f"({i}) {r:<5} {index}")
            
        valid_input = False
        while not valid_input:
            res = input("Please enter the number corresponding to the selected strcuture type: \n")
            try:
                res = int(res)
                assert res <= len(stypes)
                selected_stype = stypes.index[res-1]
                valid_input = True
            except Exception as e:
                print(e)
                print("Invalid input. Try again.")
    else:
        selected_stype = str(cif_file_data.iloc[0]['Structure type'])

    selected_entries = cif_file_data[cif_file_data['Structure type'] == selected_stype]['Filename'].to_list()
    shared = find_sites_with_same_wyckoff_symbols(os.path.join(root, selected_entries[0]))
    data0 = get_data_row(os.path.join(root, selected_entries[0]), shared)
    all_sites = [k for k in data0.keys() if k not in ['Filename', 'Formula', 'Notes', 'Num Elements']]
    all_sites = sorted(all_sites, key=lambda k: int(re.split('[a-z]', k)[0]))
    all_sites = sorted(all_sites, key=lambda k: get_ws(k))
    
    # Get user input for sites
    if len(all_sites) > 5:
        print(f"There are {len(all_sites)} sites present for this structure type.")
        print("Please select a maximum of five sites from the the list below.")
        print("Enter the numbers separated by space. e.g. 1 3 4 6")
        
        print("\nNo Site")
        for i, site in enumerate(all_sites, 1):
            print(f"({i}) {site}")

        valid_input = False
        while not valid_input:
            res = input("Please enter the numbers corresponding to the selected sites: \n")
            try:
                res = res.replace(',', ' ')
                res = [int(r) for r in res.split()]
                assert np.all(np.array(res) <= len(all_sites))
                sites = [all_sites[i-1] for i in res]
                valid_input = True
            except Exception as e:
                print(e)
                print("Invalid input. Try again.")
    else:
        sites = all_sites

    # collect data from cifs
    data = []
    pos_data = dict(zip(all_sites, [[] for _ in range(len(data0)-2)]))

    for i in selected_entries:
        row_data_w_coords = get_data_row(os.path.join(root, i), shared)
        for k in all_sites:
            pos_data[k].append(row_data_w_coords[k][1])
            row_data_w_coords[k] = row_data_w_coords[k][0]
        data.append(row_data_w_coords)
    
    pos_row = dict(zip(data0.keys(), ["" for _ in range(len(data0))]))
    for k in all_sites:
        all_positions = np.array(pos_data[k])
        avg_positions = ""
        for i in range(3):
            avg_positions += f"{all_positions[:, i].mean():.3f}"
            if all_positions[:, i].std() > 0.001:
                avg_positions += f"({int(round(all_positions[:, i].std(), 3)*1000)}) "
            else:
                avg_positions += " "
        pos_row[k] = avg_positions
    data.insert(0, pos_row)      
    data = pd.DataFrame(data)
    
    # Format site names for writing csv
    columns = [c for c in data.columns if c not in all_sites]
    formatted_site_names = {}
    for cname in all_sites:
        new_name = format_wyckoff_site_label(cname)
        formatted_site_names[cname] = new_name
        columns.append(new_name)
        
    data.rename(columns=formatted_site_names, inplace=True)
    data[columns].to_csv(f"{'-'.join(selected_stype.split(',')[:2])}.csv", index=False)
    data.rename(columns=dict(zip(formatted_site_names.values(), formatted_site_names.keys())), inplace=True)
    
    # Plot site heatmaps
    cmaps = ["Reds", "Blues", "Greens", "Purples", "Oranges"]
    for i, site in enumerate(sites):
        print(f"Processing {site}")
        ptable_heatmap_mpl(vals_dict=get_site_element_dist(data, site), site=site, stype=selected_stype, cmap=cmaps[i])
    
    # Cumulative heatmap
    ptable_heatmap_mpl(vals_dict=get_site_element_dist(data.iloc[1:], site=None), site=None, stype=selected_stype, cmap="Greys")