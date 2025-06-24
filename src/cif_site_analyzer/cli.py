from .heatmap import ptable_heatmap_mpl
from .utilities import find_sites_with_same_wyckoff_symbols
from .utilities import format_wyckoff_site_label
from .utilities import get_site_element_dist
from .utilities import get_data_row
from .utilities import get_ws
from cif_site_analyzer.cif_reader import read_cif
import pandas as pd
import numpy as np
import os
import re


def main():

    print("Welcome to CIF Site Analyzer")

    current_dir = os.getcwd()

    # Find all folders in the script directory
    # that contain at least one .cif file
    candidate_dirs = []
    for item in os.listdir(current_dir):
        full_path = os.path.join(current_dir, item)
        if os.path.isdir(full_path):
            cif_files = [
                f for f in os.listdir(full_path) if f.lower().endswith(".cif")
            ]
            if cif_files:
                candidate_dirs.append((item, cif_files))

    if not candidate_dirs:
        print(
            f"No folders containing .cif files found \
            in {current_dir}. Please add a folder with .cif files."
        )
        exit(0)

    print("\nFolders with .cif files:")
    for i, (dir_name, cif_files) in enumerate(candidate_dirs, start=1):
        print(f"{i}. {dir_name} ({len(cif_files)} files)")

    # Prompt user for sequential processing or selection
    while True:
        resp = (
            input(
                "\nWould you like to process each folder above "
                "sequentially? (Y/n): "
            )
            .strip()
            .lower()
        )
        if resp in ["", "y", "yes"]:
            selected_indices = list(range(1, len(candidate_dirs) + 1))
            break
        elif resp in ["n", "no"]:
            while True:
                folder_numbers_str = input(
                    "Enter the numbers corresponding to the "
                    "folders listed above, separated by spaces (e.g. 1 2 3): "
                )
                try:
                    selected_indices = list(
                        set(
                            int(number)
                            for number in folder_numbers_str.split()
                        )
                    )
                    if not all(
                        1 <= idx <= len(candidate_dirs)
                        for idx in selected_indices
                    ):
                        raise ValueError(
                            "One or more numbers are out of the valid range."
                        )
                    break
                except ValueError:
                    print(
                        "Please enter only valid numbers within the range, "
                        "separated by spaces."
                    )
            break
        else:
            print("Invalid input. Please enter Y or N.")

    print("\nSelected folders:")
    for idx in selected_indices:
        print(f"{idx}. {candidate_dirs[idx-1][0]}")

    # Process each selected folder
    for idx in selected_indices:
        dir_name, cif_filenames = candidate_dirs[idx - 1]
        cif_dir = os.path.join(current_dir, dir_name)
        print(f"\nProcessing folder: {dir_name}")

        # --- The rest of your original processing code goes here,
        # replacing 'cif_dir' and 'cif_filenames' as needed ---
        cif_file_data = []
        for cif_name in cif_filenames:
            cif_path = os.path.join(cif_dir, cif_name)
            try:
                cif = read_cif(cif_path)
                formula = cif.formula
                stype = cif.structure_type
                cif_file_data.append(
                    {
                        "Filename": cif_name,
                        "Formula": formula,
                        "Structure type": stype,
                    }
                )
            except Exception as e:
                print(f"Error reading {cif_path} file.")
                print(e)

        cif_file_data = pd.DataFrame(cif_file_data)
        stypes = cif_file_data["Structure type"].value_counts()
        stypes = stypes[stypes > 1]

        if not len(stypes):
            print(
                "No structure type has more than one cif in the provided list."
            )
            continue

        # Get user input for structure type
        if len(stypes) > 1:
            print("More than one structure types are found in the input list.")
            print("Please select one structure type from the list below.")

            print(f"\nNo  {'Count':<5} Structure type")
            for i, (index, r) in enumerate(stypes.items(), 1):
                print(f"({i}) {r:<5} {index}")

            valid_input = False
            while not valid_input:
                res = input(
                    "Please enter the number corresponding to the \
                     selected structure type: \n"
                )
                try:
                    res = int(res)
                    assert res <= len(stypes)
                    selected_stype = stypes.index[res - 1]
                    valid_input = True
                except Exception as e:
                    print(e)
                    print("Invalid input. Try again.")
        else:
            selected_stype = str(cif_file_data.iloc[0]["Structure type"])

        selected_entries = cif_file_data[
            cif_file_data["Structure type"] == selected_stype
        ]["Filename"].to_list()
        shared = find_sites_with_same_wyckoff_symbols(
            os.path.join(cif_dir, selected_entries[0])
        )
        data0 = get_data_row(
            os.path.join(cif_dir, selected_entries[0]), shared
        )
        all_sites = [
            k
            for k in data0.keys()
            if k not in ["Filename", "Formula", "Notes", "Num Elements"]
        ]
        all_sites = sorted(
            all_sites, key=lambda k: int(re.split("[a-z]", k)[0])
        )
        all_sites = sorted(all_sites, key=lambda k: get_ws(k))

        # Get user input for sites
        if len(all_sites) > 5:
            print(
                f"\nThere are {len(all_sites)} sites present for this \
                    structure type."
            )
            print("Please select the sites from the the list below.")
            print("Enter the numbers separated by space. e.g. 1 3 4 6")

            print("\nNo Site")
            for i, site in enumerate(all_sites, 1):
                print(f"({i}) {site}")

            valid_input = False
            while not valid_input:
                res = input(
                    "\nEnter 0 to plot all sites or enter the numbers \
                        corresponding to the selected sites: \n"
                )
                try:
                    res = res.replace(",", " ")
                    res = [int(r) for r in res.split()]
                    if 0 in res:
                        valid_input = True
                        print("Plotting all sites.")
                        sites = all_sites
                    else:
                        assert np.all(np.array(res) <= len(all_sites))
                        sites = [all_sites[i - 1] for i in res]
                        valid_input = True
                except Exception as e:
                    print(e)
                    print("Invalid input. Try again.")
        else:
            sites = all_sites

        # collect data from cifs
        data = []
        pos_data = dict(zip(all_sites, [[] for _ in range(len(data0) - 2)]))

        for i in selected_entries:
            row_data_w_coords = get_data_row(os.path.join(cif_dir, i), shared)
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
                    avg_positions += (
                        f"({int(round(all_positions[:, i].std(), 3)*1000)}) "
                    )
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
        data[columns].to_csv(
            f"{'-'.join(selected_stype.split(',')[:2])}_{dir_name}.csv",
            index=False,
        )
        data.rename(
            columns=dict(
                zip(formatted_site_names.values(), formatted_site_names.keys())
            ),
            inplace=True,
        )

        # Plot site heatmaps
        cmaps = [
            "Reds",
            "Blues",
            "Greens",
            "Purples",
            "Oranges",
            "YlOrBr",
            "OrRd",
            "PuRd",
            "RdPu",
            "BuPu",
            "GnBu",
            "PuBu",
            "YlGnBu",
            "PuBuGn",
            "BuGn",
            "YlGn",
            "viridis",
            "plasma",
            "inferno",
            "magma",
            "cividis",
        ]

        png_files = []
        for i, site in enumerate(sites):
            print(f"Processing {site}")
            ptable_heatmap_mpl(
                vals_dict=get_site_element_dist(data, site),
                site=site,
                stype=selected_stype,
                cmap=cmaps[i % len(cmaps)],
            )

            filename = "ElemDist_"
            filename += f"{'-'.join(selected_stype.split(',')[:2])}"
            filename += f"_{site.replace('$', '').replace(' ', '')}.png"
            png_files.append(filename)

        # Cumulative heatmap
        ptable_heatmap_mpl(
            vals_dict=get_site_element_dist(data, site=None),
            site=None,
            stype=selected_stype,
            cmap="Greys",
        )
        cum_filename = (
            f"ElemDist_{'-'.join(selected_stype.split(',')[:2])}.png"
        )
        png_files.append(cum_filename)

        csv_filename = (
            f"{'-'.join(selected_stype.split(',')[:2])}_{dir_name}.csv"
        )

        print(f"\nSaved {len(png_files)} PNG files:")
        for f in png_files:
            print(f"  {f}")
        print(f"Saved CSV file: {csv_filename}")
        print("Thanks for using cif-site-analyzer! teehee\n")
