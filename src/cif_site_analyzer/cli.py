from .get_data import load_cif
from .get_data import select_one_stype
from .get_data import prepare_data_for_engine
from .utils import get_ptable_vals_dict
from .utils import list_to_formula
from .utils import assign_labels_for_sites
from .utils import auto_group_sites
from .utils import get_valid_input
from .features import add_features
from .ptable_histogram import ptable_heatmap_mpl
from .plsda import run_pls_da
from .visualization import visualize_elements
from .projection import plot_elements_from_plsda_loadings
import pandas as pd
import argparse
import os


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


def main():
    description = """
        A High Throughput LMTO Calculator.

        This code automates the LMTO calculations using TB-LMTO program.

        """
    parser = argparse.ArgumentParser(
        description=description, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("input_path", help="Path to the input file / folder")

    args = parser.parse_args()

    path = args.input_path
    if not os.path.isdir(path):
        print(f"Cannot find {path}.")

    # DATA
    cif_data = load_cif(path, dev=True)
    stypes = select_one_stype(cif_data)

    if len(cif_data) == 0:
        print("No data found!")
        exit(0)

    if len(stypes) > 1:
        stypes = [[k, v] for k, v in stypes.items()]
        print(
            "\n\nMore than one structure types are found in the input list. \
                \nPlease select one structure type from the list below.\n"
        )

        for i, (stype, count) in enumerate(stypes, 1):
            print(f"({i}) {count:<5} {stype}")

        valid_input = False
        while not valid_input:
            res = input(
                "Please enter the number corresponding to the \
                    selected structure type: \n"
            )
            try:
                res = int(res)
                assert 0 < res <= len(stypes)
                selected_stype = stypes[res - 1][0]
                valid_input = True
            except Exception:
                # print(e)
                print("Invalid input. Try again.")
    else:
        selected_stype = list(stypes.keys())[0]

    data_for_engine = prepare_data_for_engine(cif_data, selected_stype)

    # Save histograms
    data_df = pd.DataFrame(data_for_engine)
    wyckoff_symbols = data_df.columns[4:]

    # RMX assignment
    wyckoff_symbol_elements = {}
    for ws in wyckoff_symbols:
        wyckoff_symbol_elements[ws] = list(data_df[ws].unique())

    auto_assignment = auto_group_sites(wyckoff_symbol_elements)
    site_assignment = assign_labels_for_sites(wyckoff_symbols, auto_assignment)

    print("\nGenerating periodic table heatmaps...")
    os.makedirs("outputs/plots", exist_ok=True)
    for i, (k, v) in enumerate(site_assignment.items()):
        v = list(v)
        print(f"Plotting periodic table heatmap for {k}: {v} site")
        ptable_heatmap_mpl(
            vals_dict=get_ptable_vals_dict(data_df[v[0]].tolist()),
            site=v,
            stype=data_df.iloc[0]["Entry prototype"],
            cmap=cmaps[i % len(cmaps)],
        )
    print("Done, plots saved inside the directory plots.")

    # df for recommendation engine
    data = []
    for cif_data in data_for_engine:
        row = {
            "Filename": cif_data["Filename"],
            "Formula": cif_data["Formula"],
        }
        for group_label, sites in site_assignment.items():
            gcomp = []
            for site in sites:
                gcomp.append(cif_data[site])
            row[group_label] = list_to_formula(gcomp)
        data.append(row)

    df_engine = pd.DataFrame(data)
    df_engine["Notes"] = ""

    # add features
    print("\nLoading features...")
    features_df = add_features(df_engine)
    print("Done")

    # pls-da; 2-components
    print("\nPerforming PLS-DA...")
    pls_loadings = run_pls_da(features_df)
    print("Done")

    # elements projection
    print("\nPlotting projections of elements and compounds...")
    coords = plot_elements_from_plsda_loadings(pls_loadings, df_engine)
    visualize_elements(coords, df_engine, compounds_markers=False)
    visualize_elements(coords, df_engine, compounds_markers=True)
    print("Done, plots saved in the directory plots.")

    print(df_engine.head())

    # get additional compounds for overlay
    def get_additional_cpds():
        res = get_valid_input(
            "\nDo you want to add compounds to compare ? ",
            ["Y", "y", "N", "n"],
        )[0]
        new_cpds = []

        if res.lower() == "y":
            print(f"The current site-group assignment is {site_assignment}.")

            all_added = False
            nc = 1
            while not all_added:
                new_cpd = get_valid_input(
                    "Enter the elements (without coefficients) at each \n"
                    "site-group\n for the new compound separated by comma : ",
                    None,
                )
                new_cpd = dict(zip(site_assignment.keys(), new_cpd))
                new_cpd["Formula"] = "".join(list(new_cpd.values()))
                new_cpd["Notes"] = "candidate"
                new_cpd["Filename"] = f"NewCand_{nc}"
                new_cpds.append(new_cpd)

                done = get_valid_input(
                    "\nDo you wish to add more compounds ? ",
                    ["Y", "y", "N", "n"],
                )[0]
                if done.lower() == "n":
                    all_added = True
        return new_cpds

    new_cpds = get_additional_cpds()
    new_cpds = pd.DataFrame(new_cpds)

    df_engine = pd.concat([df_engine, new_cpds], axis=0, ignore_index=True)
    visualize_elements(coords, df_engine, compounds_markers=True)

    # elements for recommendations
    # elements_for_screening = get_selected_elements(site_assignment)

    # projection
    # recommendation_engine(site_element_pools=elements_for_screening,
    #                       sites_df=df_engine,
    #                       coord_df=coords,
    #                       output_file="outputs/recommendations.csv")

    # Sm4Ir2InGe4, Tb4Rh2InGe4, Lu4Ni2InGe4
