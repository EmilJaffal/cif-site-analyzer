import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from collections import defaultdict
from .utils import _parse_formula

# from .utils import SiteElement


def prepare_coordinates(pls_loadings, element_properties_file):
    df_props = pd.read_csv(element_properties_file)
    df_props.dropna(axis=1, how="all", inplace=True)

    prop_columns = [col for col in df_props.columns if col != "Symbol"]
    pls_loadings["Feature_std"] = (
        pls_loadings["Feature"].str.strip().str.lower()
    )
    prop_columns_std = {col: col.strip().lower() for col in prop_columns}

    common_features = [
        orig_col
        for orig_col, std_col in prop_columns_std.items()
        if std_col in set(pls_loadings["Feature_std"])
    ]

    if not common_features:
        raise ValueError(
            "No common features found between PLS‑DA \
                loadings and element property file columns!"
        )

    X = df_props[common_features].values
    X_scaled = StandardScaler().fit_transform(X)
    X_normalized = MinMaxScaler().fit_transform(X_scaled)

    loadings = []
    for feature in common_features:
        key = feature.strip().lower()
        row = pls_loadings.loc[
            pls_loadings["Feature_std"] == key,
            ["Component_1_Loading", "Component_2_Loading"],
        ]
        loadings.append(row.iloc[0].values if not row.empty else [0, 0])
    loadings_matrix = np.array(loadings)

    coordinates = X_normalized.dot(loadings_matrix)

    df_props["PLSDA_Component_1"] = coordinates[:, 0]
    df_props["PLSDA_Component_2"] = coordinates[:, 1]

    return df_props, coordinates


def plot_elements(df, group_colors, output_path="plots/elements_plot.svg"):
    plt.figure(figsize=(8, 6), dpi=500)
    plt.style.use("ggplot")

    for group in df["Group"].unique():
        subset = df[df["Group"] == group]
        color = group_colors.get(group, "grey")
        plt.scatter(
            subset["PLSDA_Component_1"],
            subset["PLSDA_Component_2"],
            color=color,
            s=600,
            alpha=0.6,
            label=group,
        )
        for _, row in subset.iterrows():
            plt.text(
                row["PLSDA_Component_1"],
                row["PLSDA_Component_2"],
                str(row["Symbol"]),
                fontsize=18,
                horizontalalignment="center",
                verticalalignment="center",
            )

    plt.xlabel("PC1", fontsize=18)
    plt.ylabel("PC2", fontsize=18)
    plt.gca().set_aspect("equal", adjustable="box")
    plt.tick_params(axis="both", labelsize=18)
    plt.tight_layout()
    plt.savefig(output_path, dpi=500)
    # plt.show()


def plot_elements_from_plsda_loadings(
    pls_loadings,
    sites_df=None,
    element_properties_file="data/elemental-property-list.csv",
):
    df_props, coordinates = prepare_coordinates(
        pls_loadings, element_properties_file
    )

    # Grouping logic (optional)
    group_colors = [
        "tab:blue",
        "tab:orange",
        "tab:green",
        "tab:red",
        "tab:purple",
        "tab:brown",
        "tab:pink",
        "tab:gray",
        "tab:olive",
        "tab:cyan",
    ]

    unique_sites = defaultdict(set)
    if sites_df is not None:
        sites = [
            cn
            for cn in sites_df.columns
            if cn not in ["Filename", "Formula", "Notes"]
        ]

        for _, row in sites_df.iterrows():
            for site in sites:
                site_formula = row[site]
                site_formula = sorted(
                    [[k, v] for k, v in _parse_formula(site_formula).items()],
                    key=lambda s: s[-1],
                )
                unique_sites[site].add(site_formula[-1][0])

        def get_group(symbol):
            for site, elements in unique_sites.items():
                if symbol in elements:
                    return site

            return "Other"

        df_props["Group"] = df_props["Symbol"].apply(get_group)
    else:
        df_props["Group"] = "Other"

    # Plot
    site_colors = {
        s: group_colors[i % len(unique_sites)] for i, s in enumerate(sites)
    }
    plot_elements(df_props, site_colors)

    # Save coordinates
    coordinates_df = df_props[
        ["Symbol", "PLSDA_Component_1", "PLSDA_Component_2"]
    ].copy()
    coordinates_df.columns = ["Symbol", "x", "y"]
    coordinates_df.to_csv("outputs/coordinates.csv", index=False)
    print("Element coordinates saved to outputs/coordinates.csv")

    return coordinates_df
