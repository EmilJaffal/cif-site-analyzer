import pandas as pd
import numpy as np
import importlib.resources
from .utils import _parse_formula


def add_features(df, features_csv=None):
    if features_csv is None:
        with importlib.resources.as_file(
            importlib.resources.files("cif_site_analyzer.data")
            / "elemental-property-list.csv"
        ) as csv_path:
            features = pd.read_csv(csv_path)
    else:
        features = pd.read_csv(features_csv)

    features = features.dropna(axis=1)
    elements_available = features["Symbol"].tolist()

    site_labels = list(df.columns[2:])
    warned_elements = []

    def compute_features(formula, normalize=True):

        formula = _parse_formula(formula)
        if normalize:
            total = sum(formula.values())
            if total > 0:
                for k in formula.keys():
                    formula[k] = formula[k] / total

        symbols, coeffs = [], []
        for k, v in formula.items():
            if k not in elements_available:
                if k not in warned_elements:
                    print(
                        f"Warning: Element {k} not found in properties data."
                    )
                    warned_elements.append(k)
                continue

            symbols.append(k)
            coeffs.append(v)

        coeffs = np.array(coeffs)

        row_features = {}
        for feature in features.columns[1:]:
            fvals = np.array(
                [
                    float(
                        features[features["Symbol"] == el][feature].values[0]
                    )
                    for el in symbols
                ]
            )
            row_features[feature] = (fvals * coeffs).sum()
        return row_features

    output_rows = []
    for _, row in df.iterrows():
        for site_label in site_labels:
            if (
                site_label in row
                and pd.notna(row[site_label])
                and str(row[site_label]) != ""
            ):
                site_formula = row[site_label]
                out_row = {
                    "Formula": row["Formula"],
                    "Site_label": site_label,
                    "Site_formula": site_formula,
                }

                feature_dict = compute_features(site_formula, normalize=True)
                out_row.update(feature_dict)
                output_rows.append(out_row)

    out_df = pd.DataFrame(output_rows)
    ordered_cols = ["Formula", "Site_label", "Site_formula"] + list(
        features.columns[1:]
    )
    out_df = out_df[ordered_cols]

    return out_df
