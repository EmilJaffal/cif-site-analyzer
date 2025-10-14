import numpy as np
import os
from .cif_reader import read_cif
from collections import defaultdict
from collections import Counter


def load_cif(path, dev=False):

    if not os.path.isdir(path):
        print(f"Path {path} not found!")

    data = []
    cifs = [cif for cif in os.listdir(path) if cif.endswith("cif")]

    for cif_name in cifs:
        try:
            cif = read_cif(f"{path}{os.sep}{cif_name}")
            data.append(
                {
                    "Formula": cif.formula,
                    "Filename": cif_name,
                    "Entry prototype": cif.structure_type.replace("~", ""),
                    "num_elements": len(cif.elements),
                    "sites": cif.site_data,
                }
            )
        except Exception:
            print(f"Error reading {cif}")
    if dev and len(data) > 100:
        i100 = np.random.choice(list(range(len(data))), 100, replace=False)
        data = [data[i] for i in i100]

    return data


def select_one_stype(cif_data):
    stypes = defaultdict(int)
    for cif in cif_data:
        stypes[cif["Entry prototype"]] += 1

    return dict(stypes)


def prepare_data_for_engine(cif_data, selected_stype):

    selected_data = [
        d for d in cif_data if d["Entry prototype"] == selected_stype
    ]
    wyckoff_symbols = [
        f"{int(s['multiplicity'])}{s['Wyckoff_symbol']}"
        for s in selected_data[0]["sites"]
    ]
    symbol_counts = dict(Counter(wyckoff_symbols))
    sort_sites = any([v > 1 for k, v in symbol_counts.items()])

    data = []
    for cif in selected_data:
        site_data = cif.pop("sites")

        if sort_sites:
            site_data = sorted(
                site_data,
                key=lambda k: np.linalg.norm(
                    np.array([0.0, 0.0, 0.0])
                    - np.array([k["x"], k["y"], k["x"]])
                ),
            )

        wy_coordinates = {}

        wy_symbol_counts = defaultdict(int)
        added_coordinates = []
        for site in site_data:
            sc = np.array([site[s] for s in ["x", "y", "z"]])
            if len(added_coordinates):
                dists = np.linalg.norm(
                    np.vstack(added_coordinates) - sc, axis=1
                )
                if np.all(dists > 1e-3):
                    wy_symbol_counts[
                        f"{int(site['multiplicity'])}{site['Wyckoff_symbol']}"
                    ] += 1
            else:
                wy_symbol_counts[
                    f"{int(site['multiplicity'])}{site['Wyckoff_symbol']}"
                ] += 1
            added_coordinates.append(sc)

        for site in site_data:
            site_symbol = (
                f"{int(site['multiplicity'])}{site['Wyckoff_symbol']}"
            )
            site_coords = np.array([site[s] for s in ["x", "y", "z"]])

            if not len(wy_coordinates) or site_symbol not in wy_coordinates:
                dist_to_others = None
            else:
                dist_to_others = np.linalg.norm(
                    np.vstack(list(wy_coordinates[site_symbol])) - site_coords,
                    axis=1,
                )

            # modify labels, if necessary
            if wy_symbol_counts[site_symbol] > 1:
                if dist_to_others is None:
                    site_symbol = f"{site_symbol} (1)"
                    wy_coordinates[
                        f"{int(site['multiplicity'])}{site['Wyckoff_symbol']}"
                    ] = [site_coords]
                elif np.any(dist_to_others <= 1e-4):
                    sym_count = int(np.where(dist_to_others <= 1e-4)[0]) + 1
                    site_symbol = f"{site_symbol} ({sym_count})"
                else:
                    sym_count = len(dist_to_others) + 1
                    site_symbol = f"{site_symbol} ({sym_count})"
                    wy_coordinates[
                        f"{int(site['multiplicity'])}{site['Wyckoff_symbol']}"
                    ].append(site_coords)

            comp = f"{site['symbol']}"
            if site["occupancy"] != 1.0:
                comp += f"{round(site['occupancy'], 2)}"
            if site_symbol in cif:
                cif[site_symbol].append(comp)
            else:
                cif[site_symbol] = [comp]

        for k in cif.keys():
            if k not in [
                "Formula",
                "Entry prototype",
                "num_elements",
                "Filename",
            ]:
                cif[k] = "".join(cif[k])
        data.append(cif)

    return data


if __name__ == "__main__":
    path = "/home/bala/Documents/33_CRAFT/old/TiNiSi-oP12-62"
    load_cif(path)
    # read_cif(f"{path}/1221585.cif")
