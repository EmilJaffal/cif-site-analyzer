import re
import matplotlib
import matplotlib.patches
import numpy as np
import matplotlib.pyplot as plt
from .utilities import format_wyckoff_site_label
from .utilities import format_formula
from .utilities import format_mixed_formula


def ptable_heatmap_mpl(
    vals_dict: dict,
    site: str,
    stype: str,
    cmap: str,
    plot_filename: str,
    title=None,
):

    # format site
    elements_in_data = list(vals_dict.keys())
    if site is not None:
        site = format_wyckoff_site_label(site, italicize=True)

    elements = {
        "H": [1, 1],
        "He": [1, 18],
        "Li": [2, 1],
        "Be": [2, 2],
        "B": [2, 13],
        "C": [2, 14],
        "N": [2, 15],
        "O": [2, 16],
        "F": [2, 17],
        "Ne": [2, 18],
        "Na": [3, 1],
        "Mg": [3, 2],
        "Al": [3, 13],
        "Si": [3, 14],
        "P": [3, 15],
        "S": [3, 16],
        "Cl": [3, 17],
        "Ar": [3, 18],
        "K": [4, 1],
        "Ca": [4, 2],
        "Sc": [4, 3],
        "Ti": [4, 4],
        "V": [4, 5],
        "Cr": [4, 6],
        "Mn": [4, 7],
        "Fe": [4, 8],
        "Co": [4, 9],
        "Ni": [4, 10],
        "Cu": [4, 11],
        "Zn": [4, 12],
        "Ga": [4, 13],
        "Ge": [4, 14],
        "As": [4, 15],
        "Se": [4, 16],
        "Br": [4, 17],
        "Kr": [4, 18],
        "Rb": [5, 1],
        "Sr": [5, 2],
        "Y": [5, 3],
        "Zr": [5, 4],
        "Nb": [5, 5],
        "Mo": [5, 6],
        "Tc": [5, 7],
        "Ru": [5, 8],
        "Rh": [5, 9],
        "Pd": [5, 10],
        "Ag": [5, 11],
        "Cd": [5, 12],
        "In": [5, 13],
        "Sn": [5, 14],
        "Sb": [5, 15],
        "Te": [5, 16],
        "I": [5, 17],
        "Xe": [5, 18],
        "Cs": [6, 1],
        "Ba": [6, 2],
        "La": [8, 4],
        "Ce": [8, 5],
        "Pr": [8, 6],
        "Nd": [8, 7],
        "Pm": [8, 8],
        "Sm": [8, 9],
        "Eu": [8, 10],
        "Gd": [8, 11],
        "Tb": [8, 12],
        "Dy": [8, 13],
        "Ho": [8, 14],
        "Er": [8, 15],
        "Tm": [8, 16],
        "Yb": [8, 17],
        "Lu": [8, 18],
        "Hf": [6, 4],
        "Ta": [6, 5],
        "W": [6, 6],
        "Re": [6, 7],
        "Os": [6, 8],
        "Ir": [6, 9],
        "Pt": [6, 10],
        "Au": [6, 11],
        "Hg": [6, 12],
        "Tl": [6, 13],
        "Pb": [6, 14],
        "Bi": [6, 15],
        "Po": [6, 16],
        "At": [6, 17],
        "Rn": [6, 18],
        "Ac": [9, 4],
        "Th": [9, 5],
        "Pa": [9, 6],
        "U": [9, 7],
        "Np": [9, 8],
        "Pu": [9, 9],
        "Am": [9, 10],
        "Cm": [9, 11],
        "Bk": [9, 12],
        "Cf": [9, 13],
        "Es": [9, 14],
        "Fm": [9, 15],
        "Md": [9, 16],
        "No": [9, 17],
        "Lr": [9, 18],
    }

    pmask = np.zeros(shape=(9, 18), dtype=float)
    for k, v in elements.items():
        p, g = v
        pmask[p - 1, g - 1] = 1

    min_value = min(vals_dict.values())
    norm_const = max(max(vals_dict.values()) - min_value, min_value)

    heat_map = np.full(shape=(9, 18), fill_value=-(min_value), dtype=float)
    for k, v in elements.items():
        if k in vals_dict:
            heat_map[v[0] - 1, v[1] - 1] = vals_dict[k]
    heat_map *= pmask
    heat_map[heat_map == 0.0] = np.nan

    fig, ax = plt.subplots(figsize=(18, 9))

    cmap = plt.get_cmap(cmap)
    cmap.set_under("w", alpha=0.1)

    im = plt.imshow(heat_map, cmap=cmap, vmin=min_value)
    # Minor ticks
    ax.set_xticks(np.arange(-0.5, 18, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, 9, 1), minor=True)

    # Gridlines based on minor ticks
    # ax.grid(which='minor', color='w', linestyle='-', linewidth=3)

    ticks = sorted(
        set(
            [
                int(i)
                for i in np.linspace(min_value, max(vals_dict.values()), 4)
            ]
        )
    )
    im.set_clim(vmin=ticks[0], vmax=ticks[-1])
    # divider = make_axes_locatable(ax)

    cbar = plt.colorbar(
        im,
        ticks=ticks,
        location="top",
        shrink=0.25,
        cax=ax.inset_axes((0.19, 0.72, 0.4, 0.02)),
        orientation="horizontal",
    )
    cbar.ax.tick_params(labelsize=25)

    # element symbols
    kw = dict(
        horizontalalignment="center", verticalalignment="center", size=25
    )

    im.axes.text(2, 5.1, "*", alpha=0.6, **kw)
    im.axes.text(2, 7.1, "*", alpha=0.6, **kw)
    for k, v in elements.items():

        c = (
            "w"
            if ((vals_dict.get(k, 0) - min_value) / norm_const) > 0.6
            else "k"
        )
        if k in elements_in_data:
            im.axes.text(v[1] - 1, v[0] - 1, k, c=c, **kw)
        else:
            im.axes.text(v[1] - 1, v[0] - 1, k, c=c, alpha=0.6, **kw)

        rect = matplotlib.patches.Rectangle(
            xy=(v[1] - 1.5, v[0] - 1.5),
            width=1,
            height=1,
            facecolor=(0, 0, 0, 0),
            edgecolor=(0, 0, 0, 1),
            lw=1,
        )
        im.axes.add_patch(rect)

    # placeholder for title
    rectangle = matplotlib.patches.Rectangle(
        xy=(2.0, -0.59), width=9.0, height=1.7, facecolor="w", edgecolor=None
    )
    ax.add_patch(rectangle)
    rx, ry = rectangle.get_xy()
    cx = rx + rectangle.get_width() / 2.0
    cy = ry + rectangle.get_height() / 2.0

    # Compose title and wrap
    def safe_format_formula(formula):
        try:
            # Remove parentheses if present
            if formula.startswith("(") and formula.endswith(")"):
                formula = formula[1:-1]
            # If it looks like a mixed element formula
            # (e.g. Cr0.16Mo0.38Co0.46)
            if re.match(r"^([A-Z][a-z]*[0-9.]+)+$", formula):
                return format_mixed_formula(formula)
            # Only try to parse if it looks like a formula
            if re.match(r"^[A-Za-z0-9.]+$", formula):
                return format_formula(formula)
            else:
                return formula
        except Exception:
            return formula

    if site is not None:
        text_str = f"element distribution of {site} site in \
        {safe_format_formula(stype.split(',')[0])}-type compounds"
    else:
        if title is not None:
            text_str = title
        else:
            text_str = f"element distribution in \
                {safe_format_formula(stype.split(',')[0])}-type compounds"

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

    ax.annotate(
        wrapped_text,
        xy=(cx, cy),
        color="black",
        fontsize=25,
        ha="center",
        va="center",
        wrap=False,
    )

    # remove ticks
    ax.axes.get_xaxis().set_ticks([])
    ax.axes.get_yaxis().set_ticks([])
    # Remove minor ticks
    ax.tick_params(which="minor", bottom=False, left=False)

    # Remove all spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)

    if title is None:
        filename = "-".join(stype.split(",")[:2])
    else:
        filename = title
    if site is not None:
        plt.savefig(
            plot_filename,
            dpi=300,
        )
    else:
        plt.savefig(f"ElemDist_{filename}.png", dpi=300)
