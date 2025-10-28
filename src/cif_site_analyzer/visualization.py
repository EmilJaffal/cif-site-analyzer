import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from .utils import _parse_formula
from collections import defaultdict

# ---------------- Appearance Variables ----------------
COORD_MULTIPLIER = 5

# Circle settings
BASE_CIRCLE_RADIUS = 0.3
BASE_CIRCLE_LINEWIDTH = 2
PATCH_CIRCLE_RADIUS = 0.3
ALPHA = 0.15
LINE_ALPHA = 0.1

# Text settings
TEXT_SIZE = 18

# Figure size factors (to compute dimensions based on the coordinate range)
FIGURE_WIDTH_FACTOR = 0.8
FIGURE_HEIGHT_FACTOR = 0.8

# Convex hull outline settings
HULL_LINETYPE = "--"
HULL_LINEWIDTH = 2
HULL_FILL_ALPHA = 0.2

# Compound marker settings
COMPOUND_MARKER_SIZE = 15
COMPOUND_LINE_WIDTH = 3
COMPOUND_MARKER_COLOR = "black"  # default compound marker color

# --------------------------------------------------------


def parse_formula_elements(formula):
    return [(k, v) for k, v in _parse_formula(formula).items()]


def custom_log(value):
    if value < 0:
        lv = np.log10(-1 * value)
        return float(lv) * -1
    else:
        return float(np.log10(value))


def visualize_elements(
    coord_df,
    sites_df,
    site_colors,
    compounds_markers=True,
    figsize=None,
    circle_radius=None,
    patch_radius=None,
    text_size=None,
    group_data=None,
    log_scale=False,
):
    """Visualizes chemical element coordinates with site-based patches, convex
    hull outlines, and compound weighted coordinate markers with connecting
    lines.

    Parameters:
      coordinate_file : str
          Path to the Excel file with coordinates.
      sites_df : DataFrame or None
          DataFrame with site information. Expected columns:
          "Filename", "Formula", "M", "X", "R", "Notes".
      compounds_markers : bool
          When True, compute and plot compound markers and connecting lines.

    The function saves the plot with a dynamically generated
    file name and displays it.
    """

    if circle_radius is None:
        circle_radius = BASE_CIRCLE_RADIUS
    if patch_radius is None:
        patch_radius = PATCH_CIRCLE_RADIUS
    if text_size is None:
        text_size = TEXT_SIZE

    # ----- Load coordinates and apply scaling -----
    coord_df = coord_df.copy()
    coord_df["x"] = coord_df["x"] * COORD_MULTIPLIER
    coord_df["y"] = coord_df["y"] * COORD_MULTIPLIER

    if log_scale:
        coord_df["x"] = coord_df["x"].map(custom_log)
        coord_df["y"] = coord_df["y"].map(custom_log)

    # Determine plot dimensions based on the coordinate range.
    x_min, x_max = coord_df["x"].min(), coord_df["x"].max()
    y_min, y_max = coord_df["y"].min(), coord_df["y"].max()

    fig_width = (x_max - x_min) * FIGURE_WIDTH_FACTOR
    fig_height = (y_max - y_min) * FIGURE_HEIGHT_FACTOR

    if figsize is None:
        figsize = (fig_width, fig_height)
    else:
        if fig_height > fig_width:
            figsize = (figsize * fig_width / fig_height, figsize)
        else:
            figsize = (figsize, figsize * fig_height / fig_width)

    # Create figure with white background (no grid, no axes).
    fig, ax = plt.subplots(figsize=figsize, facecolor="white")
    fig.set_size_inches(figsize)
    ax.set_facecolor("white")
    ax.set_aspect("equal")
    ax.axis("off")

    # ----- Plot base circles and element text -----
    element_coords = {}  # maps element symbol to (x, y) coordinates.
    for idx, row in coord_df.iterrows():
        x, y, element = row["x"], row["y"], row["Symbol"]
        element_coords[element] = (x, y)
        base_circle = plt.Circle(
            (x, y),
            radius=circle_radius,
            edgecolor="black",
            facecolor="none",
            linewidth=BASE_CIRCLE_LINEWIDTH,
            alpha=1.0,
            zorder=3,
        )
        ax.add_patch(base_circle)
        ax.text(
            x,
            y,
            element,
            fontsize=text_size,
            ha="center",
            va="center",
            zorder=4,
        )

    # ----- Process site information to extract unique element sets -----
    unique_sites = defaultdict(set)
    if sites_df is not None:
        sites = site_colors.keys()
        for _, row in sites_df.iterrows():
            for site in sites:
                site_formula = row[site]
                site_formula = sorted(
                    [[k, v] for k, v in _parse_formula(site_formula).items()],
                    key=lambda s: s[-1],
                )
                unique_sites[site].add(site_formula[-1][0])

    # ----- Draw colored patches on top of base circles -----
    for element, (x, y) in element_coords.items():
        patch_color = None
        for i, site in enumerate(sites):
            if element in unique_sites[site]:
                patch_color = site_colors[site]
                break
        if patch_color is not None:
            patch_circle = plt.Circle(
                (x, y),
                radius=patch_radius,
                edgecolor=None,
                facecolor=patch_color,
                alpha=0.5,
                zorder=3.5,
            )
            ax.add_patch(patch_circle)

    # ----- Define helper function to plot convex hull outlines -----
    def plot_site_outline(element_set, color, label):
        points = []
        for elem in element_set:
            if elem in element_coords:
                points.append(element_coords[elem])
        points = np.array(points)
        n_points = len(points)
        if n_points < 2:
            return
        elif n_points == 2:
            ax.plot(
                points[:, 0],
                points[:, 1],
                linestyle=HULL_LINETYPE,
                color=color,
                linewidth=HULL_LINEWIDTH,
                label=label,
            )
        else:
            hull = ConvexHull(points)
            hull_points = points[hull.vertices]
            hull_points_closed = np.vstack([hull_points, hull_points[0]])
            ax.plot(
                hull_points_closed[:, 0],
                hull_points_closed[:, 1],
                linestyle=HULL_LINETYPE,
                color=color,
                linewidth=HULL_LINEWIDTH,
                label=label,
            )
            ax.fill(
                hull_points_closed[:, 0],
                hull_points_closed[:, 1],
                color=color,
                alpha=HULL_FILL_ALPHA,
            )

    # ----- Draw convex hull outlines for each site group -----
    if sites_df is not None:
        for i, site in enumerate(sites):
            # color = SITE_OUTLINE_COLORS[i % len(SITE_OUTLINE_COLORS)]
            color = site_colors[site]
            plot_site_outline(unique_sites[site], color, f"{site} Outline")

    # ----- Compute weighted compound markers from formula column -----
    candidate_found = False  # flag to track if there is any candidate entry
    if sites_df is not None and compounds_markers:
        for idx, row in sites_df.iterrows():
            formula = row["Formula"]
            items = parse_formula_elements(formula)
            if not items:
                continue
            total_count = sum(count for _, count in items)
            weighted_sum = np.array([0.0, 0.0])
            compound_elements = (
                []
            )  # list of elements with coordinates in the compound
            for el, count in items:
                if el in element_coords:
                    compound_elements.append(el)
                    weighted_sum += count * np.array(element_coords[el])
            if total_count > 0 and compound_elements:
                weighted_coord = weighted_sum / total_count
                note = str(row.get("Notes", "")).strip().lower()
                # Check if this compound row is marked as candidate.
                if note == "candidate":

                    candidate_found = True
                    marker_color = "#F60303"
                    size = COMPOUND_MARKER_SIZE
                    line_width = COMPOUND_LINE_WIDTH
                    ax.plot(
                        weighted_coord[0],
                        weighted_coord[1],
                        marker="o",
                        markersize=size,
                        markerfacecolor="#F60303",
                        zorder=4,
                        alpha=1.0,
                        markeredgecolor="none",
                    )
                    linestyle = "-"
                elif note == "unexplored":
                    # print("unexp")
                    size = COMPOUND_MARKER_SIZE
                    line_width = COMPOUND_LINE_WIDTH
                    ax.plot(
                        weighted_coord[0],
                        weighted_coord[1],
                        marker="o",
                        markersize=size,
                        markerfacecolor="#03C91D",
                        zorder=4,
                        alpha=ALPHA,
                        markeredgecolor="none",
                    )
                    linestyle = "-"
                else:
                    marker_color = COMPOUND_MARKER_COLOR
                    size = COMPOUND_MARKER_SIZE
                    line_width = COMPOUND_LINE_WIDTH
                    # For standard markers, use the filled marker version:
                    ax.plot(
                        weighted_coord[0],
                        weighted_coord[1],
                        marker="o",
                        markersize=size,
                        markerfacecolor=marker_color,
                        zorder=3,
                        alpha=ALPHA,
                        markeredgecolor="none",
                    )
                    linestyle = "-"

                # Draw connecting lines from the compound marker
                # to each element of the compound.
                for el in compound_elements:
                    el_coord = element_coords[el]
                    ax.plot(
                        [weighted_coord[0], el_coord[0]],
                        [weighted_coord[1], el_coord[1]],
                        color="grey",
                        linestyle=linestyle,
                        linewidth=line_width,
                        zorder=4,
                        alpha=LINE_ALPHA,
                    )
    # plt.plot([0, 0], [-1, 1])
    # plt.plot([-1, 1], [0, 0])
    plt.tight_layout()

    # ----- Build output file name based on conditions -----
    output_filename = "outputs/plots/elements_visualization"
    if compounds_markers:
        output_filename += "_compound"
    if sites_df is not None and candidate_found:
        output_filename += "_candidate"

    output_filepath = output_filename + ".svg"
    plt.savefig(output_filepath, dpi=500, bbox_inches="tight")
    # plt.show()
