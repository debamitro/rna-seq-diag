#!python3

import matplotlib.pyplot as plt
import matplotlib.lines as lines
from matplotlib.patches import Rectangle

configuration = {
    "exon_height": 20,
    "exon_color": "xkcd:mustard",
    "unscaled_exon_width": 1000,
    "unscaled_exon_start": 2000,
}


def make_exon_shapes(exons, y, color=configuration["exon_color"]):
    """Creates matplotlib patches representing
    a series of exons"""
    patches = []

    for exon in exons:
        start_x = exon[0]
        end_x = exon[1]
        rect = Rectangle(
            [start_x, y], end_x - start_x, configuration["exon_height"], color=color
        )
        patches.append(rect)

    return patches


def make_exon_exon_lines(
    exon_pairs, ax, y, height=5, draw_at="top", color="xkcd:light brown"
):
    """Creates matplotlib lines which may (or may not) represent
    the order in which a set of exons have been seen in a transcript.
    Mandatory arguments:
        - exon_pairs - a list of exon pairs, which indicate the lines
        to be drawn
        - ax - the matplotlib Axes object
        - y - the point on the y axis where the lines should be drawn"""
    if draw_at == "top":
        y_triplet = [
            y + configuration["exon_height"],
            y + configuration["exon_height"] + height,
            y + configuration["exon_height"],
        ]
    elif draw_at == "mid":
        y_triplet = [
            y + configuration["exon_height"] / 2,
            y + configuration["exon_height"] / 2,
            y + configuration["exon_height"] / 2,
        ]
    else:
        y_triplet = [y, y - height, y]

    for exon_pair in exon_pairs:
        # Draw a line from exon_pair[0][1] to exon_pair[1][0]
        # It can go straight, curvy on the top, or curvy on the bottom
        line = lines.Line2D(
            [exon_pair[0][1], (exon_pair[1][0] + exon_pair[0][1]) / 2, exon_pair[1][0]],
            y_triplet,
            color=color,
        )
        ax.add_line(line)


def make_exons_unscaled(exons):
    unscaled_mapping = {}
    cur_x = configuration["unscaled_exon_start"]
    for exon in exons:
        unscaled_mapping[exon] = (
            cur_x,
            cur_x + configuration["unscaled_exon_width"],
        )
        cur_x += configuration["unscaled_exon_width"] * 2
    unscaled_exons = [unscaled_mapping[x] for x in exons]
    return unscaled_mapping, unscaled_exons
