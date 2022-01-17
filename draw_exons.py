#!python3

import matplotlib.pyplot as plt
import matplotlib.lines as lines
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib import rcParams


def make_exon_shapes(exons, ax, y):
    """Creates matplotlib patches representing
    a series of exons as identified in a transcript"""
    patches = []

    previous_exon = None
    for exon in exons:
        start_x = exon[0]
        end_x = exon[1]
        points = [
            [start_x, y + 20],
            [end_x, y + 20],
            [end_x + 10, y + 10],
            [end_x, y],
            [start_x, y],
            [start_x + 10, y + 10],
        ]
        poly = Polygon(points, True)
        patches.append(poly)
        if previous_exon is not None:
            line1 = lines.Line2D(
                [previous_exon[1], (exon[0] + previous_exon[1]) / 2, exon[0]],
                [y + 20, y + 25, y + 20],
            )
            ax.add_line(line1)
        previous_exon = exon

    return patches


def draw_exons(exons, file_name=None, transcript_id=None):
    """Given an array of [start,end] offsets for exons,
    draws them in a diagram. Optionally writes out the diagram
    as a file, and also adds the transcript id"""
    _, ax = plt.subplots()

    y = 130
    patches = make_exon_shapes(exons, ax, y)

    p = PatchCollection(patches, alpha=0.3)

    colors = []
    for i in range(len(exons)):
        colors.append(80)
    p.set_array(colors)

    start_margin = end_margin = 100

    ax.set_xbound(exons[0][0] - start_margin, exons[len(exons) - 1][1] + end_margin)
    ax.set_ybound(0, 200)
    ax.add_collection(p)

    if transcript_id is not None:
        plt.text(exons[0][0], y + 20, transcript_id)

    if file_name is None:
        plt.show()
    else:
        plt.savefig(file_name)


def draw_transcripts(transcripts, file_name=None):
    """Given a dictionary where the keys are transcript IDs
    and the values are arrays of exon start and end offsets,
    draws them in a diagram. Optionally saves out the diagram
    in a file."""

    _, ax = plt.subplots()

    ymax = len(transcripts) * 40 + 20
    y = ymax
    patches = []
    colors = []
    xleft, xright = None, None
    for transcript_id in transcripts:
        y -= 40
        exons = transcripts[transcript_id]
        patches.extend(make_exon_shapes(exons, ax, y))
        for i in range(len(exons)):
            colors.append(80)
        if xleft is None or exons[0][0] < xleft:
            xleft = exons[0][0]
        if xright is None or exons[len(exons) - 1][1] > xright:
            xright = exons[len(exons) - 1][1]

    p = PatchCollection(patches, alpha=0.3)

    p.set_array(colors)

    start_margin = end_margin = 100

    ax.set_xbound(xleft - start_margin, xright + end_margin)
    ax.set_ybound(0, ymax)
    ax.add_collection(p)

    if file_name is None:
        plt.show()
    else:
        plt.savefig(file_name)


if __name__ == "__main__":
    # For better resolution, maybe?
    # rcParams['figure.dpi'] = 200.0

    # Example taken from transcript ENST00000456328.2 of the gene DDX11L1
    draw_exons([[11869, 12227], [12613, 12721], [13221, 14409]], file_name="out1.png")

    # Example taken from transcript ENST00000450305.2 of the gene DDX11L1
    draw_exons(
        [
            [12010, 12057],
            [12179, 12227],
            [12613, 12619],
            [12975, 13052],
            [13221, 13374],
            [13453, 13670],
        ],
        file_name="out2.png",
    )

    # Combining the two transcripts ENST00000456328.2 and ENST00000450305.2
    draw_transcripts(
        {
            "ENST00000456328.2": [
                [11869, 12227],
                [12613, 12721],
                [13221, 14409]
            ],
            "ENST00000450305.2": [
                [12010, 12057],
                [12179, 12227],
                [12613, 12619],
                [12975, 13052],
                [13221, 13374],
                [13453, 13670],
            ],
        }
    )
