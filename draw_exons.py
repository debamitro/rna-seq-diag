#!python3

import matplotlib.pyplot as plt
import matplotlib.lines as lines
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib import rcParams


def make_exon_shapes(exons, y, color="xkcd:mustard"):
    """Creates matplotlib patches representing
    a series of exons"""
    patches = []

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
        poly = Polygon(points, True, color=color)
        patches.append(poly)

    return patches


def make_exon_exon_lines(exons, ax, y, height=5, color="xkcd:light brown"):
    """Creates matplotlib lines representing
    the order in which a set of exons have been seen in a transcript"""
    previous_exon = None
    for exon in exons:
        if previous_exon is not None:
            line = lines.Line2D(
                [previous_exon[1], (exon[0] + previous_exon[1]) / 2, exon[0]],
                [y + 20, y + 20 + height, y + 20],
                color=color,
            )
            ax.add_line(line)
        previous_exon = exon


def draw_exons(exons, file_name=None, transcript_id=None):
    """Given an array of [start,end] offsets for exons,
    draws them in a diagram. Optionally writes out the diagram
    as a file, and also adds the transcript id"""
    y = 130
    if transcript_id is not None:
        plt.text(exons[0][0], y + 20, transcript_id)
    draw_exon_sequence_graph({"exons": exons, "sequences": [exons]}, y, file_name)


def draw_transcripts(transcripts, file_name=None):
    """Given a dictionary where the keys are transcript IDs
    and the values are arrays of exon start and end offsets,
    draws them in a diagram. Optionally saves out the diagram
    in a file."""

    _, ax = plt.subplots()

    ymax = len(transcripts) * 40 + 20
    y = ymax
    patches = []
    xleft, xright = None, None
    for transcript_id in transcripts:
        y -= 40
        exons = transcripts[transcript_id]
        patches.extend(make_exon_shapes(exons, y))
        make_exon_exon_lines(exons, ax, y)
        if xleft is None or exons[0][0] < xleft:
            xleft = exons[0][0]
        if xright is None or exons[len(exons) - 1][1] > xright:
            xright = exons[len(exons) - 1][1]

    p = PatchCollection(patches, alpha=0.3)

    start_margin = end_margin = 100

    ax.set_xbound(xleft - start_margin, xright + end_margin)
    ax.set_ybound(0, ymax)
    ax.add_collection(p)

    if file_name is None:
        plt.show()
    else:
        plt.savefig(file_name)


def draw_exon_sequence_graph(sequence_graph, y_exons=130, file_name=None):
    """Given a dictionary with two entries
     - 'exons' an array of exon start and end offsets
     - 'sequences' an array of exon sequences
    draws a graph using different colors for each sequence.
    The goal is to show different exon sequences formed from
    one universal set of exons"""
    _, ax = plt.subplots()

    exons = sequence_graph["exons"]
    patches = make_exon_shapes(exons, y_exons)
    p = PatchCollection(patches, alpha=0.3)

    colors = ["xkcd:indigo", "xkcd:forest green", "xkcd:navy blue"]
    sequence_height = 5
    sequence_index = 0
    for sequence in sequence_graph["sequences"]:
        make_exon_exon_lines(
            sequence, ax, y_exons, height=sequence_height, color=colors[sequence_index]
        )
        sequence_height += 5
        sequence_index += 1
        if sequence_index > len(colors):
            sequence_index = 0

    start_margin = end_margin = 100

    ax.set_xbound(exons[0][0] - start_margin, exons[len(exons) - 1][1] + end_margin)
    ax.set_ybound(0, 200)
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
            "ENST00000456328.2": [[11869, 12227], [12613, 12721], [13221, 14409]],
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

    # Contrived example using some exons from the above
    draw_exon_sequence_graph(
        {
            "exons": [
                [12010, 12057],
                [12179, 12227],
                [12613, 12619],
                [12975, 13052],
                [13221, 13374],
                [13453, 13670],
            ],
            "sequences": [
                [[12010, 12057], [12179, 12227], [12613, 12619], [12975, 13052]],
                [
                    [12010, 12057],
                    [12613, 12619],
                    [12975, 13052],
                    [13221, 13374],
                    [13453, 13670],
                ],
            ],
        },
        file_name="out4.png",
    )
