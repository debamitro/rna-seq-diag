#!python3

import matplotlib.pyplot as plt
import matplotlib.lines as lines
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib import rcParams
import numpy as np

configuration = {"unscaled_exon_width": 1000, "exon_height": 20}


def make_exon_shapes(exons, y, color="xkcd:mustard"):
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
    exons, ax, y, height=5, draw_at_top=True, color="xkcd:light brown"
):
    """Creates matplotlib lines representing
    the order in which a set of exons have been seen in a transcript"""
    previous_exon = None
    if draw_at_top:
        y_triplet = [
            y + configuration["exon_height"],
            y + configuration["exon_height"] + height,
            y + configuration["exon_height"],
        ]
    else:
        y_triplet = [y, y - height, y]

    for exon in exons:
        if previous_exon is not None:
            line = lines.Line2D(
                [previous_exon[1], (exon[0] + previous_exon[1]) / 2, exon[0]],
                y_triplet,
                color=color,
            )
            ax.add_line(line)
        previous_exon = exon


def draw_exons(exons, file_name=None, transcript_id=None):
    """Given an array of (start,end) offsets for exons,
    draws them in a diagram. Optionally writes out the diagram
    as a file, and also adds the transcript id"""
    y = 130
    if transcript_id is not None:
        plt.text(exons[0][0], y + 20, transcript_id)
    draw_exon_sequence_graph(
        {"exons": exons, "sequences": [exons]}, y, file_name, transcript_id
    )


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


def draw_exon_sequence_graph(
    sequence_graph, y_exons=130, file_name=None, title=None, to_scale=True
):
    """Given a dictionary with two entries
     - 'exons' an array of exon start and end offsets
     - 'sequences' an array of exon sequences
    draws a graph using different colors for each sequence.
    The goal is to show different exon sequences formed from
    one universal set of exons"""
    _, ax = plt.subplots()

    exons = sequence_graph["exons"]
    if not to_scale:
        unscaled_mapping = {}
        cur_x = 100
        for exon in exons:
            unscaled_mapping[exon] = (
                cur_x,
                cur_x + configuration["unscaled_exon_width"],
            )
            cur_x += configuration["unscaled_exon_width"] * 2
        unscaled_exons = [unscaled_mapping[x] for x in exons]
        exons = unscaled_exons

    patches = make_exon_shapes(exons, y_exons)
    p = PatchCollection(patches)

    colors = ["xkcd:indigo", "xkcd:forest green", "xkcd:navy blue"]
    sequence_height = 5
    sequence_index = 0
    at_top = True
    for sequence in sequence_graph["sequences"]:
        if not to_scale:
            unscaled_sequence = [unscaled_mapping[x] for x in sequence]
            sequence = unscaled_sequence

        make_exon_exon_lines(
            sequence,
            ax,
            y_exons,
            height=sequence_height,
            draw_at_top=at_top,
            color=colors[sequence_index],
        )
        if at_top:
            at_top = False
        else:
            at_top = True
        sequence_height += 5
        sequence_index += 1
        if sequence_index >= len(colors):
            sequence_index = 0

    start_margin = end_margin = 1000

    xmin = exons[0][0] - start_margin
    xmax = exons[len(exons) - 1][1] + end_margin
    xtick_interval = (xmax - xmin) / 10
    ax.set_xticks(np.arange(xmin, xmax, xtick_interval))
    ax.set_xbound(xmin, xmax)
    ax.set_ybound(0, 200)
    ax.add_collection(p)

    if title is not None:
        ax.set_title(title)

    if file_name is None:
        plt.show()
    else:
        plt.savefig(file_name)


if __name__ == "__main__":
    # For better resolution, maybe?
    # rcParams['figure.dpi'] = 200.0

    # Example taken from transcript ENST00000456328.2 of the gene DDX11L1
    draw_exons(
        [(11869, 12227), (12613, 12721), (13221, 14409)],
        file_name="out1.png",
        transcript_id="ENST00000456328.2",
    )

    # Example taken from transcript ENST00000450305.2 of the gene DDX11L1
    draw_exons(
        [
            (12010, 12057),
            (12179, 12227),
            (12613, 12619),
            (12975, 13052),
            (13221, 13374),
            (13453, 13670),
        ],
        file_name="out2.png",
        transcript_id="ENST00000450305.2",
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
                (12010, 12057),
                (12179, 12227),
                (12613, 12619),
                (12975, 13052),
                (13221, 13374),
                (13453, 13670),
            ],
            "sequences": [
                [(12010, 12057), (12179, 12227), (12613, 12619), (12975, 13052)],
                [
                    (12010, 12057),
                    (12613, 12619),
                    (12975, 13052),
                    (13221, 13374),
                    (13453, 13670),
                ],
            ],
        },
        file_name="out4.png",
        title="Contrived example using some exons from DDX11L1",
        to_scale=False,
    )
