#!python3

import matplotlib.pyplot as plt
import matplotlib.lines as lines
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib import rcParams
import numpy as np

configuration = {
    "unscaled_exon_width": 1000,
    "exon_height": 20,
    "unscaled_exon_start": 2000,
    "left_margin": 1000,
    "right_margin": 1000,
    "exon_color": "xkcd:mustard",
    "exon_label_color": "xkcd:white",
    "line_colors": ["xkcd:indigo", "xkcd:forest green", "xkcd:navy blue"],
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


def draw_exons(exons, file_name=None, transcript_id=None):
    """Given an array of (start,end) offsets for exons,
    draws them in a diagram. Optionally writes out the diagram
    as a file, and also adds the transcript id"""
    y = 130
    if transcript_id is not None:
        plt.text(exons[0][0], y + 20, transcript_id)
    draw_exon_sequence_graph(
        {"id": "gr1", "exons": exons, "sequences": [exons]}, y, file_name, transcript_id
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
        exons = transcripts[transcript_id]["exons"]
        patches.extend(make_exon_shapes(exons, y))
        exon_pairs = zip(exons, exons[1:])
        make_exon_exon_lines(exon_pairs, ax, y)
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
        unscaled_mapping, unscaled_exons = make_exons_unscaled(exons)
        exons = unscaled_exons

    patches = make_exon_shapes(exons, y_exons)
    p = PatchCollection(patches)

    sequence_height = 5
    sequence_index = 0
    draw_position = ["mid", "top", "bottom"]
    for sequence in sequence_graph["sequences"]:
        if not to_scale:
            unscaled_sequence = [unscaled_mapping[x] for x in sequence]
            sequence = unscaled_sequence

        exon_pairs = zip(sequence, sequence[1:])
        make_exon_exon_lines(
            exon_pairs,
            ax,
            y_exons,
            height=sequence_height,
            draw_at=draw_position[sequence_index],
            color=configuration["line_colors"][sequence_index],
        )
        sequence_height += 5
        sequence_index += 1
        if sequence_index >= len(configuration["line_colors"]):
            sequence_index = 0

    xmin = exons[0][0] - configuration["left_margin"]
    xmax = exons[len(exons) - 1][1] + configuration["right_margin"]

    if to_scale:
        xtick_interval = (xmax - xmin) / 10
        ax.set_xticks(np.arange(xmin, xmax, xtick_interval))
    else:
        ax.set_xticks([])

    ax.set_yticks([y_exons])
    if "id" in sequence_graph:
        ax.set_yticklabels([sequence_graph["id"]])

    ax.set_xbound(xmin, xmax)
    ax.set_ybound(0, 200)
    ax.add_collection(p)

    if title is not None:
        ax.set_title(title)

    if file_name is None:
        plt.show()
    else:
        plt.savefig(file_name)


def draw_exon_sequence_forest(forest, **kwargs):
    """Given a 'forest', i.e. a collection of decision trees,
    draw them in the same plot one row at a time."""
    _, ax = plt.subplots()

    ymax = len(forest["trees"]) * 40 + 20
    y = ymax
    patches = []
    xleft, xright = None, None
    yticks = []

    # Set some default values
    if "add_exon_labels" not in kwargs:
        kwargs["add_exon_labels"] = False
    if "merge_common_sequences" not in kwargs:
        kwargs["merge_common_sequences"] = False

    if kwargs["add_exon_labels"]:
        exon_labels = {}

    for tree in forest["trees"]:
        yticks.append(y)
        exons_from_tree = set()
        for sequence in tree:
            for exon in sequence:
                exons_from_tree.add(exon)
        exons = []
        for exon in exons_from_tree:
            exons.append(exon)
        exons.sort()

        unscaled_mapping, unscaled_exons = make_exons_unscaled(exons)
        if xleft is None or unscaled_exons[0][0] < xleft:
            xleft = unscaled_exons[0][0]
        if xright is None or unscaled_exons[len(unscaled_exons) - 1][1] > xright:
            xright = unscaled_exons[len(unscaled_exons) - 1][1]

        patches.extend(make_exon_shapes(unscaled_exons, y))

        if kwargs["add_exon_labels"]:
            for exon in exons:
                if exon not in exon_labels:
                    exon_labels[exon] = "e{0}".format(len(exon_labels) + 1)
                plt.text(
                    unscaled_mapping[exon][0] + 5,
                    y + 5,
                    exon_labels[exon],
                    color=configuration["exon_label_color"],
                    fontweight="bold",
                )

        sequence_height = 5
        sequence_index = 0
        draw_position = ["mid", "top", "bottom"]

        if kwargs["merge_common_sequences"]:
            unique_exon_pairs = set()
        for sequence in tree:
            unscaled_sequence = [unscaled_mapping[x] for x in sequence]

            if kwargs["merge_common_sequences"]:
                exon_pairs = [
                    p
                    for p in zip(unscaled_sequence, unscaled_sequence[1:])
                    if p not in unique_exon_pairs
                ]
                for p in exon_pairs:
                    unique_exon_pairs.add(p)
            else:
                exon_pairs = zip(unscaled_sequence, unscaled_sequence[1:])
            make_exon_exon_lines(
                exon_pairs,
                ax,
                y,
                height=sequence_height,
                draw_at=draw_position[sequence_index],
                color=configuration["line_colors"][sequence_index],
            )
            sequence_height += 5
            sequence_index += 1
            if sequence_index >= len(configuration["line_colors"]):
                sequence_index = 0
        y -= sequence_height * 2 + configuration["exon_height"]

    p = PatchCollection(patches)

    xmin = xleft - configuration["left_margin"]
    xmax = xright + configuration["right_margin"]

    ax.set_xticks([])

    ax.set_yticks(yticks)

    ax.set_yticklabels(["d{0}".format(y + 1) for y in range(len(forest["trees"]))])

    ax.set_xbound(xmin, xmax)
    ax.set_ybound(0, ymax + 50)
    ax.add_collection(p)

    if "title" in kwargs:
        ax.set_title(kwargs["title"])

    if "file_name" in kwargs:
        plt.savefig(kwargs["file_name"])
    else:
        plt.show()


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
            "ENST00000456328.2": {
                "exons": [[11869, 12227], [12613, 12721], [13221, 14409]],
            },
            "ENST00000450305.2": {
                "exons": [
                    [12010, 12057],
                    [12179, 12227],
                    [12613, 12619],
                    [12975, 13052],
                    [13221, 13374],
                    [13453, 13670],
                ],
            },
        }
    )

    # Contrived example using some exons from the above
    draw_exon_sequence_graph(
        {
            "id": "gr1",
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

    # Completely contrived example
    draw_exon_sequence_forest(
        {
            "trees": [
                [
                    [(1010, 1015), (1025, 1030), (1045, 1050), (1060, 1065)],
                    [
                        (1010, 1015),
                        (1025, 1030),
                        (1060, 1065),
                        (1070, 1075),
                        (1080, 1085),
                    ],
                ],
                [
                    [
                        (1010, 1015),
                        (1025, 1030),
                        (1045, 1050),
                        (1055, 1057),
                        (1060, 1065),
                    ],
                    [
                        (1010, 1015),
                        (1025, 1030),
                        (1060, 1065),
                        (1070, 1075),
                        (1080, 1085),
                    ],
                    [
                        (1010, 1015),
                        (1025, 1030),
                        (1060, 1065),
                        (1070, 1075),
                        (1090, 1095),
                        (1100, 1105),
                    ],
                ],
            ]
        },
        title="Contrived decision forest",
        add_exon_labels=True,
        merge_common_sequences=True,
    )
