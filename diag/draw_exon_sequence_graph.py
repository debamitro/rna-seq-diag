#!python3

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection

if __name__ == "__main__":
    from exons import make_exon_shapes, make_exons_unscaled, make_exon_exon_lines
else:
    from diag.exons import make_exon_shapes, make_exons_unscaled, make_exon_exon_lines

configuration = {
    "left_margin": 1000,
    "right_margin": 1000,
    "line_colors": ["xkcd:indigo", "xkcd:forest green", "xkcd:navy blue"],
}


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


if __name__ == "__main__":
    # Contrived example using some exons from DDX11L1
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
