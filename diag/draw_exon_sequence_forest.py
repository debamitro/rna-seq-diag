#!python3

import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection

if __name__ == "__main__":
    from exons import make_exon_shapes, make_exons_unscaled, make_exon_exon_lines
    from exons import configuration as exon_configuration
    from draw_exon_sequence_graph import configuration
else:
    from diag.exons import make_exon_shapes, make_exons_unscaled, make_exon_exon_lines
    from diag.exons import configuration as exon_configuration
    from diag.draw_exon_sequence_graph import configuration as graph_configuration


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
                    color=exon_configuration["exon_label_color"],
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
                color=graph_configuration["line_colors"][sequence_index],
            )
            sequence_height += 5
            sequence_index += 1
            if sequence_index >= len(graph_configuration["line_colors"]):
                sequence_index = 0
        y -= sequence_height * 2 + exon_configuration["exon_height"]

    p = PatchCollection(patches)

    xmin = xleft - graph_configuration["left_margin"]
    xmax = xright + graph_configuration["right_margin"]

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
