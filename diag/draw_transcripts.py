#!python3

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection

if __name__ == "__main__":
    from exons import make_exon_shapes, make_exon_exon_lines
    from draw_exon_sequence_graph import configuration as graph_configuration
else:
    from diag.exons import make_exon_shapes, make_exon_exon_lines
    from diag.draw_exon_sequence_graph import configuration as graph_configuration


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
    yticks = []
    for transcript_id in transcripts:
        y -= 40
        yticks.append(y)
        exons = transcripts[transcript_id]["exons"]
        patches.extend(make_exon_shapes(exons, y))
        exon_pairs = zip(exons, exons[1:])
        make_exon_exon_lines(exon_pairs, ax, y)
        if xleft is None or exons[0][0] < xleft:
            xleft = exons[0][0]
        if xright is None or exons[len(exons) - 1][1] > xright:
            xright = exons[len(exons) - 1][1]

    p = PatchCollection(patches)

    xmin = xleft - graph_configuration["left_margin"]
    xmax = xright + graph_configuration["right_margin"]

    # Let's put 10 ticks on the x axis
    xtick_interval = (xmax - xmin) / 10
    xticks = np.arange(xmin, xmax, xtick_interval)
    ax.set_xticks(xticks, ["{0}k bp".format(int(xt/1000)) for xt in xticks])

    # Let's put one tick per transcript on the y axis
    ax.set_yticks(yticks, [transcript_id for transcript_id in transcripts])

    ax.set_xbound(xmin, xmax)
    ax.set_ybound(0, ymax)
    ax.add_collection(p)

    if file_name is None:
        plt.show()
    else:
        plt.savefig(file_name)


if __name__ == "__main__":
    # Combining the two transcripts ENST00000456328.2 and ENST00000450305.2
    draw_transcripts(
        {
            "ENST00000456328.2": {
                "exons": [(11869, 12227), (12613, 12721), (13221, 14409)],
            },
            "ENST00000450305.2": {
                "exons": [
                    (12010, 12057),
                    (12179, 12227),
                    (12613, 12619),
                    (12975, 13052),
                    (13221, 13374),
                    (13453, 13670),
                ],
            },
        }
    )
