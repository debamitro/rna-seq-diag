#!python3

import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection

if __name__ == "__main__":
    from exons import make_exon_shapes, make_exon_exon_lines
else:
    from diag.exons import make_exon_shapes, make_exon_exon_lines


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
