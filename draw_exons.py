#!python3

import matplotlib.pyplot as plt
from matplotlib import rcParams
from diag.draw_exon_sequence_graph import draw_exon_sequence_graph


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
