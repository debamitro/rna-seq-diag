#!python3

import sys

from read_gtf import read_gtf
from draw_exons import draw_exons, draw_transcripts, draw_exon_sequence_graph
from analyze_sequences import analyze_sequences

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage is {0} <gtf-file-name> <gene-name>".format(sys.argv[0]))
        sys.exit(0)

    transcripts = read_gtf(sys.argv[1], sys.argv[2])
    draw_transcripts(transcripts)
    forest = analyze_sequences(transcripts)
    i = 1
    for tree in forest["trees"]:
        sequences = {"exons": [], "sequences": tree}
        all_exons = set()
        for sequence in tree:
            for exon in sequence:
                all_exons.add(exon)
        for exon in all_exons:
            sequences["exons"].append(exon)
        sequences["exons"].sort()

        draw_exon_sequence_graph(sequences, title="Tree {0}".format(i), to_scale=False)
        i += 1
