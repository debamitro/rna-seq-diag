#!python3

import sys

from read_gtf import read_gtf
from draw_exons import draw_transcripts, draw_exon_sequence_forest
from analyze_sequences import analyze_sequences

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage is {0} <gtf-file-name> <gene-name>".format(sys.argv[0]))
        sys.exit(0)

    transcripts = read_gtf(sys.argv[1], sys.argv[2])
    draw_transcripts(transcripts)
    forest = analyze_sequences(transcripts)
    draw_exon_sequence_forest(forest)
