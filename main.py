#!python3

import sys

from read_gtf import read_gtf
from draw_exons import draw_exons, draw_transcripts

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage is {0} <gtf-file-name> <gene-name>".format(sys.argv[0]))
        sys.exit (0)

    transcripts = read_gtf(sys.argv[1], sys.argv[2])
    draw_transcripts(transcripts)
