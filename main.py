#!python3

import sys

from read_gtf import read_gtf
from draw_exons import draw_exons

if __name__ == "__main__":
    transcripts = read_gtf(sys.argv[1], "DDX11L1")
    for transcript in transcripts:
        draw_exons(transcripts[transcript])
