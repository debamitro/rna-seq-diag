#!python3

import re
import sys


def read_gtf(file_name, query_gene_name):
    """Given a GTF file and a gene name to query for,
    return a dictionary where the keys are transcript IDs
    and the values are arrays of start and end offsets of the
    exons present in that transcript, e.g.,
     { 'tr1' : [[10,12],[17,27]] }"""

    def read_gtf_keyvalues(keyvaluestr):
        parts = keyvaluestr.split(";")
        for keyvalue in parts:
            m = re.match(r'\s*(\S+)\s*"(\S+)"', keyvalue)
            if m:
                yield (m.group(1), m.group(2))

    matching_transcripts = {}
    with open(file_name) as f:
        for line in f:
            parts = re.split(r"\s", line, maxsplit=8)
            if parts[2] == "exon":
                gene_name, transcript_id = "", ""
                for k, v in read_gtf_keyvalues(parts[8]):
                    if k == "gene_name":
                        gene_name = v
                    elif k == "transcript_id":
                        transcript_id = v
                if gene_name == query_gene_name:
                    exon = [int(parts[3]), int(parts[4])]
                    if transcript_id in matching_transcripts:
                        matching_transcripts[transcript_id].append(exon)
                    else:
                        matching_transcripts[transcript_id] = [exon]
    return matching_transcripts


if __name__ == "__main__":
    if len(sys.argv) > 2:
        transcripts = read_gtf(sys.argv[1], sys.argv[2])
        print(transcripts)
