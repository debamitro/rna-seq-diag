#!python3

import pprint

def analyze_sequences(transcripts):
    """Given a dictionary where the keys are transcript IDs
    and values are lists of (start,end) exon offsets, try to
    find out a set of 'decision trees' consisting of transcripts
    starting with the same sequence of exons, and diverging later on."""

    # Our representation of a forest of trees
    # 1. a mapping from parent to children
    child_exons = {}

    # 2. a set of root nodes
    root_exons = set()

    all_exons = set()
    for transcript_id in transcripts:
        parent_exon = None
        for exon in transcripts[transcript_id]:
            if parent_exon is not None:
                if parent_exon in child_exons:
                    child_exons[parent_exon].add(exon)
                else:
                    child_exons[parent_exon] = set([exon])
            else:
                root_exons.add(exon)
            parent_exon = exon
            all_exons.add(exon)

    sequence_forest = { "exons" : [],
                        "trees" : [] }
    for exon in all_exons:
        sequence_forest["exons"].append(exon)
    sequence_forest["exons"].sort()

    def sequences_for_tree(exon):
        """Given an exon, return a list of
        exon sequences that start with it. This function
        uses the tree generated in the parent function."""
        if exon in child_exons:
            seq = []
            for child_exon in child_exons[exon]:
                for child_sequence in sequences_for_tree(child_exon):
                    child_sequence.insert(0,exon)
                    seq.append(child_sequence)
            return seq
        else:
            return [[exon]]

    for root_exon in root_exons:
        sequence_forest["trees"].append(sequences_for_tree(root_exon))

    return sequence_forest

if __name__ == "__main__":
    transcripts = {
        "t1": [(20,25),(30,35),(50,55),(70,75)],
        "t2": [(20,25),(30,35),(40,45),(50,55),(70,75)]
    }
    sequences = analyze_sequences(transcripts)
    pp = pprint.PrettyPrinter()
    pp.pprint(sequences)
