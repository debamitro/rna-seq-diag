#!python3

import pprint


def analyze_sequences(transcripts):
    """Given a dictionary where the keys are transcript IDs
    and values are lists of (start,end) exon offsets, try to
    find out a set of 'decision trees' consisting of transcripts
    starting with the same sequence of exons, and diverging later on."""

    # Our representation of a forest of trees
    # 1. A mapping from parent tree node to children tree nodes
    #    This is implemented as a two-level dictionary where the keys are node IDs (integers)
    #    and the values are dictionaries themselves keyed by exons. The values of the
    #    values are node IDs
    child_nodes = {}

    # 2. A mapping of exons to root node IDs
    root_exons = {}

    # 3. A mapping of node IDs to exons
    node_ids = {}

    all_exons = set()

    def new_tree_node(exon):
        next_node_id = len(node_ids)
        node_ids[next_node_id] = exon
        return next_node_id

    for transcript_id in transcripts:
        parent_node = None
        for exon in transcripts[transcript_id]:
            if parent_node is not None:
                if parent_node in child_nodes:
                    if exon not in child_nodes[parent_node]:
                        child_nodes[parent_node][exon] = new_tree_node(exon)
                else:
                    child_nodes[parent_node] = {exon: new_tree_node(exon)}
                parent_node = child_nodes[parent_node][exon]
            else:
                if exon not in root_exons:
                    root_exons[exon] = new_tree_node(exon)
                parent_node = root_exons[exon]
            all_exons.add(exon)

    sequence_forest = {"exons": [], "trees": []}
    for exon in all_exons:
        sequence_forest["exons"].append(exon)
    sequence_forest["exons"].sort()

    def sequences_for_tree(node):
        """Given an exon, return a list of
        exon sequences that start with it. This function
        uses the tree generated in the parent function."""
        exon = node_ids[node]
        if node in child_nodes:
            seq = []
            for child_node in child_nodes[node]:
                for child_sequence in sequences_for_tree(child_nodes[node][child_node]):
                    child_sequence.insert(0, exon)
                    seq.append(child_sequence)
            return seq
        else:
            return [[exon]]

    for root_exon in root_exons:
        sequence_forest["trees"].append(sequences_for_tree(root_exons[root_exon]))

    return sequence_forest


if __name__ == "__main__":
    transcripts = {
        "t1": [(20, 25), (30, 35), (50, 55), (70, 75)],
        "t2": [
            (20, 25),
            (30, 35),
            (50, 55),
            (80, 85),
            (90, 95),
            (100, 105),
            (120, 125),
            (130, 135),
            (150, 155),
        ],
    }
    sequences = analyze_sequences(transcripts)
    pp = pprint.PrettyPrinter()
    pp.pprint(sequences)
