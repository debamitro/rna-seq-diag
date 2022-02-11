# rna-seq-diag

Some scripts for drawing RNA sequencing diagrams. Completely under construction.

## Dependencies

* Python 3
* Matplotlib
* Numpy

## Development dependencies

* Black (for formatting)

## How to run

The top-level script is main.py. It takes two inputs - a GTF file and the name of a gene.

## Project structure

### Read_gtf

The 'read_gtf' module provides a function for reading a GTF file and extracting
all the transcripts corresponding to a gene as a python dictionary. The keys of the dictionary
are the transcript IDs, and the values are lists of (start, end) tuples for each exon which
is part of the transcript. Something like this:

```python
{
  "tr1" : [(10,30),(50,60),(70,100)],
  "tr2": [(10,30),(40,45),(50,60),(70,100)]
}
```

### Analyze_sequences

The 'analyze_sequences' module provides some experimental code which returns a 'forest'
of 'decision trees' for the exon sequences present in some transcripts. The 'forest'
is a python dictionary with two keys -

* "exons" where the value is an ordered list of exons found in all the transcripts
* "trees" where the value is a list of list of list of exons (each tree is represented as a
collection of transcripts which use that decision tree)

A rough example of the returned 'forest' would be

```python
{
  "exons" : [(10,30),(40,45),(50,60),(70,100)],
  "trees": [
             [
               [(10,30),(50,60),(70,100)],
               [(10,30),(40,45),(50,60),(70,100)]
             ]
           ]
}
```

### Diag and draw_exons

The 'diag' package uses numpy and matplotlib to provide different functions around drawing exon sequences. There is also  the 'draw_exons' module which continues to exist, for historical reasons.

* Drawing one sequence of exons

```python
# Part of draw_exons
def draw_exons(exons, file_name=None, transcript_id=None):
    """Given an array of (start,end) offsets for exons,
    draws them in a diagram. Optionally writes out the diagram
    as a file, and also adds the transcript id"""
```

* Drawing multiple sequences of exons

```python
# Part of diag.draw_transcripts
def draw_transcripts(transcripts, file_name=None):
    """Given a dictionary where the keys are transcript IDs
    and the values are arrays of exon start and end offsets,
    draws them in a diagram. Optionally saves out the diagram
    in a file."""
```

* Drawing an arbitrary directed graph given a set of exons and a set of edges

```python
# Part of diag.draw_exon_sequence_graph
def draw_exon_sequence_graph(
    sequence_graph, y_exons=130, file_name=None, title=None, to_scale=True
):
    """Given a dictionary with two entries
     - 'exons' an array of exon start and end offsets
     - 'sequences' an array of exon sequences
    draws a graph using different colors for each sequence.
    The goal is to show different exon sequences formed from
    one universal set of exons"""
```

* Drawing a set of directed graphs, each of which tries to show the possible sequences formed by a bunch of exons

```python
# Part of diag.draw_exon_sequence_forest
def draw_exon_sequence_forest(forest, **kwargs):
    """Given a 'forest', i.e. a collection of decision trees,
    draw them in the same plot one row at a time."""
```

### main.py

This is where all the different parts are going to converge some day. As of now it does two things

* Draw all the transcripts for a given gene and annotation file
* Club the transcripts into a bunch of 'decision trees' which start with the same exons, and draw them logically
