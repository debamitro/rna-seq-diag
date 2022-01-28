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

### Draw_exons

The 'draw_exons' module provides different functions for

* Drawing one sequence of exons
* Drawing multiple sequences of exons
* Drawing an arbitrary directed graph given a set of exons and a set of edges

This module uses matplotlib and numpy

### main.py

This is where all the different parts are going to converge some day
