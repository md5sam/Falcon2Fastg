### Other tools

The ```utils/graph_of_overlaps.sh``` script is included for advanced users only, to visualize a different Falcon graph type: nodes now correspond to overlaps between reads. In the original ```Falcon2Fastg.py``` program, overlaps are represented by edges, not nodes.

Usage:

    ./graph_of_overlaps.sh edges_list reads.fasta

Caveats :

Assumed that sequence entry in FASTG file represents overlap between two reads.

Faked a constant coverage of "50" because Bandage expects coverage information 
in the FASTG record headers.

Assumed the "length" record in FASTG header refers to length of overlap. 






















