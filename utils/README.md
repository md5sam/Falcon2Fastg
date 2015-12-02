### Other tools

The ```utils/graph_of_overlaps.sh``` script is included for advanced users only, to visualize a different type o
f graph: one where the nodes correspond to overlaps between reads (in the ```Falcon2Fastg.py```, overlaps are re
presented by edges, not nodes).

Usage:

    ./graph_of_overlaps.sh edges_list reads.fasta

Caveats :

Assumed that sequence entry in FASTG file represents overlap between two reads.

Faked a constant coverage of "50" because Bandage expects coverage information 
in the FASTG record headers.

Assumed the "length" record in FASTG header refers to length of overlap. 






















