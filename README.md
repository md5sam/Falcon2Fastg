Usage : ./Falcon2FastG.sh edges_list reads.fasta

Caveats :

Assumed that the sequence entry in the FASTG file represents overlap between two reads.

Faked a constant coverage of "50" because Bandage expects coverage information in the FASTG record headers.

Assumed the "length" record in FASTG header refers to length of overlap. 









