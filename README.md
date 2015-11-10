Usage : python Falcon2FastG.py


Requires : 

preads4falcon.fasta in single_line FASTA form (available in output directory of FALCON)

sg_edges_list (available from output dir of FALCON)

seqtk (available at https://github.com/lh3/seqtk)



Usage for deprecated_overlap_converter.sh : 

./deprecated_overlap_converter.sh edges_list reads.fasta

Caveats :

Assumed that the sequence entry in the FASTG file represents overlap between two reads.

Faked a constant coverage of "50" because Bandage expects coverage information in the FASTG record headers.

Assumed the "length" record in FASTG header refers to length of overlap. 









