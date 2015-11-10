Usage : python Falcon2FastG.py


Required files : 

formatted_preads4falcon.fasta in SINGLE LINE fasta form (multi line preads4falcon.fasta is available in output directory of FALCON. This needs to be converted)

sg_edges_list (available from output dir of FALCON)


Tool requirements :

seqtk (available at https://github.com/lh3/seqtk)



Caveats : 

Faked a constant coverage of "50" because Bandage expects coverage information in the FASTG record headers.

Length record in FASTG header set to a constant "500". However, this does not seem to affect Bandage, which correctly calculates the read length














Usage for deprecated_overlap_converter.sh : 

./deprecated_overlap_converter.sh edges_list reads.fasta

Caveats :

Assumed that the sequence entry in the FASTG file represents overlap between two reads.

Faked a constant coverage of "50" because Bandage expects coverage information in the FASTG record headers.

Assumed the "length" record in FASTG header refers to length of overlap. 









