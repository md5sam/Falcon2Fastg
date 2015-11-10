# Falcon2Fastg 

This software converts the results of a PacBio assembly using Falcon, to a FASTG graph that can be visualized using Bandage.

### Usage

    python Falcon2FastG.py

Run it in the same directory as a Falcon assembly. It needs the following input files:

* formatted_preads4falcon.fasta in SINGLE LINE fasta form (multi line preads4falcon.fasta is available in output directory of Falcon. This needs to be converted)

* sg_edges_list (available from output dir of Falcon)


### Dependencies :

seqtk (available at https://github.com/lh3/seqtk)


### Caveats : 

Faked a constant coverage of "50" because Bandage expects coverage information in the FASTG record headers.

Length record in FASTG header set to a constant "500". However, this does not seem to affect Bandage, which correctly calculates the read length



### Output : 

Output of Bandage visualization of the converted FASTG file from a Falcon assembly. Only the edges used in the string graph ("G" flagged in sg_edges_list) were used for this visualization.

![Alt text](/img/Falcon2Fastg_after_bandage.png?raw=true "Falcon2Fastg after Bandage")







## Other tools


Usage for utils/graph_of_overlaps.sh : 

./graph_of_overlaps.sh edges_list reads.fasta

Caveats :

Assumed that the sequence entry in the FASTG file represents overlap between two reads.

Faked a constant coverage of "50" because Bandage expects coverage information in the FASTG record headers.

Assumed the "length" record in FASTG header refers to length of overlap. 









