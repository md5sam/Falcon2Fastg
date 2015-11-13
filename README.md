# Falcon2Fastg 

This software converts the results of a PacBio assembly using FALCON, to a FASTG graph that can be visualized using Bandage.

### Usage

    python Falcon2Fastg.py

Run it in the output directory of FALCON assembly (2-asm-falcon). It needs the following input files:

* preads4falcon.fasta

* sg_edges_list 


### Dependencies :

Biopython (available at http://biopython.org/wiki/Download)


### Output : 

The output of the tool is a FASTG file that can be opened with Bandage (https://github.com/rrwick/Bandage).


Below is a sample Bandage visualization of a FASTG file generated by Falcon2Fastg from a FALCON assembly (mitochondrial genome).

* Each node is a read (colors are random).
* Edges represent the overlaps between reads found by FALCON.
* Only the edges used in the string graph ("G" flagged in sg_edges_list) were used for this visualization.


![Alt text](/img/Falcon2Fastg_after_bandage.png?raw=true "Falcon2Fastg after Bandage")


The figure below represents a visualization of the E. coli test dataset distributed with the Falcon installation. In this case, all edges are visualized, including those involved in Repeats and edges which are Transitive Reducible.


![Alt text](/img/ecoli_Allnodes.png?raw=true "Ecoli all edges fastg after Bandage")



The figure below uses the same E. coli dataset, however the Repeat and Transitive Reducible edges are removed to give a single contiguous graph (along with multiple single edges), which Falcon outputs as a single contig. 


![Alt text](/img/ecoli_Gnodes.png?raw=true "Ecoli 'G' edges fastg after Bandage")









### Caveats : 

Each overlap between a pair of reads is currently reported as two edges.   

Faked a constant coverage of "50" because Bandage expects coverage information in the FASTG record headers.

Length record in FASTG header set to a constant "500". However, this does not seem to affect Bandage, which correctly calculates the read length




## Other tools


Usage for utils/graph_of_overlaps.sh : 

./graph_of_overlaps.sh edges_list reads.fasta

Caveats :

Assumed that the sequence entry in the FASTG file represents overlap between two reads.

Faked a constant coverage of "50" because Bandage expects coverage information in the FASTG record headers.

Assumed the "length" record in FASTG header refers to length of overlap. 









