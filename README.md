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



Falcon2Fastg was tested on the preads4falcon.fasta and sg_edges_list file produced by the E.coli test dataset provided with the Falcon install. Instructions on obtaining the dataset are here : https://github.com/PacificBiosciences/FALCON/wiki/Setup:-Complete-example  

The test was performed on a desktop system with Intel Xeon CPU W3520 @ 2.67GHz and 4 cores. Runtime for this dataset (449 Mb) of preads4falcon.fasta is :

    time python Falcon2Fastg.py

    real        4m18.640s
    user        2m16.457s
    sys         1m59.268s


The figure below represents a visualization of this E. coli data. In this case, all edges are visualized, including those involved in Repeats and edges which are Transitive Reducible.



![Alt text](/img/ecoli_Allnodes.png?raw=true "Ecoli all edges fastg after Bandage")





The figure below uses the same E. coli dataset, however the Repeat and Transitive Reducible edges are removed to give a single contiguous graph (along with many isolated nodes). Note the reduction in number of edges, from 10,552 in the previous figure to 1,457 here.


![Alt text](/img/ecoli_Gnodes.png?raw=true "Ecoli 'G' edges fastg after Bandage")









### Caveats : 

Faked a constant coverage of "50" because Bandage expects coverage information in the FASTG record headers.

Length record in FASTG header set to a constant "500". However, this does not seem to affect Bandage, which correctly calculates the read length




## Other tools


Usage for utils/graph_of_overlaps.sh : 

./graph_of_overlaps.sh edges_list reads.fasta

Caveats :

Assumed that the sequence entry in the FASTG file represents overlap between two reads.

Faked a constant coverage of "50" because Bandage expects coverage information in the FASTG record headers.

Assumed the "length" record in FASTG header refers to length of overlap. 









