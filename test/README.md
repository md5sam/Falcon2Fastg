# Test directory with example data

This small example dataset can be used to test Falcon2Fastg.py for your system.


### Data

The test data (3nodes.fasta) consists of 3 overlapping reads from chr14 of the hg38 build of Human Reference. Reads are of 5000, 4000 and 3000bp in length. 

The suffix of 5000bp read overlaps with prefix of 4000bp read, and suffix of 4000bp read overlaps with prefix of 3000 bp read. 

Each overlap is 1000 bp in length. 


### Directories

1. data/ contains the reads file 3nodes.fasta described above

2. falcon_rundir/ contains the expected output from running Falcon on this toy dataset, with the assembly parameters specified in the .cfg file

3. falcon2fastg_output/ contains the expected output of Falcon2Fastg.py after running it on the output of falcon (sg_edges_list and preads4falcon.fasta)
Optionally, this directory also contains a reference.fasta which is the perfect assembly that can be obtained by overlapping the input reads.

4. visualization/ contains a screenshot of viewing this data with Bandage

![Alt text](/visualization/3nodes.png?raw=true "Example data after Bandage")




