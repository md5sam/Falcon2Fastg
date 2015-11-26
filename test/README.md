# Test directory with example data

This small example dataset can be used to test Falcon2Fastg.py for your system.


### Data

The test data (circ_and_line.fasta) consists of 7 overlapping reads from chr14 of the hg38 build of Human Reference. Each overlap is 1000 bp in length. 

The first 4 reads form a circular ctg. and the last 3 reads form a linear ctg. 

The last 1000 bp of read 0 overlaps (100% ID) with the first 1000 bp of read 4.



### Directories

1. data/ contains the reads file circ_and_line.fasta described above

2. falcon_rundir/ contains the expected output from running Falcon on this toy dataset, with the assembly parameters specified in the .cfg file

3. falcon2fastg_output/ contains the expected output of Falcon2Fastg.py 

4. visualization/ contains two screenshots of viewing this data with Bandage


The figure below represents Bandage visualization of reads.fastg

![Alt text](/test/visualization/circ_lin_reads.png?raw=true "Example reads after Bandage")


The figure below represents Bandage visualization of contigs.fastg

![Alt text](/test/visualization/circ_lin_ctgs.png?raw=true "Example reads after Bandage")



