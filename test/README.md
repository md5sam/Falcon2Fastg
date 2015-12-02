# Test directory with example data

Two small example datasets are provided here to test Falcon2Fastg.py for your system.


### Data

- The first test dataset is in one_contig_test 

The test data (3nodes.fasta) consists of 3 overlapping reads from chr14 of the hg38 build of Human Reference. Re
ads are of 5000, 4000 and 3000bp in length. 

The suffix of 5000bp read overlaps with prefix of 4000bp read, and suffix of 4000bp read overlaps with prefix of
 3000 bp read. 

Each overlap is 1000 bp in length (100% ID).



- The second test dataset is in two_contigs_test/

The test data (circ_and_line.fasta) consists of 7 overlapping reads from chr14 of the hg38 build of Human Reference. Each overlap is 1000 bp in length. 

The first 4 reads form a circular ctg. and the last 3 reads form a linear ctg. 

The last 1000 bp of read 0 overlaps (100% ID) with the first 1000 bp of read 4.



### Directories

Each test dataset has the following folders : 

1. data/ contains the reads files described above

2. falcon_rundir/ contains the expected output from running Falcon on this toy dataset, with the assembly parameters specified in the .cfg file

3. falcon2fastg_output/ contains the expected output of Falcon2Fastg.py 

4. visualization/ contains two screenshots of viewing this data with Bandage


The figure below represents Bandage visualization of reads.fastg from two_contigs_test/

![Alt text](/test/two_contigs_test/visualization/circ_lin_reads.png?raw=true "Example reads after Bandage")


The figure below represents Bandage visualization of contigs.fastg from two_contigs_test/

![Alt text](/test/two_contigs_test/visualization/circ_lin_ctgs.png?raw=true "Example reads after Bandage")



