### This script is an earlier version of Falcon2Fastg before porting to Python.
### It represents the overlaps between reads as nodes. Hence, this is a different representation.

### Output is the SPAdes flavor of FASTG representation, for visualization by Bandage
### Dependencies : seqtk, samtools

### CAUTION : removes existing fastg and intermediate file in directory


if [ $# -ne 2 ]
then
    echo "Usage: $0 EDGES_LIST REF_READS.FASTA"
    exit 1
fi

# Create FASTG-like header from sg_edges_list found in Falcon dir 2-asm-falcon
# Sample header output looks like :
# >Node1_NDID_length_OVLPLEN_cov_OVLPCOV:Node2_NDID2_length_OVLPLEN_cov_OVLPCOV
awk 'BEGIN{OFS=""} {print ">NODE","_",$1,"_","length","_",($4-$5)<0?($4-$5)*-1:($4-$5),"_","cov","_","50",":","NODE","_",$2,"_","length","_",($4-$5)<0?($4-$5)*-1:($4-$5),"_","cov","_","50"}' $1 > part1_headers

# Save the read headers in the format R1:Ovlpstartpos-Ovlpendpos
awk 'BEGIN{OFS=""}{print $3,":",$5==0?$5:$4,"-",$5==0?$4:$5}' $1 > ctg_seqs_to_extract

# Extract the sequence from preads4falcon.fasta
for name in `cat ctg_seqs_to_extract`; do samtools faidx $2 $name; done > part2_extracted_seqs.fasta

# Convert to single line fasta (as in https://www.biostars.org/p/9262/#9264)
# Remove the header line
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < part2_extracted_seqs.fasta | tail -n +2 | grep -v ">" > part3_headless_extracted_seqs.fasta

# Paste (alternately) the headers and the extracted seqs together
paste -d"\n" part1_headers part3_headless_extracted_seqs.fasta > single_line.fastg

# Fix the formatting back to multiline, constant length mode 
seqtk seq -l 80 single_line.fastg > multi_line.fastg

# Remove :E and :B (the 3' and 5' end labels used by FALCON) 
# Bandage might complain if you do not remove these
sed 's/:E//g' multi_line.fastg | sed 's/:B//g' > FINAL.fastg

rm multi_line.fastg
rm single_line.fastg


















