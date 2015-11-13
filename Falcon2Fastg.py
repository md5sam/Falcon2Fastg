import sys
import csv
from itertools import groupby
from Bio import SeqIO
from Bio.SeqIO import FastaIO

list_of_tuples = []
all_nodes_list = []
all_source_nodes = []

sys.stdout = open('output.fastg', 'w')


def convert_multiline_to_single_line_FASTA () :
    sequences = [] 
    input_handle = open("preads4falcon.fasta", "rU")

    for record in SeqIO.parse(input_handle, "fasta"):
	sequences.append(record)
    
    output_handle = open("formatted_preads4falcon.fasta","w")
    fasta_out = FastaIO.FastaWriter(output_handle, wrap=None)
    fasta_out.write_file(sequences)
    output_handle.close()


def create_read_pair_tuples () :
    with open("sg_edges_list") as sg_entries:
        read_pairs_TMI = csv.reader(sg_entries, delimiter=' ')
        for row in read_pairs_TMI :
	    read_pair_tuple = ((row[0])[:-2],(row[1])[:-2]) 	
            list_of_tuples.append(read_pair_tuple)
    list_of_tuples.sort()	    


def collapse_ctg_list() :
    for key, sink_group in groupby(list_of_tuples, lambda x: x[0]):
	source_sink_list = []
	source_sink_list.append(key)
	for sink_node in sink_group:
	    source_sink_list.append(sink_node[1])
	    
	all_nodes_list.append(source_sink_list)


def corresponding_FASTA(record_name): 
    FASTA_flag = "not_found"
    fp = open("formatted_preads4falcon.fasta")
    for line in fp : 
	if line[0]==">" :
	    if line[1:-1] == record_name :
	        FASTA_flag = "found"
	else :
	    if FASTA_flag == "found":
	        sys.stdout.write(line)
		fp.close()
                return 


def make_header() :
    for source_to_sinks in all_nodes_list : 
	header = []
	row_node_list = []
	flag = "source"
	for node in source_to_sinks :
	    if flag == "source" :
		source = node
		all_source_nodes.append(node)
		header.append(">NODE_"+node+"_length_"+"500"+"_cov_50:")
		flag = "sink"
    	    else : 
	        header.append("NODE_"+node+"_length_"+"500"+"_cov_50,")    
	str_header = ''.join(header)
	sys.stdout.write(str_header[:-1]+';')
	print 
	corresponding_FASTA(source)


def print_non_sources() :
    fp_non = open("formatted_preads4falcon.fasta")
    NON_flag = "non_source"
    for line in fp_non :
	    if line[0] == ">":
		if line[1:-1] not in all_source_nodes :
		    sys.stdout.write(">NODE_"+line[1:-1]+"_length_"+"500"+"_cov_50,")
		    print
		    NON_flag = "non_source"
	        else :  
		    NON_flag = "source"	
	    else :
		if NON_flag == "non_source" : 
		    print line[:-1] 
    fp_non.close()     



convert_multiline_to_single_line_FASTA ()
create_read_pair_tuples()
collapse_ctg_list()
make_header()
print_non_sources()



