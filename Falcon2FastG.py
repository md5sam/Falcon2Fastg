import sys
import os
import csv
from itertools import groupby


list_of_tuples = []
all_nodes_list = []
all_source_nodes = []

sys.stdout = open('intmdt_NodeMode_out.fastg', 'w')

def create_read_pair_tuples () :
    #list_of_tuples = []
    with open("sg_edges_list") as sg_entries:
        read_pairs_TMI = csv.reader(sg_entries, delimiter=' ')
        for row in read_pairs_TMI :
	    #print (row[1])[:-2]
	    read_pair_tuple = ((row[0])[:-2],(row[1])[:-2]) 	
            #print type(read_pair_tuple)
            list_of_tuples.append(read_pair_tuple)
    list_of_tuples.sort()	    
    #print list_of_tuples[0]	



def collapse_ctg_list() :
    for key, sink_group in groupby(list_of_tuples, lambda x: x[0]):
	source_sink_list = []
	source_sink_list.append(key)
	for sink_node in sink_group:
	    #print (thing[1])
	    source_sink_list.append(sink_node[1])
	    
	all_nodes_list.append(source_sink_list)
    #print all_nodes_list	


def corresponding_FASTA(record_name): 
    FASTA_flag = "not_found"
    fp = open("formatted_preads4falcon.fasta")
    for line in fp : 
	if line[0]==">" :
            #print line[1:-1]
	    #print record_name
	    if line[1:-1] == record_name :
	        FASTA_flag = "found"
	else :
	    if FASTA_flag == "found":
	        sys.stdout.write(line)
		#print
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
	    #row_node_list.append(node)
    	    else : 
	        header.append("NODE_"+node+"_length_"+"500"+"_cov_50,")    
	for item in header :
	    sys.stdout.write(item)
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

	    


create_read_pair_tuples()
collapse_ctg_list()
make_header()
print_non_sources()

os.system("sed 's/,$/;/' intmdt_NodeMode_out.fastg > intmdt_NMout_SemiCol.fastg")

os.system("seqtk seq -l 80 intmdt_NMout_SemiCol.fastg > FINAL_nm.fastg")

os.system("rm intmdt*")


