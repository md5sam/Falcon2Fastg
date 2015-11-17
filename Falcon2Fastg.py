import sys
import csv
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqIO import FastaIO

list_of_tuples = []
all_nodes_list = []
all_source_nodes = []
direction = defaultdict(list)


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
	    read_pair_tuple = ((row[0]),(row[1])) 	
            list_of_tuples.append(read_pair_tuple)
    list_of_tuples.sort()	    


def make_direction_dict () :
    for pair in list_of_tuples :
	if pair[0][-1] == "E" and pair[1][-1] == "E" :
	    direction[pair[0][:-2]].append(pair[1][:-2])
            
	elif pair[0][-1] == "E" and pair[1][-1] == "B" :
	    direction[pair[0][:-2]].append(pair[1][:-2]+"'")
                        
	elif pair[0][-1] == "B" and pair[1][-1] == "E" :
            direction[pair[0][:-2]+"'"].append(pair[1][:-2])

	elif pair[0][-1] == "B" and pair[1][-1] == "B" :
            direction[pair[0][:-2]+"'"].append(pair[1][:-2]+"'")
	    
        else : 
	    print "ERROR : the tuples have not been parsed correctly"    


def corresponding_FASTA(record_name): 
    if record_name[-1] == "'" :
        corrected_name = record_name[:-1]
    else :
	corrected_name = record_name
    FASTA_flag = "not_found"
    fp = open("formatted_preads4falcon.fasta")
    for line in fp : 
	if line[0]==">" :
	    if line[1:-1] == corrected_name :
	        FASTA_flag = "found"
	else :
	    if FASTA_flag == "found":
		if record_name[-1] == "'" :    
	            from Bio.Seq import Seq
		    line_revcomp = Seq(line[:-1]).reverse_complement()
		    sys.stdout.write(str(line_revcomp))
		    print
	        else : 
		    sys.stdout.write(line)
		    placeholder = 0
		fp.close()
                return 


def headerify (node) :
    if node[-1] == "'" :
	return str("NODE_"+node[:-1]+"_length_"+"500"+"_cov_50'") 
    else : 
	return str("NODE_"+node+"_length_"+"500"+"_cov_50")    


def make_fastg() :
    for source, sinks in direction.items() : 
        all_source_nodes.append(source)	    
	header = []
	header.append(">"+headerify(source)+":")
	for element in sinks :
	    header.append(headerify(element)+",")
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
		sys.stdout.write(">NODE_"+line[1:-1]+"_length_"+"500"+"_cov_50;")
		print
		NON_flag = "non_source"
	    else :  
		NON_flag = "source"	    	
	else :
            if NON_flag == "non_source" : 
		print line[:-1]
		
    fp_non.close()     


def print_non_sources_complement() :
    fp_non = open("formatted_preads4falcon.fasta")
    NON_flag = "non_source" 
    for line in fp_non :
	if line[0] == ">":		
	    if str (line[1:-1]+"'") not in all_source_nodes :
	        sys.stdout.write(">NODE_"+line[1:-1]+"_length_"+"500"+"_cov_50';")								                
		print
		NON_flag = "non_source"
	    else :
		NON_flag = "source"
	else :							
            if NON_flag == "non_source" : 
		from Bio.Seq import Seq
	        line_revcomp = Seq(line[:-1]).reverse_complement()
		sys.stdout.write(str(line_revcomp))
		print

    fp_non.close() 


convert_multiline_to_single_line_FASTA ()
create_read_pair_tuples()
make_direction_dict()
make_fastg()
print_non_sources()
print_non_sources_complement()



