import sys
import csv
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqIO import FastaIO

list_of_tuples = []
all_nodes_list = []
all_source_nodes = []
follows = defaultdict(list)


sys.stdout = open('output.fastg', 'w')


# converts multiline preads4falcon.fasta into a single line fasta
def convert_multiline_to_single_line_FASTA () :
    sequences = [] 
    input_handle = open("preads4falcon.fasta", "rU")

    for record in SeqIO.parse(input_handle, "fasta"):
	sequences.append(record)
    
    output_handle = open("formatted_preads4falcon.fasta","w")
    fasta_out = FastaIO.FastaWriter(output_handle, wrap=None)
    fasta_out.write_file(sequences)
    output_handle.close()


# creates a tuple for every pair in sg_edges_list
def create_read_pair_tuples () :
    with open("sg_edges_list") as sg_entries:
        read_pairs_TMI = csv.reader(sg_entries, delimiter=' ')
        for row in read_pairs_TMI :
	    read_pair_tuple = ((row[0]),(row[1])) 	
            list_of_tuples.append(read_pair_tuple)
    list_of_tuples.sort()	    

# collapses multiple tuples into a dictionary
# key is the first entry in a tuple; each key represents "Source" node
# value(s) is(are) the second entry(entries) in tuple(s); are the "Sink" nodes
def make_follows_dict () :
    for pair in list_of_tuples :
	start_node_location = pair[0][-1] 
        end_node_location = pair[1][-1]
        
        source_node = pair[0][:-2]
        sink_node = pair[1][:-2]

        if start_node_location == "E" and end_node_location == "E" :
	    follows[source_node].append(sink_node)
            
	elif start_node_location == "E" and end_node_location == "B" :
	    follows[source_node].append(sink_node+"'")
                        
	elif start_node_location == "B" and end_node_location == "E" :
            follows[source_node+"'"].append(sink_node)

	elif start_node_location == "B" and end_node_location == "B" :
            follows[source_node+"'"].append(sink_node+"'")
	    
        else : 
	    print "ERROR : the tuples have not been parsed correctly"    


# used by make_fastg() to provide padding for the header
def headerify (node) :
    if node[-1] == "'" :
        return str("NODE_"+node[:-1]+"_length_"+"500"+"_cov_50'")
    else :
        return str("NODE_"+node+"_length_"+"500"+"_cov_50")


# used by make_fastg() to extract sequence given record name in FASTA file
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


# calls headerify() and corresponding_fasta() on 'follows' dictionary
# makes the first part of FASTG file which has nodes with overlaps 
def make_fastg() :
    for source, sinks in follows.items() : 
        all_source_nodes.append(source)	    
	header = []
	header.append(">"+headerify(source)+":")
	for element in sinks :
	    header.append(headerify(element)+",")
	str_header = ''.join(header)
	sys.stdout.write(str_header[:-1]+';')
	print
	corresponding_FASTA(source)


# makes the second part of FASTG file which has fwd nodes without overlaps
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


# makes the third part of FASTG file which has revcomp nodes without overlaps
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
make_follows_dict()
make_fastg()
print_non_sources()
print_non_sources_complement()



