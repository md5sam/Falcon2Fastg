#!/usr/bin/env python

import sys
import argparse
import csv
import os.path


from collections import defaultdict
from Bio import SeqIO
from Bio.SeqIO import FastaIO
from pyfaidx import Fasta

list_of_tuples = []
list_of_ctg_start_and_end_node_tuples = []
all_nodes_list = []
all_source_nodes = []
all_ctg_list = []
all_ctg_source_nodes = []
follows = defaultdict(list)
ctg_follows = defaultdict(list)
contigs = defaultdict(list)
circular_ctgs = defaultdict(list)
linear_ctgs = defaultdict(list)
p_ctg_names_set = set ()


parser = argparse.ArgumentParser(description='Falcon2Fastg converts FALCON output to FASTG format')
parser.add_argument('-m','--mode', help='Enter MODE as either "read", "contig" or "both"', required=False)
args = vars(parser.parse_args())

if args['mode'] == "read" :
    mode = "read"
elif args['mode'] == "contig" :
    mode = "contig"
elif args['mode'] == "both" :
    mode = "both"
else :
    mode = "both"	


# checks if sg_edges_list, p_ctg.fa, ctg_paths and preads4falcon.fasta are present 
def check_files_exist () :
    if os.path.isfile("preads4falcon.fasta") == True and os.path.isfile("sg_edges_list") == True :
	if mode == "contig" or mode == "both" :
	    if os.path.isfile("p_ctg.fa") == True and os.path.isfile("ctg_paths") : 
	        return True
	    else :
		print 
		print "ERROR! Please make sure p_ctg.fa and ctg_paths are in this directory"
		print 
		return False 
        return True
    else : 
        print
        print "ERROR! Please make sure sg_edges_list, preads4falcon.fasta, ctg_paths and p_ctg.fa are in this directory"
        print
        return False        


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


# creates a tuple for every "G" pair in sg_edges_list
def create_read_pair_tuples () :
    with open("sg_edges_list") as sg_entries:
        read_pairs_TMI = csv.reader(sg_entries, delimiter=' ')
        for row in read_pairs_TMI :
	    if row[-1] == "G" :
                read_pair_tuple = ((row[0]),(row[1]))
                list_of_tuples.append(read_pair_tuple)
    list_of_tuples.sort()	    


# these are all the ctgs for which you actually have the sequence in p_ctg.fa
def create_p_ctg_names_set () :
    p_ctg_name_list = []
    p_ctg_headers_list = []
    global p_ctg_names_set
    fp_p_ctg = open ("p_ctg.fa")
    for line in fp_p_ctg :
        if line.startswith('>'):
            p_ctg_headers_list.append(line)
    p_ctg_headers = csv.reader(p_ctg_headers_list, delimiter=' ', skipinitialspace=True)
    for row in p_ctg_headers :
        p_ctg_name = str(row[0][1:])
        p_ctg_name_list.append(p_ctg_name)
        p_ctg_name_list.append(p_ctg_name+"'")
    p_ctg_names_set = set (name for name in p_ctg_name_list)


# extract headers from ctg_paths, and converts it into a dictionary
# of the form [ctg_name] : {ctg_start_read, ctg_end_read}
def create_contig_dict () :
    ctg_headers_list = []
    fp_ctg = open ("ctg_paths")
    for line in fp_ctg :
        ctg_headers_list.append(line)
    ctg_headers = csv.reader(ctg_headers_list, delimiter=' ', skipinitialspace=True)
    for row in ctg_headers :
        ctg_name = str(row[0])
        ctg_start_node = (row[2].split("~"))[0]
        ctg_end_node = row[3]
        ctg_type = str(row[1])
        if ctg_type == "ctg_linear" :
            if ctg_name not in p_ctg_names_set : 
                if ctg_name[-1] == "R":
                    ctg_name = ctg_name[:-1]+"F'"
                else :
                    ctg_name = ctg_name[:-1]+"R'"
        else : 
            for existing_ctg in circular_ctgs :
                if ctg_start_node[:-1]+"E" in circular_ctgs[existing_ctg] or ctg_start_node[:-1]+"B" in circular_ctgs[existing_ctg] :
                    ctg_name = existing_ctg+"'" 
           
            circular_ctgs[ctg_name].append(ctg_start_node)
        all_ctg_list.append(ctg_name)
        contigs[ctg_name].append(ctg_start_node)
        contigs[ctg_name].append(ctg_end_node)         
          
              
# collapses multiple tuples into a dictionary
# key is the first entry in a tuple; each key represents "Source" node
# value(s) is(are) the second entry(entries) in tuple(s); are the "Sink" nodes
def make_read_connections () :
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


# returns True if overlap exists between a pair of reads
# used by make_contig_connections()
def ovlp_exists(curr_node, tgt_node) :
    curr_extremity = curr_node[-1]
    tgt_extremity = tgt_node[-1]  
    curr_read = curr_node[:-2]
    tgt_read = tgt_node[:-2]
   
    if curr_extremity == "E" and tgt_extremity == "B" :
        if tgt_read in follows[curr_read] :
            return True
    elif curr_extremity == "E" and tgt_extremity == "E" :
        if tgt_read+"'" in follows[curr_read] :
            return True
    elif curr_extremity == "B"and tgt_extremity == "B" :
        if tgt_read in follows[curr_read+"'"] :
            return True
    elif curr_extremity == "B" and tgt_extremity == "E" :
        if tgt_read+"'" in follows[curr_read+"'"] :
            return True
    return False


# this function creates a ctg_follows dictionary
# key is the ctg, value(s) is (are) its successor(s)
def make_contig_connections() :  
    count = 0
    for curr_ctg in contigs :
        ctgs_to_investigate = []
        for all_ctgs in all_ctg_list :    
            if curr_ctg != all_ctgs : 
                ctgs_to_investigate.append(all_ctgs)
        
        curr_ctg_end_node = contigs[curr_ctg][1]
        curr_ctg_start_node = contigs[curr_ctg][0]
      
        # check if contig is circular
        if curr_ctg_start_node == curr_ctg_end_node :
            ctg_follows[curr_ctg].append(curr_ctg)        

        for tgt_ctg in ctgs_to_investigate :
            tgt_ctg_start_node = contigs[tgt_ctg][0]
            tgt_ctg_end_node = contigs[tgt_ctg][1]

            # curr_end is the same as tgt_start
            if curr_ctg_end_node == tgt_ctg_start_node :
                ctg_follows[curr_ctg].append(tgt_ctg)

            # curr_end overlaps with tgt_start
            if ovlp_exists(curr_ctg_end_node,tgt_ctg_start_node) == True :
                ctg_follows[curr_ctg].append(tgt_ctg)


# used by make_fastg() to provide padding for the header
def headerify (node) :
    if node[-1] == "'" :
        return str("NODE_"+node[:-1]+"_length_"+"500"+"_cov_50'")
    else :
        return str("NODE_"+node+"_length_"+"500"+"_cov_50")


# calls headerify() on 'follows' dictionary
# makes the first part of FASTG file which has nodes with overlaps 
def make_fastg(node_mode) :
    if node_mode == "read" :
	dict_ = follows
        idx = Fasta('formatted_preads4falcon.fasta')
    else : 
	dict_ = ctg_follows    
        idx = Fasta('p_ctg.fa')
    for source, sinks in dict_.items() : 
	if node_mode == "read" :
	    all_source_nodes.append(source)	    
        else :
	    all_ctg_source_nodes.append(source)	
	header = []
        header.append(">"+headerify(source)+":")
        for element in sinks :
	    header.append(headerify(element)+",")
        str_header = ''.join(header)
	sys.stdout.write(str_header[:-1]+';')
        print
        if source[-1] != "'" :
            corrected_name = source
            print idx[corrected_name][:]
        else :
            corrected_name = source[:-1]
            print idx[corrected_name][:].complement.reverse
                

# makes the second part of FASTG file which has fwd nodes without overlaps
def print_non_ovlp_sources(node_mode) :
    if node_mode == "read" :
        fp_non = open("formatted_preads4falcon.fasta")
    else :
	fp_non = open("p_ctg.fa")    
    NON_flag = "non_source"
    for line in fp_non :
	if line[0] == ">" :
            if node_mode == "read" :
	        record_name = line[1:-1]
		nodes_with_overlaps = all_source_nodes
	    else :
		record_name = line.split(" ")[0][1:]
		nodes_with_overlaps = all_ctg_source_nodes
	    if record_name not in nodes_with_overlaps :
		sys.stdout.write(">NODE_"+record_name+"_length_"+"500"+"_cov_50;")
		print
		NON_flag = "non_source"
	    else :  
		NON_flag = "source"	    	
	else :
            if NON_flag == "non_source" : 
		print line[:-1]
    fp_non.close()     


# makes the third part of FASTG file which has revcomp nodes without overlaps
def print_non_ovlp_sources_complement(node_mode) :
    if node_mode == "read" :
        fp_non = open("formatted_preads4falcon.fasta")
    else :
        fp_non = open("p_ctg.fa")
    NON_flag = "non_source" 
    for line in fp_non :
	if line[0] == ">":
	    if node_mode == "read" :
		record_name = str(line[1:-1]+"'")
	        nodes_with_overlaps = all_source_nodes
	    else :
		record_name = str(line.split(" ")[0][1:]+"'")
	        nodes_with_overlaps = all_ctg_source_nodes
	    if record_name not in nodes_with_overlaps :
		sys.stdout.write(">NODE_"+record_name[:-1]+"_length_"+"500"+"_cov_50';")
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


if __name__ == "__main__": 
    #print "Mode is : ", mode 
    if check_files_exist() == True :
        convert_multiline_to_single_line_FASTA ()
        create_read_pair_tuples()
        make_read_connections()
        if mode == "read" or mode == "both" :    
            sys.stdout = open('reads.fastg', 'w')
	    make_fastg("read")
            print_non_ovlp_sources("read")
            print_non_ovlp_sources_complement("read")
	    sys.stdout.close() 
        if mode == "contig" or mode == "both":
	    sys.stdout = open('contigs.fastg','w')	
            create_p_ctg_names_set()
            create_contig_dict()
            make_contig_connections()
            make_fastg("contig")
            print_non_ovlp_sources("contig")
	    print_non_ovlp_sources_complement("contig")
	    sys.stdout.close()

