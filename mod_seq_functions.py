#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 31 16:13:36 2025

@author: jrwill
"""

import networkx as nx

from pytheas_global_vars import pgv, pgc


def parse_mod_seq(raw_seq): # parse sequence for mod synonyms in brackets ["mod"] # digest
    
    length = len(raw_seq)
    lb = [pos for pos, char in enumerate(raw_seq) if char == "["] + [length] # find bracket positions
    rb = [pos for pos, char in enumerate(raw_seq) if char == "]"] + [length]
    slist = []
    
    i = 0
    lb_flag = 0
    for i in range(length): # iterate thru sequence, capturing mods between brackets
        if i in lb:
            lb_flag = 1
            mod = []
            continue
        if i in rb:
            lb_flag = 0
            slist = slist + ["".join(mod)]
            continue
        if lb_flag == 1:
            mod.append(raw_seq[i])
        else:
            slist.append(raw_seq[i])

    slist = [s for s in slist if s != ""]
    slist = [s for s in slist if s != []]
    bogus_list = [s for s in slist if s not in pgv.nt_key_dict["Pytheas_ID"]]
    if bogus_list != []:
        s3list = "MODIFICATION INPUT ERROR " + " ".join(bogus_list)
    else:
        s3list = [pgv.nt_key_dict["Pytheas_ID"][s] for s in slist] # generate 3-base code for mods
    
    return s3list

def parse_1_letter(raw_seq):  # parse raw fasta sequence to 3-letter...all rna_mods except "fasta" 
    if pgv.rna_mods == "none" or pgv.rna_mods == "modfile":
        code_key = "Pytheas_ID"
    else:
        code_key = pgv.rna_mods

    length = len(raw_seq)
    s3list = []
    for i in range(length):
        m1 = raw_seq[i]
        b3 = pgv.nt_key_dict[code_key][m1]
        if b3 in pgv.nt_def_dict:
            s3list.append(b3)
        else:
            print ("problem with code ", pgv.rna_mods, m1, " in raw sequence ", raw_seq)
    return s3list

def generate_mod_seq(seq3): # convert 3-letter seq_list to 1-base + [mod] string 
    mod_seq_list = []
    for base in seq3:
        if base == "XXX":  # what's this for???
            mod_seq_list.append('X')
        elif pgv.nt_def_dict[base]["Type"] == "natural":
            mod_seq_list.append(pgv.nt_def_dict[base]["Pytheas_1_letter_code"])
        else:
            if type(pgv.nt_def_dict[base]["Pytheas_ID"])  != float:
                mod_seq_list.append("[" + pgv.nt_def_dict[base]["Pytheas_ID"] + "]")
            else:
                mod_seq_list.append("[" + base + "]")
    mod_seq = "".join(mod_seq_list)
    
    return mod_seq

def generate_mod_seq_ends(frag3): # convert 3-letter seq_list to 1-base + [mod] string 
    end5 = frag3[0]
    seq3 = frag3[1:-1]
    end3 = frag3[-1]
    mod_seq_list = []
    for base in seq3:
        if base == "XXX":  # what's this for???
            mod_seq_list.append('X')
        elif pgv.nt_def_dict[base]["Type"] == "natural":
            mod_seq_list.append(pgv.nt_def_dict[base]["Pytheas_1_letter_code"])
        else:
            if type(pgv.nt_def_dict[base]["Pytheas_ID"])  != float:
                mod_seq_list.append("[" + pgv.nt_def_dict[base]["Pytheas_ID"] + "]")
            else:
                mod_seq_list.append("[" + base + "]")
    mod_seq = "".join(mod_seq_list)
    e5 = [key for key, val in pgv.end_dict["end5"].items() if val == end5]
    e3 = [key for key, val in pgv.end_dict["end3"].items() if val == end3]
    mod_seq = "_".join([e5[0], mod_seq, e3[0]])
    return mod_seq

def generate_end5_mod_seq_end3(frag3): # convert 3-letter seq_list to 1-base + [mod] string 
    end5 = frag3[0]
    seq3 = frag3[1:-1]
    end3 = frag3[-1]
    mod_seq_list = []
    for base in seq3:
        if base == "XXX":  # what's this for???
            mod_seq_list.append('X')
        elif pgv.nt_def_dict[base]["Type"] == "natural":
            mod_seq_list.append(pgv.nt_def_dict[base]["Pytheas_1_letter_code"])
        else:
            if type(pgv.nt_def_dict[base]["Pytheas_ID"])  != float:
                mod_seq_list.append("[" + pgv.nt_def_dict[base]["Pytheas_ID"] + "]")
            else:
                mod_seq_list.append("[" + base + "]")
    mod_seq = "".join(mod_seq_list)
    e5 = [key for key, val in pgv.end_dict["end5"].items() if val == end5]
    e3 = [key for key, val in pgv.end_dict["end3"].items() if val == end3]
    # mod_seq = "_".join([e5[0], mod_seq, e3[0]])
    return e5[0], mod_seq, e3[0]


def generate_molecular_graph(f3, label): # new toplogy
    G = nx.Graph() # initialize molecular graph

    for resno in range(len(f3)): # add node to graph for each group
        for g in pgv.nt_fragment_dict[label][f3[resno]].groups:
            node = "_".join([f3[resno], str(resno), g])
            G.add_node(node, base = f3[resno], resno = resno, group = g)

#TODO calculate node mass here 

    for node in G.nodes:  # fill out graph nodes with data for fragmentation from topo dict
        nd = G.nodes[node]
        res =nd["resno"]
        base = nd["base"]
        # print("generate_molecular_graph", res, base,nd["group"])
        nd.update(pgv.child_dict[nd["group"]]) # add topology info
        parent_group, child_group = nd["added_edge"].split("_")
        child_res = res + nd["child_offset"]
        child_base = f3[child_res]

        if nd["parent_offset"] == 0: # standard linkage
            parent_node = "_".join([base, str(res), parent_group])
            child_node = "_".join([child_base, str(child_res), child_group])
        else: # 3'-terminal linkage
            parent_res = res + nd["parent_offset"] 
            parent_base = f3[parent_res]
            parent_node = "_".join([parent_base, str(parent_res), parent_group])
            child_node = "_".join([base, str(res), child_group])

        nd["parent_node"] = parent_node
        if child_node not in G.nodes: # patch for 3'-linkage
            nd["child_node"] = "none"
            continue
        else:
            nd["child_node"] = child_node
        
        G.add_edge(parent_node, child_node)

#TODO   reconcile with discovery and next_node    
        # nd["group_mass"] = pgv.nt_fragment_dict["light"][G.nodes[node]["base"]].mass_dict[G.nodes[node]["group"]]
        # print("generate_mol_graph: nd.keys()", nd.keys())

    return G

