#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 24 08:28:42 2023

@author: jrwill
"""

import os
import copy
from random import shuffle
import xlsxwriter


import networkx as nx
import pandas as pd

from pytheas_global_vars import pgv
from pytheas_IO import read_pytheas_file
from pytheas_objects import fragment_sequence
from worksheet_functions import format_worksheet_columns


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

def add_modifications(): # parse raw fasta seq and add modifications

    # generate seq3 from raw_seq
    for seq_id, sdict in pgv.mol_dict.items():
        if pgv.rna_mods != 'fasta':  # none, modfile, or 1-letter 
            # print(seq_id, sdict.keys())
            sdict["seq3"] = parse_1_letter(sdict["raw_seq"])
        else:
            sdict["seq3"] = parse_mod_seq(sdict["raw_seq"])  # parse mods in brackets from fasta
            
    if pgv.rna_mods == 'modfile':  # read modifications in from mod def file
    
        read_pytheas_file("mod_file")
        for key, md in pgv.mod_dict.items():
            mk = md["Molecule"]
            mp = md["Position"]
            base = md["ID_ext"].replace("[","").replace("]","")
            if base in pgv.nt_key_dict["Pytheas_ID"]:
                base3 = pgv.nt_key_dict["Pytheas_ID"][base]
            else:
                print("Problem with modification ", key, md["ID"], md["ID_ext"])
                continue

            orig_base = pgv.nt_def_dict[base3]["Originating_base"] # 3-letter code
            if mk in pgv.mol_dict:  
                base = pgv.mol_dict[mk]["seq3"][mp-1]   # check parent base is correct
                if base != orig_base:
                    print("check input: modification ", base, "has wrong parent base ", orig_base)
                else:
                     pgv.mol_dict[mk]["seq3"][mp - 1] = base3   # need to substitute mod into seq_list
 
    active_mods = []
    pgv.mod_dict = {}
    midx = 0
    for seq_id, sdict, in pgv.mol_dict.items():
        print("seq_id", seq_id, sdict["seq3"])
        sdict["modseq"] = generate_mod_seq(sdict["seq3"])
        for i in range(len(sdict["seq3"])):
            base = sdict["seq3"][i]
            if pgv.nt_def_dict[base]["Type"] == "mod":
                sdict["mods"][i+1] = base   # could expand this, but not sure why
                if base not in active_mods:
                    active_mods.append(base)
                pgv.mod_dict[midx] = {"Molecule": seq_id, "Position": i+1, "ID": base}
                midx += 1
            
    return active_mods

def enzyme_digest():
    nfrags = 0
    unique_frag_dict = {}
    for enzyme in pgv.enzymes:
        print("digesting with ", enzyme)
        for mol, mdict in pgv.mol_dict.items():    
            frag_seq_list = digest_sequence(mol, mdict, unique_frag_dict, enzyme)  # do enzyme digest -> unique_frag_seq_dict
            mdict["frag_seq_key_list"] = frag_seq_list

            mdict["digest"] = pgv.enzymes
    nfrags += len(unique_frag_dict.keys())

    return nfrags, unique_frag_dict

def count_mods(frag3):
    mod_bases = [b for b in frag3 if pgv.nt_def_dict[b]["Type"] == "mod"]
    ct = len(mod_bases)
    return ct

def unique_fragment(mol, end5, seq3, end3, fr, to, length, miss, ufdict, frag_seq_key_list):
    frag3 = [pgv.end_dict["end5"][end5]] + seq3[fr:to] + [pgv.end_dict["end3"][end3]]
    f3 = fragment_sequence(frag3) # object with 1-letter/3-letter seqs and ends
    full_mod_count = count_mods(frag3)
    
    partial_key = mol + ":" + str(fr + 1) + "_" + str(to-1) + ":" + f3.frag
    partial = 0
    if pgv.partial_mods == 'y':
        frag3_list = []
        for base in frag3:
            if pgv.nt_def_dict[base]["Precursors"] == "none":
                frag3_list.append([base])
            else:
                partial += 1
                plist = [pgv.nt_key_dict["Pytheas_ID"][p] for p in pgv.nt_def_dict[base]["Precursors"].split(":")]
                frag3_list.append([base] + plist)
        sequences = []
        expand_partial_mod_sequences(0,[],frag3_list,sequences)
    else:
        sequences = [frag3]

    for frag3 in sequences:                
        f3 = fragment_sequence(frag3) # object with 1-letter/3-letter seqs and ends
        partial_mod_count = count_mods(frag3)
        if full_mod_count == 0:
            frac_mod = -1.0
        else:
            frac_mod = round(partial_mod_count/full_mod_count, 2)
    
        frag_seq_key = ":".join([mol, "_".join([str(fr+1), str(to-1)]), f3.frag])
        ufdict[frag_seq_key] = {"mol": mol, "frag3": frag3}
        ufdict[frag_seq_key].update({"fr": fr, "to": to, "length": length, "miss": miss, 
                                     "partial_mod": frac_mod, "partial_key": partial_key})
        frag_seq_key_list.append(frag_seq_key)

def expand_partial_mod_sequences(i, current, candidates, sequences):   # recursive algorithm for expanding combinations of sequences
    # i is pointer index, current is current sequence, candidates is list of nucleotides, sequences has results
    if i == len(candidates):
        sequences.append(current.copy())
        return
     
    for j in range(len(candidates[i])):
        current.append(candidates[i][j])
        expand_partial_mod_sequences(i+1, current, candidates, sequences)
        current.pop()

def generate_cut_list(seq3, enzyme):
    
    pattern, cut_idx, end_5, end_3 = [pgv.enzyme_dict[enzyme][key] for key in ["pattern", "cut_idx", "end_5", "end_3"]]
    
    match_index_list = []
    for patt in pattern :
        for i in range(len(seq3)):
            if  patt == seq3[i:i + len(patt)]:
                match_index_list.append(i)

    length = len(seq3)
    last = length + 1
    
    cut_index_list = [0] + [m + cut_idx for m in match_index_list] + [last] # add cut offset from match
    misses = min(pgv.miss, len(cut_index_list)-2) # of matches - 2 ends
    cut_list = [[cut_index_list[i], cut_index_list[i+miss + 1], miss] # loop thru cut_index_list for each miss 
                for miss in range(misses+1) 
                for i in range(len(cut_index_list) - miss - 1)]

    return cut_list

def generate_custom_cut_list(seq3, mol, ufdict, frag_seq_key_list):
    
    #TODO check for sensible range with missed cleavages
    cut_list = []
    length = len(seq3)
    last = length + 1
    
    for key, cd in pgv.custom_cleavage_dict.items():
        pattern, cut_idx, end_5, end_3 = [cd[key] for key in ["pattern", "cut_idx", "end_5", "end_3"]]
        
        match_index_list = []
        for patt in pattern :
            for i in range(len(seq3)):
                if  patt == seq3[i:i + len(patt)]:
                    match_index_list.append(i)
                    
        cut_index_list = [0] + [m + cut_idx for m in match_index_list] + [last] # add cut offset from match

        misses = min(pgv.miss, len(cut_index_list)-2) # of matches - 2 ends
        cut_list = [[cut_index_list[i], cut_index_list[i+miss + 1], miss] # loop thru cut_index_list for each miss 
                    for miss in range(misses+1) 
                    for i in range(len(cut_index_list) - miss - 1)]
        
        process_cut_list(seq3, cut_list, mol, ufdict, frag_seq_key_list)

def generate_random_cut_list(seq3):
    
    min_length = pgv.nonspecific_min_length
    max_length = pgv.nonspecific_max_length
    
    length = len(seq3)
    cut_list = []
    for fr in range(0, length - min_length + 1):
        max_to = min(length + 1, fr + max_length)
        for to in range(fr + min_length, max_to):
             if 0 < to <= length:
                cut_list.append([fr, to, -1])
                    
    return cut_list

def process_cut_list(seq3, cut_list, mol, ufdict, frag_seq_key_list):

     for fr,to,miss in cut_list:     
       if fr == 0:
           end5_list = pgv.mol_end5
       else:
           end5_list = pgv.frag_end5
       if pgv.use_iTRAQ == 'y':
           end3_list = ['tag']
       else:
           if to == len(seq3) + 1:
               end3_list = pgv.mol_end3
           else:
               end3_list = pgv.frag_end3
       
       for end5 in end5_list:
          for end3 in end3_list:
               if end3 == "tag":
                   s3 = seq3[fr:to] + ["CTG"]
               else:
                   s3 = seq3[fr:to]
               
               length = len(s3)
               if pgv.min_length <= length <= pgv.max_length:
                   print(mol, end5, end3, fr, to, length, miss)
                   unique_fragment(mol, end5, seq3, end3, fr, to, length, miss, ufdict, frag_seq_key_list)


def digest_sequence(mol, mdict, ufdict, enzyme): # process one sequence # digest
     
     seq3 = mdict["seq3"]
     frag_seq_key_list = []
     miss = 0

     if enzyme == "none":
         fr = 0
         to = len(seq3)
         length = len(seq3)
         end5_list = mdict["mol_5_end"]
         end3_list = mdict["mol_3_end"]
          
         for end5 in end5_list:
             for end3 in end3_list:
                 unique_fragment(mol, end5, seq3, end3, fr, to, length, miss, ufdict, frag_seq_key_list)
         return frag_seq_key_list

     elif enzyme == "nonspecific":
        cut_list = generate_random_cut_list(seq3)
        
     elif enzyme == "custom":
        read_pytheas_file("custom_cleavage_file")
        generate_custom_cut_list(seq3, mol, ufdict, frag_seq_key_list)
        return frag_seq_key_list
        
     else:        
          cut_list = generate_cut_list(seq3, enzyme)
         
     process_cut_list(seq3, cut_list, mol, ufdict, frag_seq_key_list)
         
     return frag_seq_key_list

def add_precursor_ions(frag_dict, z_limit): # zlimit is specific z or "all" to take from charge tables # digest
    nprec = 0
    idx = 0
    for seq, sdict in frag_dict.items(): # seq is in the form end5_seq_end3
        ion_frag_dict = {}
        for label in pgv.isotopic_species: # loop thru labels starting here
            G = generate_molecular_graph(sdict["frag3"], label) # needs frag3 argument
            ion_fragment_dict = calculate_m0_mz1(G, label, sdict["frag3"], z_limit)
            ion_frag_dict[label] = ion_fragment_dict
 
        sdict["ion_frag_dict"] = ion_frag_dict

        for label in pgv.isotopic_species:
            nprec += len(frag_dict[seq]["ion_frag_dict"][label]["mz1"].keys())
        idx += 1
               
    return nprec


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

    return G


def calculate_m0_mz1(G, label, seq3n, z_limit): 

    m0 = 0.0
    for node in G.nodes:  #   calculate monoisotopic mass
        nd = G.nodes[node]
        base =nd["base"]
        grp = nd["group"]
#   TODO get mass from graph node_dict
        m0 += pgv.nt_fragment_dict[label][base].mass_dict[grp] # this is a bit wasteful
   
    if pgv.ion_mode == "-":  
        zsign = -1
    else:
        zsign = 1
    length = min(len(seq3n) - 2 , pgv.max_length) # subtract 2 for ends
    if z_limit == "all":
        z1_list = [z * zsign for z in pgv.MS_charge_dict["MS1_charges"][length]]
    else:
        z1_list = [int(z_limit)]

    z_dict = {z: (m0 + z * pgv.hmass)/abs(z) for z in z1_list} 
    ion_fragment_dict = {"M0": m0, "mz1": z_dict}

    return ion_fragment_dict


def build_frag_dict(unique_frag_dict): # fragments as dictionary
    frag_dict = {}
    for fkey,fdict  in unique_frag_dict.items(): 
            frag3 = fdict["frag3"]
            f3 = fragment_sequence(frag3)
            if f3.frag not in frag_dict:
                frag_dict[f3.frag] = {}
                frag_dict[f3.frag]["seq_list"] = [fkey]
                for key,value in fdict.items(): # copy ufd
                    frag_dict[f3.frag][key] = value
                frag_dict[f3.frag]["frag_type"] = "target"
            else:
                  frag_dict[f3.frag]["seq_list"].append(fkey)
                  
    if pgv.add_decoys == 'y':
        decoy_dict = {}
        didx = 0
        for frag, fdict in frag_dict.items():
            f3 = fragment_sequence(fdict["frag3"])
            
            flag = "n"
            for i in range(len(f3.seq3) ** 4):
                shuffle(f3.seq3)
                f3.generate_mod_seq()
                decoy_frag = "_".join([f3.end5, f3.modseq, f3.end3])
                if decoy_frag not in frag_dict:
                    decoy_dict[decoy_frag] = {}
                    for key,value in fdict.items(): # copy dict
                        decoy_dict[decoy_frag][key] = value

                    decoy_dict[decoy_frag]["seq_list"] = ["decoy_" + str(didx) + ":1_" + str(len(f3.seq3)) + ":" + decoy_frag]
                    decoy_dict[decoy_frag]["mol"] = "decoy_" + str(didx)
                    decoy_dict[decoy_frag]["frag3"] = [f3.end5n] + f3.seq3 + [f3.end3n]
                    decoy_dict[decoy_frag]["frag_type"] = "decoy"
                    
                    flag = "y"
                    didx += 1
                    break
            if flag == "n":
                print(" no decoy generated for ", frag)
        frag_dict.update(decoy_dict)
          
    return frag_dict

def build_precursor_dict():   # expand frag_dict by z and label for complete list of precursors
    unique_precursor_dict = {}  # key is updkey = end5_seq_end3_z_label
    prec_idx = 0
    for fkey, fdict in pgv.frag_dict.items():
        for label,ldict in fdict["ion_frag_dict"].items():
            m0 = ldict["M0"]
            for z,mz in ldict["mz1"].items():
                updkey = "_".join([fkey, str(z), label])
                if updkey in unique_precursor_dict:
                    print("duplicate updkey ", updkey)
                unique_precursor_dict[updkey] = copy.deepcopy(fdict)
                unique_precursor_dict[updkey].pop("ion_frag_dict") # each precursor has its own charge
                unique_precursor_dict[updkey].update({"m0": m0, "mz1": mz, "z": z, "label":label})
                prec_idx += 1

    return unique_precursor_dict

def output_digest_file(output_file): 
    
    digest_df = pd.DataFrame.from_dict(pgv.unique_precursor_dict, orient="index").reset_index()
    
    fd = pgv.output_format_dict["digest_output"]  # fd = "format_dict" is formatting information
    # reorder columns according to fd
    column_order = [col for col in sorted(list(fd.keys()), key = lambda x: fd[x]["order"]) if col in digest_df.columns] # sort order
    column_order = [col for col in column_order if fd[col]["order"] < 100] # eliminate those with 999 as "not shown"
    digest_df = digest_df[column_order]
    
    digest_df.frag3 = digest_df.frag3.apply(generate_mod_seq_ends) # change list to mod_seq

    digest_output_file = os.path.join(pgv.job_dir, output_file + ".xlsx")
    workbook = xlsxwriter.Workbook(digest_output_file,{"nan_inf_to_errors": True})
    worksheet = workbook.add_worksheet(output_file)
    
    for okey,odict in fd.items(): # set up  column formats for xlsxwriter, list of two for green/white alternation

    # gridlines hidden if bg color
       if "0" in odict["format"]:
           odict["xformat"] =workbook.add_format({'num_format': odict["format"]})
       else:
           odict["xformat"] = workbook.add_format({'num_format': "0"})


    bold = workbook.add_format()
    bold.set_bold()
    
    format_worksheet_columns(worksheet, digest_df, fd)   # sets the column widths

    for col in digest_df.columns:  # output header
        worksheet.write(0, digest_df.columns.get_loc(col), fd[col]["label"], bold)  # fd has user-defined label column if desired
    
    # write out each cell
    rowidx = 1
    for row in digest_df.iterrows():
        colidx = 0
        for col in row[1].keys():
            if type(row[1][col]) == list:
                row[1][col] = ",".join(row[1][col])
            worksheet.write(rowidx, colidx, row[1][col], fd[col]["xformat"])
            colidx +=1
        rowidx +=1
   
    workbook.close()


