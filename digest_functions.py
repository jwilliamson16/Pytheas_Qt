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
import pickle
import math


import pandas as pd
import numpy as np

from pytheas_global_vars import pgv, pgc
from pytheas_IO import read_pytheas_file
from pytheas_objects import fragment_sequence
from worksheet_functions import format_worksheet_columns
from matrix_plot_functions import (Labeled_Matrix, make_mod_text_dict, sequence_color, matrix_plot, 
                                   make_seq_text_dict, make_long_text_dict, make_text_dict,
                                   matrix_plot_new, pad_end3_gray_matrix, color_mods_matrix,
                                   color_long_matrix_by_seq, pad_end3_gray, color_mods,
                                   color_matrix_by_seq)
from mod_seq_functions import (parse_mod_seq, generate_mod_seq, generate_mod_seq_ends,
                               parse_1_letter, generate_molecular_graph)


def add_modifications(): # parse raw fasta seq and add modifications

    # generate seq3 from raw_seq
    for seq_id, sdict in pgv.mol_dict.items():
        # mol = pgv.molecule_dict[seq_id]
        # mol = pgv.mol_dict[seq_id]
        if pgv.rna_mods != 'fasta':  # none, modfile, or 1-letter 
            # print(seq_id, sdict.keys())
            sdict["seq3"] = parse_1_letter(sdict["raw_seq"])
            # mol.seq3 = parse_1_letter(mol.raw_seq)
        else:
            sdict["seq3"] = parse_mod_seq(sdict["raw_seq"])  # parse mods in brackets from fasta
            # mol.seq3 = parse_mod_seq(mol.seq3)
    # pickle.dump(pgv.molecule_dict, open(os.path.join(pgv.job_dir, "molecule_2+dict.pkl"), "wb" ))
            
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
 
    # mod_dict primarily used in match sequence plotting
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
        if enzyme == "custom":
            read_pytheas_file("custom_cleavage_file")
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
    
    partial_key = mol + ":" + str(fr + 1) + "_" + str(to) + ":" + f3.frag
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
    
        frag_seq_key = ":".join([mol, "_".join([str(fr+1), str(to)]), f3.frag])
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
    last = length
    
    cut_index_list = [m + cut_idx for m in match_index_list] # add cut offset from match
    if cut_index_list[0] > 0:
        cut_index_list = [0] + cut_index_list
    if cut_index_list[-1] < last:
        cut_index_list = cut_index_list + [last]
    
    misses = min(pgv.miss, len(cut_index_list)-2) # of matches - 2 ends
    cut_list = [[cut_index_list[i], cut_index_list[i+miss + 1], miss] # loop thru cut_index_list for each miss 
                for miss in range(misses+1) 
                for i in range(len(cut_index_list) - miss - 1)]

    return cut_list  # [fr:to]

def generate_custom_cut_list(seq3, mol, ufdict, frag_seq_key_list):
    
    #TODO check for sensible range with missed cleavages
    cut_list = []
    template_cut_index_list = []
    length = len(seq3)
    last = length
    # last = length + 1  # changed for digest plot bug
    
    for key, cd in pgv.custom_cleavage_dict.items():
        for cut_idx in cd["cut_idx"]:
            pattern, end_5, end_3 = [cd[key] for key in ["pattern", "end_5", "end_3"]]
            match_index_list = []
            for patt in pattern :
                for i in range(len(seq3)):
                    if  patt == seq3[i:i + len(patt)]:
                        match_index_list.append(i)
            print("match_index_list", mol, key, match_index_list)          
                    
        # cut_index_list = [0] + [m + cut_idx for m in match_index_list] + [last] # add cut offset from match
            template_cut_index_list += [m + cut_idx for m in match_index_list] # add cut offset from match
            print("template_cut_index_list", template_cut_index_list)
    cut_index_list = sorted(list(set([0] + template_cut_index_list + [last])))  # sorted master list
    # cut_index_list = list(set([0] + template_cut_index_list + [last])) # sorted master list
    print("cut_index_list", cut_index_list)
    misses = min(pgv.miss, len(cut_index_list)-2) # of matches - 2 ends
    cut_list = [[cut_index_list[i], cut_index_list[i+miss + 1], miss] # loop thru cut_index_list for each miss 
                for miss in range(misses+1) 
                for i in range(len(cut_index_list) - miss - 1)]
    
    print("cut list", mol)
    print(cut_list)
    print()
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
     # print("cut_list", cut_list)
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
                   # print("process_cut_list: ", mol, end5, end3, fr, to, length, miss)
                   unique_fragment(mol, end5, seq3, end3, fr, to, length, miss, ufdict, frag_seq_key_list)


def digest_sequence(mol, mdict, ufdict, enzyme): # process one sequence # digest
     
     seq3 = mdict["seq3"]
     frag_seq_key_list = []
     miss = 0

     if enzyme == "none":
         fr = 0
         # to = len(seq3) + 1
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
        # read_pytheas_file("custom_cleavage_file")
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
    # length = min(len(seq3n) - 2 , pgv.max_length) # subtract 2 for ends
    length = min(len(seq3n) - 2 , 20) # subtract 2 for ends

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

def output_digest_file(output_file): # move to worksheet functions
    
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

def make_digest_plot(output_file):
    
    # determine matrix dimensions nr = n_seq, nc = max_seq_len
    max_seq_len = 0
    for mol, mdict in pgv.mol_dict.items():
        slen = len(mdict["seq3"])
        if slen > max_seq_len:
            max_seq_len = slen
    # n_seq = len(list(pgv.mol_dict.keys()))

    row_labels = [mol for mol in pgv.mol_dict.keys()]
    col_labels = [str(i+1) for i in range(max_seq_len)]
    
    lm = Labeled_Matrix(row_labels, col_labels)
    pad_end3_gray(lm)
    lm.text_dict = make_text_dict(row_labels, col_labels, 100)
    color_mods(lm, pgc.red) # default color is red
    cleavage_box = True

    for f, fdict in pgv.frag_dict.items():
        if fdict["frag_type"] != "target": # skip decoys
            continue
        seq_list = fdict["seq_list"]
        n_frags = len(seq_list)
        seq3 = fdict["frag3"][1:-1]
        # print("frag", f, seq_list)
        for seq in seq_list:
            mol, r, _ = seq.split(":")
            length = len(pgv.mol_dict[mol]["raw_seq"])
            row_idx = row_labels.index(mol)
            fr, to = map(int,r.split("_"))   # these are sequence indices
            
            # print(row_idx, fr, to, seq3, n_frags )
            color_matrix_by_seq(lm, row_idx, fr, to, seq3, n_frags, length, cleavage_box)
    
    lm.output_file = output_file
    lm.title = "Digest Plot for " + pgv.fasta_file
    matrix_plot_new(lm)       
 
def make_long_digest_plot(output_file):  # plot long sequences as blocks of 100, one plot for each sequence
    
    nc = 100
    for molecule, mdict in pgv.mol_dict.items():
        slen = len(mdict["seq3"])
        nr = math.floor(slen/nc)
        if slen%100 != 0:
            nr += 1
            
        row_labels = [molecule + "(" + str(100*row + 1) + ":" + str(100*(row + 1)) + ")"  for row in range(nr)]
        col_labels = [str(i+1) for i in range(nc)]
       
        lm = Labeled_Matrix(row_labels, col_labels) # matrix plot object
        pad_end3_gray_matrix(lm, slen) # gray out boxes beyond 3'-end
        color_mods_matrix(lm, molecule, pgc.red) # default color is red
        lm.text_dict = make_long_text_dict(row_labels, col_labels, molecule, nc) # text dict with mods and optionally all bases
        
        cleavage_box = True  # make optional?
        
        for f, fdict in pgv.frag_dict.items(): # go thru frag_dict
            if fdict["frag_type"] != "target": # skip decoys
                continue
            seq_list = fdict["seq_list"]
            n_frags = len(seq_list)
            seq3 = fdict["frag3"][1:-1]
            for seq in seq_list:
                mol, r, _ = seq.split(":")
                length = len(pgv.mol_dict[mol]["raw_seq"])

                if mol != molecule:
                    continue
                
                fr, to = map(int,r.split("_"))   # these are sequence indices
                color_long_matrix_by_seq(lm, fr, to, seq3, n_frags, length, cleavage_box)

        lm.output_file = output_file + "_" + molecule
        lm.title = "Digest Plot for " + molecule + " in " + pgv.fasta_file
        matrix_plot_new(lm)
