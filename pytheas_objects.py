#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 10:05:36 2024

@author: jrwill
"""

import math
import numpy as np
import pandas as pd

from pytheas_global_vars import pgv, pgc

# class nt_fragment: # class for nucleotide group fragments
#     def __init__(self, base, df, label):
#         df_dict = df.to_dict(orient = "index")
#         self.base = base
#         self.label = label
#         self.atom_names = df.iloc[:,0].tolist()   # don't need to keep
#         self.atom_symbols = df.iloc[:,1].tolist()  # don't need to keep
#         self.atom_stois = df.iloc[:,2].tolist()
#         self.atom_groups = df.iloc[:,3].tolist() # don't need to keep
#         glist = []
#         for g in self.atom_groups:
#             if g not in glist:
#                 glist.append(g)
#         self.groups = glist
#         self.group_dict = {g:[] for g in self.groups}
#         self.mass_dict = {}
#         for g,s,stoi in zip(self.atom_groups,self.atom_symbols, self.atom_stois):
#               if g not in self.mass_dict.keys():
#                 self.mass_dict[g] = stoi * pgv.atomic_dict[s]["am"]
#               else:
#                 self.mass_dict[g] += stoi * pgv.atomic_dict[s]["am"]
#         self.mass = sum(self.mass_dict.values())
        
class nt_fragment: # class for nucleotide group fragments
     def __init__(self, base, df, label):
         self.base = base
         self.label = label
         atom_symbols = df["Atom_symbol"] # don't need to keep
         atom_stois = df["N_atoms"]
         atom_groups = df["Atom_group"] # don't need to keep
         self.groups = list(atom_groups.unique())
         self.mass_dict = {g:0.0 for g in self.groups}
         for g,s,stoi in zip(atom_groups, atom_symbols, atom_stois):
                 self.mass_dict[g] += stoi * pgv.atomic_dict[s]["am"]
         self.mass = sum(self.mass_dict.values())
 
class MS2_spectrum:
    def __init__(self, mz1, rt, z, mzarray, intarray):  # is mz1 needed???
        self.mz1 = mz1
        self.rt = rt
        self.z = z
        self.raw_mz = mzarray 
        self.raw_int = intarray
        self.ms2 = {mz: ia for mz, ia in zip(self.raw_mz, self.raw_int)}
    
    def initialize_ms2(self):  # reset ms2 dictionary from original input
        self.ms2 = {mz: ia for mz, ia in zip(self.raw_mz, self.raw_int)}

    def find_peaks(self):  # crude peak picking, implemented for mzML files that are not stick spectra
        idx_list = [i for i in range(1,self.raw_int.size-1) 
                    if self.raw_int[i] > self.raw_int[i+1] and self.raw_int[i] > self.raw_int[i-1] ]
        self.ms2 = {self.raw_mz[i]:self.raw_int[i] for i in idx_list}
        
    def remove_precursor(self, m, window): #window in mz (not ppm)
        pre_mz_list = []
        for mz in self.ms2.keys():
            if abs(mz - m) <= window:
                pre_mz_list.append(mz)
        for mz in pre_mz_list:
            self.ms2.pop(mz,None)
            
    def normalize_intensities(self, scale): 
        self.max_int = float(max(self.ms2.values())) # have to float to avoid "float32" type ???
        factor = scale/self.max_int
        for mz,it in self.ms2.items():
            self.ms2[mz] = it * factor
            
    def threshold(self,thresh): # absolute threshold
        t_list = []
        for mz, it in self.ms2.items():
            if it <= thresh:
                t_list.append(mz)
        for mz in t_list:
            self.ms2.pop(mz,None)
    
    def squash_spectrum(self):
        for mz, it in self.ms2.items():
            self.ms2[mz] = math.sqrt(it)
            
    def bin_spectrum(self, nbins):
        spec_width = pgv.MS2_mzhigh - pgv.MS2_mzlow  #noqa
        delta = spec_width/nbins
        # print("bin_spectrum: ", pgv.MS2_mzhigh, pgv.MS2_mzlow, spec_width, delta)
        self.bspec = [0.0] * nbins
        # bins = [pgv.mz_low + i  * delta for i in range(nbins)]
        for mz, ia in self.ms2.items():  
            idx = int((mz - pgv.MS2_mzlow)/delta)-1 #noqa
            if idx >= nbins:
                continue
                # print(idx)
            else:
                self.bspec[idx] += ia
    
    def fft_spectrum(self):
        self.ft = np.fft.fft(self.bspec)
        # self.ft = np.fft.fftshift(self.ft)
                                                       
    def complex_conjugate(self):
        self.ft = np.conjugate(self.ft)
        
    def fft_MS2(self):
        self.bin_spectrum(4096)
        self.fft_spectrum()
        self.complex_conjugate()
        
    def gaussian_filter(self,width):
        gf = np.array([math.exp(-i**2/(width**2 * len(self.ft)**2)) for i in range(len(self.ft))])
        self.ft = np.multiply(self.ft,gf)

class ms2_ion: # class for theoretical ms2 ions
    def __init__(self, mz2, intensity, series, index, z):
        self.mz2 = mz2
        self.intensity = intensity
        self.series = series
        self.index = index
        self.z = z

class ms2_match: # class for match between theo ion and exptl ion
    def __init__(self, ms2_ion, mz2, intensity, ppmo): # theo ion and observed parameters
        self.theo_mz2 = ms2_ion.mz2
        self.theo_int = ms2_ion.intensity
        self.series = ms2_ion.series
        self.index = ms2_ion.index
        self.z = ms2_ion.z
        self.obs_mz2 = mz2
        self.obs_int = intensity
        self.ppmo = ppmo # ppm_offset


class group: # class for individual fragments in context of oligomer 
    def __init__(self, base, resno, group):
        self.base = base
        self.group = group
        self.resno = resno
        self.parent = pgv.topo_dict[group]['parent']
        self.parent_resno = resno + pgv.topo_dict[group]['parent_resno']
        self.node = "_".join([base,str(resno),group])
        self.frag_left = pgv.topo_dict[group]['fragment_left']
        self.frag_right = pgv.topo_dict[group]['fragment_right']
        self.hcorr_right = pgv.topo_dict[group]['H_corr_right']
        self.hcorr_left = pgv.topo_dict[group]['H_corr_left']
        self.frag_type = pgv.topo_dict[group]['main_side']
        self.frag_offset = pgv.topo_dict[group]['fragment_number_offset']
        self.props = [self.resno, [self.frag_left, self.frag_right], [self.hcorr_left,self.hcorr_right], self.frag_type, self.frag_offset]
    
    def add_node(self,G,node_dict): # each group is a node in the molecular graph
        G.add_node(self.node)
        node_dict[self.node] = self
    
    def add_edge(self,G,parent, edge_dict): # each edge is a fragmentation point in the molecular graph
        if parent in G.nodes:
            v1 = self.node
            v2 = parent
            self.parent_node = parent
            G.add_edge(v1,v2)
            edge_dict[(v1,v2)] = (v1,v2) # alas, but have to keep track of ordered nodes in edge separately
            edge_dict[(v2,v1)] = (v1,v2)


class fragment_sequence: # class to hold various sequence representations of a fragment
    #convention:
        # frag3 is 3-letter code including end groups as list:  primary internal representation
        # seq3  is 3-letter code without end groups as list
        # frag is 1-letter code with brackets for mods, an terminal groups as prefix_ _suffix for output
        # seq is 1-letter code with brackets for mods, no terminal groups for output
        # end5n is 3-letter code for 5'end
        # end5 is 1-letter code for 5'end
        # end3n is 3-letter code for 3'end
        # end3 is 1-letter code for 3'end
       
    def __init__(self,frag3):
        self.frag3 = frag3
        mod_seq_list = []
        fr = 0
        to = len(self.frag3)
        self.end5n, self.end5 = "", ""
        prefix, suffix = "", ""
        self.end3n, self.end3 = "", ""

        for base in self.frag3:
            base_id = pgv.nt_def_dict[base]["Pytheas_ID"]
            base_type = pgv.nt_def_dict[base]["Type"]
            if base == "XXX":
                mod_seq_list.append('X')
            elif base_type == "natural":
                mod_seq_list.append(base_id)
            else:
                if type(base_id)  == float or base_id =="":
                    mod_seq_list.append("[" + base + "]")
                else:
                    mod_seq_list.append("[" + base_id + "]")
                    
            if base_type == "end5":
                self.end5n = base
                self.end5 = base_id
                prefix = self.end5 + "_"
                fr = 1
            if base_type == "end3":
                self.end3n = base
                self.end3 = base_id
                suffix = "_" + self.end3
                to = len(self.frag3) -1
        
        self.modseq = mod_seq_list[fr:to]
        self.seq3 = frag3[fr:to]
        self.seq= "".join(mod_seq_list[fr:to])
        self.frag = prefix + self.seq + suffix

#TODO resolve two different defs of modseq
    def generate_mod_seq(self): # convert 3-letter seq_list to 1-base + [mod] string # digest
        mod_seq_list = []
        # for base in seq_list:
        for base in self.seq3:
            if base == "XXX":
                mod_seq_list.append('X')
            elif pgv.nt_def_dict[base]["Type"] == "natural":
                mod_seq_list.append(pgv.nt_def_dict[base]["Pytheas_1_letter_code"])
            else:
                if type(pgv.nt_def_dict[base]["Pytheas_ID"])  != float:
                    mod_seq_list.append("[" + pgv.nt_def_dict[base]["Pytheas_ID"] + "]")
                else:
                    mod_seq_list.append("[" + base + "]")
        self.modseq = "".join(mod_seq_list)
        
    
class molecule:
    def __init__(self, mdict):
        for key, value in mdict.items():
            setattr(self, key, value)
        
class atom:
    def __init__(self, adict):
        for key, value in adict.items():
            setattr(self, key, value)
            
class residue:
    def __init__(self, nt_dict):
        for key, value in nt_dict.items():
            setattr(self, key, value)