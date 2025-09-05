#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 14:18:29 2024

@author: jrwill
"""

# from itertools import permutations, combinations, product

from pytheas_global_vars import pgv, pgc
from pytheas_objects import fragment_sequence
from digest_functions import generate_mod_seq
from scoring_functions import ppm_range, ppm_offset
from match_functions import scoring, rank_matches
from pytheas_IO import save_pickle
# import digest_functions as dg
# import match_functions as ma
# from pytheas_modules.Combination_sum_class_work import Combination_Sum
# from pytheas_modules.Combination_Product_work import Combination_Product

import networkx as nx
import bisect
import math
import json
import numpy as np
import xlsxwriter
import copy
import pandas as pd
# import gc
from pathlib import Path
from copy import deepcopy

import os
from datetime import datetime
from collections import Counter
import statsmodels.api as sm # recommended import according to the docs
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from collections import deque


#TODO move to pytheas_objects
class MS2_ion_match:
    __slots__ = ("theo_mz2", "theo_int", "series", "index", "z", "obs_mz2", "obs_int", "ppmo")
                   # matches.append({"theo_mz2": mz2,"theo_int": cid["intensity"], "series": cid["series"], 
                   #                     "index": cid["index"], "z": cid["z"], "obs_mz2": match_mass, "obs_int": match_int, 
                   #                     "ppmo": calc_offset})

    def __init__(self, theo_mz2, theo_int, series, index, z, obs_mz2, obs_int, ppmo):
        self.theo_mz2 = theo_mz2
        self.theo_int = theo_int
        self.series = series
        self.index = index
        self.z = z
        self.obs_mz2 = obs_mz2
        self.obs_int = obs_int
        self.ppmo = ppmo
        
        # for key, val in ion_dict:
        #     setattr(self, key, val)

class Combination_Sum:
    def __init__(self, candidates, target, tol):
        self.candidates = candidates
        self.target =  target # mz1
        self.tol = tol  # pgv.MS1_ppm
        self.results = []  # combinations to return
        self.ctr = 0  # counter for recursion
        self.bctr = 0 # counter for bailing
        
    def depth_first_search(self, i, current, total):   # recursive algorithm for combination sum
        # i is pointer index, current is current sequence, total is sum of current, candidates is list of masses, 
        self.ctr +=1
        count = sum([1 for b in current if b in pgv.mod_mass_list]) # bail if too many mods in candidates
        if count > pgv.max_mods:
            self.bctr += 1
            return
        if abs(total - self.target) < self.tol:
            self.results.append(current.copy())
            return
        if i >= len(self.candidates) or total > self.target:
            return
        
        current.append(self.candidates[i])
        self.depth_first_search(i,current, total + self.candidates[i])
        current.pop()
        self.depth_first_search( i+1, current, total)
    
    def combination_sum(self):
        self.depth_first_search(0, [], 0)
        return self.results, self.ctr, self.bctr
# %% Combination_Product


class Combination_Product:
    def __init__(self, node_dict):
        for key, val in node_dict.items():
            setattr(self,key,val)
        # globals:
        #   candidate_factors = target composition of sequence given as prime factors
        #   prime_multiple = target product of candidate_factors
        #   node_node_list = list of all possible node values for each node in the sequence graph
        #   z2_list = list of allowed charge states for MS2_ions
        #   m0 = monoisotopic mass of target
        #   ms2_vals = sorted_list of experimentl mz2 values for matching
        #   ms2_intensities = sorted list of experimental intensities
        #   length
        self.ctr = 0
        self.bctr = 0
        self.nctr = 0
        self.ectr = 0
        self.complete_list = []
        self.parent_dict = {}
        self.depth_score_dict = {}
        # complete_dict
        #   factors_list = current list of prime factors for completed sequences
        #   ms2_ions_list = current list of ms2_ions for completed sequences
        #   matched_ion_list - current list of matched_ions for completed sequences

    def matching_ms2_dfs(self, cid_ions):  #cid ions is a list of ion dicts
     
        matches = []
    
        for cid in cid_ions:
            mz2 = cid["mz2"]
            ppmo = ppm_range(mz2, pgv.MS2_ppm_offset)
            ppm = ppm_range(mz2 + ppmo , pgv.MS2_ppm) 
            lower = mz2 + ppmo - ppm
            upper = mz2 + ppmo + ppm
            idxlo = max(bisect.bisect_left(self.ms2_vals,lower) - 2,0)
            idxhi = min( bisect.bisect_left(self.ms2_vals,upper) + 2, len(self.ms2_vals))
            
            for idx in range(idxlo,idxhi):  # check for ppm tolerance
                match_mass = self.ms2_vals[idx]
                match_int = self.ms2_ints[idx]
                calc_offset = ppm_offset(mz2, match_mass)
                if abs(calc_offset) < pgv.MS2_ppm:
                    matches.append(MS2_ion_match(mz2, cid["intensity"], cid["series"], cid["index"], cid["z"], match_mass, match_int, calc_offset))
                    # matches.append({"theo_mz2": mz2,"theo_int": cid["intensity"], "series": cid["series"], 
                    #                     "index": cid["index"], "z": cid["z"], "obs_mz2": match_mass, "obs_int": match_int, 
                    #                     "ppmo": calc_offset})
        return matches

    def dfs_product_ms2(self, i, current_dict, parent): 
        # i is current pointer
        # current_dict:
        #       current_factors = current list of prime factors for each 
        #       ms2_ions = current list of ms2_ions up to current node
        #       m_list = current list of cumulative masses up to current node
        #       matches = current list of matches thru current node
        
        # self.parent = parent
        cf = current_dict["current_factors"]
        lcf = len(cf)
        ncf = len([c for c in cf if c > 1])
        prod = simple_product(cf)
        
        if prod == self.prime_multiple and lcf == len(self.candidate_factors): # solution!
            self.complete_list.append(copy.deepcopy(current_dict))
            # print("     ", "prime multiple_match")
            return
        
        
            
        if (self.prime_multiple/simple_product(cf)%1 != 0.0 and lcf) != 0:
            self.bctr += 1# terminating by not a factor
            # print("     ", "bad composition")
            return

# either track score by depth, or do breadth first search

        # print()
        # print("current_dict_matches")
        # for m in current_dict["matches"]:
        #     for ml in m:
        #         print(ml)
            
        # print(current_dict["matches"])
        n_matches = len([m for ml in current_dict["matches"] for m in ml if m != []])
        score = pgv.match_intcpt + pgv.match_slope * n_matches
        
        if len(current_dict["current_factors"]) > 0: # branch point
            factor = current_dict["current_factors"][-1]
            # if lcf in self.depth_score_dict:
            #     if score < 0.5 * self.depth_score_dict[lcf]:
            #         # pass
            #         return
            #     if score > self.depth_score_dict[lcf]:
            #         self.depth_score_dict[lcf] = score
            # else:
            #     self.depth_score_dict[lcf] = score
                                         
        else:
            factor = 0

        # check progress of matching
        # if pgv.check_dfs_matching == 'y':
        #     # if i > 7:
        #     #     msum = sum([len(m) for m in current_dict["matches"]])/i
        #     #     if msum < pgv.matching_cutoff:
        #     #         print("hitting eject", i,msum, pgv.matching_cutoff, current_dict["current_factors"])
        #     #         return

        #      if lcf > 0:
        #         # print(current_dict["matches"])
        #         # print("n_matches", i, n_matches, ncf, score)

        #         if n_matches < pgv.match_intcpt +  pgv.match_slope * ncf:
        #             # n_matches = len([m for m in current_dict["matches"] if m != []])
        #             print(" terminating poor match", i, ncf, n_matches)
        #             self.ectr += 1
        #             return
        
        self.ctr += 1
        self.parent_dict[self.ctr] = {"parent": parent, "index": i, "n_matches": n_matches, "depth": ncf, "score": score, "factor": factor}
        parent = self.ctr

    
        for j in range(len(self.candidate_factors[i])):
            # if len(self.candidate_factors[i]) > 1:
            #        print(i,  sum([len(m) for m in current_dict["matches"]]), self.candidate_factors[i][j])
            #        print(current_dict["current_factors"])
            # # add candidate_node
            current_dict["current_factors"].append(self.candidate_factors[i][j])
            mtot, cid_ions = self.next_node(current_dict["m_list"][i], self.node_list_list[i][j])
            matched_ions = self.matching_ms2_dfs(cid_ions)
            current_dict["ms2_ions"].append(cid_ions)
            current_dict["m_list"].append(mtot)
            current_dict["matches"].append(matched_ions)
            # self.parent_dict[self.ctr] = self.parent
            #recursion...next node
            self.dfs_product_ms2(i+1, current_dict, parent)
            # reset to try next candidate at position i
            current_dict["current_factors"].pop()
            current_dict["ms2_ions"].pop()
            current_dict["m_list"].pop()
            current_dict["matches"].pop()
            
    def next_node(self, mtot, node):
        nd = self.G.nodes[node] # node dict 
        mtot += nd["group_mass"] # add mass of new node
        cid_ions = []
        fl = nd["fragment_left"]
        fr = nd["fragment_right"]
        ml = mtot + nd["H_corr_left"] * pgv.hmass
        mr = self.m0 - ml + nd["H_corr_right"] * pgv.hmass

    #   left fragment
        if fl != "none":
            if nd["main_side"] == "S": # side-chain does not split backbone
                left = nd["resno"]
            else:
                left =  nd["resno"]+ nd["fragment_left_offset"]
            for zs in self.z2_list:
                mz2 = (ml + zs * pgv.hmass)/abs(zs)
                if pgv.MS2_mzlow < mz2 < pgv.MS2_mzhigh: 
                    cid_ions.append({"mz2": mz2, "intensity": 1.0, "series": fl, "index": left, "z": zs})
            
    #   right_fragment
        if fr != "none":
            if nd["main_side"] == "S": # side-chain does not split backbone
                right = nd["resno"]
            else:
                right = self.length - nd["resno"] + nd["fragment_right_offset"]  # right ions are n-resno
            for zs in self.z2_list:
                mz2 = (mr + zs * pgv.hmass)/abs(zs)
                if pgv.MS2_mzlow < mz2 < pgv.MS2_mzhigh: 
                    cid_ions.append({"mz2": mz2, "intensity": 1.0, "series": fr, "index": right, "z": zs})
                    if fr == 'y':  # add y-P loss
                        mz2 = (mr - pgv.losses_dict['y-P']["loss_mass"] + zs * pgv.hmass)/abs(zs)
                        cid_ions.append({"mz2": mz2, "intensity": 1.0, "series": 'y-P', "index": right, "z": zs})

                    if fr == 'z':  # add z-P loss
                        mz2 = (mr - pgv.losses_dict['z-P']["loss_mass"] + zs * pgv.hmass)/abs(zs)
                        cid_ions.append({"mz2": mz2, "intensity": 1.0, "series": 'z-P', "index": right, "z": zs})

        self.nctr += 1
        return mtot, cid_ions

    def combination_product(self):
    
        self.dfs_product_ms2(0, {"current_factors": [], "ms2_ions": [], "m_list": [0], "matches": []}, 0)
        
        # print("parent_dict", self.parent_dict)
        
        return self.complete_list, self.ctr, self.bctr, self.nctr, self.ectr, self.parent_dict
# %% Combination_Product_BFS

class Combination_Product_BFS:
    def __init__(self, node_dict):
        for key, val in node_dict.items():
            setattr(self,key,val)
        # globals:
        #   candidate_factors = target composition of sequence given as prime factors
        #   prime_multiple = target product of candidate_factors
        #   node_node_list = list of all possible node values for each node in the sequence graph
        #   z2_list = list of allowed charge states for MS2_ions
        #   m0 = monoisotopic mass of target
        #   ms2_vals = sorted_list of experimentl mz2 values for matching
        #   ms2_intensities = sorted list of experimental intensities
        #   length
        self.ctr = 0
        self.bctr = 0
        self.nctr = 0
        self.ectr = 0
        self.complete_list = []
        self.parent_dict = {}
        self.depth_score_dict = {}
        self.queue = []
        # complete_dict
        #   factors_list = current list of prime factors for completed sequences
        #   ms2_ions_list = current list of ms2_ions for completed sequences
        #   matched_ion_list - current list of matched_ions for completed sequences

    def matching_ms2_dfs(self, cid_ions):  #cid ions is a list of ion dicts
     
        matches = []
    
        for cid in cid_ions:
            mz2 = cid["mz2"]
            ppmo = ppm_range(mz2, pgv.MS2_ppm_offset)
            ppm = ppm_range(mz2 + ppmo , pgv.MS2_ppm) 
            lower = mz2 + ppmo - ppm
            upper = mz2 + ppmo + ppm
            idxlo = max(bisect.bisect_left(self.ms2_vals,lower) - 2,0)
            idxhi = min( bisect.bisect_left(self.ms2_vals,upper) + 2, len(self.ms2_vals))
            
            for idx in range(idxlo,idxhi):  # check for ppm tolerance
                match_mass = self.ms2_vals[idx]
                match_int = self.ms2_ints[idx]
                calc_offset = ppm_offset(mz2, match_mass)
                if abs(calc_offset) < pgv.MS2_ppm:
                    matches.append(MS2_ion_match(mz2, cid["intensity"], cid["series"], cid["index"], cid["z"], match_mass, match_int, calc_offset))
                    # matches.append({"theo_mz2": mz2,"theo_int": cid["intensity"], "series": cid["series"], 
                    #                     "index": cid["index"], "z": cid["z"], "obs_mz2": match_mass, "obs_int": match_int, 
                    #                     "ppmo": calc_offset})
        return matches

    def bfs_product_ms2(self, i): 
        # i is current pointer = level pointer
        
        # current dict set is in self.queue
        # current_dict:
        #       current_factors = current list of prime factors for each 
        #       ms2_ions = current list of ms2_ions up to current node
        #       m_list = current list of cumulative masses up to current node
        #       matches = current list of matches thru current node
        
        # self.parent = parent
        level_factors = self.candidate_factors[i]
        
        level_queue = []
        
            
        if len(level_factors) == 1:
            for current_queue_dict in self.queue:

                self.ctr += 1
                current_dict = deepcopy(current_queue_dict)
                current_dict["parent"] = current_dict["node_idx"]
                current_dict["current_factors"].append(level_factors[0])
                cf = current_dict["current_factors"]
                lcf = len(cf)
                mtot, cid_ions = self.next_node(current_dict["m_list"][i], self.node_list_list[i][0])
                matched_ions = self.matching_ms2_dfs(cid_ions)
                current_dict["ms2_ions"].append(cid_ions)
                current_dict["m_list"].append(mtot)
                current_dict["matches"].append(matched_ions)
                n_matches = len([m for ml in current_dict["matches"] for m in ml if m != []])
                score = pgv.match_intcpt + pgv.match_slope * n_matches
                current_dict["node_idx"] = self.ctr
                
                # score_list.append(score)
                # current_dict_list.append(deepcopy(current_dict))
                # current_dict = deepcopy(current_layer_dict)
                level_queue.append(deepcopy(current_dict))
                factor = level_factors[0]
                self.parent_dict[self.ctr] = {"parent": current_dict["parent"], "index": i, "n_matches": n_matches, "depth": i, "score": score, "factor": factor}
                    # parent = self.ctr

            
        else:
            score_list = []
            branch_list = []
            match_list = []
            fac_list = []
            for current_queue_dict in self.queue:
                for j in range(len(level_factors)):
                    current_dict = deepcopy(current_queue_dict)
                    current_dict["parent"] = current_dict["node_idx"]
                    # print("branching ", j)
                    current_dict["current_factors"].append(self.candidate_factors[i][j])
                    cf = current_dict["current_factors"]
                    lcf = len(cf)
                    
                    # BAIL!
                    if (self.prime_multiple/simple_product(cf)%1 != 0.0 and lcf) != 0:
                        self.bctr += 1# terminating by not a factor
                        # print("     ", "BAIL: bad composition layer ", i, "factor ", j, " = ", level_factors[j])
                        # print(cf)
                        continue
                    
                    # # DONE!
                    # if simple_product(cf) == self.prime_multiple and lcf == len(self.candidate_factors): # solution!
                    #      self.complete_list.append(copy.deepcopy(current_dict))
                    #      print("     ", "DONE: prime multiple_match")
                    #      # continue
                
    
                    mtot, cid_ions = self.next_node(current_dict["m_list"][i], self.node_list_list[i][j])
                    matched_ions = self.matching_ms2_dfs(cid_ions)
                    current_dict["ms2_ions"].append(cid_ions)
                    current_dict["m_list"].append(mtot)
                    current_dict["matches"].append(matched_ions)
                    n_matches = len([m for ml in current_dict["matches"] for m in ml if m != []])
                    score = pgv.match_intcpt + pgv.match_slope * n_matches
                    score_list.append(score)
                    match_list.append(n_matches)
                    fac_list.append(level_factors[j])
                    # print("scores ", score, score_list)
                    branch_list.append(deepcopy(current_dict))
            if len(score_list) != 0:
                max_score = max(score_list)
                scores_sorted = sorted(score_list, reverse = True)
                if len(score_list) > pgv.max_branches:
                    score_cut = scores_sorted[pgv.max_branches]
                else:
                    score_cut = pgv.match_factor * max_score
                for score, match, factor, branch in zip(score_list, match_list, fac_list, branch_list):
                    if score >= score_cut:
                        self.ctr += 1
                        branch["node_idx"] = self.ctr
                        # factor = self.candidate_factors[i][j]
                        self.parent_dict[self.ctr] = {"parent": branch["parent"], "index": i, "n_matches": match, "depth": i, "score": score, "factor": factor}
                        level_queue.append(branch)

        
        self.queue = deepcopy(level_queue)
        
       
        print("level ", i, " has ", len(self.queue), " branches")

                
                    
                    
            
#         cf = current_dict["current_factors"]
#         lcf = len(cf)
#         ncf = len([c for c in cf if c > 1])
#         prod = simple_product(cf)
        
#         # print("bfs_product_ms2: i, current_dict,parent", i, current_dict, parent)
#         # print("prod", prod, self.prime_multiple)
#         # DONE!
#         if prod == self.prime_multiple and lcf == len(self.candidate_factors): # solution!
#             self.complete_list.append(copy.deepcopy(current_dict))
#             # print("     ", "DONE: prime multiple_match")
#             return
        
        
        
#         # BAIL!
#         if (self.prime_multiple/simple_product(cf)%1 != 0.0 and lcf) != 0:
#             self.bctr += 1# terminating by not a factor
#             # print("     ", "BAIL: bad composition layer ", i)
#             return

# # either track score by depth, or do breadth first search

#         # print()
#         # print("current_dict_matches")
#         # for m in current_dict["matches"]:
#         #     for ml in m:
#         #         print(ml)
            
#         # print(current_dict["matches"])
#         n_matches = len([m for ml in current_dict["matches"] for m in ml if m != []])
#         score = pgv.match_intcpt + pgv.match_slope * n_matches
        
#         if len(current_dict["current_factors"]) > 0: # branch point
#             factor = current_dict["current_factors"][-1]
#             # if lcf in self.depth_score_dict:
#             #     if score < 0.5 * self.depth_score_dict[lcf]:
#             #         # pass
#             #         return
#             #     if score > self.depth_score_dict[lcf]:
#             #         self.depth_score_dict[lcf] = score
#             # else:
#             #     self.depth_score_dict[lcf] = score
                                         
#         else:
#             factor = 0
        
        
        
#         self.ctr += 1
#         self.parent_dict[self.ctr] = {"parent": parent, "index": i, "n_matches": n_matches, "depth": ncf, "score": score, "factor": factor}
#         parent = self.ctr

#         current_layer_dict = deepcopy(current_dict)
#         # current_factor_list = current_factors
#         # initial_factor_list = self.candidate_factors[i]
#         score_list = []
#         current_dict_list = []
#         n_matches = len([m for ml in current_dict["matches"] for m in ml if m != []])
#         print("before branching, n_matches = ", n_matches)

#         for j in range(len(self.candidate_factors[i])): # go thru all possible nodes
#             print("layer iteration ", i, j, self.candidate_factors[i][j])
#             # if len(self.candidate_factors[i]) > 1:
#             #        print(i,  sum([len(m) for m in current_dict["matches"]]), self.candidate_factors[i][j])
#             #        print(current_dict["current_factors"])
#             # # add candidate_node
#             # print("layer current_factors before append", j, current_dict["current_factors"])
#             current_dict["current_factors"].append(self.candidate_factors[i][j])
#             cf = current_dict["current_factors"]
#             lcf = len(cf)
            
#             # BAIL!
#             if (self.prime_multiple/simple_product(cf)%1 != 0.0 and lcf) != 0:
#                 self.bctr += 1# terminating by not a factor
#                 # print("     ", "BAIL: bad composition layer ", i)
#                 continue

#             # print("layer current_factors after append", j, current_dict["current_factors"])
#             # current_factors.append(self.candidate_factors[i][j])
#             mtot, cid_ions = self.next_node(current_dict["m_list"][i], self.node_list_list[i][j])
#             matched_ions = self.matching_ms2_dfs(cid_ions)
#             current_dict["ms2_ions"].append(cid_ions)
#             current_dict["m_list"].append(mtot)
#             current_dict["matches"].append(matched_ions)
#             n_matches = len([m for ml in current_dict["matches"] for m in ml if m != []])
#             score = pgv.match_intcpt + pgv.match_slope * n_matches
#             score_list.append(score)
#             current_dict_list.append(deepcopy(current_dict))
#             current_dict = deepcopy(current_layer_dict)
#             # print("layer current_factors after restore", j, current_dict["current_factors"])
#             # print("layer current_factors", j, current_dict["current_factors"])
#         if len(score_list) == 0:
#             return
#             # current_factors = current_factor_list
#         max_score = max(score_list)
#         keep_dict_list = [c for s,c in zip(score_list, current_dict_list) if s >= max_score/2]
#         if len(score_list) > 1: 
#             print("score_list", score_list)
#             print("keeping ", len(keep_dict_list), " of ", len(score_list))
#         for c in keep_dict_list:

#             # self.parent_dict[self.ctr] = self.parent
#             #recursion...next node
#             # self.ctr += 1
#             # n_matches = len([m for ml in c["matches"] for m in ml if m != []])
#             # score = pgv.match_intcpt + pgv.match_slope * n_matches
#             # factor = c["current_factors"][-1]
#             # cf = c["current_factors"]
#             # ncf = len([c for c in cf if c > 1])
#             # self.parent_dict[self.ctr] = {"parent": parent, "index": i, "n_matches": n_matches, "depth": ncf, "score": score, "factor": factor}
#             # parent = self.ctr
#             # print(i,c["current_factors"], parent, score, factor)
#             print("going to next level", i, parent)
#             self.bfs_product_ms2(i + 1, c, parent)
#             # # reset to try next candidate at position i
#             # current_dict["current_factors"].pop()
#             # current_dict["ms2_ions"].pop()
#             # current_dict["m_list"].pop()
#             # current_dict["matches"].pop()
            
    def next_node(self, mtot, node):
        nd = self.G.nodes[node] # node dict 
        mtot += nd["group_mass"] # add mass of new node
        cid_ions = []
        fl = nd["fragment_left"]
        fr = nd["fragment_right"]
        ml = mtot + nd["H_corr_left"] * pgv.hmass
        mr = self.m0 - ml + nd["H_corr_right"] * pgv.hmass

    #   left fragment
        if fl != "none":
            if nd["main_side"] == "S": # side-chain does not split backbone
                left = nd["resno"]
            else:
                left =  nd["resno"]+ nd["fragment_left_offset"]
            for zs in self.z2_list:
                mz2 = (ml + zs * pgv.hmass)/abs(zs)
                if pgv.MS2_mzlow < mz2 < pgv.MS2_mzhigh: 
                    cid_ions.append({"mz2": mz2, "intensity": 1.0, "series": fl, "index": left, "z": zs})
            
    #   right_fragment
        if fr != "none":
            if nd["main_side"] == "S": # side-chain does not split backbone
                right = nd["resno"]
            else:
                right = self.length - nd["resno"] + nd["fragment_right_offset"]  # right ions are n-resno
            for zs in self.z2_list:
                mz2 = (mr + zs * pgv.hmass)/abs(zs)
                if pgv.MS2_mzlow < mz2 < pgv.MS2_mzhigh: 
                    cid_ions.append({"mz2": mz2, "intensity": 1.0, "series": fr, "index": right, "z": zs})
                    if fr == 'y':  # add y-P loss
                        mz2 = (mr - pgv.losses_dict['y-P']["loss_mass"] + zs * pgv.hmass)/abs(zs)
                        cid_ions.append({"mz2": mz2, "intensity": 1.0, "series": 'y-P', "index": right, "z": zs})

                    if fr == 'z':  # add z-P loss
                        mz2 = (mr - pgv.losses_dict['z-P']["loss_mass"] + zs * pgv.hmass)/abs(zs)
                        cid_ions.append({"mz2": mz2, "intensity": 1.0, "series": 'z-P', "index": right, "z": zs})

        self.nctr += 1
        return mtot, cid_ions

    def combination_product(self):
        
        current_dict = {"current_factors": [], "ms2_ions": [], "m_list": [0], "matches": [], "node_idx":0, "parent": 0}
        self.queue.append(current_dict)
        parent = 0

        for level in range(len(self.candidate_factors)):
    
            self.bfs_product_ms2(level)
            parent += 1
        
        # print("final queue")
        # for current_dict in self.queue:
        #     print(current_dict["current_factors"])
            
        # print("complete_list")
        # for cl in self.complete_list:
        #     print(cl)
        # print("parent_dict", self.parent_dict)
        
        print()
        print("***********done with bfs******************")
        print()
        

        
        return self.queue, self.ctr, self.bctr, self.nctr, self.ectr, self.parent_dict
# %%


def find_closest_index(a, x): # a is sorted low to hi   NOT NEEDED?
    i = bisect.bisect_left(a, x)
    if i >= len(a):
        i = len(a) - 1
    elif i and a[i] - x > x - a[i - 1]:
        i = i - 1
    return (i, a[i])

def build_mass_dict():
# def build_mass_dict(mass_dict, base_list):
   
    print("discovery: modification_set", pgv.modification_set, type(pgv.modification_set))
    base_list = [b for s in  pgv.modification_set for b in pgv.set_dict[s] ]
    print("mod list", base_list)
    
    # #TODO    # need to loop thru labels
    mass_dict = {key:val.mass for key,val in pgv.nt_fragment_dict["light"].items() if "end" not in pgv.nt_def_dict[key]["Type"]}

    
    sub_mass_dict = {base:mass for base, mass in mass_dict.items() if base in base_list}

    tol = .01
    for bi, mi in sub_mass_dict.items():
        for bj,mj in sub_mass_dict.items():
            if bi == bj:
                continue
            if abs(mi - mj) < tol:
                print(bi,bj, abs(mi-mj))
            
    pgv.iso_mass_list = sorted(list(set([round(m,3) for b,m in sub_mass_dict.items()]))   )       
    pgv.iso_mass_dict = {m:[] for m in pgv.iso_mass_list}

    for bi,mi in sub_mass_dict.items():
        mkey = round(mi,3)
        if mkey in pgv.iso_mass_dict.keys():
            pgv.iso_mass_dict[mkey].append(bi)
    
    pgv.mod_mass_list = []
    for key, blist in pgv.iso_mass_dict.items():
        for b in blist:
            print(key, b)
            if b in pgc.natural:
                break
            if key not in pgv.mod_mass_list:
                print("addin", key)
                pgv.mod_mass_list.append(key)


    print("iso_mass_dict", pgv.iso_mass_dict)
    print("iso_mass_list", pgv.iso_mass_list)
    print("mod_mass_list", pgv.mod_mass_list)

def expand_seq_dfs(i, current, candidates, sequences):   # recursive algorithm for expanding combinations of sequences
    # i is pointer index, current is current sequence, candidates is list of nucleotides, sequences has results
    if i == len(candidates):
        sequences.append(current.copy())
        return
    if i > len(candidates): #just in case
        print("too long")
        return
    
    for j in range(len(candidates[i])):
        current.append(candidates[i][j])
        expand_seq_dfs(i+1, current, candidates, sequences)
        current.pop()

def find_compositions(mobs, mobs_adj): # find base compositions within tol of mobs

#TODO bail if pgv.max_mods exceeded
    tc = datetime.now()

    ppm_tol = pgv.MS1_ppm
    tol = mobs * ppm_tol/1000000.0
    ctr = 0
    nt_mass_list = list(pgv.iso_mass_dict.keys())
    
    combo_sums, ctr, bctr = Combination_Sum(nt_mass_list, mobs_adj, tol).combination_sum()

    seq_comps = []
    for c in combo_sums:
        mseq = [pgv.iso_mass_dict[m] for m in c] # turn numbers into sequences
        expanded_sequences = []
        expand_seq_dfs(0,[],mseq,expanded_sequences)
        seq_comps.extend(expanded_sequences)
  
    print("find_compositions took " , ctr, "iterations", datetime.now() - tc, "with mod_bails ", bctr)
    
    return seq_comps, ppm_tol, tol

def generate_molecular_graph(f3, label):   # same as digest_functions
    G = nx.Graph() # initialize molecular graph
    
    for resno in range(len(f3.frag3)): # add node to graph for each group
        for g in pgv.nt_fragment_dict[label][f3.frag3[resno]].groups:
            node = "_".join([f3.frag3[resno], str(resno), g])
            G.add_node(node, base = f3.frag3[resno], resno = resno, group = g)
            
    for node in G.nodes:  # fill out graph nodes with data for fragmentation from topo dict
        nd = G.nodes[node]
        res =nd["resno"]
        base = nd["base"]
        nd.update(pgv.child_dict[nd["group"]]) # add topology info
        parent_group, child_group = nd["added_edge"].split("_")
        child_res = res + nd["child_offset"]
        child_base = f3.frag3[child_res]

        if nd["parent_offset"] == 0: # standard linkage
            parent_node = "_".join([base, str(res), parent_group])
            child_node = "_".join([child_base, str(child_res), child_group])
        else: # 3'-terminal linkage
            parent_res = res + nd["parent_offset"] 
            parent_base = f3.frag3[parent_res]
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


def build_node_list(f3, label): #build list of lists for nodes  ## NOT NEEDED??

    base_list, base_node_list, base_factors_list = [], [], []
    
    G = generate_molecular_graph(f3, label)
    nb = 0
    for node in G.nodes:
        G.nodes[node]["group_mass"] = pgv.nt_fragment_dict["light"][G.nodes[node]["base"]].mass_dict[G.nodes[node]["group"]]
        if G.nodes[node]["group"] == 'B':
            if G.nodes[node]["base"] not in base_list:
                base_list.append(G.nodes[node]["base"])
                base_node_list.append(node)
                base_factors_list.append(pgc.primes[nb])
                nb +=1      
        
    node_list_list = []
    factors_list = []
    prime_multiple = 1
    
    for node in G.nodes:
        if G.nodes[node]["group"] == "B":
            node_list_list.append(base_node_list)
            prime_multiple *= pgc.primes[base_list.index(G.nodes[node]["base"])]
            factors_list.append(base_factors_list)
        else:
            node_list_list.append([node])
            factors_list.append([1])
    
    return G, base_list, base_factors_list, node_list_list, factors_list, prime_multiple

def build_node_dict(f3, label): #build list of lists for nodes
    
    base_list, base_node_list, base_factors_list = [], [], []
    
    G = generate_molecular_graph(f3, label)
    nb = 0
    for node in G.nodes:
        G.nodes[node]["group_mass"] = pgv.nt_fragment_dict["light"][G.nodes[node]["base"]].mass_dict[G.nodes[node]["group"]]
        if G.nodes[node]["group"] == 'B':
            if G.nodes[node]["base"] not in base_list:
                base_list.append(G.nodes[node]["base"])
                base_node_list.append(node)
                base_factors_list.append(pgc.primes[nb])
                nb +=1      
        
    node_list_list = []
    factors_list = []
    prime_multiple = 1
    
    for node in G.nodes:
        if G.nodes[node]["group"] == "B":
            node_list_list.append(base_node_list)
            prime_multiple *= pgc.primes[base_list.index(G.nodes[node]["base"])]
            factors_list.append(base_factors_list)
        else:
            node_list_list.append([node])
            factors_list.append([1])
            
    node_dict = {"G": G, "base_list": base_list, "base_factors_list": base_factors_list,
                 "node_list_list": node_list_list, "candidate_factors": factors_list, 
                 "prime_multiple": prime_multiple}
    
    # return G, base_list, base_factors_list, node_list_list, factors_list, prime_multiple
    return node_dict


def calculate_mz1(G, label, z):

    m0 = 0.0
    for node in G.nodes:
        base = G.nodes[node]["base"]
        grp = G.nodes[node]["group"]
        m0 += pgv.nt_fragment_dict[label][base].mass_dict[grp] 
   
    if pgv.ion_mode == "-":  
        zsign = -1
    else:
        zsign = 1
        
    zs = abs(z)*zsign
    mz1 = (m0 + zs * pgv.hmass)/abs(z)
 
    return m0, mz1

def simple_product(vals):
    p = 1
    for v in vals:
        p *= v
    return p

def next_node(G, mtot, node, z2_list, m0, length): # not needed, but keep to incorporate into digest/match
    # global nctr
    nd = G.nodes[node] # node dict 
    # base = nd["base"]
    # grp = nd["group"]
    # mgrp = pgv.nt_fragment_dict["light"][base].mass_dict[grp]
    # mtot += mgrp
    # mgrp = pgv.nt_fragment_dict["light"][base].mass_dict[grp]
    mtot += nd["group_mass"] # add mass of new node
    cid_ions = []
    fl = nd["fragment_left"]
    fr = nd["fragment_right"]
    ml = mtot + nd["H_corr_left"] * pgv.hmass
    mr = m0 - ml + nd["H_corr_right"] * pgv.hmass

#   left fragment
    if fl != "none":
        if nd["main_side"] == "S": # side-chain does not split backbone
            left = nd["resno"]
        else:
            left =  nd["resno"]+ nd["fragment_left_offset"]
        for zs in z2_list:
            mz2 = (ml + zs * pgv.hmass)/abs(zs)
            if pgv.MS2_mzlow < mz2 < pgv.MS2_mzhigh: 
                cid_ions.append({"mz2": mz2, "intensity": 1.0, "series": fl, "index": left, "z": zs})
        
#   right_fragment
    if fr != "none":
        if nd["main_side"] == "S": # side-chain does not split backbone
            right = nd["resno"]
        else:
            right = length - nd["resno"] + nd["fragment_right_offset"]  # right ions are n-resno
        for zs in z2_list:
            mz2 = (mr + zs * pgv.hmass)/abs(zs)
            if pgv.MS2_mzlow < mz2 < pgv.MS2_mzhigh: 
                cid_ions.append({"mz2": mz2, "intensity": 1.0, "series": fr, "index": right, "z": zs})
                if fr == 'y':  # add y-P loss
                    mz2 = (mr - pgv.losses_dict['y-P']["loss_mass"] + zs * pgv.hmass)/abs(zs)
                    cid_ions.append({"mz2": mz2, "intensity": 1.0, "series": 'y-P', "index": right, "z": zs})

                if fr == 'z':  # add z-P loss
                    mz2 = (mr - pgv.losses_dict['z-P']["loss_mass"] + zs * pgv.hmass)/abs(zs)
                    cid_ions.append({"mz2": mz2, "intensity": 1.0, "series": 'z-P', "index": right, "z": zs})

    # nctr += 1
    return mtot, cid_ions


def repack_current_factors(factors_list, base_list, base_factor_list):
    seq = []
    for f in factors_list:
        if f == 1:
            continue
        seq.append(base_list[base_factor_list.index(f)])
        
    return seq

#TODO redo with list comprehension
def repack_ion_list(ms2_ions_list):  #
    ilist = []
    idict = {}
    idx = 0
    for i in ms2_ions_list:
        ilist.extend(i)
    for i in ilist:
        idict[idx] = i
        idx +=1
    idict_sorted = dict(sorted(idict.items(), key=lambda item: item[1]["mz2"]))
    return idict_sorted

def repack_match_list(ms2_ions_list):
    # ilist = []
    # idict = {}
    # idx = 0
    # for i in ms2_ions_list:
    #     ilist.extend(i)
    # for i in ilist:
    #     idict[idx] = i
    #     idx +=1
    # idict_sorted = dict(sorted(idict.items(), key=lambda item: item[1]["obs_mz2"]))
    alist = ["theo_mz2", "theo_int", "series", "index", "z", "obs_mz2", "obs_int", "ppmo"]
    ilist = [ion for ion_list in ms2_ions_list for ion in ion_list if ion_list != []]
    idict = {idx:ilist[idx] for idx in range(len(ilist)) }
    idict_ms2_sorted = dict(sorted(idict.items(), key=lambda item: item[1].obs_mz2))
    idict_sorted = {}
    for key, val in idict_ms2_sorted.items():
        mdict = {a: getattr(val, a) for a in alist}
        # print(key, mdict)
        idict_sorted[key] = mdict
        
    
    return idict_sorted

def score_discovery(ms2_key, prec_dict, ion_series): # Add the Sp score to the precursor ion matches
#TODO  got rid of "all" option  fix!
    if ion_series == ['all']:
         ion_series_list = pgv.all_CID_series
    else:
        ion_series_list = ion_series
    thresh = float(pgv.MS2_peak_threshold)
    sumi_all, max_int = sumI_all(pgv.ms2_dict[ms2_key], thresh) # should this be the theo digest???? NO!

    for key, precursor in prec_dict.items():
        
        # pdict = 
        precursor["sumi"] = sumI(precursor["matched_ions"])
        precursor["sumi_all"] = sumi_all
        precursor["n"] = n_calc(precursor["matched_ions"], [], thresh)
        precursor["L"] = L_calc(precursor["CID_ions"], [], thresh)
        precursor["beta_list"] = [round(consecutive_series(precursor["matched_ions"],s),3) for s in ion_series_list]
        precursor["beta"] = 1.0 + sum(precursor["beta_list"])
#TODO put alpha in this eqn
        precursor["alpha"] = pgv.alpha
        if precursor["sumi_all"] == 0 or  precursor["L"] == 0:
            calc_score = 0
        else:
            calc_score = precursor["sumi"] * precursor["n"] * precursor["beta"]/(precursor["sumi_all"] * precursor["L"])
        precursor["Sp"] = calc_score
        precursor["score_details"] = [":".join([sk, str(precursor[sk])]) for sk in pgv.score_keys]

        precursor["max_int"] = max_int

        # calculate Xcorr score
        tspec = calc_theo_spectrum(precursor)
        corr_fft = np.multiply(tspec, pgv.ms2_dict[ms2_key].ft)
        corr_ift = np.fft.ifft(corr_fft)
        corr = np.real(corr_ift)
        corr = np.fft.fftshift(corr)
        xcorr = max(corr)/np.mean(corr[pgc.np2 - pgv.xcorr_avg_width: pgc.np2 - pgv.xcorr_excl_width] + corr[pgc.np2 + pgv.xcorr_excl_width: pgc.np2 + pgv.xcorr_avg_width])
        precursor["Xc"] = xcorr
        precursor["Xc_sc"] = xcorr * pgv.xcorr_length_scale * precursor["length"]


def score_rank_discovery_dfs(ms2_key, precursor_dict): # precursor_dict built from combos that match precursor mz1
    
    #no matching...aleady done by matching_ms2_dfs
    # t1 = datetime.now()
    scoring(ms2_key, precursor_dict, pgv.CID_series)  # add Sp scores to ms2_match_dict
    # print("score took " ,datetime.now() - t1)

    # t1 = datetime.now()
    rank_matches(precursor_dict, pgv.ntop) # all matches
    # print("rank took ",datetime.now() - t1)

    top_keys = [key for key in list(precursor_dict.keys()) if "Sp_rank" in precursor_dict[key]]
    top_sequences = [ [key, precursor_dict[key]["frag"], round(precursor_dict[key]["Sp"],3), precursor_dict[key]["Sp_rank"]] for key in top_keys]

    return top_sequences

def score_distribution_plot(sp, score, file):

    hfont = {'fontname':'Helvetica',
            'size'   : 10}
    
    histopts = {"edgecolor": "black", "linewidth": 0.25, "alpha": 0.2}
                
    if len(sp) < 5:
        return
    spa = np.asarray(sorted(sp))
    ecdf = sm.distributions.ECDF(spa)

    x = np.linspace(0,max(spa),num = 100)
    y = ecdf(spa)
    
    try:
        fit_pars, covar = curve_fit(evd_cdf,spa,y, p0 = [0.2,0.2,0.2])
    except:
        print("ECDF curve fit failed, skipping ", file)
        return
    
    mu, sig, eps = fit_pars
    print("fit parameters ", fit_pars)
    
    cdf = evd_cdf(x,mu,sig, eps)
    pdf = evd_pdf(x,mu,sig, eps)
    
    top_p = 1.0 - evd_cdf(spa[-1], mu, sig, eps)
    top_p2 = 1.0 - evd_cdf(spa[-2], mu, sig, eps)

    nSp = len(spa)
    
    # cdf_n = evd_cdf(x,mu,sig, eps)**nSp
    # cdf_n1 = nSp * evd_cdf(x,mu,sig, eps)**(nSp-1)*(1-evd_cdf(x,mu,sig, eps)) + evd_cdf(x,mu,sig,eps)**nSp
    
    fig, (ax1,ax2) = plt.subplots(1, 2)
    fig.set_size_inches(6,3)
    
    #PDF subplot
    ax1.hist(spa, density=True, bins=int(math.sqrt(len(spa))), color = "blue", label = "histogram", **histopts)
    ax1.set_xlim([spa[0], spa[-1]])
    ax1.legend(loc='best', frameon=False)

    ax1.plot(x,pdf, 'k-', lw=1, label='fitted pdf')
    ax1.set_xlabel(score, **hfont)
    ax1.set_ylabel('EPDF', **hfont)
    ax1.legend()
    
    #CDF_subplot
    ax2.plot(spa, y, 'r-', lw = 1, label = "empirical distribution")
    ax2.plot(x, cdf, 'ko', lw=1, markersize = 0.2, label='fitted cdf')
    ax2.legend(loc='best', frameon=False)
    ax2.plot(spa[-1],evd_cdf(spa[-1],mu,sig, eps), 'ro', label = "top Sp = " + str(round(spa[-1],3)) + " p* = " +  str(round(top_p, 4)))
    ax2.plot(spa[-2],evd_cdf(spa[-2],mu,sig, eps), 'kv', label = "2nd Sp = " + str(round(spa[-2],3)) + " p* = " + str(round(top_p2, 4)))
    # ax2.plot(x, cdf_n, 'g-', lw = 1, label = "nth order statistic cdf")
    # ax2.plot(x, cdf_n1, 'g--', lw = 1, label = "n-1th order statistic cdf")
    
    ax2.set_xlabel(score, **hfont)
    ax2.set_ylabel('CDF', **hfont)
    ax2.legend(loc="center right")
    ms2_key, seq, Sp, score, _ = os.path.basename(file).split("_")
    fig.suptitle("#" + ms2_key + " " + seq + " " + score + ' dist: mu = ' + str(round(mu,2)) + " sig = " + str(round(sig,2)) + " eps = " + str(round(eps,2)) + " n = " + str(nSp), **hfont)
    print(" saving to ", file)
    plt.savefig(file)
    plt.close(fig)

# TODO either skip or incorporate into score distribution plot...
def rank_score_plot(sp, score, file):

    hfont = {'fontname':'Helvetica',
            'size'   : 10}
    
    if len(sp) < 5:
        return
    sps = sorted(sp)
    rank = [1 - (i+1)/len(sp) for i in range(len(sp))]
    
    top = sps[-1]
    second = sps[-2]
    delta = top - second
    pe = 1.0/(len(sp) + 1)
    fig, ax = plt.subplots()
    fig.set_size_inches(6,3)
    
    ax.plot(sps, rank, 'k-', lw = 1, label = "rank")
    ax.plot(sps[-1],rank[-1], 'ro', label = "best")
    ax.plot(sps[-2],rank[-2], 'kv', label = "2nd best")
    
    ax.set_xlabel(score, **hfont)
    ax.set_ylabel(score + ' rank/N', **hfont)
    ax.legend(loc="upper right")
    ms2_key, seq, _, _, _ = os.path.basename(file).split("_")
    fig.suptitle("ms2_key = " + ms2_key + " " + seq +  " N = " + str(len(sp)) + " delta = " + str(round(delta,3)) + " p(e) = "  + str(round(pe,3)), **hfont)
    plt.savefig(file)
    plt.close(fig)


def match_permutations_dfs(f3, ms2_key, label, cidx, n_comps):
    #node_list_list and candidate_factors are trees
    #prime_multiple and m0 are precursor proprties
    #ms2_vals and ms2_ints are spectrum properties
    t2 = datetime.now()
    
    node_dict = build_node_dict(f3, label)
    # print("node_dict", f3, node_dict)
 
    spec = pgv.ms2_dict[ms2_key]
    ms2_vals = [mz2 for mz2, intensity in spec.ms2.items()]
    ms2_ints = [intensity for mz2, intensity in spec.ms2.items()]
    if pgv.ion_mode == "-":
        zsign = -1
    else:
        zsign = 1
    z2_list = [z * zsign for z in pgv.MS_charge_dict["MS2_charges"][abs(spec.z)]]
    m0 = spec.mz1*abs(spec.z) - abs(spec.z)*zsign*pgv.hmass
    
    precursor_dict = {}
    pidx = 0

    m0, mz1 = calculate_mz1(node_dict["G"], label, spec.z)
    length = len(f3.seq3)
    
    node_dict["ms2_vals"] = ms2_vals
    node_dict["ms2_ints"] = ms2_ints
    node_dict["m0"] = m0
    node_dict["z2_list"] = z2_list
    node_dict["length"] = length
    
    start = datetime.now()
    # complete_list, ctr, bctr, nctr, ectr, pgv.parent_dict= Combination_Product(node_dict).combination_product()
    complete_list, ctr, bctr, nctr, ectr, pgv.parent_dict= Combination_Product_BFS(node_dict).combination_product()
    print(" done with new Combination Product")
    print("Combination_Product took", datetime.now() - start, "for ", len(complete_list), "matches")

    
    print(" number of dfs evals = ", ctr)
    print(" number of node evals = ", nctr)
    print(" number of mod bails = ", bctr)
    print("number of n_match bails = ", ectr)

    frag_seq_key_list = []

    offset = 1000000.0 * (spec.mz1 - mz1)/spec.mz1

    for complete_dict in complete_list:
        # print("current_factors", complete_dict["current_factors"])
        seq3 = repack_current_factors(complete_dict["current_factors"], node_dict["base_list"], node_dict["base_factors_list"])
        ms2_ions = repack_ion_list(complete_dict["ms2_ions"])
        matches = repack_match_list(complete_dict["matches"])
        match_length_list = [len(m) for m in complete_dict["matches"]]

#TODO generalize ends
        fr3 = fragment_sequence([f3.end5n] + seq3 + [f3.end3n])  # just have seq3 be frag3
        # fr3 = fragment_sequence(["5OH"] + seq3 + ["3PO"])  # just have seq3 be frag3
        name = "permutation_" + str(pidx)
        fr = 1
        to = len(seq3)
        frag_seq_key = ":".join([name, str(fr),str(to), fr3.frag])
        
        precursor_dict[pidx] = fr3.__dict__
        precursor_dict[pidx].update({"seq_list": [name], "fr": fr, "to": to, "length": to, "miss": 0, "mod_onoff": 0, 
                                              "frag_type": "seed", "mz1": mz1, "offset": offset, 
                                              "mz_exp": spec.mz1, "z": spec.z, "label": label, "rt": spec.rt,
                                              "CID_ions": ms2_ions, "matched_ions": matches, "match_length_list": match_length_list})
           
        frag_seq_key_list.append(frag_seq_key)  #not used??
        pidx += 1

    print("ion generation and matching took ", datetime.now() - t2)

    # t3 = datetime.now()
    top = score_rank_discovery_dfs(ms2_key, precursor_dict) # top = [key, frag, Sp, rank]
    # print("scor + rank took ", datetime.now() - t3)
    top_sort = sorted(top, key=lambda x: x[2], reverse = True)
    top_keys = [t[0] for t in top_sort]
    top_frags = [t[1] for t in top_sort]
 
    if len(top) > 0:
        top_match = top_sort[0]
    else:
        top_match = "none"

    if pgv.Sp_stats_plots == 'y':
        sp_list = [precursor_dict[key]['Sp'] for key in precursor_dict.keys()]
        if len(sp_list) != 0:
            top_score = round(max(sp_list),2)
        else:
            top_score = 0.0
        spfile = str(ms2_key) + "_" + f3.seq + "_" + str(top_score) + "_Sp_plot.pdf"
        score_distribution_plot(sp_list, "Sp", pgv.plot_dir + "/" + spfile)
        rankfile = str(ms2_key) + "_" + f3.seq + "_" + str(top_score) + "_rank_plot.pdf"
        rank_score_plot(sp_list, "Sp", pgv.plot_dir + "/" + rankfile)
    
    nperms = number_permutations(f3.seq3)

    print()
    print("discovery dfs match for mz1 = ", round(spec.mz1,3), "nperms = ", nperms, "best match of ",len(top_sort), "matches is ", top_match)
    print("composition ", cidx + 1, " of ", n_comps, f3.seq3, " took ", datetime.now() - t2)
       
    return precursor_dict, top_sort, top_keys, top_frags, top_match

def number_permutations(composition):
    bases = {}
    for base in composition:
        if base in bases:
            bases[base] += 1
        else:
            bases[base] = 1
    
    denom = 1
    for base, count in bases.items():
        denom *= math.factorial(count)
    
    num_perms = int(math.factorial(len(composition))/denom)
    return num_perms



def discover_spectra():
    
    build_mass_dict()
    
    ms2_ctr = 0
    # master_match_dict = {}
    ms2_file_list = []
    
    
    # for ms2_key, spec in pgv.ms2_dict.items():  # discover_spectra() return master_match_dict
    for ms2_key in pgv.ms2_key_list:  # discover_spectra() return master_match_dict
        spec = pgv.ms2_dict[ms2_key]
        # if pgv.max_spectra != "all":
        #    if ms2_ctr > int(pgv.max_spectra):
        #        break
        
        pgv.spec_dir = os.path.join(pgv.job_dir, "discovery_data_" + str(ms2_key))
        pgv.plot_dir = pgv.spec_dir
        Path(pgv.spec_dir).mkdir(parents=True, exist_ok=True)

        if ms2_ctr > 300:
            break
        ms2_ctr += 1
        t1 = datetime.now()
        cidx = 0 # composition index
        master_composition_dict = {}
        composition_dict = {}
        if pgv.ion_mode == "-":  
            zsign = -1
        else:
            zsign = 1
             
        z = abs(spec.z)

        mobs = spec.mz1*z - zsign * z * pgv.hmass  # m0 for mz1
        print()
        print("*******")
        print("starting: ", ms2_key, spec.mz1, spec.z, zsign, z, mobs)
        for end5 in pgv.frag_end5:
            for end3 in pgv.frag_end3:
                end5n = pgv.end_dict["end5"][end5]
                end3n = pgv.end_dict["end3"][end3]
                for label in pgv.isotopic_species:

                    mobs_adj = mobs - (pgv.nt_fragment_dict[label][end5n].mass + pgv.nt_fragment_dict[label][end3n].mass)
        
                    seq_comps, ppm_tol, tol = find_compositions(mobs, mobs_adj) # need to add label 
                    print("Compositions:")
                    for comp in seq_comps:
                        print("   ", comp)
                    # print("m0 tolerance: ", ms2_key, ppm_tol, round(tol, 5), "# comps = ", len(seq_comps))
                    if len(seq_comps) > 25:
                        print("TOO MANY COMPOSITIONS...Skip")
                        continue                    # print(seq_comps)
                    composition_dict[ms2_key] = {"n_comps": len(seq_comps), "seq_comps":seq_comps, "ppm_tol": ppm_tol, "tol":tol} # needed??? output separately?
                    composition_dict[ms2_key]["prec_dict"] = {"mz1": spec.mz1, "z": spec.z, "m0": mobs, "end5": end5, "end3": end3, "end5_3": end5n, "end3_3": end3n, 
                                  "label": label}
  
                    for seq_comp in seq_comps:  # match_composition() function
                        nperms = number_permutations(seq_comp)
                        # if nperms > pgv.max_permutations:
                        if len(seq_comp) > 5:
                            pgv.check_dfs_matching = 'y'
                            print()
                            print("using check_dfs on ms2_key ", ms2_key, " of length ", str(len(seq_comp)), "with", nperms, "permutations")
                            print()
                            # continue
                        else:
                            pgv.check_dfs_matching = 'y'
                        f3 = fragment_sequence([end5n] + seq_comp + [end3n])
                        precursor_dict, top_sort, top_keys, top_frags, top_match = match_permutations_dfs(f3, ms2_key, label, cidx, len(seq_comps))
                        unpacked_precursor_dict = prune_precursor_dict(precursor_dict)
                        # insert here:  polish permutations
                        
                        # maybe do Sp histogram first
                        master_composition_dict[cidx] = unpacked_precursor_dict
                        cidx += 1
        
        pickle_file = os.path.join(pgv.spec_dir, "discovery_" + str(ms2_key) + ".pkl")
        ms2_file_list.append(pickle_file)
        print("Saved: ", pickle_file)

        save_pickle(master_composition_dict, pickle_file)
        print("spectrum ", ms2_key, " took ", datetime.now() - t1 )
        print("*******")
        print()
    return ms2_file_list


def prune_precursor_dict(precursor_dict):
    pruned_dict = {}
             
    pkeys_sorted = sorted(list(precursor_dict.keys()), key = lambda x: precursor_dict[x]['Sp'], reverse = True)
    
    ntop = min(len(pkeys_sorted), pgv.ntop)
    pruned_dict = {i:precursor_dict[pkeys_sorted[i]] for i in range(ntop)}

    return pruned_dict

def unpack_precursor_dict(precursor_dict):
    unpacked_comp_dict = {}
    ckey = 0
    print("precursor_dict keys", precursor_dict.keys())
    for cidx, cdict in precursor_dict.items(): # go thru each composition
        print("cdict keys", cdict.keys())
        # print("cdict", cdict)
        for midx, match_dict in cdict.items(): # go thru each match
            if len(match_dict) > 0:
                print("match_dict", match_dict)
                unpacked_comp_dict[ckey] = {"comp_index": cidx, "match_index": midx} # concatentate
                unpacked_comp_dict[ckey].update(match_dict)
                ckey += 1
            else:
                print("empty match_dict", cidx, midx)
            
    ckeys_sorted = sorted(list(unpacked_comp_dict.keys()), key = lambda x: unpacked_comp_dict[x]['Sp'], reverse = True)
    
    ntop = min(len(ckeys_sorted), pgv.ntop)
    unpacked_precursor_dict = {i:unpacked_comp_dict[ckeys_sorted[i]] for i in range(ntop)}
 
    return unpacked_precursor_dict
    


def unpack_master_discovery_dict(master_match_dict):
    unpacked_dict = {}
    ukey = 0
    for ms2_key, mdict in master_match_dict.items(): 
        ms2_comp_dict = {}
        ckey = 0
        for cidx, cdict in mdict.items(): # go thru each composition
            for midx, match_dict in cdict.items(): # go thru each match
                ms2_comp_dict[ckey] = {"ms2_key": ms2_key, "comp_index": cidx, "match_index": midx} # concatentate
                ms2_comp_dict[ckey].update(match_dict)
                ckey += 1
                
        ckeys_sorted = sorted(list(ms2_comp_dict.keys()), key = lambda x: ms2_comp_dict[x]['Sp'], reverse = True)
        for ckey in ckeys_sorted[0:pgv.ntop]:
            unpacked_dict[ukey] = copy.deepcopy(ms2_comp_dict[ckey])
            ukey += 1
    
    return unpacked_dict
 
# def evd_pdf(x,mu,sig):  # old Mathematica def....may be wrong
#     return np.exp((mu-x)/sig) * np.exp(-np.exp((mu-x)/sig))/sig

def evd_pdf(x, mu, sig, eps): # Generalized Extreme Value (GEV) distribution -- Wikipedia
    # mu is location param
    # sig is scale param
    # eps is shape param
    # t_of_x = (1 + eps*(x-mu)/sig)**(-1/eps))
    pdf = (1 + eps*(x-mu)/sig)**(-(eps + 1)/eps) * np.exp( -(1 + eps*(x-mu)/sig)**(-1/eps) )/sig
    return pdf


# def evd_cdf(x,mu,sig):
#     return np.exp(-np.exp((mu-x)/sig))
#     # formulas from Mathematica for extreme value distribution
#     # evdpdf[mu_, sig_, x_] := Exp[(mu - x)/sig] Exp[-Exp[(mu - x)/sig]]/sig
#     # evdcdf[mu_, sig_, x_] := Exp[-Exp[(mu - x)/sig]]

def evd_cdf(x, mu, sig, eps): # Generalized Extreme Value (GEV) distribution -- Wikipedia
    cdf = np.exp( -(1 + eps*(x-mu)/sig)**(-1/eps) )
    return cdf
    # formulas from Mathematica for extreme value distribution
    # evdpdf[mu_, sig_, x_] := Exp[(mu - x)/sig] Exp[-Exp[(mu - x)/sig]]/sig
    # evdcdf[mu_, sig_, x_] := Exp[-Exp[(mu - x)/sig]]

# def inv_evd_cdf(u, mu, sig, eps): # Inverse Generalized Extreme Value (GEV) distribution -- Wikipedia ### DOES NOT WORK
#     # icdf = -(1 + np.log(u)**(-eps))*(sig/eps) + mu
#     icdf =sig * ((-np.log(u)**(-eps)) - (1- eps * mu/sig))/eps
#     return icdf

def validate_discovery():
    upm_df = pd.DataFrame.from_dict(pgv.unpacked_match_dict, orient = "index") # match output
    dis_df = pd.DataFrame.from_dict(pgv.unpacked_discovery_dict, orient = "index") # discovery output

    ms2_keys = list(set(list(upm_df["ms2_key"].unique()) + list(dis_df["ms2_key"].unique())))

    ms2_keys = pgv.ms2_key_list
    n_matches = 0
    n_sub_matches = 0

    for key in ms2_keys:
        m_df = upm_df[upm_df["ms2_key"] == key]
        d_df = dis_df[dis_df["ms2_key"] == key]
        
        if len(d_df) != 0:
            d_top = d_df.iloc[0]
            d_seq = generate_mod_seq(d_top["frag3"])
         
        if len(m_df) != 0:
            m_top = m_df.iloc[0]
            m_seq = generate_mod_seq(m_top["frag3"])

        if len(d_df) == 0:
            print("ms2_key ", key, "no discovery:  match = ", m_seq)
            continue
        
        if len(m_df) == 0:
            print("ms_key ", key, "no match: discovery = ", d_seq)
            continue
        
        
        if m_seq == d_seq:
            print("ms2_key", key, "match = discovery", m_seq)
            n_matches += 1
        else:
            print("ms2_key", key, "match != discovery ", m_seq, d_seq)
            ctr = 0
            for idx, row in d_df.iterrows():
                d_seq = generate_mod_seq(row["frag3"])
                ctr += 1
                if d_seq == m_seq:
                    print("found match at row ", ctr)
                    n_sub_matches += 1
                    

    print()
    print("total matches = ", n_matches)
    print("suboptimal matches ", n_sub_matches)
          

# def get_validated_target(ms2_key):
#     if ms2_key in pgv.valid_dict.keys():
#         target =  pgv.valid_dict[ms2_key]["frag"]
#         end3, tseq, end5 = target.split("_")
#         target3 = dg.parse_mod_seq(tseq)
#     else:
#         target = "I_unknown_am"
#         target3 = ["IAM", "UNK", "KNO", "WN "]
#     return target, target3

