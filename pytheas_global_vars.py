#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pytheas_global_vars.py
global variables defined for all pytheas modules

Created on Thu Sep 30 06:54:36 2021

@author: jrwill
"""

import os
import pandas as pd


class PGV:  # create namespace for pytheas global variables to be imported as pgv
    def __init__(self):

        # print("pytheas_root",  os.path.split(os.path.dirname(__file__)))
        # print(os.path.dirname(__file__))
        # pytheas_root = os.path.split(os.path.dirname(__file__))[0] # directory for pytheas.py
        pytheas_root = os.path.dirname(__file__)
        pytheas_dir = os.path.join(pytheas_root, "pytheas_data")
        print("pytheas_dir", pytheas_dir)
        working_dir = os.getcwd()
        # #TODO don't override pytheas_dir when loading parameters -- allow sharing of parameter files
        # #TODO eliminate user directory -- needlessly complex

        # # Build pytheas global variable dictionary
        # # dictionary has info for building Qt interface, with default values and list of options
        pgvfile = os.path.join(pytheas_root, "pytheas_data/pytheas_global_parameters.xlsx")
        pgvdf = pd.read_excel(pgvfile).fillna("")
        self.pgvdict = pgvdf.set_index('variable_name').T.to_dict() # this holds the variables for GUI setup

        # "value" key not loaded directly from xlsx file due to multiple data types (int, float, char, list)
        # # load current values into pgvdict
        for key in self.pgvdict.keys():
            dtype = self.pgvdict[key]["data_type"]
            def_val = self.pgvdict[key]["default_value"]
            self.pgvdict[key]["value"] = PGVStrToType(dtype, def_val)
            self.pgvdict[key]["pgv"] = key # store variable name....for use by widgets

        # # set directories based on execution environment
        self.pgvdict["pytheas_dir"]["value"] = pytheas_dir
        self.pgvdict["working_dir"]["value" ] = working_dir
        self.pgvdict["pytheas_dir"]["default_value"] = pytheas_dir
        self.pgvdict["working_dir"]["default_value" ] = working_dir

        for key, val in self.pgvdict.items():  # make values available in pgv namespace
            if "value" in val:
                value = val["value"]
            else:
                value = None
            setattr(self,key, value)
 
            
class PGC: # class for global constants
    def __init__(self):
        
        neutron_mass = 1.008665
        
        #worksheet colors
        white = "white"
        green = 'D0FFD0'
  
        # colors for sequence_map plots
        white_hsv = [0.0,0.0,1.0]
        dark_gray = [0.0, 0.0, 0.5]
        red = [0.0, 0.75, 1.0]
        dark_blue = [0.62, 0.75, 1.0]
        light_blue = [0.62, 0.25, 1.0]
        dark_green = [0.3, 0.75, 1.0]
        light_green = [0.3, 0.25, 1.0]

        primes = primes_list(25)
 
        ions = pgvdict["CID_series"]["value"]
        wgts = [.5, .2, 2, .5, 1, .2, 2, 1, 1, 1, 1,1] # combo 2:  
        iw_dict = {ion:wgt for ion, wgt in zip(ions,wgts)}
 
        np2 = 2048 # for fft

        natural = ['ADE', 'CYT', 'GUA', 'URI']

        # floating option panel placement
        panel_y_pos = 50
        panel_x_pos = 1075
        panel_y_delta = 50
        
        # pytheas dicts to save to json files
        digest_json = ["mol", "mod", "unique_precursor",
                                                "unique_frag", "frag"]
        match_json =  ["unpacked_match", "ms2_match", "match", "top_match"]
        discovery_json = ["master_discovery", "unpacked_discovery", "discovery",
                             "iso_mass"]
  
        # ms2_plot constants
        ion_series = ["ag", "Rep", "a","a-B","b","c","d","w","x","y","y-P","z","z-P"] # color is first match in list, default = black
        ion_colors = ["gray", "gray", "green", "green","blue", "magenta", "red", "green", "blue", "magenta", "magenta", "red", "red"]
        ion_color_dict =  {s:c for s,c in zip(ion_series,ion_colors)}

        for var, val in locals().items():  # shortcut to avoid self.var for each var above
            if var == "self":
                continue
            setattr(self, var, val)


def PGVStrToType(dtype, sval):  # convert strings from widgets or pytheas_global_variables.xlsx to data types
    value = ""
    if dtype == "integer":
        try:
            value = int(sval)
        except:
            value = sval
    elif dtype == "float":
        try:
            value = float(sval)
        except:
            value = sval
    elif dtype == "string" or dtype == "directory" or dtype == "file":
        sval = str(sval)
        value = sval.replace("'","").strip()  # strip out quotes
    elif dtype == "list":
        if isinstance(sval,list):
            value = sval
        else:
            value = str(sval).split(",") # lists saved as comma delimited string
    return(value)

def PGVTypeToStr(dtype, value):  # convert data variables to strings for saving parameters
    sval = ""
    if dtype == "integer":
        sval = str(value)
    elif dtype == "float":
        sval = str(value)
    elif dtype == "string" or dtype == "directory" or dtype == "file":
        sval = str(value)
        sval = sval.replace("'","").strip()  # strip out quotes
    elif dtype == "list":
        if isinstance(sval,list):
            sval = ",".join(value)
        else:
            sval =",".join(value)
    return(sval)


def primes_list(n):
    factor = 10
    max_n = factor * n  # works up to n = 1000
    plist = []
    sieve = [True] * (max_n+1)
    for p in range(2, max_n+1):
        if sieve[p] and sieve[p]%2 == 1:
            plist.append(p)
            for i in range(p, max_n+1, p):
                sieve[i] = False
    if len(plist) < n:
        print("increase factor = ", factor, "to get ", n, " primes")
        return []
    return plist[0:n]

pgv = PGV() 
pgvdict = pgv.pgvdict
pgc = PGC()

