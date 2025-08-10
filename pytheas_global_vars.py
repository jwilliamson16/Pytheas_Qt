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
    def __init__(self, pdict):
        
        print("SETTING UP PYTHEAS GLOBAL VARIABLES")

        # # pytheas_root = "/Users/jrwill/prog/pytheas_tk_interface/pytheas_root"
        # pytheas_root = os.path.split(os.path.dirname(__file__))[0] # directory for pytheas.py
        # pytheas_dir = os.path.join(pytheas_root, "pytheas_data")
        # working_dir = os.getcwd()
        # user_dir = ""

        # #TODO don't override pytheas_dir when loading parameters -- allow sharing of parameter files
        # #TODO eliminate user directory -- needlessly complex

        # # Build pytheas global variable dictionary
        # # dictionary has info for building Qt interface, with default values and list of options
        # pgvfile = os.path.join(pytheas_root, "pytheas_data/pytheas_global_parameters.xlsx")
        # pgvdf = pd.read_excel(pgvfile).fillna("")
        # pgvdict = pgvdf.set_index('variable_name').T.to_dict() # this holds the variables for GUI setup

        # # load current values into pgvdict

        # for key in pgvdict.keys():
        #     vtype,vval, group, wtype = list( map(pgvdict.get(key).get, ['data_type', 'default_value', "widget_group", "widget_type"] ))
        #     pgvdict[key]["value"] = PGVStrToType(vtype,vval)
        #     # pgvdict[key]["value"] = str(vval)
        #     pgvdict[key]["default_value"] = str(vval)
        #     pgvdict[key]["pgv"] = key


        # # set defaults based on execution
        # pgvdict["pytheas_dir"]["value"] = pytheas_dir
        # pgvdict["working_dir"]["value" ] = working_dir
        # pgvdict["pytheas_dir"]["default_value"] = pytheas_dir
        # pgvdict["working_dir"]["default_value" ] = working_dir

        
        
        
        for key, val in pdict.items():
            if "value" in val:
                value = val["value"]
            else:
                value = None
            setattr(self,key, value)
    
    def update_widget(self, var, val): # primarily for Load_Global_Vars to set widgets
        setattr(self,var, val)
        if var in self.widget_dict and self.run == "CL":  # This updates the widget
            print("updating widget in PGV class", var, vval)
            self.widget_dict[var].load_value(str(val))
            print("just updated widget in PGV class", var, vval)

            
class PGC: # class for global constants
    def __init__(self):
        
        
        # pgcdict = {}
        # pgcdict["neutron_mass"] = {"value": 1.008665}
        neutron_mass = 1.008665
        white = "white"
        green = 'D0FFD0'
        
        # pgcdict["white"] = {"value": "white"}
        # pgcdict["green"] = {"value": 'D0FFD0'}

        # colors for sequence_map plots
        white_hsv = [0.0,0.0,1.0]
        dark_gray = [0.0, 0.0, 0.5]
        red = [0.0, 0.75, 1.0]
        dark_blue = [0.62, 0.75, 1.0]
        light_blue = [0.62, 0.25, 1.0]
        dark_green = [0.3, 0.75, 1.0]
        light_green = [0.3, 0.25, 1.0]
        # pgcdict["white_hsv"] = {"value":[0.0,0.0,1.0]} # gotta be a float @#$%^&^*
        # pgcdict["dark_gray"] = {"value": [0.0, 0.0, 0.5]}
        # pgcdict["red"] = {"value":[0.0, 0.75, 1.0]}
        # pgcdict["dark_blue"] = {"value":[0.62, 0.75, 1.0]}
        # pgcdict["light_blue"] = {"value":[0.62, 0.25, 1.0]}
        # pgcdict["dark_green"] = {"value":[0.3, 0.75, 1.0]}
        # pgcdict["light_green"] = {"value":[0.3, 0.25, 1.0]}


        primes = primes_list(25)
        # pgcdict["primes"] = {"value": primes}

        ions = pgvdict["CID_series"]["value"]
        wgts = [.5, .2, 2, .5, 1, .2, 2, 1, 1, 1, 1,1] # combo 2:  

        iw_dict = {ion:wgt for ion, wgt in zip(ions,wgts)}
        # pgcdict["iw_dict"] = {"value": iw_dict}
        # weight_beta = 'y'
        np2 = 2048 
        # pgcdict["np2"] = {"value": np2}

        natural = ['ADE', 'CYT', 'GUA', 'URI']
        # pgcdict["natural"] = {"value": natural}

        panel_y_pos = 50
        panel_x_pos = 1075
        panel_y_delta = 50
        
        # pgcdict["panel_y_pos"] = {"value": 50}  # ypos of top floating widget
        # pgcdict["panel_x_pos"] = {"value": 1075} # xpos of floating widgets
        # pgcdict["panel_y_delta"] = {"value": 50} # padding between floating widgets

        digest_json = ["mol", "mod", "unique_precursor",
                                                "unique_frag", "frag"]
        match_json =  ["unpacked_match", "ms2_match", "match", "top_match"]
        # pgcdict["digest_json"] = {"value": ["mol", "mod", "unique_precursor",
        #                                         "unique_frag", "frag"]}
        # pgcdict["match_json"] = {"value": ["unpacked_match", "ms2_match", "match", "top_match"]}

        # ms2_plot constants
        ion_series = ["ag", "Rep", "a","a-B","b","c","d","w","x","y","y-P","z","z-P"] # color is first match in list, default = black
        ion_colors = ["gray", "gray", "green", "green","blue", "magenta", "red", "green", "blue", "magenta", "magenta", "red", "red"]
        ion_color_dict =  {s:c for s,c in zip(ion_series,ion_colors)}
        # pgcdict["ion_color_dict"] = {"value": {s:c for s,c in zip(ion_series,ion_colors)}}
        # pgcdict["ion_series"] = {"value": ion_series}

        for var, val in locals().items():  # shortcut to avoid self.var for each var above
            if var == "self":
                continue
            setattr(self, var, val)


        # for key, val in cdict.items():
        #     setattr(self, key, val["value"])

def PGVStrToType(vtype, vval):  # convert strings from widgets or pytheas_global_variables.xlsx to data types
    val = ""
    if vtype == "integer":
        try:
            val = int(vval)
        except:
            val = vval
    elif vtype == "float":
        try:
            val = float(vval)
        except:
            val = vval
    elif vtype == "string" or vtype == "directory" or vtype == "file":
        vval = str(vval)
        val = vval.replace("'","").strip()  # strip out quotes
    elif vtype == "list":
        if isinstance(vval,list):
            val = vval
        else:
            # val = str(vval).replace(" ",",").split(",") # have to save lists as space delimited
            val = str(vval).split(",") # lists saved as comma delimited string
    return(val)

def PGVTypeToStr(vtype, vval):  # convert data variables to strings for saving parameters
    val = ""
    if vtype == "integer":
        val = str(vval)
    elif vtype == "float":
        val = str(vval)
    elif vtype == "string" or vtype == "directory" or vtype == "file":
        vval = str(vval)
        val = vval.replace("'","").strip()  # strip out quotes
    elif vtype == "list":
        if isinstance(vval,list):
            val = ",".join(vval)
        else:
            # val = str(vval)
            val = ",".join(vval)
    return(val)


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

    

print("SETTING UP PYTHEAS GLOBAL VARIABLES")

# pytheas_root = "/Users/jrwill/prog/pytheas_tk_interface/pytheas_root"
pytheas_root = os.path.split(os.path.dirname(__file__))[0] # directory for pytheas.py
pytheas_dir = os.path.join(pytheas_root, "pytheas_data")
working_dir = os.getcwd()
user_dir = ""

#TODO don't override pytheas_dir when loading parameters -- allow sharing of parameter files
#TODO eliminate user directory -- needlessly complex

# Build pytheas global variable dictionary
# dictionary has info for building Qt interface, with default values and list of options
pgvfile = os.path.join(pytheas_root, "pytheas_data/pytheas_global_parameters.xlsx")
pgvdf = pd.read_excel(pgvfile).fillna("")
pgvdict = pgvdf.set_index('variable_name').T.to_dict()

# load current values into pgvdict

for key in pgvdict.keys():
    vtype,vval, group, wtype = list( map(pgvdict.get(key).get, ['data_type', 'default_value', "widget_group", "widget_type"] ))
    pgvdict[key]["value"] = PGVStrToType(vtype,vval)
    # pgvdict[key]["value"] = str(vval)
    pgvdict[key]["default_value"] = str(vval)
    pgvdict[key]["pgv"] = key


# set defaults based on execution
pgvdict["pytheas_dir"]["value"] = pytheas_dir
pgvdict["working_dir"]["value" ] = working_dir
pgvdict["pytheas_dir"]["default_value"] = pytheas_dir
pgvdict["working_dir"]["default_value" ] = working_dir



# add miscellaneous constants to namespace and pgvdict


# pgcdict = {}
# pgcdict["neutron_mass"] = {"value": 1.008665}
# pgcdict["white"] = {"value": "white"}
# pgcdict["green"] = {"value": 'D0FFD0'}

# # colors for sequence_map plots
# pgcdict["white_hsv"] = {"value":[0.0,0.0,1.0]} # gotta be a float @#$%^&^*
# pgcdict["dark_gray"] = {"value": [0.0, 0.0, 0.5]}
# pgcdict["red"] = {"value":[0.0, 0.75, 1.0]}
# pgcdict["dark_blue"] = {"value":[0.62, 0.75, 1.0]}
# pgcdict["light_blue"] = {"value":[0.62, 0.25, 1.0]}
# pgcdict["dark_green"] = {"value":[0.3, 0.75, 1.0]}
# pgcdict["light_green"] = {"value":[0.3, 0.25, 1.0]}


# primes = primes_list(25)
# pgcdict["primes"] = {"value": primes}

# ions = pgvdict["CID_series"]["value"]
# wgts = [.5, .2, 2, .5, 1, .2, 2, 1, 1, 1, 1,1] # combo 2:  

# iw_dict = {ion:wgt for ion, wgt in zip(ions,wgts)}
# pgcdict["iw_dict"] = {"value": iw_dict}
# # weight_beta = 'y'
# np2 = 2048 
# pgcdict["np2"] = {"value": np2}

# natural = ['ADE', 'CYT', 'GUA', 'URI']
# pgcdict["natural"] = {"value": natural}

# pgcdict["panel_y_pos"] = {"value": 50}  # ypos of top floating widget
# pgcdict["panel_x_pos"] = {"value": 1075} # xpos of floating widgets
# pgcdict["panel_y_delta"] = {"value": 50} # padding between floating widgets

# pgcdict["digest_json"] = {"value": ["mol", "mod", "unique_precursor",
#                                         "unique_frag", "frag"]}
# pgcdict["match_json"] = {"value": ["unpacked_match", "ms2_match", "match", "top_match"]}

# # ms2_plot constants
# ion_series = ["ag", "Rep", "a","a-B","b","c","d","w","x","y","y-P","z","z-P"] # color is first match in list, default = black
# ion_colors = ["gray", "gray", "green", "green","blue", "magenta", "red", "green", "blue", "magenta", "magenta", "red", "red"]
# pgcdict["ion_color_dict"] = {"value": {s:c for s,c in zip(ion_series,ion_colors)}}
# pgcdict["ion_series"] = {"value": ion_series}

# initialize namespace with dictionary
pgv = PGV(pgvdict) 
# pgc = PGC(pgcdict)
pgc = PGC()

