#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pytheas.py
Main module for pytheas
PyQt widgets for GUI

Created on Tue Sep 28 06:51:01 2021

@author: jrwill
"""

import argparse
import os
import sys
from datetime import datetime
from pathlib import Path
import glob

# print("sys.argv at start: ", sys.argv)
import pandas as pd
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QDesktopWidget, QFileDialog, QMessageBox
# from IPython import get_ipython
import xlsxwriter



# pytheas modules
from pytheas_global_vars import pgv, pgc, pgvdict, PGVStrToType, PGVTypeToStr
from digest import digest
from match import match
from match_functions import match_output_keys
from pytheas_IO import read_pytheas_file, load_pickle_set
from Pytheas_Qt_widgets import (PytheasRadioButtonBar, PytheasCheckBoxBar, PytheasEntry,
                                                PytheasDirectory, PytheasFile, PytheasFileOutput,
                                                PytheasLabel, PytheasPanel, PytheasFixedPanel, 
                                                PytheasOptionPanel, PytheasButton, PytheasButtonPanel)
from worksheet_functions import write_parameter_worksheet
from isodist import isodist
from discovery import discovery


class Logger(object):
    # duplicates sys.stdout to a log file 
    #source: https://stackoverflow.com/q/616645
    def __init__(self, filename):
        self.stdout = sys.stdout
        self.file = open(filename, mode="a")
        sys.stdout = self

    def __del__(self):
        self.close()

    def __enter__(self):
        pass

    def __exit__(self, *args):
        self.close()

    def write(self, message):
        self.stdout.write(message)
        self.file.write(message)

    def flush(self):
        self.stdout.flush()
        self.file.flush()
        os.fsync(self.file.fileno())

    def close(self):
        if self.stdout != None:
            sys.stdout = self.stdout
            self.stdout = None

        if self.file != None:
            self.file.close()
            self.file = None
                   
# sys.stdout = Logger(pgv.working_dir + "/pytheas.log")

##########  BUTTON FUNCTIONS

def Save_Globals_File(path_list):
        
    name_dict = {var: getattr(pgv,var) for var in pgvdict.keys()}
    
    for key, vval in name_dict.items():
        vtype = pgvdict[key]["data_type"]
        val = PGVTypeToStr(vtype, vval)
        if type(val) != str:
            val = str(val)
        name_dict[key] = val

    par_df = pd.DataFrame(name_dict.items(), columns=["names", "values"])
    ppar_file = os.path.join(pgv.working_dir, *path_list)

    write_parameter_worksheet(par_df, ppar_file)
    print("saving parameters to ", ppar_file)
    
def Load_Globals_File(file_path):
    default = pgv.pytheas_parameters_save
    pgv.pytheas_parameters_save = file_path
    Load_Global_Vars()
    pgv.pytheas_parameters_save = default

def Load_Global_Vars(): # button fuction

    # recent_parameters = "pytheas_recent_parameters.xlsx"
    ppar_file = os.path.join(pgv.working_dir, pgv.pytheas_parameters_save)

    print("Loading Parameters from ", ppar_file)
    if not os.path.exists(ppar_file):
        Error("ERROR: pytheas.py: Load_Global_Vars():" + ppar_file + " does not exist!", 
              "Try Load_Previous Globals or enter parameters manually ")
        return
    
    df = pd.read_excel(ppar_file)
    
    # global vars in three places
        # pgvdict:  this is used for convenient Save_Globals_File
        # pgv namespace for convenient access in code
        # pgv.widget_dict to update Qt widgets 
    for index,row in df.iterrows():
        var = row["names"]
        val = row["values"]
        if var in pgvdict:
            vtype = pgvdict[var]["data_type"]
            vval = PGVStrToType(vtype,val)
        try:
             setattr(pgv, var, vval)   # update global var
             if pgv.run == "CL" and var in pgv.widget_dict:
                 pgv.widget_dict[var].load_value(str(val))  # update the widget
        except:
            print("No global variable for ", var,  "...skipping")
                    
def Load_Previous_Globals(): # button function
    file_dialog = QFileDialog()
    file_dialog.setWindowTitle("Select File")
    file_path, file = QFileDialog.getOpenFileName(None, pgv.working_dir, "Select Excel File", "")
    Load_Globals_File(file_path)
    # default = pgv.pytheas_parameters_save
    # pgv.pytheas_parameters_save = file_path
    # Load_Global_Vars()
    # pgv.pytheas_parameters_save = default

def Save_Global_Variables(): # button function
    # path_list = [pgv.pytheas_parameters_save]
    path_list = ["pytheas_recent_parameters.xlsx"]
    print("save file", path_list)
    Save_Globals_File(path_list)

def Print_Global_Vars(): # button function
    print("PYTHEAS GLOBAL VARIABLES:")
    for var, widget in pgv.widget_dict.items():
        try:
            print(var, widget.value, getattr(pgv,var))
        except:
            pass

def LoadDigest():  # button function
    
    load_dir = QFileDialog.getExistingDirectory(None,"Select Job Directory", pgv.working_dir)
    if load_dir == None:
        return
    print("loading previous digest from ", load_dir)
    
    load_pickle_set(pgv.digest_pickle, load_dir)
    # for obj in pgc.digest_pickle:
    #     pickle_file = os.path.join(load_dir, obj + ".pickle")
    #     load_pickle(obj, pickle_file)


    # load_pickle_files(pgc.digest_pickle, load_dir)
    # load_json_files(pgc.digest_json, load_dir)    
    par_file = glob.glob(os.path.join(load_dir,'*parameters*.xlsx'))[0]
    Load_Globals_File(par_file)
    
def LoadMatch(): # button function
    
    load_dir = QFileDialog.getExistingDirectory(None,"Select Job Directory", pgv.working_dir)
    if load_dir == None:
        return

    # load_json_files(pgc.match_json, load_dir)
    par_file = glob.glob(os.path.join(load_dir,'*parameters*.xlsx'))[0]
    # print("par file", par_file)
    Load_Globals_File(par_file)
    load_dir = pgv.match_job
    print("loading previous match from ", load_dir)
    load_pickle_set(pgc.match_pickle, load_dir)
    # for obj in pgc.match_pickle:
    #     pickle_file = os.path.join(load_dir, obj + ".pickle")
    #     load_pickle(obj, pickle_file)

    
    load_dir = pgv.digest_job
    print("loading previous digest from ", load_dir)
    load_pickle_set(pgc.digest_pickle, load_dir)
    # load_json_files(pgc.digest_json, load_dir)
    # for obj in pgc.digest_pickle:
    #     pickle_file = os.path.join(load_dir, obj + ".pickle")
    #     load_pickle(obj, pickle_file)

        
def Quit(): # button function
    sys.exit()
     
########## FUNCTIONS FOR PYTHEAS MODULES
 
def read_standard_files():
    for w in pgv.std_file_list:
        print("reading standard file ", w)
        if "light" in w and "light" not in pgv.isotopic_species:
            continue
        if "heavy" in w and "heavy" not in pgv.isotopic_species:
            continue        
        error = read_pytheas_file(w)
        if error != None:
            print("error", error)
            error_msg, info = error
            Error(error_msg, info)
            return error
    return None

def setup_job_dir(module):
    dir_list = [d.split("/")[-2] for d in glob.glob(pgv.working_dir + '/*Job*/')]
    if dir_list == []:
        last_dir = module + "_Job_000"
    else:
        last_dir = sorted(dir_list, key=lambda x: int(x.split("_")[-1]))[-1]

    m, j, n = last_dir.split("_")  
    next_dir = "_".join([module, j, format(int(n)+1, '03d')])
    pgv.job_dir = os.path.join(pgv.working_dir, next_dir)
    Path(pgv.job_dir).mkdir(parents=True, exist_ok=True)
    logfile = os.path.join(pgv.job_dir, next_dir + ".log")
    logger = Logger(logfile)
    return logger, next_dir

def Error(error_msg, info):
    # msg = QMessageBox()
    # msg.setIcon(QMessageBox.Critical)
    # msg.setText("\n" + error_msg + "\n")
    # msg.setInformativeText("\n" + info + "\n")
    # msg.setWindowTitle("Error")
    # msg.exec_()
    print("ERROR: ", error_msg, info)

def SaveRunParameters(module, next_dir):

    if not os.path.isdir(pgv.job_dir):
        print ("time stamp making directory....should not happen")
        os.makedirs(pgv.job_dir)

    file_base = module + "_parameters"
    job = next_dir.split("_")[-1]
    par_file =  file_base + "_" + job + ".xlsx"
    path_list = [pgv.job_dir, par_file]
    print("parameter_file: ", path_list)
    if pgv.run == "CL":
        pgv.widget_dict["pytheas_parameters_save"].load_value(par_file) # apparently save file name to itself
    Save_Globals_File(path_list)
    
    par_file = "pytheas_recent_parameters.xlsx"
    path_list = [pgv.working_dir, par_file]
    Save_Globals_File(path_list)


def inSilicoDigest():  # main pytheas module
    logger, next_dir = setup_job_dir("Digest")
    pgv.digest_job = next_dir
    if pgv.run == "CL":
        pgv.widget_dict["digest_job"].load_value(next_dir)

    print("Performing In Silico Digestion:")
    print("isotopic species:", pgv.isotopic_species)
    
    error = read_standard_files()
    
    if error != None:
        return
    
    read_pytheas_file("rna_mod_defs_light")
    read_pytheas_file("rna_mod_defs_heavy")

    error = digest()
    if error != None:
        Error("Error in " + error, "check fasta and mod input files!")
    SaveRunParameters("digest", next_dir)
    print("Done with In Silico Digestion!")
    print()
    print()
    
    logger.close()
        
def matchSpectra(): # main pytheas module
    logger, next_dir = setup_job_dir("Match")
    pgv.match_job = next_dir
    if pgv.run == "CL":
        pgv.widget_dict["match_job"].load_value(next_dir)

    print('Matching MS2 spectra')
    
    error = read_standard_files()
    if error != None:
        return

    try:
        if pgv.unique_precursor_dict:
            pass
    except:
        Error("No Digest dictionaries present", "Either Run In_Silico_Digest or Load Digest Files")
        return
    match()
    SaveRunParameters("matching", next_dir)
    print("Done with Matching MS2 spectra")
    print()
    print()
    logger.close()

def Isodist():
    logger, next_dir = setup_job_dir("Isodist")
    pgv.isodist_job = next_dir
    print("Quantifying fragments by Isodist")

    error = read_standard_files()
    if error != None:
        return

    # print("initial mod set", pgv.modification_set)
    isodist()
    # print("after discovery", pgv.modification_set)
    SaveRunParameters("isodist", next_dir)
   
    print("Done with Isodist")


def Discovery():
    logger, next_dir = setup_job_dir("Discovery")
    pgv.discovery_job = next_dir
    print("Sequencing fragments by discovery")

    error = read_standard_files()
    if error != None:
        return

    print("initial mod set", pgv.modification_set)
    discovery()
    print("after discovery", pgv.modification_set)
    SaveRunParameters("discovery", next_dir)
    print("after timestamp", pgv.modification_set)
    print("Done with Discovery")

# def SetUp():
#     print("Reading standard files ")
#     error = read_standard_files()
#     if error != None:
#         return
    
#     pgv.read_pytheas_file(pgv.widget_dict["modification_set_file"])
#     pgv.read_pytheas_file(pgv.widget_dict["MS_data_file"])  # read in MS data to pgv.ms2_dict

#     for ms2_key, mdict in pgv.ms2_dict.items():
#         print(ms2_key, "max_int", mdict.max_int, "# ms2", len(list(mdict.ms2.keys())))


      
# def Print_Match_Keys():
#     try:
#         print("ms2_match_dict keys: ", pgv.match_dict_keys)
#     except:
#         Error("No match dict!", "Either Load_Match_Files or run matching")

# pytheas_root = os.path.split(os.path.dirname(__file__))[0] # directory for pytheas.py

# def Print_System_Info():
#     # print("screen: ", screen_width, screen_height, width_scale, height_scale, screen_scale )
#     print("****** SYSTEM and SCREEN DIMENSIONS ******")
#     print("native screen width: ", screen_width)
#     print("native screen height: ", screen_height)
#     print("calculated screen scale: ", screen_scale)
#     print("Platform detected: ", platform_system)
#     print("default font size:", def_font_size)
#     print("font scale: ", font_scale)
#     print("display font used: ", entry_font)
#     print("scaled Tk window dimensions: ", xs,"x",ys)

#     print("*****************************(************")

# Handle command line arguments

parser = argparse.ArgumentParser("Pytheas_Qt")
# parser.add_argument("--screen_scale", dest = "screen_scale", type=float, help="enforce a relative screen scale as a refuge for those with a truly small monitor. range (0.5 to 2.0)")
parser.add_argument("--load_vars",  action='store_true')
# parser.add_argument("--print_sys_info", action = 'store_true', help = "print info about os and screen dimensions")
parser.add_argument("--digest", action = 'store_true')
parser.add_argument("--match", action = 'store_true')
parser.add_argument("--nogui", action = 'store_true')

clargs = parser.parse_args()
print("arguments = ", clargs)

# valid_args = ['nogui', 'load_vars', 'screen_scale', 'all', 'in_silico_digest', 'matching', 'statistics', 'final_report', 'visualization']

def build_Pytheas_gui():
    
    # pgvdict has all of the info for building the GUI
    
    fixed_panels = ["input_files_dir", "input_files_req", "input_files_opt"] # fixed panels at top of layout

    pgv.option_dict = {}    # defines the widgets for each option panel, key = group
    option_group_dict = {}   # defines the option groups for the master option panel  key = widget group
    std_file_list = [] # files to be read on startup
    label_file_list = [] # heavy/light files

    for key, pdict in pgvdict.items(): # extract info for each panel
    
        if pdict["hide"] == "hide":
            continue
    
        if "option_group" in pdict:
            og = pdict["option_group"]
            wg = pdict["widget_group"]
            if og not in pgv.option_dict.keys():
                pgv.option_dict[og] = [] # preserves input order
            pgv.option_dict[og].append(key)
    
            if pgvdict[key]["option_group"] == "required_input_files" and pgvdict[key]["widget_type"] != "PytheasLabel":
                std_file_list.append(key)
            if "rna_mod_defs" in key:   # label files have to be read last
                label_file_list.append(key)
            
            if og in fixed_panels:
                continue
            if wg in option_group_dict.keys():
                if og not in option_group_dict[wg]:
                    option_group_dict[wg].append(og)
            else:
                option_group_dict[wg] = [og]
    
     
    pgv.std_file_list = std_file_list + label_file_list # label has to be read last
    
    
    print("INITIALIZING MAIN WINDOW")
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    main_window = QWidget()
    main_window.setWindowTitle("Pytheas 2.0")
    
    # screen_geo = QDesktopWidget().screenGeometry()
    # widget_geo = main_window.geometry()
    # print("screen geo", screen_geo.width(), screen_geo.height())
    # print("widget geo", widget_geo.width(), widget_geo.height())
    # x = screen_geo.width() - widget_geo.width()
    # y = screen_geo.height() - widget_geo.height()
    # print("xy", x, y)
    # main_window.move(x, y)
    main_window.move(0,0)
    main_window_layout = QVBoxLayout()
    main_window.setLayout(main_window_layout)
    main_window.show()
    
    floating_windows = True
    
    #   build all panel widgets
    
    option_panel_list, fixed_panel_list = [], []
    for group, glist in pgv.option_dict.items():
        
        if group not in fixed_panels:
            panel = PytheasPanel(group, glist)
            option_panel_list.append(panel)
            panel.hide_panel()
        else:
            panel = PytheasFixedPanel(group, glist)
            fixed_panel_list.append(panel)
            main_window_layout.addWidget(panel)
            
    # add module buttons
    
    button_dict = {"Global Variables:": {"Load Last Globals": Load_Global_Vars,
                                         "Choose Previous Globals": Load_Previous_Globals,
                                         "List Globals": Print_Global_Vars,
                                         "Save Globals": Save_Global_Variables                                     
                                          },
                   "Pytheas Modules: ": {"In Silico Digest": inSilicoDigest,
                                         "Match Spectra": matchSpectra, 
                                         "Isodist": Isodist,
                                         "Discovery": Discovery,
                                         # "Dev": Development,
                                         "Quit": Quit
                                         },
                   "Load Pytheas Files": {"Load Previous Digest": LoadDigest,
                                          "Load Previous Match and Digest": LoadMatch}
                   }
    
    for label, bdict in button_dict.items():
        button_panel = PytheasButtonPanel(label, bdict)
        main_window_layout.addWidget(button_panel)
    
    # add master option panel
    
    pgv.option_panel_dict = {panel.group: panel for panel in option_panel_list}
    option_panel = PytheasOptionPanel(option_group_dict, pgv.option_panel_dict, floating_windows)    # option panel widget
    main_window_layout.addWidget(option_panel)
    
    if not floating_windows:    # turn this off to get floating windows
        for panel in option_panel_list:        
            main_window_layout.addWidget(panel)      
        
    all_panel_list = fixed_panel_list + option_panel_list
    master_widget_list = [w for panel in all_panel_list for w in panel.widget_list]
    
    pgv.widget_dict = {w.pgv: w for w in master_widget_list}   # for save/load global variables use
    
    return main_window, app

# screen_geo = QDesktopWidget().screenGeometry()
# print("Screen geo", screen_geo)
# widget_geo = main_window.geometry()
# print("widget geo", widget_geo)
# print("available geo ", QDesktopWidget().availableGeometry())
# print("screen geo", screen_geo.width(), screen_geo.height())
# print("widget geo", widget_geo.width(), widget_geo.height())
# x = screen_geo.width() - widget_geo.width()
# y = screen_geo.height() - widget_geo.height()
# print("xy", x, y)
# main_window.move(x, y)


main_window, app = build_Pytheas_gui()

# crude hack to allow development in Spyder

try:
    run_mode = get_ipython().__class__.__name__
except:
    run_mode = "command_line"
    print("command_line_arguments: ", sys.argv)

if run_mode == 'SpyderShell':
    print("This program is running inside Spyder") 
    pgv.run = "Spyder"
else:
    print("This program is running on command line")
    pgv.run = "CL"


# handle command line args and launch app
      
print("launching Pytheas...")

if clargs.load_vars:
    Load_Global_Vars()
 
if clargs.digest:
    inSilicoDigest()
    
if clargs.match:
    matchSpectra()
    
if clargs.nogui:
    pass
else:
    if pgv.run == "CL":
        app.exec()
        print("Pytheas ready to Go!")
        
        






