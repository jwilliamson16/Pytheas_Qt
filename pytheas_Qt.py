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
from pytheas_IO import read_pytheas_file, read_json, load_json_files
from Pytheas_Qt_widgets import (PytheasRadioButtonBar, PytheasCheckBoxBar, PytheasEntry,
                                                PytheasDirectory, PytheasFile, PytheasFileOutput,
                                                PytheasLabel, PytheasPanel, PytheasFixedPanel, 
                                                PytheasOptionPanel, PytheasButton, PytheasButtonPanel)




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

# # print("sys.argv after pytheas_module import: ", sys.argv)

# import pytheas_final_report as fr
# import pytheas_mapping as ma
# import pytheas_visualization_html as vi
# import digest_new_8 as dn

def Load_Global_Vars():

    print("Load_Global_Vars")
    print("working_dir", pgv.working_dir)
    print("save_value", pgv.pytheas_parameters_save)
    ppar_file = os.path.join(pgv.working_dir, pgv.pytheas_parameters_save)

    print("loading from ", ppar_file)
    if not os.path.exists(ppar_file):
        Error("ERROR: pytheas.py: Load_Global_Vars():" + ppar_file + " does not exist!", 
              "Try Load_Previous Globals or enter parameters manually ")
        return
    
    df = pd.read_excel(ppar_file)
             
    # for index,row in df.iterrows():
    #     if row['names'] in pgv.widget_dict:
    #         pgv.widget_dict[row['names']].load_value(str(row['values']))
    #         # print(row['names'],row['values'])
    #     else:
    #         print("No widget for ", row['names'], "...skipping")

    for index,row in df.iterrows():
        var = row["names"]
        val = row["values"]
        if var in pgvdict:
            vtype = pgvdict[var]["data_type"]
            vval = PGVStrToType(vtype,val)
        try:
            print("updating", var, vval )
            if var == "enzymes":
                continue
            setattr(pgv, var, vval)   # HOW DOES THIS UPDATE WIDGET????
            print("just updated", var, vval)
            # print(var, vtype, val, vval)
        except:
            print("No global variable for ", var,  "...skipping")
            
        if pgv.run == "CL":
            if var in pgv.widget_dict:
                pgv.widget_dict[var].load_value(str(val))
        # try:
        #    w = pgv.widget_dict[var]
        #    w.value = val
        # except:
        #    print("Couldn't update widget for ", var,  "...skipping")
            

def Load_Previous_Globals():
    file_dialog = QFileDialog()
    file_dialog.setWindowTitle("Select File")
    file_path, file = QFileDialog.getOpenFileName(None, pgv.working_dir, "Select Excel File", "")
    default = pgv.pytheas_parameters_save
    pgv.pytheas_parameters_save = file_path
    Load_Global_Vars()
    pgv.pytheas_parameters_save = default

def Load_Globals_File(file_path):
    default = pgv.pytheas_parameters_save
    pgv.pytheas_parameters_save = file_path
    Load_Global_Vars()
    pgv.pytheas_parameters_save = default

def Error(error_msg, info):
    # msg = QMessageBox()
    # msg.setIcon(QMessageBox.Critical)
    # msg.setText("\n" + error_msg + "\n")
    # msg.setInformativeText("\n" + info + "\n")
    # msg.setWindowTitle("Error")
    # msg.exec_()
    print("ERROR: ", error_msg, info)

def Save_Global_Vars(path_list):
        
    # name_dict = {var: pdict["value"] for var, pdict in pgvdict.items()}
    name_dict = {var: getattr(pgv,var) for var in pgvdict.keys()}
    
    for key, vval in name_dict.items():
        vtype = pgvdict[key]["data_type"]
        val = PGVTypeToStr(vtype, vval)
        # print(key, type(vval), val)
        if type(val) != str:
            val = str(val)
        name_dict[key] = val

    par_df = pd.DataFrame(name_dict.items(), columns=["names", "values"])
    
    ppar_file = os.path.join(pgv.working_dir, *path_list)
    
    workbook = xlsxwriter.Workbook(ppar_file,{"nan_inf_to_errors": True})
    worksheet = workbook.add_worksheet(pgv.job_dir.split("/")[-1])
    
    bold = workbook.add_format()
    bold.set_bold()
    worksheet.set_column(0,0,25)
    worksheet.set_column(1,1,100)
    
    # format_worksheet_columns(worksheet, par_df, fd)   # sets the column widths

    for col in par_df.columns:  # output header
        worksheet.write(0, par_df.columns.get_loc(col), col, bold)  # fd has user-defined label column if desired
    
    # write out each cell
    rowidx = 1
    for row in par_df.iterrows():
        colidx = 0
        for col in row[1].keys():
            if type(row[1][col]) == list:
                row[1][col] = ",".join(row[1][col])
            worksheet.write(rowidx, colidx, row[1][col])
            colidx +=1
        rowidx +=1

    # worksheet.show_gridlines = True
    
    workbook.close()
    
    # writer = pd.ExcelWriter(ppar_file) 
    # df.style.set_properties(**{'text-align': 'left'})
    # df.to_excel(writer, index=False, sheet_name='pytheas_parameters', na_rep='NaN')
    # workbook  = writer.book
    # worksheet = writer.sheets['pytheas_parameters']

    # header_fmt = workbook.add_format({'bold': True})
    # worksheet.set_row(0, None, header_fmt)
    
    # left_format = workbook.add_format({'align': 'left'})
    
    # col_idx = df.columns.get_loc("names")
    # worksheet.set_column(col_idx, col_idx, 25)
    # col_idx = df.columns.get_loc("values")
    # worksheet.set_column(col_idx, col_idx, 100, left_format)

    # writer.save()
    print("saving parameters to ", os.path.basename(ppar_file))
    print("saving parameters to ", ppar_file)

def TimeStampParameters(module, next_dir):

    # parameter_folder = os.path.join(pgv.working_dir, "parameter_history")      
    # if not os.path.isdir(parameter_folder):
    #     os.makedirs(parameter_folder)

    if not os.path.isdir(pgv.job_dir):
        print ("time stamp making directory....should not happen")
        os.makedirs(pgv.job_dir)

    file_base = module + "_parameters"
    # timestamp = datetime.now().strftime("%Y_%m_%d_%H_%M")
    job = next_dir.split("_")[-1]
    # par_file =  file_base + "_" + timestamp + ".xlsx"
    par_file =  file_base + "_" + job + ".xlsx"

    path_list = [pgv.job_dir, par_file]
    print("parameter_file: ", path_list)
    if pgv.run == "CL":
        pgv.widget_dict["pytheas_parameters_save"].load_value(par_file)
    Save_Global_Vars(path_list)
 
def Save_Global_Variables():
    path_list = [pgv.pytheas_parameters_save]
    Save_Global_Vars(path_list)

def Print_Global_Vars():
    print("PYTHEAS GLOBAL VARIABLES:")
    for var, widget in pgv.widget_dict.items():
        try:
            print(var, widget.value, getattr(pgv,var))
        except:
            pass

def Quit():
    # sys.exit(app.exec())
    sys.exit()
     
 #  functions for pytheas modules
 
def read_standard_files():
    for w in std_file_list:
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

        # last_dir = sorted(dir_list)[-1]
    m, j, n = last_dir.split("_")  
    next_dir = "_".join([module, j, format(int(n)+1, '03d')])
    pgv.job_dir = os.path.join(pgv.working_dir, next_dir)
    Path(pgv.job_dir).mkdir(parents=True, exist_ok=True)
    logfile = os.path.join(pgv.job_dir, next_dir + ".log")
    logger = Logger(logfile)
    return logger, next_dir


def inSilicoDigest():
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
    TimeStampParameters("digest", next_dir)
    print("Done with In Silico Digestion!")
    logger.close()
        
def matchSpectra():
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
    TimeStampParameters("matching", next_dir)
    print("Done with Matching MS2 spectra")
    logger.close()

# def Discovery():
#     print("Sequencing fragments by discovery")

#     error = read_standard_files()
#     if error != None:
#         return

#     print("initial mod set", pgv.modification_set)
#     dis.discovery()
#     print("after discovery", pgv.modification_set)
#     TimeStampParameters("discovery")
#     print("after timestamp", pgv.modification_set)
#     print("Done with Discovery")

# def SetUp():
#     print("Reading standard files ")
#     error = read_standard_files()
#     if error != None:
#         return
    
#     pgv.read_pytheas_file(pgv.widget_dict["modification_set_file"])
#     pgv.read_pytheas_file(pgv.widget_dict["MS_data_file"])  # read in MS data to pgv.ms2_dict

#     for ms2_key, mdict in pgv.ms2_dict.items():
#         print(ms2_key, "max_int", mdict.max_int, "# ms2", len(list(mdict.ms2.keys())))


# def Development():
#     print("Under construction")
#     pass

def LoadDigest():
    
    load_dir = QFileDialog.getExistingDirectory(None,"Select Job Directory", pgv.working_dir)
    if load_dir == None:
        return
    
    print("loading previous digest from ", load_dir)

    load_json_files(pgc.digest_json, load_dir)    
    par_file = glob.glob(os.path.join(load_dir,'*parameters*.xlsx'))[0]
    print("par file", par_file)
    Load_Globals_File(par_file)
    

    
def LoadMatch():
    
    load_dir = QFileDialog.getExistingDirectory(None,"Select Job Directory", pgv.working_dir)
    if load_dir == None:
        return
    
    print("loading previous match from ", load_dir)

    load_json_files(pgc.match_json, load_dir)
    par_file = glob.glob(os.path.join(load_dir,'*parameters*.xlsx'))[0]
    print("par file", par_file)
    Load_Globals_File(par_file)
    
    load_dir = pgv.digest_job
    print("loading previous digest from ", load_dir)
    
    load_json_files(pgc.digest_json, load_dir)
 

# def statistics():
#     print('Statistical analysis...')
#     ps.statistics()
     
# def report():
#     # Update_Global_Vars()
#     print(" Generating final report...")
#     fr.final_report()
#     # printlog("Generating sequence mapping ")
#     ma.mapping()
      
# def visualization():
#     # Update_Global_Vars()
#     print("Generating Visualization...")
#     vi.visualization()
    
# def CommandLineRun(args):
#     Load_Global_Vars()
#     if "in_silico_digest" in args or "all" in args:
#         inSilicoDigest()
#     if "matching" in args or "all" in args:
#         matchSpectra()
#     if "statistics" in args or "all" in args:
#         statistics()
#     if "final_report" in args or "all" in args:
#         report()
#     if "visualizationt" in args or "all" in args:
#         visualization()
#     sys.exit()

#TODO update Loads to look for latest Job directory

# def Load_Digest_Files():
#     try:
#         pgv.mol_dict = pgv.read_json(os.path.join(pgv.pytheas_data_folder, "mol.json"))
#         pgv.mod_dict = pgv.read_json(os.path.join(pgv.pytheas_data_folder, "mod.json"))
#         pgv.unique_frag_dict = pgv.read_json(os.path.join(pgv.pytheas_data_folder, "unique_frag.json"))
#         pgv.frag_dict = pgv.read_json(os.path.join(pgv.pytheas_data_folder, "frag.json"))
#         pgv.unique_precursor_dict = pgv.read_json(os.path.join(pgv.pytheas_data_folder, "unique_precursor.json"))
#         pgv.precursor_dict = pgv.read_json(os.path.join(pgv.pytheas_data_folder, "precursor.json"))
#     except:
#         Error("Problem reading pytheas digest files", "Check " + pgv.pytheas_data_folder + " or re-run digest")
        
# def Load_Match_Files():
#     try:       
#         pgv.unpacked_match_dict = pgv.read_json(os.path.join(pgv.pytheas_data_folder, "unpacked_match.json"))
#         pgv.ms2_match_dict = pgv.read_json(os.path.join(pgv.pytheas_data_folder, "ms2_match.json"))
#         pgv.match_dict = pgv.read_json(os.path.join(pgv.pytheas_data_folder, "match.json"))
#         pgv.top_match_dict = pgv.read_json(os.path.join(pgv.pytheas_data_folder, "top_match.json"))

#         match_output_keys()
#     except:
#         Error("Problem reading pytheas match files", "Check " + pgv.pytheas_data_folder + " or re-run matching")

      
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

# parser = argparse.ArgumentParser("Pytheas_tk")
# parser.add_argument("--screen_scale", dest = "screen_scale", type=float, help="enforce a relative screen scale as a refuge for those with a truly small monitor. range (0.5 to 2.0)")
# parser.add_argument("--load_vars", action = 'store_true',  help = "load saved parameters in pytheas_parameters.xlsx upon execution")
# parser.add_argument("--print_sys_info", action = 'store_true', help = "print info about os and screen dimensions")

# clargs = parser.parse_args()
# # print("arguments = ", clargs)
# # print("args.screen_scale", clargs.screen_scale)
# # print("args.load_vars", clargs.load_vars)

# valid_args = ['nogui', 'load_vars', 'screen_scale', 'all', 'in_silico_digest', 'matching', 'statistics', 'final_report', 'visualization']

fixed_panels = ["input_files_dir", "input_files_req", "input_files_opt"] # fixed panels at top of layout

option_dict = {}    # defines the widgets for each option panel
for key, pdict in pgvdict.items():
    if "option_group" in pdict:
        group = pdict["option_group"]
        if group not in option_dict.keys():
            option_dict[group] = [] # preserves input order
        option_dict[group].append(key)

option_group_dict = {}   # defines the option groups for the master option widget
for key, odict in pgvdict.items():
    if "option_group" in odict:
        wg = odict["widget_group"]
        og = odict["option_group"]
        if og in fixed_panels:
            continue
        if wg in option_group_dict.keys():
            if og not in option_group_dict[wg]:
                option_group_dict[wg].append(og)
        else:
            option_group_dict[wg] = [og]

# standard files 
    std_file_list = []
    label_file_list = []
    for w, wdict in pgvdict.items():
        if "option_group" in wdict:
            if pgvdict[w]["option_group"] == "input_files" and pgvdict[w]["widget_type"] != "PytheasLabel":
                std_file_list.append(w)
    # std_file_list = [w for w in pgvdict.keys() if pgvdict[w]["option_group"] == "input_files" and pgvdict[w]["widget_type"] != "PytheasLabel"]
            if "rna_mod_defs" in w:
                label_file_list.append(w)
    # label_file_widgets = [w for w in master_widget_list if "rna_mod_defs" in w.pgv]
    # std_file_widgets = std_file_widgets + label_file_widgets
    std_file_list = std_file_list + label_file_list


def build_Pytheas_gui():
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
    for group, glist in option_dict.items():
        
        if group not in fixed_panels:
            panel = PytheasPanel(group, glist)
            option_panel_list.append(panel)
            panel.hide_panel()
        else:
            panel = PytheasFixedPanel(group, glist)
            fixed_panel_list.append(panel)
            # panel.show_panel()
            # print("panel geo", group, panel.geometry())
            
    #   Assemble GUI layout
    
    # add file panels at top
    # print("adding panels")
    # print()
    # print("***** IGNORE QT WARNINGS *****")
    for panel in fixed_panel_list:
    #     print("panel ", panel.group)
         main_window_layout.addWidget(panel)
    # print("***** END BOGUS QT WARNINGS *****")
    # print()
    
    
        # add module buttons
    
    button_dict = {"Global Variables:": {"Load Global Variables": Load_Global_Vars,
                                         "Load Previous Globals": Load_Previous_Globals,
                                         "List Global Variables": Print_Global_Vars,
                                         "Save Global Variables": Save_Global_Variables                                     
                                          },
                   # "Load Previous Files:": {"Load Digest Files": Load_Digest_Files,
                   #                          "Load Match Files": Load_Match_Files,
                   #                          "Print match_dict keys": Print_Match_Keys
                   #                          },
                   "Pytheas Modules: ": {"In Silico Digest": inSilicoDigest,
                                         "Match Spectra": matchSpectra, 
                                         # "Setup: Read Std Files": SetUp,
                                         # "Discovery": Discovery,
                                         # "Dev": Development,
                                         "Quit": Quit
                                         # "Statistics": statistics, 
                                         # "Final Report": report,
                                         # "Visualization": visualization
                                         },
                   "Load Pytheas Files": {"Load Previous Digest": LoadDigest,
                                          "Load Previous Match and Digest": LoadMatch}
                   }
    
    for label, bdict in button_dict.items():
        button_panel = PytheasButtonPanel(label, bdict)
        main_window_layout.addWidget(button_panel)
    
    # add master option panel
    
    option_panel_dict = {panel.group: panel for panel in option_panel_list}
    option_panel = PytheasOptionPanel(option_group_dict, option_panel_dict, floating_windows)    # option panel widget
    
    main_window_layout.addWidget(option_panel)
    
    if not floating_windows:    # turn this off to get floating windows
        for panel in option_panel_list:        
            main_window_layout.addWidget(panel)      
        
    all_panel_list = fixed_panel_list + option_panel_list
    master_widget_list = [w for panel in all_panel_list for w in panel.widget_list]
    
    pgv.widget_dict = {w.pgv: w for w in master_widget_list}   # for save/load global variables use
    
    # read standard definition files
    std_file_widgets = [w for w in master_widget_list if w.option_group == "input_files" and w.widget_type != "PytheasLabel"]
    label_file_widgets = [w for w in master_widget_list if "rna_mod_defs" in w.pgv]
    std_file_widgets = std_file_widgets + label_file_widgets

    return main_window, app
# pgv.pytheas_data_folder = os.path.join(pgv.working_dir,"pytheas_data_files")  # folder for json outputs
# if not os.path.isdir(pgv.pytheas_data_folder):
#     os.makedirs(pgv.pytheas_data_folder)

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


# if clargs.load_vars:
#     Load_Global_Vars()
print("Pytheas ready to Go!")
# if get_ipython().__class__.__name__ == 'SpyderShell':

#     print("This program is running inside Spyder") 
#     pgv.run = "Spyder"
# else:
#     print("This program is running on command line")
#     pgv.run = "CL"

pgv.run = "CL"
# if __name__ == "__main__":
print("launching Pytheas app...")
main_window, app = build_Pytheas_gui()
main_window.show()
app.exec()

# CL setup and run   

# pgv.pytheas_root = '/Users/jrwill/prog/pytheas_tk_interface/pytheas_root'
# # pgv.working_dir = '/Users/jrwill/prog/pytheas_tk_interface/pytheas_root/../Met_tRNA'
# # pgv.working_dir = '/Users/jrwill/prog/pytheas_tk_interface/pytheas_root/../Example_data'
# Load_Global_Vars()

# pgv.plot_sequence_map = 'y'
# inSilicoDigest()
# matchSpectra()





