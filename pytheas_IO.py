#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 10:08:53 2024

@author: jrwill
"""
import os
import json

from Bio import SeqIO
# from iteration_utilities import duplicates
from pyteomics import mgf,mzml
import pandas as pd
import pickle

from pytheas_global_vars import pgv, pgc, pgvdict
from pytheas_objects import nt_fragment, MS2_spectrum, molecule, atom, residue
from Pytheas_Qt_widgets import PytheasPanel

def read_pytheas_file(file):
    error = None
    base = getattr(pgv, file)
    print("working_dir", pgv.working_dir)
    for parent_dir in [pgv.working_dir, pgv.user_dir, pgv.pytheas_dir]:
        infile = os.path.join(parent_dir, base)
        print(infile)
        if not os.path.isfile(infile):
            continue
        else:
            ffunc = pgvdict[file]["file_read_function"]
            method = globals()[ffunc]
            method(infile)
            print("read ", infile)
            return error
    
    error = [" unable to open file!", "please find file: "  + infile]
    return error

    
def read_atomic_mass_file(file):
    with open(file,'r') as jfile:
        pgv.atomic_dict = json.load(jfile)
    pgv.hmass = pgv.atomic_dict["H1"]["am"]
    # pgv.atomic_mass_dict = {key: atom(adict) for key, adict in pgv.atomic_dict.items()}

"""
Masses from:
(1) http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=html&isotype=some
(2) http://www.sisweb.com/referenc/source/exactmaa.htm
"""

def read_nt_def_file(file):
    nt_def_df = pd.read_excel(file)
    pgv.nt_key_dict = {"Pytheas_ID":{}, "Modomics_1":{}, "Pytheas_1":{}}  # dict to translate 3 letter code to pytheas 1 letter code
    pgv.end_dict = {"end5":{}, "end3":{}}
    for i, row in enumerate(nt_def_df.itertuples(), 1):
        if pd.isna(row.AMBER_3_letter_code):
            continue
        code3 = row.AMBER_3_letter_code
        code1 = row.Pytheas_1_letter_code
        synonym = row.Pytheas_ID
        mod_1 = row.Modomics_1
        
        pgv.nt_key_dict["Pytheas_1"][code1] = code3   
        pgv.nt_key_dict["Pytheas_ID"][synonym] = code3
        # pgv.nt_key_dict[code3] = code3
        pgv.nt_key_dict["Modomics_1"][mod_1] = code3
        
        if "end" in row.Type:
            pgv.end_dict[row.Type][row.Pytheas_ID] = row.AMBER_3_letter_code
    
    pgv.nt_def_dict = nt_def_df.set_index('AMBER_3_letter_code').to_dict(orient="index")   # has all nucleotide names/mods and atom stoichiometries
    
    pgv.residue_dict = {key: residue(nt_dict) for key, nt_dict in pgv.nt_def_dict.items()}

def read_nt_mod_def_file(file):
    try:
        pgv.nt_fragment_dict
    except:
        pgv.nt_fragment_dict = {}
    if "light" in file:
        label = "light"
        print("label is ", label)
    if "heavy" in file:
        label = "heavy"
        print("label is ", label)

    nt_sheet_dict = pd.read_excel(file, None)
    pgv.nt_fragment_dict[label] = {nt: nt_fragment(nt,nt_df,label) 
                        for nt,nt_df in nt_sheet_dict.items()}

# def read_topo_def_file(file):
#     topo_sheet_dict = pd.read_excel(file, None)
#     topo_df = topo_sheet_dict["topology"]
#     pgv.topo_dict = topo_df.set_index('group').to_dict('index')
#     losses_df = topo_sheet_dict["secondary_losses"]
#     pgv.losses_dict = losses_df.set_index('fragment').to_dict('index')
#     for loss,ldict in pgv.losses_dict.items():
#         if ldict["loss_mass_list"] == "none":
#             ldict["loss_mass"] = 0
#         else:
#             atom_list = ldict["loss_mass_list"].split(",")
#             ldict["loss_mass"] = sum([pgv.atomic_dict[a]["am"]for a in atom_list])

def read_topo_def_file(file): # from discovery_functions.py
    topo_sheet_dict = pd.read_excel(file, None)
    topo_df = topo_sheet_dict["normal_topology"]
    pgv.topo_dict = topo_df.set_index('group').to_dict('index')
    
    patch_df = topo_sheet_dict["patch_topology"]
    pgv.patch_dict = patch_df.set_index('added_edge').to_dict('index') # not used??

    child_df = topo_sheet_dict["child_topology"]
    pgv.child_dict = child_df.set_index('group').to_dict('index') # not used??
 
    losses_df = topo_sheet_dict["secondary_losses"]
    pgv.losses_dict = losses_df.set_index('fragment').to_dict('index')
    for loss,ldict in pgv.losses_dict.items():
        if ldict["loss_mass_list"] == "none":
            ldict["loss_mass"] = 0
        else:
            atom_list = ldict["loss_mass_list"].split(",")
            ldict["loss_mass"] = sum([pgv.atomic_dict[a]["am"]for a in atom_list])


def parse_enzyme_pattern(pattern): # return list of lists to iterate over
    pattern_list = [s.replace("]","").split(",") for s in pattern.replace(" ","").replace('[','').split('],')]
    return pattern_list

def read_enzyme_def_file(file):
    enzyme_sheet_dict = pd.read_excel(file,None)
    enzyme_df = enzyme_sheet_dict["enzymes"]
    pgv.enzyme_dict = enzyme_df.set_index("name").to_dict('index')
    
    for valdict in pgv.enzyme_dict.values(): # split elements into lists
        valdict["pattern"] = parse_enzyme_pattern(valdict["pattern"])
        valdict["end_3"] = valdict["end_3"].split(",")
        valdict["end_5"] = valdict["end_5"].split(",")
        valdict["cut_idx"] = int(valdict["cut_idx"])

def read_custom_cleavage(file):
    cleavage_sheet_dict = pd.read_excel(file,None)
    cleavage_df = cleavage_sheet_dict["enzymes"]
    pgv.custom_cleavage_dict = cleavage_df.set_index("name").to_dict('index')
    
    for valdict in pgv.custom_cleavage_dict.values(): # split elements into lists
        valdict["pattern"] = parse_enzyme_pattern(valdict["pattern"])
        valdict["end_3"] = valdict["end_3"].split(",")
        valdict["end_5"] = valdict["end_5"].split(",")
        valdict["cut_idx"] = int(valdict["cut_idx"])

def read_charge_file(file):
    chg_sheet_dict = pd.read_excel(file, None)
    pgv.MS_charge_dict = {chg: {int(l):[int(zz) for zz in str(z).split(",")] 
                for l,z in zip(chg_df.Length, chg_df.Charges)} 
                for chg, chg_df in chg_sheet_dict.items()}

# def read_modification_sets(pgvname, setfile):
    
#TODO  read this at startup to avoid rebuilding the panel
def read_modification_sets(file):
    # global set_dict
    # print("initial mod set", modification_set)
    set_sheet_dict = pd.read_excel(file, None)
    set_df = set_sheet_dict["modification_sets"]
    pgv.set_dict = {"natural": ["ADE", "CYT", "GUA", "URI"]}
    sets = [set_name for set_name in set_df.columns if "set" in set_name]
    print("sets in df: ", sets)
    for set_name in sets:
        # pgv.set_dict[set_name] = [base for base in list(set_df[set_name]) if str(base) != 'nan' and base in pgv.nt_key_dict]
        pgv.set_dict[set_name] = [base for base in list(set_df[set_name]) if str(base) != 'nan']
        print("   set name ", set_name, pgv.set_dict[set_name])
    # set_widget = pgv.widget_dict["modification_set"]  #noqa
    # cb_set = [cb.text() for cb in set_widget.cb_list]
    # pgvdict["modification_set"]["option_list"] = sets
    sets = pgv.set_dict.keys()
    pgvdict["modification_set"]["option_list"] = ",".join(sets)
    # pgv.option_panel_dict["Modifications"] = PytheasPanel("Modifications", pgv.option_dict["Modifications"])
    print("read_modification_set: pgv.set_dict", pgv.set_dict)


def write_json(j_dict, j_file):
    j = json.dumps(j_dict, indent=4)
    with open(j_file, 'w') as f:
        print(j, file=f)

def read_json(j_file):
    print("reading json file: ", j_file)
    with open(os.path.normpath(j_file), 'r') as f:
        j_dict = json.load(f)
    return j_dict
  
def read_output_format(file):
    # global output_format_dict
    output_dict = pd.read_excel(file,None)
    pgv.output_format_dict = {okey: odf.set_index("variable").to_dict("index") for okey, odf in output_dict.items()}   

def read_fasta_file(fasta_file): # digest
    # global molecules, mol_dict
    print("read_fasta_file", fasta_file)
    sequences = SeqIO.parse(fasta_file,"fasta")
    # print("read_fasta_file: sequences", [seq for seq in sequences])
    sequence_list = [[seq.id, str(seq.seq).replace("_","")] for seq in sequences]
    pgv.mol_dict = {}
    for seq_id, seq in sequence_list:
        if seq_id not in pgv.mol_dict:
            pgv.mol_dict[seq_id] = {
                "raw_seq": seq,  # this may be unmodified seq or seq with [mod]   
                # "seq_list": parse_mod_seq(seq),
                "mol_5_end": pgv.mol_end5, 
                "mol_3_end": pgv.mol_end3,
                "mods": {}}
        else:
            print("Duplicate sequence ID: ", seq_id)
            
    pgv.molecule_dict = {key: molecule(mdict) for key, mdict in pgv.mol_dict.items()}
    pickle.dump(pgv.molecule_dict, open(os.path.join(pgv.job_dir, "molecule_dict.pkl"), "wb" ))

     
def read_mod_file(mod_file):
    pgv.mod_dict = pd.read_excel(mod_file).to_dict("index")

def read_MS_data_file(file):
    ext = os.path.splitext(file)[-1]
    
    if ext == ".mzML":
        pgv.ms2_dict = read_mzML_file(file)
    elif ext == ".mgf":
        pgv.ms2_dict = read_mgf_file(file)
    else:
        print("unsupported data format: ", ext)
        print("current support for .mzML and .mgf")
 
    if pgv.Xcorr =='y':  #noqa
        print("Using Xcorr")

    for key, spectrum in pgv.ms2_dict.items():
        
        if pgv.precursor_window_removal > 0.0:  # in ppm  #noqa
            pre_win = spectrum.mz1 * pgv.precursor_window_removal/1000000.0  #noqa
            spectrum.remove_precursor(spectrum.mz1, pre_win)

        spectrum.normalize_intensities(pgv.normalize_MS2)  #noqa
        
        if pgv.MS2_peak_threshold > 0.0:  #noqa
            spectrum.threshold(pgv.MS2_peak_threshold)  #noqa
            
        if pgv.squash_spectrum == 'y':
            spectrum.squash_spectrum()
            
        if pgv.Xcorr =='y':  #noqa
            # print("Using Xcorr")
            spectrum.fft_MS2()

    # make a text copy for saving to .json 
    att_list = ["mz1", "rt", "z", "ms2", "max_int"]  # don't include raw ms2 data, saving 100's of MB
    pgv.ms2_key_dict = {}
    for key, spec in pgv.ms2_dict.items():
        pgv.ms2_key_dict[key] = {att:getattr(spec,att)for att in att_list} 
    
def read_mzML_file(file):
    
    data_dict = {}
    key = 0
    eof = 0
    with mzml.read(file) as reader:
       while eof == 0:
            
            try:
                scan = next(reader)
            except:
                eof = 1
                break
            
            try:
                selectedIon = scan["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]
                mz1 = selectedIon["selected ion m/z"]
                z = selectedIon["charge state"]
                peak_in = selectedIon["peak intensity"]
                rt = scan["scanList"]["scan"][0]["scan start time"]
                mz = scan["m/z array"]
                ia = scan["intensity array"]
                print(key,mz1, peak_in, rt)
                data_dict[key] = MS2_spectrum(mz1, rt, z, mz, ia)
                key += 1
            except:
                continue
   
    print("MZML file: ", len(data_dict), " precursor ions read")
        
    return data_dict

def read_mgf_file(file):
    mgf_file=open(file)
    mgf_data_dict = mgf.read(mgf_file)  # indexed MGF objectMS2_max read using pyteomics
    print("mgf_data_dict type: ", type(mgf_data_dict))
    data_dict = {}
    key = 0
    for entry in mgf_data_dict:
        if pgv.max_spectra != "all":  #noqa
            if key > int(pgv.max_spectra): #noqa
                break

        scan = entry["params"]["title"]
        # print(entry)
        pm, slug = entry["params"]["pepmass"]
        mz1 = float(pm)
        rt = float(entry["params"]["rtinseconds"])
        try:
            z = int(entry["params"]["charge"][0])
        except:
            print("no charge for scan", scan)
            z =0
        mz = entry["m/z array"] # numpy array
        ia = entry["intensity array"] # numpy array
        data_dict[key] = MS2_spectrum(mz1, rt, z, mz, ia)
        key += 1
    print("MGF file: ", len(data_dict), " precursor ions read")
        
    return data_dict

# TODO fix this hot mess
# def read_modification_sets(setfile):
#     global set_dict
#     print("initial mod set", modification_set)
#     set_df = pd.read_excel(setfile)
#     set_dict = {"natural": ["ADE", "CYT", "GUA", "URI"]}
#     sets = [set_name for set_name in set_df.columns if "set" in set_name]
#     for set_name in sets:
#         set_dict[set_name] = [base for base in list(set_df[set_name]) if str(base) != 'nan' and base in pgv.nt_key_dict["Pytheas_ID"]]
#     set_widget = widget_dict["modification_set"]  #noqa
#     cb_set = [cb.text() for cb in set_widget.cb_list]
#     print("initial cb set ", cb_set)
#     pgvdict["modification_set"]["option_list"] = sets
#     ctr = len(set_widget.cb_list) + 1
#     for set_name in sets:
#         if set_name not in cb_set:
#             cb = QCheckBox(set_name)
#             set_widget.cb_list.append(cb)
#             set_widget.layout.addWidget(cb, 0, ctr)
#             cb.stateChanged.connect(set_widget.set_value)
#             ctr += 1
#     print(" after read set", modification_set)

#     print("loading json files from : ", load_dir)

def load_json_files(module, load_dir):
    # print("digest_json_dict", pgc.digest_json_dict)
    
    for jroot in module:
        jfile = jroot + ".json"
        jvar = jroot + "_dict"
        setattr(pgv, jvar, read_json(os.path.join(load_dir,"pytheas_json_files", jfile)))

def save_json_files(module, json_dir):
    print("json files for ", module)
    for jroot in module:
        jfile = jroot + ".json"
        jvar = getattr(pgv, jroot + "_dict")
        write_json(jvar, os.path.join(json_dir, jfile))
        
        
        
