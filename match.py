#!/usr/bin/python3

"""
Last update: April 2021
Author: Luigi D'Ascenzo, PhD - The Scripps Research Institute, La Jolla (CA)
Contact info: dascenzoluigi@gmail.com
GitHub project repository: https://github.com/ldascenzo/pytheas

"""
from datetime import datetime
import os
from pathlib import Path


from pytheas_global_vars import pgv, pgc
from pytheas_IO import read_pytheas_file, save_json_files
from match_functions import (match_spectra, unpack_match_dict, iTRAQ_quantitation, 
                                 output_match_dict_file, match_output_keys, 
                                 consolidated_match_output, make_sequence_plot, 
                                 make_long_sequence_plot, plot_ms2_spectra)

from match_functions import ppm_offset_plot,top_Sp_histogram, match_output_for_massacre

# nctr = 0
def match():    
 
    t1 = datetime.now()

    print()
    print("STEP 1:  READ MS DATA")
    read_pytheas_file("MS_data_file")  # read in MS data to pgv.ms2_dict
    pgv.n_ms2_keys = len(list(pgv.ms2_dict.keys()))

    print()
    print("STEP 2:  MATCH SPECTRA")
    pgv.ms2_match_dict = match_spectra() # primary matching
    
    print()
    print("STEP 3:  GENERATE CONSOLIDATED MATCH REPORT")
    pgv.unpacked_match_dict, pgv.top_match_dict = unpack_match_dict(pgv.ms2_match_dict)  # unpack ms2_key level to make flat dict with one entry per match
 
    if pgv.use_iTRAQ == 'y':
        iTRAQ_quantitation() # adds iTRAQ to unpacked dict

    match_job = "_" + pgv.match_job.split("_")[-1]
    output_match_dict_file(pgv.unpacked_match_dict, "match_output" + match_job)  # output the complete matching report
    if pgv.run == "CL":
        pgv.widget_dict["match_output"].load_value("match_output" + match_job + ".xlsx")
    match_output_keys(pgv.unpacked_match_dict) # makes list of match_dict keys to print out....for convenience

#TODO which of these are needed...        
    top_match_dict, pgv.seq_match_dict, pgv.match_dict = consolidated_match_output(pgv.unpacked_match_dict, "consolidated_match_output" + match_job)
 
    if pgv.massacre_input == "y":
        match_output_for_massacre()
    print()
    print("STEP 4: GENERATE SEQUENCE MAP AND PLOT SPECTRA")
    if pgv.run == "CL":
        pgv.widget_dict["consolidated_match_output"].load_value("consolidated_match" + match_job + ".xlsx")

    if pgv.plot_sequence_map == 'y':
        max_seq_len = max([len(mdict["seq3"]) for mdict in pgv.mol_dict.values()])
        if max_seq_len < 100:
            make_sequence_plot("match_sequence_map" + match_job) 
        else:
            make_long_sequence_plot("match_sequence_map" + match_job)
        
    if pgv.output_match_json == "y":  # write out json files 
        json_dir = os.path.join(pgv.job_dir, "pytheas_json_files")
        Path(json_dir).mkdir(parents=True, exist_ok=True)
        save_json_files(pgc.match_json, json_dir)
    
    print("pgv.plot_MS2_spectra: ", pgv.plot_MS2_spectra)
    if pgv.plot_MS2_spectra == 'y':
        plot_ms2_spectra()
 
    ppm_offset_plot()
    top_Sp_histogram()
 
    ctr, cid = 0, 0
    for ukey, udict in pgv.unique_precursor_dict.items():
        ctr +=1
        if "CID_ions" in udict:
            cid +=1
        
    print("# unique precursors = ", ctr, "# of matched precursors = ", cid)
    print("matching took :", datetime.now() - t1)
         
    return
