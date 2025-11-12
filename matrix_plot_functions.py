#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 31 15:36:58 2025

@author: jrwill
"""

import os
import math

import matplotlib.pyplot as plt

from pytheas_global_vars import pgv, pgc
from ms2_plot import labeled_matrix_plot


def matrix_plot(color_matrix, row_labels, col_labels, text_dict, output_file):
    
    fs = 36  # fontsize for labels
    pfont = 'Arial'
    plt.rcParams.update({'font.size': fs, 'font.family': "sans-serif", "font.sans-serif": pfont})
    plt.rcParams.update({'font.size': fs, 'font.family': "monospace", "font.sans-serif": pfont})
    
    fig, ax = plt.subplots(layout="constrained")  # plot internal coords are 1 unit/row-col

    nr, nc, _ = color_matrix.shape
    xsize = nc + 2*max([len(s) for s in row_labels]) *fs/72 + 1
    ysize = nr + 2*max([len(s) for s in col_labels]) *fs/72 + 1
    xscale = xsize
    yscale = ysize
    
    if pgv.scale_matrix_plots == 'y':
        xtarget = 8
        ytarget = 11
        xscale = xsize/xtarget
        yscale = ysize/ytarget
        xyscale = max(xscale, yscale, 1.0)
    else:
        xyscale = 1.0
    
    print("plot size: ", xsize/xyscale, ysize/xyscale, "scale", xyscale)
    fig.set_size_inches(xsize/xyscale, ysize/xyscale) # absolute plot size for display in inches
    plt.rcParams["figure.autolayout"] = True
    plt.axis('off')
    ax.set_aspect('equal')

    labeled_matrix_plot(color_matrix, row_labels, col_labels, text_dict, fs, pgv.match_sequence_labels, xyscale, ax)
   
    fig_width, fig_height = plt.gcf().get_size_inches()
    pdffile = os.path.join(pgv.job_dir, output_file + ".pdf")
    fig.savefig(pdffile, format='pdf', dpi=300)

    plt.close()


def make_mod_text_dict(row_labels, col_labels, find_row, molecule, nc):
    mod_text_dict = {} # text for modified positions
    idx = 0 
    for mod, mdict in pgv.mod_dict.items():
        mpos = str(mdict["Position"])
        mol = mdict["Molecule"]
        if find_row:
            col_idx = col_labels.index(mpos)
            midx = col_idx
            row_idx = row_labels.index(mol)
        else:
            if mol == molecule:
                midx = mdict["Position"]
                row_idx = math.floor((midx-1)/nc)
                col_idx = midx - 100 * row_idx - 1
            else:
                continue

        mbase = pgv.nt_def_dict[mdict["ID"]]["Pytheas_ID"]
        mod_text_dict[midx] = {"row": row_idx, "col": col_idx, "text": mbase}
        idx += 1
    return mod_text_dict

def make_seq_text_dict(row_labels, col_labels, find_row, molecule, nc):
    seq_text_dict = {} # text for modified positions
    idx = 1 
    
    for base in pgv.mol_dict[molecule]["raw_seq"]:
        if find_row:
            col_idx = idx -1
            row_idx = row_labels.index(molecule)
        else:
            # if idx == 0:
            #     print("first base = ", base)
        # for mod, mdict in pgv.mod_dict.items():
                # midx = mdict["Position"]
            row_idx = math.floor((idx-1)/nc)
            col_idx = idx - 100 * row_idx - 1
        # else:
        #     continue

        # mbase = pgv.nt_def_dict[mdict["ID"]]["Pytheas_ID"]
        seq_text_dict[idx] = {"row": row_idx, "col": col_idx, "text": base}
        idx += 1
    return seq_text_dict

def make_text_dict(row_labels, col_labels, nc):
    seq_text_dict = {} # text for modified positions
    idx = 0
    
    rc_index_dict = {}

    if pgv.base_labels_seq_map == "y":
        for mol, mdict in pgv.molecule_dict.items():
            row_idx = row_labels.index(mol)
            midx = 1 

            for base in pgv.mol_dict[mol]["raw_seq"]:
                col_idx = midx - 1
                seq_text_dict[idx] = {"row": row_idx, "col": col_idx, "text": base}
                rc_index_dict[str(col_idx) + "_" + str(row_idx)] = idx  # index to find mod positions for over_write
                midx += 1
                idx += 1
                
    for mod, mdict in pgv.mod_dict.items():
        mpos = str(mdict["Position"])
        mol = mdict["Molecule"]

        col_idx = col_labels.index(mpos)
        row_idx = row_labels.index(mol)
        if str(col_idx) + "_" + str(row_idx) in rc_index_dict:
            midx = rc_index_dict[str(col_idx) + "_" + str(row_idx)]
        else:
            midx = idx
        mbase = pgv.nt_def_dict[mdict["ID"]]["Pytheas_ID"]
        seq_text_dict[midx] = {"row": row_idx, "col": col_idx, "text": mbase}
        idx += 1

    return seq_text_dict


def make_long_text_dict(row_labels, col_labels, molecule, nc):
    seq_text_dict = {} # text for modified positions
    idx = 1 
    
    if pgv.base_labels_seq_map == "y":
        for base in pgv.mol_dict[molecule]["raw_seq"]:
            row_idx = math.floor((idx-1)/nc)
            col_idx = idx - 100 * row_idx - 1
            seq_text_dict[idx] = {"row": row_idx, "col": col_idx, "text": base}
            idx += 1
            
    for mod, mdict in pgv.mod_dict.items():
        mol = mdict["Molecule"]
        if mol == molecule:
            midx = mdict["Position"]
            row_idx = math.floor((midx-1)/nc)
            col_idx = midx - 100 * row_idx - 1
        else:
            continue
     
        mbase = pgv.nt_def_dict[mdict["ID"]]["Pytheas_ID"]
        seq_text_dict[midx] = {"row": row_idx, "col": col_idx, "text": mbase}


    return seq_text_dict

    
def sequence_color(n_matches, seq):
    col = pgc.dark_gray
    if n_matches == 1 and seq in pgc.natural:
        col = pgc.dark_blue
    if n_matches > 1 and seq in pgc.natural:
        col = pgc.light_blue
    if n_matches == 1 and seq not in pgc.natural:
        col = pgc.dark_green
    if n_matches > 1 and seq not in pgc.natural:
        col = pgc.light_green
    return col
