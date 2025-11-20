#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 31 15:36:58 2025

@author: jrwill
"""

import os
import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import hsv_to_rgb

from pytheas_global_vars import pgv, pgc
# from ms2_plot import labeled_matrix_plot

class Labeled_Matrix:
    def __init__(self, row_labels, col_labels):
        self.row_labels = row_labels
        self.nr = len(row_labels)
        self.col_labels = col_labels
        self.nc = len(col_labels)
        self.color_matrix = np.asarray([[pgc.white_hsv for i in range(self.nc)] for j in range(self.nr)])
        self.lw_matrix = np.asarray([[1 for i in range(self.nc)] for j in range(self.nr)])
        self.text_dict = {}
        self.title = ""
        self.output_file = None
  

def labeled_matrix_plot_new(lm, fs, labels, scale, ax):
    
    # nr, nc, _ = hsv_data.shape
    # colors = hsv_to_rgb(hsv_data)
    
    scale = 1/scale
    col_sc = scale # column scale (width)
    row_sc = scale        # row scale (height)
    nrs = lm.nr * row_sc
    ncs = lm.nc * col_sc
          
    rctr = 0.5 * row_sc # box centers
    cctr = 0.5 * col_sc
    roff = 0.1 * row_sc # label offsets
    coff = 0.1 * col_sc
    xbox = col_sc
    ybox = row_sc
    fss = fs * 0.75 * scale
    lws = 1 * scale
    
    # draw color grid
    for row in range(lm.nr):
        rs = (lm.nr-row-1) * row_sc # make rows go from top to bottom
        for col in range(lm.nc):
            cs = col * col_sc
            ax.plot(cs, rs,'k ') # have to plot points, or ax has no dimension scaled #$^&*#$
            color = hsv_to_rgb(lm.color_matrix[row][col])
            square = plt.Rectangle((cs, rs), xbox, ybox, facecolor = color,linewidth  = lws, edgecolor='black')
            ax.add_patch(square)

    # row and column labels
    for row in range(lm.nr):
        rs = (lm.nr-row-1) * row_sc   # make rows go from top to bottom
        if "left" in labels:
            ax.text(-coff ,rs + rctr, lm.row_labels[row], size = fss, ha = 'right', va = 'center') # left side
        if "right" in labels:
            ax.text(ncs + coff ,rs + rctr, lm.row_labels[row], size = fss, ha ='left', va = 'center') # right side
            
    for col in range(lm.nc):
        cs = col * col_sc
        if "bottom" in labels:
            ax.text(cs + cctr, -roff,  lm.col_labels[col], size = fss, ha = 'center', va = 'top', rotation = 'vertical') # bottom side
        if "top" in labels:
            ax.text(cs + cctr, nrs + roff,  lm.col_labels[col], size = fss, ha = 'center', va = 'bottom', rotation = 'vertical') # top side

    # text in boxes
    for key, td in lm.text_dict.items():
        rs = (lm.nr - td["row"] - 1) * row_sc
        cs = td["col"] * col_sc
        if len(td["text"]) > 3:
            fs_box = scale * fs/2
        else:
            fs_box = scale * fs
        
        ax.text(cs + cctr,rs +rctr,td["text"], size = fs_box, ha = 'center', va = 'center')

    # bold rectangles for cleavages
    for row in range(lm.nr):
        rs = (lm.nr-row-1) * row_sc # make rows go from top to bottom
        for col in range(lm.nc):
            if lm.lw_matrix[row,col] != 1:
                cs = col * col_sc
                lwx = lm.lw_matrix[row,col] * scale
                square = plt.Rectangle((cs, rs), xbox, ybox, fill = False,linewidth  = lwx, edgecolor='black')
                ax.add_patch(square)



        
def matrix_plot_new(lm):
    
    fs = 36  # fontsize for labels
    pfont = 'Arial'
    plt.rcParams.update({'font.size': fs, 'font.family': "sans-serif", "font.sans-serif": pfont})
    plt.rcParams.update({'font.size': fs, 'font.family': "monospace", "font.sans-serif": pfont})
    
    fig, ax = plt.subplots(layout="constrained")  # plot internal coords are 1 unit/row-col

    # nr, nc, _ = color_matrix.shape
    nr = lm.nr
    nc = lm.nc
    xsize = nc + 2*max([len(s) for s in lm.row_labels]) *fs/72 + 1
    ysize = nr + 2*max([len(s) for s in lm.col_labels]) *fs/72 + 1
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

    labeled_matrix_plot_new(lm, fs, pgv.match_sequence_labels, xyscale, ax)
    
    plt.title(lm.title, fontsize = fs / xyscale)
   
    fig_width, fig_height = plt.gcf().get_size_inches()
    pdffile = os.path.join(pgv.job_dir, lm.output_file + ".pdf")
    fig.savefig(pdffile, format='pdf', dpi=300)
    # pdffile = os.path.join(pgv.job_dir, lm.output_file + ".png")
    # fig.savefig(pdffile, format='png', dpi=300)

    plt.close()

def matrix_plot(color_matrix, row_labels, col_labels, text_dict, output_file, lw_matrix):
    
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

    labeled_matrix_plot(color_matrix, row_labels, col_labels, text_dict, fs, pgv.match_sequence_labels, xyscale, ax, lw_matrix)
   
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
    # print(n_matches, seq, col)
    return col

def pad_end3_gray_matrix(lm, slen):
             
    for row in range(lm.nr): # gray out cells > length of seq
        for col in range(lm.nc):
            idx = 100*row + col + 1
            if idx > slen:
                lm.color_matrix[row, col] = pgc.dark_gray # dark gray
                
def pad_end3_gray(lm):
    for mol, mdict in pgv.mol_dict.items(): # gray out entries with no sequence
        slen = len(mdict["seq3"])
        row_idx = lm.row_labels.index(mol)
        for i in range(slen, lm.nc): # gray
            lm.color_matrix[row_idx,i] = pgc.dark_gray


def color_mods_matrix(lm, molecule, color):
       for key, mdict in pgv.mod_dict.items():
           if mdict["Molecule"] != molecule:
               continue
           idx = mdict["Position"]
           row_idx = math.floor((idx-1)/lm.nc)
           col_idx = idx - 100 * row_idx - 1
           lm.color_matrix[row_idx,col_idx] = color
           
def color_mods(lm, color):
    for key, mdict in pgv.mod_dict.items():
        molecule = mdict["Molecule"]
        idx = mdict["Position"]
        if molecule in lm.row_labels:
            row_idx = lm.row_labels.index(molecule)
            col_idx = idx - 1
            lm.color_matrix[row_idx, col_idx] = color

def color_matrix_by_seq(lm, row_idx, fr, to, seq3, n_frags, length, cleavage_box):
    seq3_idx = 0     # sequence index within fragment
    for idx in range(fr - 1, to ):    # index in molecule                
        col_idx = idx
        # print(idx)
        lm.color_matrix[row_idx, col_idx] = sequence_color(n_frags, seq3[seq3_idx])
        seq3_idx += 1
        if cleavage_box and idx == to - 1 and idx != length - 1:
            lm.lw_matrix[row_idx, col_idx] = 7


def color_long_matrix_by_seq(lm, fr, to, seq3, n_frags, length, cleavage_box):
    seq3_idx = 0     
    for idx in range(fr, to + 1):                    
        row_idx = math.floor((idx-1)/lm.nc)
        col_idx = idx - 100 * row_idx - 1
        lm.color_matrix[row_idx, col_idx] = sequence_color(n_frags, seq3[seq3_idx])
        seq3_idx += 1
        if cleavage_box and idx == to:
            lm.lw_matrix[row_idx,col_idx] = 7

