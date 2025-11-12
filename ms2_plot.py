#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 30 10:27:59 2023

@author: jrwill
"""

import os
import io
import networkx as nx
import numpy as np

import math
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb
import matplotlib.style as mplstyle
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rcParams
from PIL import Image

from PyQt5.QtWidgets import (QVBoxLayout, QWidget, QLabel, QFrame, QHBoxLayout)
from PyQt5.QtGui import QPixmap, QImage
from PyQt5.QtCore import Qt


from pytheas_global_vars import  pgv, pgc
from mod_seq_functions import generate_mod_seq

class ion_text_label:
    def __init__(self, text, text_object, x, y):
        self.text = text
        self.to = text_object
        self.sx = x  # coords of top of peak
        self.sy = y
        self.tx = x # position to draw text
        self.ty = y

    def get_bb(self,r,ax):
        self.bb = self.to.get_window_extent(renderer=r).transformed(ax.transData.inverted())
        # coords are in plot coords: x ~[100, 1600], y ~[0, 6]
        # bb x0, y0 is lower left corner...
        self.height = self.bb.ymax - self.bb.y0
        self.width = self.bb.xmax - self.bb.x0
        self.x0 = self.bb.x0 # center label over peak
        self.xmax = self.bb.xmax
        self.y0 = self.bb.y0
        self.ymax = self.bb.ymax
        self.xc = self.x0 + self.width/2 # center for drawing connector lines
        self.yc = self.y0 + self.height/2
        
    def offset_text(self, delx,dely): # think tx and x0 are the same...
        # print("offset_text", delx,dely,type(delx),type(dely))
        # print(self.tx,type(self.tx))
        self.tx = self.tx + delx
        self.ty = self.ty + dely
        self.x0 = self.x0 + delx
        self.y0 = self.y0 + dely
        self.xc = self.xc + delx
        self.yc = self.yc + dely
        
    def center_x(self):
        self.offset_text(-self.width/2,0.0)


class ImageWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.image_label = QLabel()
        layout = QVBoxLayout()
        layout.addWidget(self.image_label)
        self.setLayout(layout)


    def set_image(self, image_data):
            if isinstance(image_data, np.ndarray):
                if image_data.dtype == np.uint8:
                    height, width, channel = image_data.shape
                    bytes_per_line = 3 * width
                    q_image = QImage(image_data.data, width, height, bytes_per_line, QImage.Format_RGB888).rgbSwapped()
                else:
                     raise ValueError("Numpy array must have dtype uint8")
            elif isinstance(image_data, QImage):
                q_image = image_data
            else:
                raise TypeError("Image data must be a numpy array or a QImage object")
            pixmap = QPixmap.fromImage(q_image)
            self.image_label.setPixmap(pixmap)
            self.image_label.setScaledContents(True)

class PytheasImagePanel(QWidget):
    def __init__(self, label, image_data):
        super().__init__()
        self.panel_layout = QVBoxLayout()   # overall layout for widget set
        
        self.frame = QFrame()  # frame widget to hold widget set
        self.frame_layout = QVBoxLayout(self.frame)
        self.frame.setFrameShape(QFrame.Panel)
        self.label = QLabel(label)
        self.label.setFixedWidth(200)
        self.label.setStyleSheet("font-weight: bold")
        self.frame_layout.addWidget(self.label)

        # if isinstance(image_data, np.ndarray):
        #     if image_data.dtype == np.uint8:
        #         height, width, channel = image_data.shape
        #         bytes_per_line = 3 * width 
        #         print("image_data h, w, bpl ", height, width, channel, bytes_per_line)
        #         q_image = QImage(image_data.data, width, height, bytes_per_line, QImage.Format_RGB888).rgbSwapped()
        #     else:
        #          raise ValueError("Numpy array must have dtype uint8")
        # elif isinstance(image_data,PIL.PngImagePlugin.PngImageFile):
        print("image panel size: ", image_data.size)
        data = image_data.tobytes("raw", "RGBA")
        q_image = QImage(data, image_data.size[0], image_data.size[1], QImage.Format_RGBA8888)

        # elif isinstance(image_data, QImage):
        #     q_image = image_data
        # else:
        #     raise TypeError("Image data must be a numpy array or a QImage object")
        pixmap = QPixmap.fromImage(q_image)
        # scaled_pixmap = pixmap.scaled(800, 600, Qt.KeepAspectRatio)

        self.image_label = QLabel()
        self.image_label.setPixmap(pixmap)
        self.image_label.setScaledContents(True)
        # self.image_label.resize(800,800)
        self.frame_layout.addWidget(self.image_label)

        self.panel_layout.addWidget(self.frame)
        self.setLayout(self.panel_layout)


def labeled_matrix_plot(hsv_data, row_labels, col_labels, text_dict, fs, labels, scale, ax):
    
    nr, nc, _ = hsv_data.shape
    colors = hsv_to_rgb(hsv_data)
    
    col_sc = scale # column scale (width)
    row_sc = scale        # row scale (height)
    nrs = nr * row_sc
    ncs = nc * col_sc
          
    rctr = 0.5 * row_sc # box centers
    cctr = 0.5 * col_sc
    roff = 0.1 * row_sc # label offsets
    coff = 0.1 * col_sc
    xbox = col_sc
    ybox = row_sc
    fss = fs * 0.75
    lws = 1 * scale
    
    # draw color grid
    for row in range(nr):
        rs = (nr-row-1) * row_sc # make rows go from top to bottom
        for col in range(nc):
            cs = col * col_sc
            ax.plot(cs, rs,'k ') # have to plot points, or ax has no dimension scaled #$^&*#$
            color = colors[row][col]
            square = plt.Rectangle((cs, rs), xbox, ybox, facecolor = color,linewidth  = lws, edgecolor='black')
            ax.add_patch(square)

    # row and column labels
    for row in range(nr):
        rs = (nr-row-1) * row_sc   # make rows go from top to bottom
        if "left" in labels:
            ax.text(-coff ,rs + rctr, row_labels[row], size = fss, ha = 'right', va = 'center') # left side
        if "right" in labels:
            ax.text(ncs + coff ,rs + rctr, row_labels[row], size = fss, ha ='left', va = 'center') # right side
            
    for col in range(nc):
        cs = col * col_sc
        if "bottom" in labels:
            ax.text(cs + cctr, -roff,  col_labels[col], size = fss, ha = 'center', va = 'top', rotation = 'vertical') # bottom side
        if "top" in labels:
            ax.text(cs + cctr, nrs + roff,  col_labels[col], size = fss, ha = 'center', va = 'bottom', rotation = 'vertical') # top side

    # text in boxes
    for key, td in text_dict.items():
        rs = (nr - td["row"] - 1) * row_sc
        cs = td["col"] * col_sc
        if len(td["text"]) > 3:
            fs_box = fs/2
        else:
            fs_box = fs
        ax.text(cs + cctr,rs +rctr,td["text"], size = fs_box, ha = 'center', va = 'center')



def series_name(m):
    sn = m["series"] + "(" + pgv.ion_mode + str(abs(int(m["z"]))) + ")"
    return sn

def ion_series_matrix(ukey, ax, tfs):
    udict = pgv.unpacked_match_dict[ukey]
    ax.axis('off')
    pfont = 'Arial'
        
    ion_series = ["a","a-B","b","c","d","w","x","y","y-P","z","z-P"]
    a_series = ["a","a-B","b","c","d"]
    w_series = ["w","x","y","y-P","z","z-P"]
    ion_hues = [0.33, 0.33, 0.66, 0.84, 0.0, 0.33, 0.66, 0.84, 0.84, 0.0, 0.0]
    ion_hue_dict = {s:c for s,c in zip(ion_series,ion_hues)}

    f3 = udict["frag3"]
    seq3 =f3[1:-1]
    seq3_len = len(seq3)
 
    matches = udict["matched_ions"]
    cid = udict["CID_ions"] # list of theoretical 

    # list of all series_charge states
    szlist = sorted(list(set([series_name(m) for m in matches.values() if m["series"] in ion_series])))

    # left and right labels
    abcd_idx = [i + 1 for i in range(seq3_len)] # pos 0 + [1,2,3,4...n] + (n+1)
    wxyz_idx = [seq3_len + 1 - i for i in abcd_idx]
    bilist = [pgv.nt_def_dict[b]["Pytheas_ID"] if b in pgv.nt_def_dict else "" for b in seq3]
    alist = [str(a) for a in abcd_idx] # indices for abcd
    wlist = [str(w) for w in wxyz_idx] # indices for wxyz
    
    nr = len(abcd_idx)
    nc = len(szlist)
    iarray = -np.ones((nr,nc))  # initialize intensity array

    # make text dict with ms2 values
    text_dict = {}
    idx = 0
    for mkey, m in cid.items():
        sn = series_name(m)
        if sn not in szlist:
            continue
        if m["series"] in a_series and m["index"] in abcd_idx:
            ri = abcd_idx.index(m["index"])
        elif m["series"] in w_series and m["index"] in wxyz_idx:
            ri = wxyz_idx.index(m["index"])
        else:
            continue
        ci = szlist.index(sn)
        mz2 = round(m["mz2"],2)
        text_dict[idx] = {"row": ri, "col": ci, "text": str(mz2)}
        iarray[ri, ci] = 0.0 # set to indicate theo ion present

        idx += 1
        
    # find max intensity for normalization
    match_int_list = [mdict["obs_int"] for m, mdict in matches.items() if mdict["series"] in ion_series]
    if len(match_int_list) > 0:
        match_max_int = max(match_int_list)   
    else:
        match_max_int = 1.0
  
    for mkey, m in matches.items():  # set intensity array for matches
        sn = series_name(m)
        if sn not in szlist:
            continue
        if m["series"] in a_series and m["index"] in abcd_idx:
            ri = abcd_idx.index(m["index"])
        elif m["series"] in w_series and m["index"] in wxyz_idx:
            ri = wxyz_idx.index(m["index"])
        else:
            continue
        ci = szlist.index(sn)
        iarray[ri, ci] = m["obs_int"]/match_max_int

    color_matrix = np.zeros((nr,nc, 3))  # set up color matrix
            
    for row in range(nr): # color image array
        irow = nr - row -1 # first row at top, last row at bottom
        for col in range(nc):
            sat = iarray[irow,col]
            if sat == 0.0:  # theo ion not matched
                color_matrix[irow, col] = [0.0, sat, 0.75] #light gray
            elif sat == -1:  # no theo ion
                color_matrix[irow, col] = [0.0,0.0,0.5] # dark gray,
            else: # matched ion
                color_matrix[irow, col] = [ion_hue_dict[szlist[col][0]], sat, 1.0]
 
    if nc == 0 or nr == 0:
        return 

    labeled_matrix_plot(color_matrix, [], [], text_dict, tfs, [], 1.0, ax) # was fsbox
    # plot the labels
    
    rctr = 0.5 # box centers
    cctr = 0.5
    rq = rctr/2.0 # quarter of a box
    cq = cctr/2.0
 
    for row in range(nr): # row labels have base-idx and abcd idx
        irow = nr - row - 1 # rows start at top, going down to 0
        plt.text(-2*cctr, irow + rctr , bilist[row], size = tfs, ha='right', va = "center", font = pfont, fontweight = 'bold')
        plt.text(-cq, irow + rctr , alist[row], size = tfs, ha='right', va = "center", font = pfont, style = 'italic')
        plt.text(nc + cq, irow + rctr, wlist[row], size = tfs, ha='left', va = "center", font = pfont, style = 'italic') 
    for col in range(nc): # column labels 
        plt.text(col + cctr, nr + rq , szlist[col] + " ", size = tfs * 0.75, ha='center', va = "bottom" , rotation = 90, font = pfont)

    plt.text(-cctr, nr + rq , "abcd", size = tfs,  ha='center', va = "bottom", rotation = 90, font = pfont, style = 'italic')
    plt.text(nc + cctr, nr + rq , "wxyz", size = tfs,  ha='center', va = "bottom", rotation = 90, font = pfont, style = 'italic')

def ms2_plot_title(ukey):
    udict = pgv.unpacked_match_dict[ukey]
    pre_mz = str(round(udict["mz1"],3))
    rt = str(round(udict["rt"],1))
    pre_z = str(udict["z"])
    seq = generate_mod_seq(udict["frag3"][1:-1])
    ms2_key = udict["ms2_key"]
    Sp = str(round(udict["Sp"],3))
    title = "Precursor M/z: %s, spec# %i, seq =  %s, z = %s, RT = %s, Sp = %s" %(pre_mz, ms2_key, seq, pre_z, rt, Sp)
    plot_filename = "_".join([seq,pre_mz,pre_z,rt]) + ".pdf"

    return title, plot_filename

def plot_ms2_spectrum(ukey):
    udict = pgv.unpacked_match_dict[ukey]  
    ms2_key =udict["ms2_key"]
    ion_dict = udict["matched_ions"]

    # constants for line annotations
    del0 = 0.1 # space between (spectrum or text) and connector
    del1 = 0.1 # length of terminal connector
    del2 = 0.3
    delx = 5
    
    unmatched_linewidth = 0.25
    matched_linewidth = 0.5
    conn_linewidth = 0.5 # width of label connector lines
    dely = 0.5 # increment for shifting overlapped labels
    
    # font setup
    pfont = 'Arial'
    fs = 9    #fontsize for peak labels in points
    axis_fs = 12
    rcParams.update({'font.size': fs, 'font.family': "sans-serif", "font.sans-serif": pfont})
    label_font_size = 10 # tick label size
    plt.rcParams['xtick.labelsize'] = label_font_size
    plt.rcParams['ytick.labelsize'] = label_font_size
    mplstyle.use('fast')
    
    title, plot_filename = ms2_plot_title(ukey)
     
    # experimental spectrum
    mz2 = list(pgv.ms2_dict[ms2_key].ms2.keys()) # entire MS2 spectrum
    if pgv.log_ms2_plot == 'y':
        yplot_max = 1.2*math.log(pgv.normalize_MS2) + 1.0
        yplot_min = math.log(pgv.MS2_peak_threshold)    
        it2 = [math.log(x) for x in list(pgv.ms2_dict[ms2_key].ms2.values())] # intensities
        y_label = "Ln[Intensity]"
        dely = 0.5 # increment for shifting overlapped labels
    
    else:
        yplot_max = 1.2 * pgv.normalize_MS2
        yplot_min = 0
        it2 = [x for x in list(pgv.ms2_dict[ms2_key].ms2.values())]
        y_label = "Intensity"
        dely = 5
    
    # matched spectrum
    mzm = [ikey["obs_mz2"] for ikey in ion_dict.values()] # mz2 list
    if pgv.log_ms2_plot == 'y':
        itm = [math.log(ikey["obs_int"]) for ikey in ion_dict.values()] # intensity list
    else:
        itm = [ikey["obs_int"] for ikey in ion_dict.values()]
    ilab = [ikey["series"] + r'$%s_{%s}$' %("",str(ikey["index"])) + " (" + pgv.ion_mode + str(abs(ikey["z"])) + ")"  for ikey in ion_dict.values()]
    xyz = sorted(zip(mzm,itm,ilab))
    mzms = [x for x,y,z in xyz]
    itms = [y for x,y,z in xyz]
    ilabs = [z for x,y,z in xyz]
    
    xsize = 12  # spectrum plot size in inches
    ysize = 6
    
    #   initial dummy plot to get position of text labels.  Matplotlib needs to render before bounding box can be determined
    tdict = initialize_label_positions(mz2, it2, mzms, itms, ilabs, yplot_max, del2, pfont, fs, xsize, ysize)
             
    # new plot for adjusted labels 
    fig1, ax1 = plt.subplots(figsize=(xsize, ysize)) # first plot is MS2 spectrum
      
    ax1.set_ylim([yplot_min, yplot_max])
    ax1.set_xlim([0, pgv.MS2_mzhigh])
    ax1.set_xlabel("M/z", fontsize = axis_fs, font = pfont)
    ax1.set_ylabel(y_label, fontsize = axis_fs, font = pfont)
    ax1.set_title(title, font = pfont, fontsize = 12)
    
    stem_all = ax1.stem(mz2, it2, 'black', markerfmt=' ')
    stem_all[1].set_linewidth(unmatched_linewidth)
    stem_all[1].set_linestyles("dashed")
        
    for x,y,lab in xyz: # instead of stem plot, to allow different colors
        ax1.plot([x,x],[0,y],color_ion(lab), linewidth = matched_linewidth)
    
    # initial label positions are centered at top of each peak
    osets = [[o for o in tdict.keys()]] # initialize overlap sets to entire list of labels
    adjust_labels_new(tdict, osets, delx, dely, yplot_max) # two rounds for good measure
    osets, _ = find_label_overlaps(tdict,osets)
    adjust_labels_new(tdict, osets, delx, dely, yplot_max)
    osets, _ = find_label_overlaps(tdict,osets)
    
    # at this point, labels are separated along y.  
    for tk, t in tdict.items():
        ax1.text(t.tx, t.ty, t.text, color = color_ion(t.text), font = pfont, fontsize=fs, rotation = 90)
        xcon = [t.xc, t.xc, t.sx, t.sx]
        ycon = [t.ty - del0, t.ty - del0 - del1, t.sy + del0 + del1, t.sy + del0]
        ax1.plot(xcon,ycon, 'gray', linewidth = conn_linewidth )
    
    return fig1

def plot_ion_series_matrix(ukey):
    fig_x = 60
    fig_y = 20
    pfont = 'Arial'
    fs_box = 36  # font size for ion matrix plot
    title, plot_filename = ms2_plot_title(ukey)
    nr = len(pgv.unpacked_match_dict[ukey]["frag3"]) - 2
    fig2, ax2 = plt.subplots(figsize=(fig_x, fig_y)) # second plot is ion_series matrix
    ax2.set_xlim([-3, fig_x - 3])  # space to left and right for labels
    ax2.set_ylim([0, fig_y]) # space at top for title and labels
    ax2.axis('off')
    ax2.text(0,nr + 2, title, font = pfont, fontsize = fs_box)

    ion_series_matrix(ukey, ax2, fs_box)  

    return fig2

def ms2_plot(plot_folder, ukey):
    
    udict = pgv.unpacked_match_dict[ukey]
    min_ions = 5
    ion_dict = udict["matched_ions"]
    if len(list(ion_dict.keys())) < min_ions:
        return
    
    # fig1 = plot_ms2_spectrum(ukey)    
    # buf = io.BytesIO()
    # fig1.savefig(buf, format='png')
    # buf.seek(0)
    # image = Image.open(buf)
 
    fig1 = plot_ms2_spectrum(ukey)    
    fig2 = plot_ion_series_matrix(ukey)
    title, plot_filename = ms2_plot_title(ukey)
    
    buf = io.BytesIO()
    fig2.savefig(buf, format='png')   # for some reason this works with fig2 but not fig1???
    buf.seek(0)
    image = Image.open(buf)

    
    # width_inches, height_inches = fig2.get_size_inches()
    # dpi = fig2.get_dpi()

    # width_pixels = width_inches * dpi
    # height_pixels = height_inches * dpi

    
    # buf = io.BytesIO()
    # fig1.savefig(buf, format='png')
    # buf.seek(0)
    # image = Image.open(buf)
    # image_array = np.array(image, dtype = np.uint8)

    # image = np.array(buf.read(), dtype = np.uint8)
    # image = np.frombuffer(buf.read(), dtype=np.uint8).reshape(width_pixels, height_pixels)
    print("image type ", type(image))
    print("image size ", image.size)
    image_panel = PytheasImagePanel("title",image)

    # image_widget = ImageWidget()
    # image_widget.set_image(image_array) # Use image_qimage in the second case
   
    # pgv.main_window_layout.addWidget(image_panel)
    image_panel.show()

    
    with PdfPages(os.path.join(plot_folder, plot_filename)) as pdf: # multipage PDF
            plt.figure(fig1)
            pdf.savefig()
            plt.figure(fig2)
            pdf.savefig()

    plt.close(fig1)
    plt.close(fig2)

   
def bbox_overlap(bbox_a, bbox_b):
    a_bot_x, a_bot_y, a_width, a_height = [getattr(bbox_a,x) for x in ["x0", "y0", "width", "height"]]
    b_bot_x, b_bot_y, b_width, b_height = [getattr(bbox_b,x) for x in ["x0", "y0", "width", "height"]]
    
    a_top_x = a_bot_x + a_width
    a_top_y = a_bot_y + a_height
    b_top_x = b_bot_x + b_width
    b_top_y = b_bot_y + b_height

    if a_bot_x > b_top_x or b_bot_x > a_top_x:
       return False
   
    if a_bot_y > b_top_y or b_bot_y > a_top_y:
       return False
   
    return True

def format_ion_key(key,seq):    # not used?  keep for iTRAQ
    series, idx, z = key.split(":")
    if series in ["Rep", "BalTag"]:
        idx = ""
    if idx == "0" and series in "wxyz":
        idx = "tagC"
    if idx == "-1" and series in "wxyz":
        idx = "tag"
    fkey = r'$%s^{%s}_{%s}$' % (series, z, idx)
    return fkey

def color_ion(key): # key is full text of ion label
    color = 'black'
    for series, col in pgc.ion_color_dict.items():
        if series in key:
            return col
    
    return color

def adjust_labels_new(td, osets, delx, dely, ymax): # recursive adjustment of overlaps by y-increment
  drawn_keys = []
  for key in td.keys():
      for dkey in drawn_keys:
          while bbox_overlap(td[key],td[dkey]):
              # print(td[key].text, td[dkey].text)
              if td[key].ty + 5 * dely < ymax:
                  td[key].offset_text(0.0, dely) # adjust y-coord of label
              else:
                  td[key].offset_text(delx,0.0)
                  
      drawn_keys.append(key)
          
def find_label_overlaps(tdict, osets): # find overlapping bounding boxes for set of labels
    odict = {}
    for oset in osets:
        for tk1 in oset:
            olist = []
            for tk2 in oset:
                if tk1 == tk2:
                    continue
                if bbox_overlap(tdict[tk1],tdict[tk2]):
                    olist.append(tk2)
            odict[tk1] = olist

    # use graph to find sets of overlapping peaks.  key is index val is <ion_text_label>
    OG = nx.Graph() # overlap graph
    for o in odict.keys():
        OG.add_node(o)
  
    for o, ol in odict.items():
        for ool in ol:
            if (o,ool) not in OG.edges():
                OG.add_edge(o,ool)
    
    cc = nx.connected_components(OG)  # groups of overlapped keys
    cc_sets = [list(c) for c in cc]
    if len(cc_sets) == len(list(tdict.keys())): # termination condition for recursion
        done = True
    else:
        done =  False

    return cc_sets, done
        
def initialize_label_positions(mz2, it2, mzms, itms, ilabs, yplot_max, del0, pfont, fs, xsize, ysize):
    # dummy plot to get label positions...have to plot in order to extract them

    fig = plt.figure(num=1, clear=True,figsize=(xsize, ysize))
    ax = fig.add_subplot()
    fig.canvas.draw()
    ax.set_ylim([0, yplot_max])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    stem_all = ax.stem(mz2, it2, 'black', markerfmt=' ')
    stem_all[1].set_linewidth(0.25)
    stem_all[1].set_linestyles("dashed")

#   matched spectrum plot    
    stem_matched = ax.stem(mzms, itms, 'green', markerfmt=' ')  # stick plot
    stem_matched[1].set_linewidth(0.5)
    
    tdict = {}
    tidx = 0
    for x, y, lab in zip(mzms, itms, ilabs):
        t = plt.text(float(x), float(y), lab, font = pfont, fontsize =fs, rotation = 90)
        tdict[tidx] = ion_text_label(lab, t, x, y )
        tidx += 1
    
    fig.canvas.draw()
    r=fig.canvas.renderer
    for tk, t in tdict.items():
        t.get_bb(r,ax)  # get text bounding box.  Do this only once.  seems to be for horizontal text so hgt/wid are swapped?
        t.center_x()    # center text box over peak
        t.offset_text(0.0,del0) # shift label up

    draft_buf = io.BytesIO()  # object to hold image in bytes format
    plt.tight_layout()
    fig.savefig(draft_buf, format = 'png', dpi = 100)
    draft_buf.seek(0) # rewind to start (?)     # plt.savefig(wdir + 'test.png', dpi = 300) # save temp file
    plt.close()

    return tdict # return objects with label positions and width/height
