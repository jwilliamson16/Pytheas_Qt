#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 13:59:54 2025

@author: jrwill
"""

#isodist_global_variables

class IGV:
    def __init__(self):
        self.npt = 131072	     # np is the number of points calculated
        self.ncp=int(self.npt/2)+1		     # ncp is the number of complex points				
        self.scale_mz=1000.0	     # scale_mz determines the number of points calculated per dalto for oversampling				
        self.outlabels=["file","mol_id","seq","mw","z_charge","chisq","mz"]  # headers for csv
        self.sig_global = 100.0
        self.rtpad = 60.0   # RT window for peak fitting
        self.mzpad = 2.0
        
        
        
        
igv = IGV()
