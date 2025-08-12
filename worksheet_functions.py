#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 13:40:23 2024

@author: jrwill
"""

import os
import pandas as pd
import xlsxwriter

from pytheas_global_vars import pgv


def format_worksheet_columns(worksheet, df, fd):
    max_col_length = 100
    for col in df:   # set column widths
        # print("worksheet column ", col)
        if "matched_" in col:
            column_length = max(df[col].astype(str).map(len).max(), len(col))
            column_length = min(column_length,max_col_length)        
            col_idx = df.columns.get_loc(col)
            worksheet.set_column(col_idx, col_idx, column_length)
            
        if col in fd:
            if "matched_" in col or fd[col]["col_width"] == "var":
                column_length = max(df[col].astype(str).map(len).max(), len(col))
            else: 
                column_length = fd[col]["col_width"]
            column_length = min(column_length,max_col_length)        
            col_idx = df.columns.get_loc(col)
            worksheet.set_column(col_idx, col_idx, column_length)

def write_worksheet_rows(worksheet, df, fd, bold):
    for col in df.columns:  # output header
        if col in fd:
            worksheet.write(0, df.columns.get_loc(col), col, bold)
    rowidx = 1
    for row in df.iterrows():
        colidx = 0
        for col in row[1].keys():
            if col in fd:
                if "CID" not in col:
                    # print(rowidx, colidx,row[1][col], type(row[1][col]), fd[col]["xformat"] )
                    cell = row[1][col]
                    if type(cell) == list:
                        cell = " ".join(cell)
                    worksheet.write(rowidx, colidx, cell, fd[col]["xformat"])
                else:
                    worksheet.write(rowidx, colidx, row[1][col])
                colidx +=1
            rowidx +=1
            
def write_parameter_worksheet(par_df, ppar_file):          

    print("write_parameter_worksheet", ppar_file)
    workbook = xlsxwriter.Workbook(ppar_file,{"nan_inf_to_errors": True})
    # worksheet = workbook.add_worksheet(pgv.job_dir.split("/")[-1])
    worksheet = workbook.add_worksheet(ppar_file.split("/")[-1])

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

