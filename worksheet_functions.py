#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 13:40:23 2024

@author: jrwill
"""

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
