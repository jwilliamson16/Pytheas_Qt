#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 13:50:41 2024

@author: jrwill
"""

from PyQt5.QtWidgets import (QWidget, QLabel, QPushButton, QGridLayout, 
                              QLineEdit, QRadioButton, QButtonGroup, QFrame, 
                              QHBoxLayout, QCheckBox, QFileDialog, QVBoxLayout)
from PyQt5.QtCore import pyqtSignal
import os

from pytheas_global_vars import pgv, pgc, pgvdict, PGVStrToType


# PytheasWidgets 

class PytheasRadioButtonBar(QWidget):
    def __init__(self, pgvd):
        super().__init__()
        for key,val in pgvd.items(): # set up widget attributes
            setattr(self, key, val)
            # print(key, val)

        self.layout = QGridLayout()
        self.label = QLabel(self.widget_text)
        self.label.setFixedWidth(200)
        self.layout.addWidget(self.label, 0, 0)
        self.rb_list = []
        ctr = 1
        self.group = QButtonGroup()  
        for b in self.option_list.split(","): # always str
            rb = QRadioButton(b)
            self.rb_list.append(rb)
            self.layout.addWidget(rb, 0, ctr)
            ctr += 1
            rb.clicked.connect(self.set_value)
            self.group.addButton(rb)
        n_spacers = 12 - len(self.rb_list) # crude hack to get alignment in panels
        for i in range(n_spacers):
            self.layout.addWidget(QWidget(), 0, ctr)            
            ctr += 1

        self.setLayout(self.layout)
        self.layout.setContentsMargins(0, 0, 0, 0)
        
        self.load_value(self.default_value)
 
    def set_value(self): # monitor signal from widget
        # print("RBB set_value", self.pgv, self.sender().text())
        # print(self.sender.__dict__.keys())
        # print(self.sender.__dir__)
        self.load_value(self.sender().text())
        
    def load_value(self, value):
        # print("RBB load_value before", self.pgv, self.value)
        for rb in self.rb_list: 
            if rb.text() == str(value):
                rb.setChecked(True)
            else:
                rb.setChecked(False)
                
        self.value = value
        # print("RBB load_vaue after", self.pgv, self.value)

        setattr(pgv, self.pgv, PGVStrToType(self.data_type, self.value))

class PytheasCheckBoxBar(QWidget):
    def __init__(self, pgvd):
        super().__init__()
        for key,val in pgvd.items():  # set up widget attributes
            setattr(self,key,val)
 
        self.layout = QGridLayout()
        self.label = QLabel(self.widget_text)
        self.label.setFixedWidth(200)
        self.layout.addWidget(self.label, 0, 0)

        self.cb_list = []
        ctr = 1
        for b in self.option_list.split(","): # always str
            cb = QCheckBox(b)
            self.cb_list.append(cb)
            self.layout.addWidget(cb, 0, ctr)
            # cb.stateChanged.connect(self.set_value)
            cb.stateChanged.connect(self.set_value)
            ctr += 1
            
        n_spacers = 12 - len(self.cb_list)
        for i in range(n_spacers):
            self.layout.addWidget(QWidget(), 0, ctr)            
            ctr += 1

        self.setLayout(self.layout)
        self.layout.setContentsMargins(0, 0, 0, 0)
        
        self.load_value(self.default_value)
        # print("PYTHEAS_CHECKBOXBAR: init", self.default_value)

    def set_value(self):  # monitor signal from widget
        # print("start PYTHEAS_CHECKBOXBAR: set_value", self.pgv)
        # self.load_value([cb.text() for cb in self.cb_list if cb.isChecked()])
        value = [cb.text() for cb in self.cb_list if cb.isChecked()]
        value_list = PGVStrToType(self.data_type, value)
        setattr(pgv, self.pgv, value_list)
        # self.load_value([cb.text() for cb in self.cb_list if cb.isChecked()])
        # print(" done with PYTHEAS_CHECKBOXBAR: set_value", [cb.text() for cb in self.cb_list if cb.isChecked()])
 
    def load_value(self, value):
        # print("start PYTHEAS_CHECKBOXBAR: load_value", self.pgv, self.value)

        value_list = PGVStrToType(self.data_type, value)
        for cb in self.cb_list:
            if cb.text() in value_list:
                cb.setChecked(True)
            else:
                cb.setChecked(False)
                
        self.value = value_list

        setattr(pgv, self.pgv, self.value)
        # print("done with PYTHEAS_CHECKBOXBAR: load_value", self.pgv, self.value)


class PytheasEntry(QWidget):
    def __init__(self, pgvd):
        super().__init__()
        for key,val in pgvd.items():  # set up widget attributes
            setattr(self,key,val)
        self.layout = QHBoxLayout()
        self.label = QLabel(self.widget_text)
        self.label.setFixedWidth(200)
        self.layout.addWidget(self.label)
        self.entry = QLineEdit()
        self.layout.addWidget(self.entry)
        self.setLayout(self.layout)
        self.entry.textEdited.connect(self.set_value)
        self.entry.textChanged.connect(self.set_value)
        self.entry.editingFinished.connect(self.set_value)
        self.layout.setContentsMargins(0, 0, 0, 0)
        
        self.load_value(self.default_value)

    def set_value(self): # monitor signal from widget
        self.load_value(self.sender().text())
        
    def load_value(self, value):
        self.entry.setText(str(value))
        
        self.value = PGVStrToType(self.data_type, value)
        setattr(pgv, self.pgv, self.value)

 
class PytheasDirectory(QWidget):
    def __init__(self, pgvd):
        super().__init__()
        for key,val in pgvd.items():  # set up widget attributes
            setattr(self,key,val)

        self.dir_layout = QHBoxLayout()
        self.label = QLabel(self.widget_text)
        self.label.setFixedWidth(120)
        self.button = QPushButton("Choose Directory")
        self.button.clicked.connect(self.choose_directory)
        self.button.setFixedSize(120,20)
        self.dir_label = QLabel(self.default_value)
        self.dir_label.setFixedWidth(600)
        self.dir_layout.addWidget(self.label)
        self.dir_layout.addWidget(self.button)
        self.dir_layout.addWidget(self.dir_label)
        self.dir_layout.setContentsMargins(0, 0, 0, 0)

        self.setLayout(self.dir_layout)
        
        self.load_value(self.default_value)

    def choose_directory(self):
        self.directory = QFileDialog.getExistingDirectory(self, 'Choose Directory', '/Users/jrwill/prog/pytheas_tk_interface')

        if self.directory:
            self.load_value(self.directory)
            self.dir_label.setText(self.directory)
        
    def load_value(self, value):
        self.value = value
        self.dir_label.setText(value)
        setattr(pgv, self.pgv, self.value)

class PytheasFile(QWidget):
    def __init__(self, pgvd):
        super().__init__()
        for key,val in pgvd.items():  # set up widget attributes
            setattr(self,key,val)

        self.layout = QHBoxLayout()        
        self.label = QLabel(self.widget_text)
        self.label.setFixedWidth(200)
        self.file_widget = QLineEdit()
        self.file_widget.setFixedWidth(450)
        self.file_widget.textEdited.connect(self.set_value)
        self.file_widget.textChanged.connect(self.set_value)
        self.file_widget.editingFinished.connect(self.set_value)

        self.file_button = QPushButton("Browse")
        self.file_button.setFixedWidth(75)

        self.file_button.clicked.connect(self.show_file_dialog)
        self.layout.addWidget(self.label)
        self.layout.addWidget(self.file_widget)
        self.layout.addWidget(self.file_button)
        self.setLayout(self.layout)
        self.layout.setContentsMargins(0, 0, 0, 0)
        filter_list = self.option_list.split(",")
        self.filter = "Choose File: (" + " ".join(["*." + f for f in filter_list]) + ")"
        
        self.load_value(self.default_value)

    def show_file_dialog(self):
        file_dialog = QFileDialog(self)
        file_dialog.setWindowTitle("Select File")
        self.file_path, self.file = QFileDialog.getOpenFileName(self, "Select Excel File", "", self.filter)
        if self.file_path ==  "": # in case of cancel
            return
        self.load_value(self.file_path)

    def set_value(self): # monitor signal from widget
        self.load_value(self.sender().text())
             
    def load_value(self, file_path):
        file_dir = os.path.dirname(file_path)
        # trim file_path if file is in work/user/pytheas dirs
        if file_dir in [pgv.working_dir, pgv.user_dir, pgv.pytheas_dir]:
             value = os.path.basename(file_path)
        else:
            value = file_path
        
        self.value = value
        self.file_widget.setText(value)
        setattr(pgv, self.pgv, self.value)

class PytheasFileOutput(QWidget):    # identical to PytheasEntry
    def __init__(self, pgvd):
        super().__init__()
        for key,val in pgvd.items():  # set up widget attributes
            setattr(self,key,val)
        self.layout = QHBoxLayout()
        self.label = QLabel(self.widget_text)
        self.label.setFixedWidth(200)
        self.layout.addWidget(self.label)
        self.entry = QLineEdit()
        self.entry.setFixedWidth(450)

        self.layout.addWidget(self.entry)
        self.entry.setText(str(self.default_value))
        self.setLayout(self.layout)
        self.entry.textEdited.connect(self.set_value)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.load_value(self.default_value)

    def set_value(self): # monitor signal from widget
        self.load_value(self.sender().text())
 
    def load_value(self, value):
        # print("PFO load_value", value)
        self.value = value
        self.entry.setText(str(value))
        setattr(pgv, self.pgv, self.value)

class PytheasLabel(QWidget):
    def __init__(self,pgvd):
        super().__init__()
        for key,val in pgvd.items():  # set up widget attributes
            setattr(self,key,val)
        self.label_layout = QHBoxLayout()
        self.label = QLabel(self.widget_text)
        self.label.setStyleSheet("font-weight: bold")
        self.label_layout.addWidget(self.label)
        self.setLayout(self.label_layout)
        self.label_layout.setContentsMargins(0, 0, 0, 0)
        
    def load_value(self,value):
        pass
    
class PytheasButton(QWidget):
    def __init__(self, label, command):
        super().__init__()
        self.label = label
        self.command = command           
        self.widget = QWidget()
        self.layout = QHBoxLayout()
        self.button = QPushButton(self.label)
        self.button.clicked.connect(self.command)
        self.layout.addWidget(self.button)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(self.layout)
        
class PytheasPanel(QWidget): # Floating panel for options
    closed = pyqtSignal(str) # pass widget group with signal to PytheasOptionPanel

    def __init__(self, group, pgv_list):
        super().__init__()
 
        self.group = group
        self.pgv_list = pgv_list
        self.panel_layout = QVBoxLayout()   # overall layout for widget set
        
        self.frame = QFrame()  # frame widtet to hold widget set
        self.frame_layout = QVBoxLayout(self.frame)
        self.frame.setFrameShape(QFrame.Panel)

        self.widget_list = []
        for p in pgv_list:
            # print("PytheasPanel", p)
            pwidget = pgvdict[p]["widget_type"]
            # print("widget_type", pwidget)
            if pwidget == "PytheasLabel":
                self.setWindowTitle(pgvdict[p]["widget_text"])
                continue
            # print(pgvdict[p])
            w = globals()[pwidget](pgvd=pgvdict[p])
            self.frame_layout.addWidget(w)
            self.widget_list.append(w)

        self.panel_layout.addWidget(self.frame)
        self.setLayout(self.panel_layout)
        self.rect = self.frame.geometry()
        
    def show_panel(self):
        self.show()
        self.shown =  True
        
    def hide_panel(self):
        self.hide()
        self.shown = False
        
    def closeEvent(self, event): # this is called when widget is closed from upper left
        self.closed.emit(self.group) # send signal to PytheasOptionPanel
        super().closeEvent(event) # signal monitored by PytheasOptionPanel

class PytheasFixedPanel(QWidget):

    def __init__(self, group, pgv_list):
        super().__init__()
        self.group = group
        self.pgv_list = pgv_list
        self.panel_layout = QVBoxLayout()   # overall layout for widget set
        
        self.frame = QFrame()  # frame widtet to hold widget set
        self.frame_layout = QVBoxLayout(self.frame)
        self.frame.setFrameShape(QFrame.Panel)

        self.widget_list = []
        for p in pgv_list:
            pwidget = pgvdict[p]["widget_type"]
            w = globals()[pwidget](pgvd=pgvdict[p])
            self.frame_layout.addWidget(w)
            self.widget_list.append(w)


        self.panel_layout.addWidget(self.frame)
        self.setLayout(self.panel_layout)
        self.rect = self.frame.geometry()
        

class PytheasOptionPanel(QWidget): # CheckBox widget to control display of option panels
    def __init__(self, option_group_dict, option_panel_dict, floating_windows):
        super().__init__()
        self.option_group_dict = option_group_dict
        self.option_panel_dict = option_panel_dict
        self.floating_windows = floating_windows
        self.panel_layout = QVBoxLayout()   # overall layout for widget set
        
        self.frame = QFrame()  # frame widget to hold widget set
        self.frame_layout = QGridLayout(self.frame)
        self.frame.setFrameShape(QFrame.Panel)

        self.label = QLabel("Pytheas Module Options")
        self.label.setStyleSheet("font-weight: bold")
        self.label.setFixedWidth(200)
        self.frame_layout.addWidget(self.label, 0, 0)

        self.group = QButtonGroup()  
        self.cb_dict = {}# one big grid of CheckBoxes
        row = 1
        for group, group_list in self.option_group_dict.items():
            col = 0
            label = QLabel(group)
            label.setFixedWidth(200)
            self.frame_layout.addWidget(label, row, col)

            for b in group_list: # always str
                col += 1
                cb = QCheckBox(b)
                self.cb_dict[b] = cb
                self.frame_layout.addWidget(cb, row, col)
                cb.stateChanged.connect(self.update_panels) # this handles closing panel from checkbox
            row += 1

        self.panel_layout.addWidget(self.frame)
        self.setLayout(self.panel_layout)
        
        for cb, w in option_panel_dict.items(): # connect to closeEvent defined in Panel
            # print("Option_Panel: cb", cb)
            w.closed.connect(self.close_panel) # this handles closing panel from panel upper corner

    def update_panels(self):  # monitor signal from widgets, and redraw floating panels
        panel_y = pgc.panel_y_pos # starting pos for stacked floating widgets
        for cb, w in self.cb_dict.items():
            if w.isChecked():
                self.option_panel_dict[cb].show_panel()
                self.option_panel_dict[cb].move(pgc.panel_x_pos, panel_y)
                # get size of floating panel widget AFTER it is shown
                panel_y += self.option_panel_dict[cb].size().height() + pgc.panel_y_delta
            else:
                 self.option_panel_dict[cb].hide_panel()
               
    def close_panel(self, group): #upon closeEvent signal, change cb state to trigger update and close
         self.cb_dict[group].setChecked(False)


class PytheasButtonPanel(QWidget):
    def __init__(self, label, button_dict):
        super().__init__()
        self.panel_layout = QVBoxLayout()   # overall layout for widget set
        
        self.frame = QFrame()  # frame widtet to hold widget set
        self.frame_layout = QHBoxLayout(self.frame)
        self.frame.setFrameShape(QFrame.Panel)
        self.label = QLabel(label)
        self.label.setFixedWidth(200)
        self.label.setStyleSheet("font-weight: bold")
        self.frame_layout.addWidget(self.label)

        self.widget_list = []
        for button, command in button_dict.items():
            b = PytheasButton(button, command)
            self.frame_layout.addWidget(b)
            self.widget_list.append(b)

        self.panel_layout.addWidget(self.frame)
        self.setLayout(self.panel_layout)

