#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 14:00:34 2025

@author: jrwill
"""

import copy

from scipy.fft import rfft, irfft,rfftfreq
import numpy as np
import matplotlib.pyplot as plt


from isodist_global_variables import igv
from pytheas_global_vars import pgv


class ISODIST_MODEL:
    def __init__(self, species_dict):
        for key, val in species_dict.items():
            setattr(self, key, val)
            

# def read_control_data(infile):
    
#     with open(infile) as file:
#         lines = [line.split(" ") for line in file]
    
#     igv.opt = lines[0][0]       # program option = fitit, tryit  ( not used)
#     igv.batchfile = lines[1][0] # batch input file of peaks to be fit
#     igv.atomfile = lines[2][0]  # atom definition file  = atom_definitions.txt
#     igv.resfile = lines[3][0]   # residue definition file
#     igv.niter = int(lines[4][0])   # number of iterations for least squares (not used)
#     igv.sig_global = float(lines[5][0])  # global error ( not used)
    
#     igv.x_init = [float(lines[6][0]),float(lines[7][0]),float(lines[8][0])]
    
#     igv.auto_baseline = lines[6][1]
#     print("auto_baseline = ", igv.auto_baseline)
     
#     # return(opt,batchfile,atomfile,resfile,niter,sig_global,x_init,auto_baseline)


def create_atom_dict():
    
    with open(igv.atomfile) as file:
            lines = [line.strip().split() for line in file]

    nl = 0
    igv.atom_dict = {}
    while nl < len(lines):
        # nat = int(lines[nl][0])
        atom_symb = lines[nl][0]
        niso = int(lines[nl][1])
        # var = lines[nl][3]
        atom_comment = " ".join(lines[nl][2:])
        nl = nl + 1
        miso_list = []
        fiso_list = []
        for i in range(niso):
            miso_list.append(float(lines[nl][0]))
            fiso_list.append(float(lines[nl][1]))
            nl = nl + 1
        atom_entry = {"niso": niso, "miso": miso_list,"fiso": fiso_list, "comment": atom_comment}
        igv.atom_dict[atom_symb] = atom_entry
        
    # return(atom_dict)

def create_model():

    n_var_atom = 0
    with open(igv.model_file) as file:
            lines = [line.strip().split() for line in file]
    
    reslib_comment = " ".join(lines[0]).split("=")[0]
    print("MODEL: ",reslib_comment)
    
    # SPECIES DICT
    n_species = int(lines[1][0])
    species_lab = []
    amp = []
    ptr = 2
    igv.species_dict = {}
    for species in range(n_species):
        species_lab.append(lines[ptr][0])
        amp.append(float(lines[ptr][1]))
        ptr = ptr + 1
    natom_types = int(lines[ptr][0])
    ptr = ptr + 1
    
    
    igv.species_dict["n_species"] = n_species
    igv.species_dict["species_lab"] = species_lab
    igv.species_dict["amp"] = amp
    igv.species_dict["natom_types"] = natom_types
    
    igv.n_species = n_species
    
    # ATOM ORDER DICT
    igv.atom_order_dict = {}
    igv.atom_fix_dict = {}
    atom_idx = 0
    var_atom = []
    frac = []
    for atoms in range(natom_types):
        atom_symb = lines[ptr][0]
        atom_fix = lines[ptr][1]
        if atom_fix == "default":
            fix = 0
        elif atom_fix == "fixed":
            fix = 0
            frac.append(float(lines[ptr][2]))
        elif atom_fix == "variable":
            fix = 1
            frac.append(float(lines[ptr][2]))
            var_atom.append(atom_symb)
            n_var_atom = n_var_atom+1
        igv.atom_order_dict[atom_symb] = atom_idx
        igv.atom_fix_dict[atom_symb]=[atom_fix,fix,frac]
        atom_idx = atom_idx + 1
        ptr = ptr + 1
    igv.species_dict["frac"] = frac
    
    # print("constructing model")
    igv.model = ISODIST_MODEL(igv.species_dict)
    

    igv.model.n_var_atom = n_var_atom
    igv.model.var_atom = var_atom
    # omitting variable residue labeling for now
    # igv.model.n_var_res = n_var_res
    # igv.model.n_fix_res = n_fix_res
    # igv.model.var_res = var_res
    # igv.model.phi = phi
    # return(species_dict,atom_order_dict,atom_fix_dict,residue_dict)    

# def create_residue_dicts():

#     n_var_atom = 0
#     with open(igv.resfile) as file:
#             lines = [line.strip().split() for line in file]
    
#     reslib_comment = " ".join(lines[0]).split("=")[0]
#     print("RESIDUE LIBRARY: ",reslib_comment)
    
#     # SPECIES DICT
#     n_species = int(lines[1][0])
#     species_lab = []
#     amp = []
#     ptr = 2
#     igv.species_dict = {}
#     for species in range(n_species):
#         species_lab.append(lines[ptr][0])
#         amp.append(float(lines[ptr][1]))
#         ptr = ptr + 1
#     natom_types = int(lines[ptr][0])
#     ptr = ptr + 1
    
    
#     #TODO create model object here
    
#     igv.species_dict["n_species"] = n_species
#     igv.species_dict["species_lab"] = species_lab
#     igv.species_dict["amp"] = amp
#     igv.species_dict["natom_types"] = natom_types
    
#     igv.n_species = n_species
    
#     # ATOM ORDER DICT
#     igv.atom_order_dict = {}
#     igv.atom_fix_dict = {}
#     atom_idx = 0
#     var_atom = []
#     frac = []
#     for atoms in range(natom_types):
#         atom_symb = lines[ptr][0]
#         atom_fix = lines[ptr][1]
#         if atom_fix == "default":
#             fix = 0
#         elif atom_fix == "fixed":
#             fix = 0
#             frac.append(float(lines[ptr][2]))
#         elif atom_fix == "variable":
#             fix = 1
#             frac.append(float(lines[ptr][2]))
#             var_atom.append(atom_symb)
#             n_var_atom = n_var_atom+1
#         igv.atom_order_dict[atom_symb] = atom_idx
#         igv.atom_fix_dict[atom_symb]=[atom_fix,fix,frac]
#         atom_idx = atom_idx + 1
#         ptr = ptr + 1
#     igv.species_dict["frac"] = frac
    
#     # print("constructing model")
#     igv.model = ISODIST_MODEL(igv.species_dict)
    
#     # RESIDUE DICT
    
#     nr=0
#     n_var_res=0
#     n_fix_res=0
#     igv.residue_dict = {}
#     phi_fix = []
#     phi = []
#     phi_res = []
#     var_res = []
    
    
#     while ptr < len(lines):
#         resi = lines[ptr][0]
#         res_fix = lines[ptr][1]
#         resform = []
#         phi_res = 1.0
#         if res_fix == "default":
#             fixed_res = 0
#         elif res_fix == "fixed":
#             fixed_res = 1
#             n_fix_res = n_fix_res  + 1
#             phi_res = float(lines[ptr][2])
#             phi_fix.append(phi_res)
#         elif res_fix == "variable":
#             fixed_res = 2
#             n_var_res = n_var_res + 1
#             phi_res = float(lines[ptr][2])
#             phi.append(phi_res) 
#             var_res.append(resi)
#         ptr = ptr + 1
     
#         for species in range(n_species):
#             resform.append(list(map(int,lines[ptr][:natom_types])))
#             ptr = ptr + 1
#         entry = {"res_idx":nr, "res_fix":res_fix,"fixed_res":fixed_res,"phi_res":phi_res,
#                  "resform":resform,"phi_fix":phi_fix,"phi":phi,"var_res":var_res}
        
#         igv.residue_dict[resi] = entry
#         nr = nr + 1
#     # igv.residue_dict["n_var_atom"] = n_var_atom
#     # igv.residue_dict["var_atom"] = var_atom
#     # igv.residue_dict["n_var_res"] = n_var_res
#     # igv.residue_dict["n_fix_res"] = n_fix_res
#     # igv.residue_dict["var_res"] = var_res
#     # igv.residue_dict["phi"] = phi   # print("var_atom",var_atom)
    
#     igv.model.n_var_atom = n_var_atom
#     igv.model.var_atom = var_atom
#     igv.model.n_var_res = n_var_res
#     igv.model.n_fix_res = n_fix_res
#     igv.model.var_res = var_res
#     igv.model.phi = phi
#     # return(species_dict,atom_order_dict,atom_fix_dict,residue_dict)

def create_atom_mu():
# generate mu domain spectra for elements in order for residue definitions
    atom_mu = []   # build array as list
    igv.atom_keys = list(igv.atom_order_dict.keys())
    for key in igv.atom_keys:        
        niso,miso_list,fiso_list = [igv.atom_dict[key][k] for k in ["niso","miso","fiso"]]
        el_mz = [0] * igv.npt
        for i in range(niso):
            n=int(miso_list[i]*igv.scale_mz+0.5)
            el_mz[n] = fiso_list[i]
        el_mu = rfft(el_mz)
        atom_mu.append(el_mu)
    igv.atom_mu = np.array(atom_mu)  # convert to numpy array
    # return(atom_mu_array)

def initialize_residue_mu(active_residues):
    res_keys = [key for key in active_residues if len(key) == 1]
    res_keys = active_residues
    res_fixed_mu = []
    idx = 0
    igv.res_order_dict = {}
    for reskey in res_keys:
        igv.res_order_dict[reskey] = idx
        idx += 1
        # resform = igv.residue_dict[reskey]["resform"]
        resform = [pgv.nt_fragment_dict[label][reskey].mol_form for label in pgv.isotopic_species]
        species_mu = []
        species_var_mu = []
        for j in range(igv.n_species):
            res_mu = [1.0 + 0.0j] * igv.ncp  # initialize complex mu
            species_var_mu.append(res_mu)
            # for atomkey in igv.atom_keys:
            for atom, atidx in igv.atom_order_dict.items():
                # atidx = igv.atom_order_dict[atomkey]
                atom_fix,fix,frac = igv.atom_fix_dict[atom]
                print(reskey, j, atom)
                if fix == 0:
                    if atom not in resform[j]:
                        continue
                    n = resform[j][atom]
                    amu = list(map(lambda x:pow(x,n),igv.atom_mu[atidx]))
                    res_mu = np.multiply(res_mu,amu)
            species_mu.append(res_mu)
        res_fixed_mu.append(species_mu)
    igv.res_fixed_mu = np.array(res_fixed_mu)    # convert to numpy array
    igv.res_tot_mu = copy.deepcopy(igv.res_fixed_mu)     # copy for multiplying fixed and var parts       

def calculate_mol_form(frag3, label, z):
    mol_form = {}
    for res in frag3:
        res_form = pgv.nt_fragment_dict[label][res].mol_form
        for atom, stoi in res_form.items():
            if atom in mol_form:
                mol_form[atom] += stoi
            else:
                mol_form[atom] = stoi
        # print("mol_form: ", mol_form)
    mol_form["H1"] -= z
    
    return mol_form


def calc_isodist(x, return_opt, seq, yobs, mz_hd, z):
     
     #  fit function to take fit parameters x and return residuals eft_obs - sp_mz  
 
    # generate mu domain spectra for variable atoms
   
    xround = []
    for i in range(len(x)):  # for printing 
        xround.append(round(x[i],3))
 
    # calculate mu domain spectra for variable atoms
    for j in range(igv.model.n_var_atom):
        atom = igv.model.var_atom[j]   # N15
        miso_list = igv.atom_dict[atom]["miso"]
        vidx = igv.atom_order_dict[atom]
        fr = x[vidx]
        el_mz = [0] * igv.npt
        nu = int(miso_list[0]*igv.scale_mz+0.5)
        nf = int(miso_list[1]*igv.scale_mz+0.5)
        el_mz[nu] = 1 - fr
        el_mz[nf] = fr
        igv.atom_mu[vidx] = rfft(el_mz)
     
    active_res = set(list(seq))  # list of active residues to recalc variable part

    # calculate variable part of residue spectra 
    # indices are residue,species 
     
    for reskey in active_res:
        res_idx = igv.res_order_dict[reskey]
        
        resform = []
        for label in pgv.isotopic_species:
            rf = pgv.nt_fragment_dict[label][reskey].mol_form
            resform.append(rf)
            
        for j in range(igv.model.n_species):
            res_mu = np.ones(igv.ncp, dtype = complex)
            for k in range(igv.model.n_var_atom):
                atom = igv.model.var_atom[k]
                atidx = igv.atom_order_dict[atom]
                atom_fix,fix,frc = igv.atom_fix_dict[atom]
                if fix == 1:
                    if atom not in resform[j]:
                        continue
                    n = resform[j][atom]
                    amu = np.ones(igv.ncp,dtype = complex)
                    for nn in range(n):
                        amu = np.multiply(amu, igv.atom_mu[atidx])

                    res_mu = np.multiply(res_mu,amu)
 
                igv.res_tot_mu[res_idx][j] = np.multiply(igv.res_fixed_mu[res_idx][j],res_mu)         
    
    baseline = x[0]
    off = x[1]
    gw = x[2]
 
    cgw = np.sqrt(2)*(gw + 0.0j)
    gw2p=gw*np.sqrt(2.0*np.pi)
    cfac_const = 2.0*np.pi*((off+mz_hd)*igv.scale_mz)

    cpft = np.exp((0.0 + igv.cmz*cfac_const*1j))
    cgft = np.exp(-((igv.cmz +0j)/cgw)**2)/gw2p
    
    tot_mu = np.zeros(igv.ncp,dtype = complex)   #initialize total spectrum
    residues = list(seq)
    for j in range(igv.model.n_species):
        sp_mu = np.ones(igv.ncp, dtype = complex)  #initialize
        for reskey in residues:
            # res_idx = igv.residue_dict[reskey]["res_idx"]
            res_idx = igv.res_order_dict[reskey]
            sp_mu = np.multiply(sp_mu, igv.res_tot_mu[res_idx][j])  # accumulate product of residue spectra
        for k in range(z):
            sp_mu = np.multiply(sp_mu, np.conj(igv.atom_mu[igv.atom_order_dict['H1']])) # subtract z protons
        sp_mu = np.multiply(sp_mu,cgft)  # gaussian convolution
        sp_mu = np.multiply(sp_mu,cpft)  # heterodyne frequency shift
        tot_mu = tot_mu + x[igv.model.amp_idx[j]] * sp_mu   # accumulate species
    
        
# back transform and add baseline

    sp_mz = irfft(tot_mu) * igv.ncp/2 + baseline # normalization applied for consistency with isodist.f
   
    if return_opt == "residuals":
    
    # calculate residuals yobs - sp_mz using mz_ptr
        residuals = []
        mz_residuals = []

        for i in range(len(yobs)):
            residuals.append(yobs[i]-sp_mz[igv.mz_ptr[i]])
            mz_residuals.append(igv.mz[i])
        return(residuals)
    
    elif return_opt == "spectrum":
    
        spectrum =[]
        mz_spectrum = []
        for i in range(len(igv.mz_ptr)):
            spectrum.append(sp_mz[igv.mz_ptr[i]])
            mz_spectrum.append(igv.mz[igv.mz_ptr[i]])

        mz_spectrum = np.array(mz_spectrum) /z + mz_hd
        spectrum = np.array(spectrum)
        
        return mz_spectrum, spectrum
     

def calculate_residuals(x_fit, fit, resid):
#   RESIDUAL EXTRAVAGANZA
#   used outside of isodist for neural net analysis of goodness of fit

    # numpy ARRAYS FOR RESIDUAL CALCULATION:  resid, netfit, resid_m, resid_mm, netfit_m, netfit_mm
    # resid is residuals = "dy"
    
    resid = - np.array(resid)  # old resid = fit - dat
    fit = np.array(fit)
    
    # if igv.auto_baseline == "auto":
    #     netfit = fit - b
    # else:
    netfit = fit - x_fit[0]
              
    max_fit = float(np.amax(netfit))
    min_fit = float(np.amin(netfit))
    onepmax = 0.01*max_fit
        
    min_rsd = float(np.amin(resid))
    max_rsd = float(np.amax(resid))
    
# masked arrays 
    resid_m = np.where(netfit>= 10, resid, 0)  #masked at 10
    resid_mm = np.where(netfit >= onepmax, resid, 0)  #masked at 10
    netfit_m = np.where(netfit >= 10, netfit, 0)  #masked at 10
    netfit_mm = np.where(netfit >= onepmax,netfit,0)  #masked at 1% 
    
    nmz = len(resid)
    nmz_m = np.count_nonzero(netfit_m)
    nmz_mm = np.count_nonzero(netfit_mm)
  
# comments are from fortran version
      # 		
    # c Calculate the residuals
    # c These are mostly the old residuals by Mike Sykes.
    # c In this version:
    # c Removed duplicates of the large resids matrix with different normalizations
    # c Reports relevant fit sums for proper normalization by user 
    # c The only normalization here is by the number of summed data points
    # c Some signed meausures (not abs or **2) are also included -FA
    
 # intermediate numpy arrays
   
    abs_resid = np.abs(resid)
    sqrt_abs_resid = np.sqrt(abs_resid)
    abs_netfit = np.abs(netfit)
    sqrt_abs_netfit = np.sqrt(abs_netfit)
    sqrt_abs_resid = np.sqrt(abs_resid)
    weight = netfit/max_fit
    
    
    avg_fit = []   # list for fit quantities output
    # 	avg_fit1=avg_fit1+netfit			!avg fit
    # 	avg_fit1=avg_fit1/nmz			!avg fit
    fit_avg = np.sum(netfit)/nmz                            #fit1 - checked
    avg_fit.append(fit_avg)

    # 	avg_fit2=avg_fit2+netfit*netfit		!avg square fit
    # 	avg_fit2=avg_fit2/nmz			!avg square fit
    fit_sq = np.sum(netfit*netfit)/nmz                      #fit2 - checked
    avg_fit.append(fit_sq)

    # 	avg_fit3=avg_fit3+sqrt(abs(netfit))	!avg root fit
    # 	avg_fit3=avg_fit3/nmz			!avg root fit
    fit_rt = np.sum(sqrt_abs_netfit)/nmz                    #fit3 - checked
    avg_fit.append(fit_rt)

    # 	avg_fit4=avg_fit4+fit_m				!avg masked fit, mask >10
    # 	avg_fit4=avg_fit4/nmz_m			!avg masked fit, mask >10
    fit_m = np.sum(netfit_m)/nmz_m                          #fit4 - checked
    avg_fit.append(fit_m)

    # 	avg_fit5=avg_fit5+fit_sq_m			!avg masked square fit, mask >10
    # 	avg_fit5=avg_fit5/nmz_m			!avg masked square fit, mask >10
    fit_sq_m = np.sum(netfit_m*netfit_m)/nmz_m              #fit5 - checked
    avg_fit.append(fit_sq_m)

    # 	avg_fit6=avg_fit6+fit_rt_m			!avg masked root fit, mask >10
    # 	avg_fit6=avg_fit6/nmz_m			!avg masked root fit, mask >10
    fit_rt_m = np.sum(np.sqrt(np.abs(netfit_m)))/nmz_m      #fit6 - checked
    avg_fit.append(fit_rt_m)

    # 	avg_fit7=avg_fit7+fit_mm			!avg masked fit, mask > 0.01*max_fit
    # 	avg_fit7=avg_fit7/nmz_mm		!avg masked fit, mask > 0.01*max_fit
    fit_mm = np.sum(netfit_mm)/nmz_mm                       #fit7 - checked
    avg_fit.append(fit_mm)

    # 	avg_fit8=avg_fit8+fit_sq_mm			!avg masked square fit, mask > 0.01*max_fit
    # 	avg_fit8=avg_fit8/nmz_mm		!avg masked square fit, mask > 0.01*max_fit
    fit_sq_mm = np.sum(netfit_mm*netfit_mm)/nmz_mm          #fit8 - checked
    avg_fit.append(fit_sq_mm)

    # 	avg_fit9=avg_fit9+fit_rt_mm			!avg masked root fit, mask > 0.01*max_fit
    # 	avg_fit9=avg_fit9/nmz_mm		!avg masked root fit, mask > 0.01*max_fit
    fit_rt_mm = np.sum(np.sqrt(np.abs(netfit_mm)))/nmz_mm   #fit9 - checked
    avg_fit.append(fit_rt_mm)

    # 	avg_fit10=avg_fit10+netfit*netfit*netfit !avg cube fit
    # 	avg_fit10=avg_fit10/nmz			!avg cube fit
    fit_cu = np.sum(netfit*netfit*netfit)/nmz               #fit10 - checked
    avg_fit.append(fit_cu)

    # 	avg_fit11=avg_fit11+abs(netfit)		!avg abs fit	
    # 	avg_fit11=avg_fit11/nmz			!avg abs fit
    fit_abs = np.sum(abs_netfit)/nmz                        #fit11  - checked
    avg_fit.append(fit_abs)

    avg_wr = []   # list for weighted residual output

    # C	(1)
    # 	rsd(i)=dy !Keep this one in an array for output
    # 	avg_wr1=avg_wr1+abs(dy)			!avg abs resids
    # 	avg_wr1=avg_wr1/nmz				!avg abs resids
    avgabsres = np.sum(abs_resid)/nmz                                #1  - checked
    avg_wr.append(avgabsres)

    # C	(2)
    # 	rsd_m=dy
    # 	avg_wr2=avg_wr2+abs(rsd_m)		!avg masked abs resids, mask >10
    # 	avg_wr2=avg_wr2/nmz_m			!avg masked abs resids, mask >10
    rsd_m = np.sum(np.abs(resid_m))/nmz_m                       #2 - checked
    avg_wr.append(rsd_m)

     # C	 (3)
     # 	rsd_w=dy*abs(netfit)
     # 	avg_wr3=avg_wr3+abs(rsd_w)		!avg weighted abs resids, weight = abs(netfit)
     # 	avg_wr3=avg_wr3/nmz				!avg weighted abs resids, weight = abs(netfit)
    rsd_w = np.sum(np.abs(resid*abs_netfit))/nmz                    #3 - checked
    avg_wr.append(rsd_w)

    # C	(4)
    # 	rsd_sq=dy*dy
    # 	avg_wr4=avg_wr4+rsd_sq			!avg square resids
    # 	avg_wr4=avg_wr4/nmz				!avg square resids
    rsd_sq = np.sum(resid * resid)/nmz                       #4  - checked
    avg_wr.append(rsd_sq)

    # C	(5)
    # 		rsd_sq_m=dy*dy    
    # 	avg_wr5=avg_wr5+rsd_sq_m		!avg masked square resids, mask: >10
    # 	avg_wr5=avg_wr5/nmz_m			!avg masked square resids, mask: >10
    rsd_sq_m = np.sum(resid_m*resid_m)/nmz_m                    #5 - checked
    avg_wr.append(rsd_sq_m)

    # C	(6)
    # 	rsd_sq_w=dy*dy*abs(netfit)
    # 	avg_wr6=avg_wr6+rsd_sq_w		!avg weighted square resids, weight = abs(netfit)
    # 	avg_wr6=avg_wr6/nmz				!avg weighted square resids, weight = abs(netfit)
    rsd_sq_w = np.sum(resid*resid*abs_netfit)/nmz           #6 - checked
    avg_wr.append(rsd_sq_w)

    # C	(7)
    # 	rsd_rt=sqrt(abs(dy))
    # 	avg_wr7=avg_wr7+rsd_rt			!avg root resids
    # 	avg_wr7=avg_wr7/nmz				!avg root resids
    rsd_rt = np.sum(sqrt_abs_resid)/nmz                      #7 - checked
    avg_wr.append(rsd_rt)

    # C	(8)
    # 	rsd_rt_m=sqrt(abs(dy))
    # 	avg_wr8=avg_wr8+rsd_rt_m		!avg masked root resids, mask >10
    # 	avg_wr8=avg_wr8/nmz_m			!avg masked root resids, mask >10
    rsd_rt_m = np.sum(np.sqrt(np.abs(resid_m)))/nmz_m         #8 - checked
    avg_wr.append(rsd_rt_m)

    # C	(9)
    # 	rsd_rt_w=sqrt(abs(dy))*abs(netfit)
    # 	avg_wr9=avg_wr9+rsd_rt_w		!avg weighted root resids, weight = abs(netfit)
    # 	avg_wr9=avg_wr9/nmz				!avg weighted root resids
    rsd_rt_w = np.sum(sqrt_abs_resid*abs_netfit)/nmz       #9 - checked
    avg_wr.append(rsd_rt_w)
  
    # C	(10)
    # 	rsd_g=dy*sqrt(abs(netfit))
    # 	avg_wr10=avg_wr10+abs(rsd_g)	!avg gently weighted abs resids, weight = sqrt(abs(netfit))
    # 	avg_wr10=avg_wr10/nmz			!avg gently weighted abs resids, weight = sqrt(abs(netfit))
    rsd_g = np.sum(np.abs(resid*sqrt_abs_netfit))/nmz               #10 - checked
    avg_wr.append(rsd_g)

    # C	(11)
    # 	rsd_sq_g=dy*dy*sqrt(abs(netfit))
    # 	avg_wr11=avg_wr11+rsd_sq_g		!avg gently wieghted square resids, weight = sqrt(abs(netfit))
    # 	avg_wr11=avg_wr11/nmz			!avg gently wieghted square resids, weight = sqrt(abs(netfit))
    rsd_sq_g = np.sum(resid*resid*sqrt_abs_netfit)/nmz      #11 - checked
    avg_wr.append(rsd_sq_g)

    # C	(12)
    # 	rsd_rt_g=sqrt(abs(dy))*sqrt(abs(netfit))
    # 	avg_wr12=avg_wr12+rsd_rt_g		!avg gently weighted root resids, weight = sqrt(abs(netfit))
    # 	avg_wr12=avg_wr12/nmz			!avg gently weighted root resids, weight = sqrt(abs(netfit))
    rsd_rt_g = np.sum(sqrt_abs_resid*sqrt_abs_netfit)/nmz  #12 - checked
    avg_wr.append(rsd_rt_g)
  
    # C	(13)
    # 	rsd_mm=dy    
    # 	avg_wr13=avg_wr13+abs(rsd_mm)	!avg masked abs resids, mask > 0.01*max_fit
    # 	avg_wr13=avg_wr13/nmz_mm		!avg masked abs resids, mask > 0.01*max_fit
    rsd_mm = np.sum(np.abs(resid_mm))/nmz_mm                            #13 - checked
    avg_wr.append(rsd_mm)

    # C	(14)
    # 	rsd_sq_mm=dy*dy
    # 	avg_wr14=avg_wr14+rsd_sq_mm		!avg masked square resids, mask > 0.01*max_fit
    # 	avg_wr14=avg_wr14/nmz_mm		!avg masked square resids, mask > 0.01*max_fit
    rsd_sq_mm = np.sum(resid_mm*resid_mm)/nmz_mm                        #14 - checked
    avg_wr.append(rsd_sq_mm)

    # C	(15)
    # 	rsd_rt_mm=sqrt(abs(dy))
    # 	avg_wr15=avg_wr15+rsd_rt_mm		!avg masked root resids, mask > 0.01*max_fit
    # 	avg_wr15=avg_wr15/nmz_mm		!avg masked root resids, mask > 0.01*max_fit
    rsd_rt_mm = np.sum(np.sqrt(np.abs(resid_mm)))/nmz_mm                #15 - checked
    avg_wr.append(rsd_rt_mm)
   
    # C	(16)
    # 	rsd_ww=dy*abs(weight)   
    # 	avg_wr16=avg_wr16+abs(rsd_ww)	!avg weighted abs resids, weight = abs(netfit/max_fit)
    # 	avg_wr16=avg_wr16/nmz			!avg weighted abs resids, weight = abs(netfit/max_fit)
    rsd_ww = np.sum(np.abs(resid * np.abs(weight)))/nmz                 #16 - checked
    avg_wr.append(rsd_ww)

    # C	(17)
    # 	rsd_sq_ww=dy*dy*abs(weight)
     # 	avg_wr17=avg_wr17+rsd_sq_ww		!avg weighted square resids, weight = abs(netfit/max_fit)
    # 	avg_wr17=avg_wr17/nmz			!avg weighted square resids, weight = abs(netfit/max_fit)
    rsd_sq_ww = np.sum(resid*resid*np.abs(weight))/nmz                  #17  - checked
    avg_wr.append(rsd_sq_ww)

    # C	(18)
    # 	rsd_rt_ww=sqrt(abs(dy))*abs(weight)
    # 	avg_wr18=avg_wr18+rsd_rt_ww		!avg weighted root resids, weight = abs(netfit/max_fit)
    # 	avg_wr18=avg_wr18/nmz			!avg weighted root resids, weight = abs(netfit/max_fit)
    rsd_rt_ww = np.sum(np.sqrt(np.abs(resid))*np.abs(weight))/nmz       #18 - checked
    avg_wr.append(rsd_rt_ww)

    # C	(19)
    # 	rsd_mmww=dy*abs(weight)    
    # 	avg_wr19=avg_wr19+abs(rsd_mmww)	!avg masked weighted abs resids, mask > 0.01*max_fit, weight = abs(netfit/max_fit)
    # 	avg_wr19=avg_wr19/nmz_mm		!avg masked weighted abs resids, mask > 0.01*max_fit, weight = abs(netfit/max_fit)
    rsd_mmww = np.sum(np.abs(resid*np.abs(weight)))/nmz_mm              #19 - checked
    avg_wr.append(rsd_mmww)

    # C	(20)
    # 	rsd_sq_mmww=dy*dy*abs(weight)
    # 	avg_wr20=avg_wr20+rsd_sq_mmww	!avg masked weighted square resids, mask > 0.01*max_fit, weight = abs(netfit/max_fit)
    # 	avg_wr20=avg_wr20/nmz_mm		!avg masked weighted square resids, mask > 0.01*max_fit, weight = abs(netfit/max_fit)
    rsd_sq_mmww = np.sum(resid*resid*np.abs(weight))/nmz_mm             #20 - checked
    avg_wr.append(rsd_sq_mmww)

    # C	(21)
    # 	rsd_rt_mmww=sqrt(abs(dy))*abs(weight)
    # 	avg_wr21=avg_wr21+rsd_rt_mmww	!avg masked weighted root resids, mask > 0.01*max_fit, weight = abs(netfit/max_fit)
    # 	avg_wr21=avg_wr21/nmz_mm		!avg masked weighted root resids, mask > 0.01*max_fit, weight = abs(netfit/max_fit)
    rsd_rt_mmww = np.sum(np.sqrt(np.abs(resid))*np.abs(weight))/nmz_mm  #21 - checked
    avg_wr.append(rsd_rt_mmww)

     # C (22) 
     # 	avg_wr22=avg_wr22+dy			!avg resids
     # 	avg_wr22=avg_wr22/nmz			!avg resids
    avgres = np.sum(resid)/nmz                                 #22 - checked
    avg_wr.append(avgres)

    # C (23)  
    # 	avg_wr23=avg_wr23+dy*dy*dy		!avg cube resids
    # 	avg_wr23=avg_wr23/nmz			!avg cube resids
    rsd_cu = np.sum(resid * resid * resid)/nmz               #23 - checked
    avg_wr.append(rsd_cu)

    return max_fit, min_fit, max_rsd, min_rsd, avg_fit, avg_wr


def plot_obs_calc(obs_x, obs_y, calc_x, calc_y, title):
        fix, ax = plt.subplots()
        ax.plot(obs_x, obs_y, ".k")
        plt.ylabel("Amplitude")
        plt.xlabel("m/z")
        plt.title (title)
        plt.plot(calc_x, calc_y,"-b")
        plt.show()
