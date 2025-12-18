#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 14:00:34 2025

@author: jrwill
"""

import copy
from datetime import datetime
import os

from scipy.fft import rfft, irfft,rfftfreq
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares


# from isodist_global_variables import pgv
from pytheas_global_vars import pgv, pgc
from pytheas_IO import read_pytheas_file, load_pickle_set

from mod_seq_functions import generate_mod_seq, generate_mod_seq_ends
from minispectrum_functions import plot_minispectrum
from PyQt5.QtWidgets import  QFileDialog


class ISODIST_MODEL:
    def __init__(self, species_dict):
        for key, val in species_dict.items():
            setattr(self, key, val)
            

def find_active_residues():
    pgv.active_residues = []
    for idx, row in pgv.top_match_df.iterrows():
        frag3 = row["frag3"]
        for res in frag3:
            if res not in pgv.active_residues:
                pgv.active_residues.append(res)
    
    print("active residues: ", pgv.active_residues)


#TODO integrate this with pgv.atomic_dict
def create_atom_dict():
    
    read_pytheas_file("isodist_atomdef_file")
    lines = pgv.text_lines.copy()

    nl = 0
    pgv.atom_dict = {}
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
        pgv.atom_dict[atom_symb] = atom_entry
        
def create_model():

    read_pytheas_file("isodist_modeldef_file")
    lines = pgv.text_lines.copy()

    n_var_atom = 0

    reslib_comment = " ".join(lines[0]).split("=")[0]
    print("MODEL: ",reslib_comment)
    
    # SPECIES DICT
    n_species = int(lines[1][0])
    species_lab = []
    amp = []
    ptr = 2
    pgv.species_dict = {}
    for species in range(n_species):
        species_lab.append(lines[ptr][0])
        amp.append(float(lines[ptr][1]))
        ptr = ptr + 1
    natom_types = int(lines[ptr][0])
    ptr = ptr + 1
    
    pgv.species_dict["n_species"] = n_species
    pgv.species_dict["species_lab"] = species_lab
    pgv.species_dict["amp"] = amp
    pgv.species_dict["natom_types"] = natom_types
    
    pgv.n_species = n_species
    
    # ATOM ORDER DICT
    pgv.atom_order_dict = {}
    pgv.atom_fix_dict = {}
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
        pgv.atom_order_dict[atom_symb] = atom_idx
        pgv.atom_fix_dict[atom_symb]=[atom_fix,fix,frac]
        atom_idx = atom_idx + 1
        ptr = ptr + 1
    pgv.species_dict["frac"] = frac
    
    pgv.model = ISODIST_MODEL(pgv.species_dict)
    

    pgv.model.n_var_atom = n_var_atom
    pgv.model.var_atom = var_atom
    # omitting variable residue labeling for now
    # pgv.model.n_var_res = n_var_res
    # pgv.model.n_fix_res = n_fix_res
    # pgv.model.var_res = var_res
    # pgv.model.phi = phi
    # return(species_dict,atom_order_dict,atom_fix_dict,residue_dict)    



def create_atom_mu():
# generate mu domain spectra for elements in order for residue definitions
    atom_mu = []   # build array as list
    pgv.atom_keys = list(pgv.atom_order_dict.keys())
    for key in pgv.atom_keys:        
        niso,miso_list,fiso_list = [pgv.atom_dict[key][k] for k in ["niso","miso","fiso"]]
        el_mz = [0] * pgc.npt
        for i in range(niso):
            n=int(miso_list[i]*pgc.scale_mz+0.5)
            el_mz[n] = fiso_list[i]
        el_mu = rfft(el_mz)
        atom_mu.append(el_mu)
    pgv.atom_mu = np.array(atom_mu)  # convert to numpy array

def initialize_residue_mu(active_residues):
    res_keys = [key for key in active_residues if len(key) == 1]
    res_keys = active_residues
    res_fixed_mu = []
    idx = 0
    pgv.res_order_dict = {}
    for reskey in res_keys:
        pgv.res_order_dict[reskey] = idx
        idx += 1
        # resform = pgv.residue_dict[reskey]["resform"]
        resform = [pgv.nt_fragment_dict[label][reskey].mol_form for label in pgv.isotopic_species]
        species_mu = []
        species_var_mu = []
        for j in range(pgv.n_species):
            res_mu = [1.0 + 0.0j] * pgc.ncp  # initialize complex mu
            species_var_mu.append(res_mu)
            # for atomkey in pgv.atom_keys:
            for atom, atidx in pgv.atom_order_dict.items():
                # atidx = pgv.atom_order_dict[atomkey]
                atom_fix,fix,frac = pgv.atom_fix_dict[atom]
                # print(reskey, j, atom)
                if fix == 0:
                    if atom not in resform[j]:
                        continue
                    n = resform[j][atom]
                    amu = list(map(lambda x:pow(x,n),pgv.atom_mu[atidx]))
                    res_mu = np.multiply(res_mu,amu)
            species_mu.append(res_mu)
        res_fixed_mu.append(species_mu)
    pgv.res_fixed_mu = np.array(res_fixed_mu)    # convert to numpy array
    pgv.res_tot_mu = copy.deepcopy(pgv.res_fixed_mu)     # copy for multiplying fixed and var parts       

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


def isodist_setup():
    find_active_residues()  
    
    
    # mz = np.zeros(0)
    # for i in range(igc.npt):            # mz axis
    #     mz = np.append(mz,float(i+1)/igc.scale_mz)
    
    #TODO check if 1st point should be zero or, as above
    
    pgv.mz = np.linspace(0, pgc.npt/pgc.scale_mz, pgc.npt)  # calc mz axis
    pgv.cmz = rfftfreq(pgc.npt) # complex calc mz axis
    create_atom_dict()  # input atom parameters   jgv.atom_dict
    create_model() # pgv.species_dict  pgv.atom_order_dict pgv.atom_fix_dict pgv.model
    create_atom_mu()  # generate mu domain spectra for atoms  (natom_types, ncp) 
    initialize_residue_mu(pgv.active_residues)  # generate mu domain spectra for residues
    
    initialize_fit_model()


def calc_isodist(x, return_opt, seq, yobs, mz_hd, z):
     
     #  fit function to take fit parameters x and return residuals eft_obs - sp_mz  
 
    # generate mu domain spectra for variable atoms
   
    xround = []
    for i in range(len(x)):  # for printing 
        xround.append(round(x[i],3))
 
    # calculate mu domain spectra for variable atoms
    for j in range(pgv.model.n_var_atom):
        atom = pgv.model.var_atom[j]   # N15
        miso_list = pgv.atom_dict[atom]["miso"]
        vidx = pgv.atom_order_dict[atom]
        fr = x[vidx]
        el_mz = [0] * pgc.npt
        nu = int(miso_list[0]*pgc.scale_mz+0.5)
        nf = int(miso_list[1]*pgc.scale_mz+0.5)
        el_mz[nu] = 1 - fr
        el_mz[nf] = fr
        pgv.atom_mu[vidx] = rfft(el_mz)
     
    active_res = set(list(seq))  # list of active residues to recalc variable part

    # calculate variable part of residue spectra 
    # indices are residue,species 
     
    for reskey in active_res:
        res_idx = pgv.res_order_dict[reskey]
        
        resform = []
        for label in pgv.isotopic_species:
            rf = pgv.nt_fragment_dict[label][reskey].mol_form
            resform.append(rf)
            
        for j in range(pgv.model.n_species):
            res_mu = np.ones(pgc.ncp, dtype = complex)
            for k in range(pgv.model.n_var_atom):
                atom = pgv.model.var_atom[k]
                atidx = pgv.atom_order_dict[atom]
                atom_fix,fix,frc = pgv.atom_fix_dict[atom]
                if fix == 1:
                    if atom not in resform[j]:
                        continue
                    n = resform[j][atom]
                    amu = np.ones(pgc.ncp,dtype = complex)
                    for nn in range(n):
                        amu = np.multiply(amu, pgv.atom_mu[atidx])

                    res_mu = np.multiply(res_mu,amu)
 
                pgv.res_tot_mu[res_idx][j] = np.multiply(pgv.res_fixed_mu[res_idx][j],res_mu)         
    
    baseline = x[0]
    off = x[1]
    gw = x[2]
 
    cgw = np.sqrt(2)*(gw + 0.0j)
    gw2p=gw*np.sqrt(2.0*np.pi)
    cfac_const = 2.0*np.pi*((off+mz_hd)*pgc.scale_mz)

    cpft = np.exp((0.0 + pgv.cmz*cfac_const*1j))
    cgft = np.exp(-((pgv.cmz +0j)/cgw)**2)/gw2p
    
    tot_mu = np.zeros(pgc.ncp,dtype = complex)   #initialize total spectrum
    residues = list(seq)
    for j in range(pgv.model.n_species):
        sp_mu = np.ones(pgc.ncp, dtype = complex)  #initialize
        for reskey in residues:
            # res_idx = pgv.residue_dict[reskey]["res_idx"]
            res_idx = pgv.res_order_dict[reskey]
            sp_mu = np.multiply(sp_mu, pgv.res_tot_mu[res_idx][j])  # accumulate product of residue spectra
        for k in range(z):
            sp_mu = np.multiply(sp_mu, np.conj(pgv.atom_mu[pgv.atom_order_dict['H1']])) # subtract z protons
        sp_mu = np.multiply(sp_mu,cgft)  # gaussian convolution
        sp_mu = np.multiply(sp_mu,cpft)  # heterodyne frequency shift
        tot_mu = tot_mu + x[pgv.model.amp_idx[j]] * sp_mu   # accumulate species
    
        
# back transform and add baseline

    sp_mz = irfft(tot_mu) * pgc.ncp/2 + baseline # normalization applied for consistency with isodist.f
   
    if return_opt == "residuals":
    
    # calculate residuals yobs - sp_mz using mz_ptr
        residuals = []
        mz_residuals = []

        for i in range(len(yobs)):
            residuals.append(yobs[i]-sp_mz[pgv.mz_ptr[i]])
            mz_residuals.append(pgv.mz[i])
        return(residuals)
    
    elif return_opt == "spectrum":
    
        spectrum =[]
        mz_spectrum = []
        for i in range(len(pgv.mz_ptr)):
            spectrum.append(sp_mz[pgv.mz_ptr[i]])
            mz_spectrum.append(pgv.mz[pgv.mz_ptr[i]])

        mz_spectrum = np.array(mz_spectrum) /z + mz_hd
        spectrum = np.array(spectrum)
        
        return mz_spectrum, spectrum
     

def initialize_fit_model():
    pgv.parlabel = ["B","OFF","GW"]
    pgv.x_init = [1.0, 0.01, 0.003]  # initial values for B, OFF, GW
    lower_bounds = [-1000.0,-1.0,0.001]
    upper_bounds = [np.inf,1.0,0.03]
    fitpar_idx = 2
    
    # amplitude parameters
    alab = "AMP_"
    pgv.model.amp_idx = []
    for i in range(pgv.model.n_species):
        fitpar_idx += 1
        pgv.model.amp_idx.append(fitpar_idx)
        pgv.parlabel.append(alab + pgv.model.species_lab[i])
        pgv.x_init.append(pgv.model.amp[i])
        lower_bounds.append(0.0)
        upper_bounds.append(np.inf)
    
    # fractional atom parameters
    flab = "FRC_"
    pgv.model.frac_idx = []
    for i in range(pgv.model.n_var_atom):
        fitpar_idx += 1
        pgv.model.frac_idx.append(fitpar_idx)
        pgv.parlabel.append(flab + pgv.model.var_atom[i])
        pgv.x_init.append(pgv.model.frac[i])
        lower_bounds.append(0.0)
        upper_bounds.append(1.0)
    
    pgv.par_err_labels = [l + "_err" for l in pgv.parlabel]  # parameter errors from fit
    pgv.par_relerr_labels = [l + "_re" for l in pgv.parlabel] # relative errors
    pgv.frac_lab_labels = ["FRAC_" + l.split("_")[-1] for l in pgv.parlabel if "AMP" in l] # frac_lab values
    pgv.xbounds = (lower_bounds,upper_bounds)
 
    print("number of params to fit :",len(pgv.x_init),pgv.x_init)
     
    # column labels for isodist output

    pgv.column_labels = pgc.outlabels + pgv.parlabel + pgv.frac_lab_labels + pgv.par_err_labels + pgv.par_relerr_labels
    pgv.column_labels = pgv.column_labels + ["max_fit","min_fit","max_rsd","min_rsd"]
    nfit = 11  
    nwr = 23
    for nf in range(nfit):
        pgv.column_labels.append("avg_fit" + str(nf+1))
    for nw in range(nwr):
        pgv.column_labels.append("avg_wr" + str(nw+1))
        

def fit_isotope_distribution(match_df_idx):
    match_row = pgv.top_match_df.loc[match_df_idx]
    ms = pgv.minispec[match_df_idx]
    
    frag3 = match_row["frag3"]
    seq = generate_mod_seq_ends(frag3)
    z = abs(ms.z)
    moz = match_row["mz1"]

    print("SEQUENCE, z = ", seq, z)
    
    peak_root = "_".join([str(match_row["ms2_key"]), ms.mol, seq, str(z)])
    seqlab = ms.mol
    
    ri_array = ms.rt_slice  # experimental spectrum
    rmz_array = ms.mini_mz

    nmz = 0 # counter to fill experimental arrays
    pgv.mz_ptr = np.zeros(len(rmz_array), dtype = int)
    xobs = np.zeros(len(rmz_array))
    yobs = np.zeros(len(rmz_array))
    mz_hd = rmz_array[0] * z  # heterodyne mass set to lower bound of mz
    
    # build arrays for fitting, mz_ptr has indices from oversampled theo distribution
    for rmz, ri in zip(rmz_array, ri_array):
        n = int((rmz * z - mz_hd)*pgc.scale_mz + 0.5)

        pgv.mz_ptr[nmz] = n
        xobs[nmz] = (rmz * z - mz_hd) * pgc.scale_mz
        yobs[nmz] = ri
        nmz += 1
        
    x_init = pgv.x_init.copy()
    x_init[3] = max(yobs)
    x_init[4] = max(yobs)
    
#TODO are GW and AMP anticorrelated?  
    
    if "showguess" in pgv.isodist_plot_options:     # calc initial guess
        
        mz_spec, spec = calc_isodist(x_init, "spectrum", frag3, yobs, mz_hd, z)
        resid = calc_isodist(pgv.x_init, "residuals", frag3, 12*yobs, mz_hd, z)
        plot_obs_calc(mz_spec, yobs, mz_spec, spec, seq + " Initial guess")
     
    prelim_end = datetime.now()
    fit_start = prelim_end

    lsq_soln = least_squares(calc_isodist, pgv.x_init, verbose = 1, bounds = pgv.xbounds, 
                             x_scale = 'jac', max_nfev=100,
                             args =("residuals", frag3, yobs, mz_hd, z))
    x_fit = lsq_soln.x
    chisq = 2.0 *lsq_soln.cost/(pgc.sig_global*pgc.sig_global*nmz)  # factor of 2 to make cost correpsond to chisquared
        
    xfit = [float(round(x,6)) for x in x_fit]
    
    print("fit time = ", datetime.now()- fit_start)
    print(" solution = ", xfit)
           
    mz_fit, fit = calc_isodist(x_fit, "spectrum", frag3, yobs, mz_hd, z)
    resid = calc_isodist(pgv.x_init, "residuals", frag3, yobs, mz_hd, z)  # can get this from lsq_soln

    if "showfit" in pgv.isodist_plot_options:
        plot_obs_calc(mz_fit, yobs, mz_fit, fit, seq + " Final Fit")
   
    # precalculate frac_lab
    tot_amp = sum([x for x, l in zip(x_fit,pgv.parlabel) if "AMP" in l])
    if tot_amp != 0:
        frac_lab = [x/tot_amp for x, l in zip(x_fit,pgv.parlabel) if "AMP" in l]
    else:
        frac_lab = [0 for x, l in zip(x_fit,pgv.parlabel) if "AMP" in l]
    
   # parameter errors
   
   # Extract the Jacobian and calculate covariance and standard errors
    J = lsq_soln.jac
    m = len(rmz_array)  # Number of data points
    n = len(lsq_soln.x)   # Number of parameters
    s_squared = 2 * lsq_soln.cost / (m - n) # Estimated variance of residuals
    
    # Calculate the covariance matrix
    try:
        cov_matrix = np.linalg.inv(J.T @ J) * s_squared
        
        perr = np.sqrt(np.diag(cov_matrix))    # Standard errors
        print(f"Fitted parameters: {lsq_soln.x}")
        print(f"Standard errors of parameters: {perr}")
    except np.linalg.LinAlgError:
        print("Could not compute covariance matrix. Check for singular Jacobian.")
        perr = np.zeros(n)
   
    rel_err = np.zeros(n)
    for i in range(n):
        if x_fit[i] != 0:
            rel_err[i] = perr[i]/x_fit[i]
        
    max_fit, min_fit, max_rsd, min_rsd, avg_fit, avg_wr = calculate_residuals(x_fit, fit, resid)
    
    pars = pgv.parlabel
    par_labels = pars
    # par_labels = ["_".join(["rt", p, "fit"]) for p in pars]
    par_err_labels = ["_".join([p, "err"]) for p in pars]
    # par_err_labels = [l + "_err" for l in parlabels]
    par_relerr_labels = ["_".join([p, "re"]) for p in pars]
    # par_relerr_labels = [l + "re" for l in parlabels]
    
    for p, v in zip(par_labels, x_fit):
        setattr(ms, p, v)
    for p, v in zip(par_err_labels, perr):
        setattr(ms, p, v)
    for p, v in zip(par_relerr_labels, rel_err):
       setattr(ms, p, v)

    
    
    
    # assemble output columns
    column_data =[peak_root, seqlab, seq, moz, z, chisq, mz_hd]
    column_data = column_data + x_fit.tolist() + frac_lab + perr.tolist() + rel_err.tolist() + [max_fit, min_fit, max_rsd, min_rsd] + avg_fit + avg_wr
    
    column_dict = {col:dat for col, dat in zip(pgv.column_labels, column_data)}

# todo store this in minispectrum

    ms.mz_fit = fit

#TODO make this a function
    if "plotfit" in pgv.isodist_plot_options:
# write out fit plot
        plotfile = os.path.join(pgv.isodist_plot_dir, peak_root + ".png")
        plot_minispectrum(ms)
        # plt.plot(ms.mini_mz, yobs,".k")
        # plt.ylabel("amplitude")
        # plt.xlabel("m/z")
        # plt.title(seq + " m0 = " + str(round(mz_hd,3)))
        # plt.plot(ms.mini_mz, fit, "-r")
        plt.savefig(plotfile, dpi = 300)
        

    return column_dict



def calculate_residuals(x_fit, fit, resid):
#   RESIDUAL EXTRAVAGANZA
#   used outside of isodist for neural net analysis of goodness of fit

    # numpy ARRAYS FOR RESIDUAL CALCULATION:  resid, netfit, resid_m, resid_mm, netfit_m, netfit_mm
    # resid is residuals = "dy"
    
    resid = - np.array(resid)  # old resid = fit - dat
    fit = np.array(fit)
    
    # if pgv.auto_baseline == "auto":
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

def Load_thru_RT():
    
    load_dir = QFileDialog.getExistingDirectory(None,"Select Job Directory", pgv.working_dir)
    if load_dir == None:
        return
    print("loading previous RT files from ", load_dir)
    load_pickle_set(pgc.isodist_rt_pickle, load_dir)

def Load_thru_Minispec():
    
    load_dir = QFileDialog.getExistingDirectory(None,"Select Job Directory", pgv.working_dir)
    if load_dir == None:
        return
    print("loading previous RT files from ", load_dir)
    load_pickle_set(pgc.isodist_rt_pickle, load_dir)
    print("loading previous minispec from ", load_dir)
    load_pickle_set(pgc.isodist_minispec_pickle, load_dir)
    
def Load_thru_RT_fit():
    
    load_dir = QFileDialog.getExistingDirectory(None,"Select Job Directory", pgv.working_dir)
    if load_dir == None:
        return
    print("loading previous RT files from ", load_dir)
    load_pickle_set(pgc.isodist_rt_pickle, load_dir)
    print("loading previous minispec from ", load_dir)
    load_pickle_set(pgc.isodist_minispec_pickle, load_dir)
    print("loading previous RT_fit from ", load_dir)
    load_pickle_set(pgc.isodist_rtfit_pickle, load_dir)
    
    pgv.minispec = pgv.minispec_rt.copy()
    
def Load_thru_Isodist():
    
    load_dir = QFileDialog.getExistingDirectory(None,"Select Job Directory", pgv.working_dir)
    if load_dir == None:
        return
    print("loading previous RT files from ", load_dir)
    load_pickle_set(pgc.isodist_rt_pickle, load_dir)
    print("loading previous minispec from ", load_dir)
    load_pickle_set(pgc.isodist_minispec_pickle, load_dir)
    # print("loading previous RT_fit from ", load_dir)
    # load_pickle_set(pgc.isodist_rtfit_pickle, load_dir)
    print("loading previous isodist fit from ", load_dir)
    load_pickle_set(pgc.isodist_fit_pickle, load_dir)
    
    pgv.minispec = pgv.minispec_isodist.copy()
   
    
    
