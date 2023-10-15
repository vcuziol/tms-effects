# -*- coding: utf-8 -*-

import numpy as np
import pickle
import sys

from neuron import h

from ..parameters import *

class TMSSimulation:
    # -----------------------------------
    def __init__(self, simnibs_session_name=None):
        
        if simnibs_session_name is not None:
            self.load_electric_field(simnibs_session_name)
    # -----------------------------------
    def load_electric_field(self, simnibs_session_name):
        # gmmcE: gray matter motor cortex, electric field
        # wmmcE: white matter motor cortex, electric field
        
        # the electric field vectors are loaded assuming they are in V/m
        # (volts per meter, output unit of SimNIBS), which is equivalent
        # to mV/mm (mV is the unit of quasipotentials and extracellular
        # potential in NEURON).
        with open(MAIN_DATA_FOLDER + "/" + simnibs_session_name + \
                  "/fields/gm_files/gm_motor_E.dat") as file_to_load:
            gmmcE = pickle.load(file_to_load)
        with open(MAIN_DATA_FOLDER + "/" + simnibs_session_name + \
                  "/fields/wm_files/wm_motor_E.dat") as file_to_load:
            wmmcE = pickle.load(file_to_load)
        with open(MAIN_DATA_FOLDER + "/" + simnibs_session_name + \
                  "/meshes/motor_cortex_meshes/gm_motor_geom.dat") \
                  as file_to_load:
            gmmcgeom = pickle.load(file_to_load)
        with open(MAIN_DATA_FOLDER + "/" + simnibs_session_name + \
                  "/meshes/motor_cortex_meshes/wm_motor_geom.dat") \
                  as file_to_load:
            wmmcgeom = pickle.load(file_to_load)
        
        self.E_field_vectors = gmmcE
        
##        if INTERPOLATE:
##            from scipy.interpolate import interpn, LinearNDInterpolator
##            points = np.concatenate((gmmcigeom['centers'], wmmcgeom['centers']), axis=0)
##            values = np.concatenate((gmmcEi, wmmcE), axis=0)
##            self.E_field_interpolator = LinearNDInterpolator(points, values)
        
        # mcfi = motor cortex for interpolation (region bigger than
        # motor cortex (mc) used for creating the neuron populations)
##        self.mcfi_points = points #remov
        
    # -----------------------------------
    def generate_rlc_waveform(self, dt, tstop, DEL, DUR, AMP, Cc, Rc, Lc):
        
        h('dt = ' + str(dt))
        h('tstop = ' + str(tstop))
        
        # loads the "stim_waveform" function
        h.load_file(ROOT_DIR + '/core/codes/RLC_waveform.hoc')
        
        # create waveform: the result is stored in the "stim_amp" NEURON Vector
        h.stim_waveform(DEL, DUR, AMP, Cc, Rc, Lc)
        wf = np.array(h.stim_amp.x)
        
        # normalization
        wf = wf/np.max(wf)
        
        # damping ratio
        #zeta = lambda R,L,C: (R/2.0) * np.sqrt(C/L)
        
        self.stim_time_course = wf
        
        return self.stim_time_course
        
    # -----------------------------------
    def build_v_ext(self, v_seg_values, v_time_course):
        
        v_ext = np.zeros((len(v_seg_values), len(v_time_course)))

        for i in range(len(v_seg_values)):
            v_ext[i,:] = v_seg_values[i] * v_time_course

        return v_ext
