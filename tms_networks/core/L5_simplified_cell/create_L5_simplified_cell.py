import matplotlib.pyplot as plt
import numpy as np 

import neuron
from neuron import h

from .params import rusu as rusu_params
#from tms_networks.core.L5_simplified_cell.params import rusu as rusu_params
from .biophysics import load_membrane_mechanisms, load_membrane_mechanisms_axon

# -----
def create_cell(mechanisms_path=''):
    # I removed some mechanisms from 'L5_simplified_cell/mechanisms' folder
    # so that only the mechanisms that are not already 
    # in 'AberraEtAl2018_edit/mechanisms' are loaded.
    
    if len(mechanisms_path) > 0: 
        neuron.load_mechanisms(mechanisms_path)
    
    h.load_file('nrngui.hoc')
    h.load_file('import3d.hoc')
    h.load_file('template.hoc')
    
    global cell
    
    cell = h.cADpyr232_L5_TTPC2_8052133265(0)
    load_membrane_mechanisms(cell, rusu_params)
    
    return cell

def load_axon_biophysics(cell):
    load_membrane_mechanisms_axon(cell, rusu_params)

# ---------------

def test_current_inj():
    from neuron import gui
    shape_window = h.PlotShape()
    
    # -----
    v_vec = h.Vector()        # Membrane potential vector
    t_vec = h.Vector()        # Time stamp vector

    seg_to_record = cell.axon[0](0.5)
    v_vec.record(seg_to_record._ref_v)
    t_vec.record(h._ref_t)
    print "seg_to_record: ", seg_to_record.sec.name(), seg_to_record.x

    # -----

    current_inj = 15.0 # nA

    stim = h.IClamp(cell.soma[0](0.5))
    stim.delay = 50.0
    stim.dur = 50.0
    stim.amp = current_inj

    h.tstop = 150.0
    h.dt = 0.025
    print "running simulation: tstop: {} \t / \t dt: {}".format(h.tstop, h.dt)
    h.run()

    plt.plot(np.array(list(t_vec)), np.array(list(v_vec)))
    plt.axvline(100.0,c='r'); plt.axvline(50.0,c='r')
    plot_title = str(cell)+" Vm during \n current injection of {} nA".format(str(current_inj))
    plt.title(plot_title)
    fig_filename = str(cell)+"_current_"+str(current_inj)+"nA.png"

    plt.show()

    ##plt.savefig(fig_filename)
    ##print "saved fig: ", fig_filename
