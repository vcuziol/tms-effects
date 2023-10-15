from neuron import h

def load_membrane_mechanisms(cell, params):
    """
    function for setting up the membrane mechanisms and their parameters
    
    input:
        dends (list): list of all dendrite sections
        axon1 (neuron section): the long axon
    """
    
    dends = list(cell.dend) + list(cell.apic)
    
    h.tauh_hh3 = params.tauh_hh3
    h.sN_hh3 = params.sN_hh3
    
    for sec in dends:               # add HH and leakage to all dendrites
        sec.nseg = params.dend_nseg
        sec.insert('pas')
        sec.e_pas = params.e_pas
        sec.g_pas = params.g_pas
        sec.cm = params.cm
        sec.insert('hh3')
        sec.gkbar_hh3 = params.apic_k1
        sec.gkbar2_hh3 = params.apic_k2
        sec.gnabar_hh3 = params.apic_na
        sec.gl_hh3 = params.apic_gl
        sec.insert('kdr')
        sec.gbar_kdr = params.apic_kdr
        
    for sec in cell.soma:
        sec.insert('hh3')            # add HH and leakage to the soma
        sec.gnabar_hh3 = params.soma_na
        sec.gkbar_hh3 = params.soma_k1
        sec.gkbar2_hh3 = params.soma_k2
        sec.gl_hh3 = params.soma_gl
        sec.nseg = params.soma_nseg
        sec.insert('pas')
        sec.e_pas = params.e_pas
        sec.g_pas = params.g_pas
        sec.cm = params.cm
        sec.insert('kdr')
        sec.gbar_kdr = params.soma_kdr


def load_membrane_mechanisms_axon(cell, params):
    
    print "running load_membrane_mechanisms_axon. "
    
    axonal_default = [cell.axon[0]] + list(cell.unmyelin)
    
    for sec in axonal_default:
        sec.insert('hh3')            # add HH, dr type K current, and leakage to the original axon
        sec.gnabar_hh3 = params.axon_na
        sec.gkbar_hh3 = params.axon_k1
        sec.gkbar2_hh3 = params.axon_k2
        sec.gl_hh3 = params.axon_gl
        sec.insert('kdr')
        sec.gbar_kdr = params.axon_kdr
        sec.nseg = params.axon_nseg
        sec.insert('pas')
        sec.e_pas = params.e_pas
        sec.g_pas = params.g_pas
        sec.cm = params.cm
    
    for sec in cell.myelin:
        sec.insert('pas')
        sec.cm = params.myelin_cm
        sec.e_pas = params.e_pas
        sec.g_pas = params.g_pas
        sec.insert('hh3')
        sec.gnabar_hh3 = params.axon_na
        sec.gkbar_hh3 = params.axon_k1
        sec.gkbar2_hh3 = params.axon_k2
        sec.insert('kdr')
        sec.gbar_kdr = params.axon_kdr
    
    for sec in cell.node:
        sec.insert('pas')
        sec.cm = params.cm
        sec.g_pas = params.g_pas_node
        sec.e_pas = params.e_pas
        sec.insert('hh3')
        sec.gnabar_hh3 = params.axon_na
        sec.gkbar_hh3 = params.axon_k1
        sec.gkbar2_hh3 = params.axon_k2
        sec.gl_hh3 = params.axon_gl
        sec.insert('kdr')
        sec.gbar_kdr = params.axon_kdr

##    all_axon = list(cell.node) + list(cell.axon) + list(cell.myelin) + list(cell.unmyelin)
##
##    for sec in cell.axon:
##        sec.insert('hh3')            # add HH, dr type K current, and leakage to the original axon
##        sec.gnabar_hh3 = params.axon_na
##        sec.gkbar_hh3 = params.axon_k1
##        sec.gkbar2_hh3 = params.axon_k2
##        sec.gl_hh3 = params.axon_gl
##        sec.insert('kdr')
##        sec.gbar_kdr = params.axon_kdr
##        sec.nseg = params.axon_nseg
##        sec.insert('pas')
##        sec.e_pas = params.e_pas
##        sec.g_pas = params.g_pas
##        sec.cm = params.cm
