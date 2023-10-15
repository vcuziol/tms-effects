# Parameters used in the simulation

output_filename = 'single_pulse'
onset = 2
measures = ['L5_spikes', 'L5_mem']

# simparams
nr_trials = 10     # number of L5 neurons
steps_per_ms = 40  # simulation time steps
tstop = 15         # simulation stopping time
celsius = 36       # simulation temperature
parallel = 0       # indicating use of parallelization


# the layer 5 cell (using the same topology as Larkum et al. 2009)
l5_file = '070603c2.cll'

# L2/3 neurons
nr_neurons = {'AMPA' : 190,    # number of neurons with AMPA-type synapses
             'NMDA' : 50,    # number of neurons with NMDA-type synapses
             'GABAa': 40,    # number of neurons with GABAa-type synapses
             'GABAb': 20}    # number of neurons with GABAb-type synapses

syn_placement_uniform = [0.001, 0.999]

# axon length and diameter
axon_L = 10     # for inital part of axon
axon_diam = 5

# number of segments
dend_nseg = 7
soma_nseg = 5
axon_nseg = 11
axon1_nseg = 11

# passive cable properties
e_pas = -70
g_pas = 1.0/20000
cm = 1

tauh_hh3 = 0.1
sN_hh3 = 0

# Biophysical parameters, see paper, p.414
apic_k1  = 0.0                         # dendritic potassium
apic_k2  = 0.0001                      # and sodium
apic_na  = 0.004                       # conductances
apic_kdr = 0
apic_gl  = 0.0001

soma_k1  = 0.06                        # somatic potassium
soma_k2  = 0.3                         # and sodium
soma_na  = 0.008                       # conductances
soma_kdr = 0
soma_gl  = 0

myelin_k1 = 0.06
myelin_k2 = 0.0001
myelin_na = 0.1
myelin_cm = 0.04

g_pas_node = 0.02
axon_k1  = 1                           # axon potassium
axon_k2  = 0                           # and sodium
axon_na  = 1.9                         # conductances
axon_kdr = 10                          # delayed-rectifier potassium current
axon_gl  = 0


# Synapses, see paper, p.403
# rise times, decay times, peak conductances, and reversal potentials
#  for each receptor type

syn_params = {}

syn_params['AMPA'] = {'TAU1': 0.2,
                    'TAU2': 1.7,
                    'GPEAK': 0.2,
                    'E': 0 }

syn_params['NMDA'] = {'TAU1': 2,
                    'TAU2': 26,
                    'GPEAK': 0.03,
                    'E': 0 }

syn_params['GABAa'] = {'TAU1': 0.3,
                    'TAU2': 2.5,
                    'GPEAK': 0.5,
                    'E': -70 }

syn_params['GABAb'] = {'TAU1': 45.2,
                     'TAU2': 175.16,
                     'GPEAK': 0.05,
                     'E': -70   }

syn_params['STD'] = {'tau_D1': 200,
                     'd1': .5}

#rates of the poisson spiketrains
l23_rates= {'AMPA'  :  6.,
            'NMDA'  :  6.,
            'GABAa' :  12.,
            'GABAb' :  12.}

l23_shift = {'AMPA'  :  0.,
             'NMDA'  :  0.,
             'GABAa' :  0.,
             'GABAb' :  0.}

mode = 'mu-rhythm-shift'

l23_5_weights_std = 0.5
l23_5_weights_i = -1.5          # excitatory synapse weights
l23_5_weights_e = -1.5          # inhibitory synapse weights
l23_5_delay = 1.0               # synaptic transmission delay

# PARAM_TMS_ONSET
TMS_I_l5_scale = 6
TMS_I_l23_scale_inh = 116
TMS_I_l23_scale_exc = 74
