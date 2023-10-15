# -*- coding: utf-8 -*-
import os
import sys
import numpy as np
from mpi4py import MPI




# ---- paths
ROOT_DIR = os.path.abspath('./tms_networks')
MAIN_DATA_FOLDER = ROOT_DIR + '/core/data'
BBP_MODELS_FOLDER = ROOT_DIR + '/core/AberraEtAl2018_edit/cells'
BBP_MECHS_FOLDER = ROOT_DIR + '/core/AberraEtAl2018_edit/mechanisms'

BMTK_DIR = ROOT_DIR + '/bmtk_env'
POOL_DIR = ROOT_DIR + '/pool'
POP_FOLDERS_PATH = ROOT_DIR + '/simulation_temp/pop_folders.txt'
TMS_NETWORK_DIR = ROOT_DIR

MSH_INPUTS_FOLDER = ROOT_DIR + '/preprocessing/msh_inputs'
msh_folders = os.listdir(MSH_INPUTS_FOLDER)

COMPILE_MECHS = True

# ---- parameters for preprocessing and simnibs
show_coil_basis = False
number_of_coil_orientations = 7
simnibs_mat_fn = MSH_INPUTS_FOLDER \
    + '/simnibs_simulation/simnibs_simulation_20200117-043700.mat' # basefile
simnibs_bin_path = os.environ['HOME'] + '/simnibs_2.1.2/bin'
simnibs_command = simnibs_bin_path + '/simnibs'
get_fields_command = simnibs_bin_path + '/get_fields_at_coordinates'

# ---- simulation preferences
RUN_SIMULATION = True
run_parallel = 0
n_procs = 7
plotFilesPrefix = 'test_cell'
DUMP_NEURON_3D_DATA = 'obj' # 'pts', 'obj' or ''
DUMP_VMEM = False
DUMP_AP_SITE_MEASURES = True

# where to store results for the 'tms_isolated_neurons' simulations.
RESULTS_FOLDER = ROOT_DIR + '/results_isolated_neurons/'

# segment of the L5 neuron which will be used to compose the epidural response.
L5_AXON_SEG_INDEX = 1000
# old: recorded_segment_index = 1200

# ---- model parameters
# depth of the neurons in the gray matter, given in percentage of
# approximated cortical depth at a point in the gray matter surface.
# numbers based loosely on the ones used by Aberra et al 2019:
# L1: 0.06, L2/3: 0.4, L4: 0.55, L5: 0.65, L6: 0.85.
relative_depths = {
##'L1': 0.06,
'L23': 0.42,
##'L4': 0.55,
'L5': 0.61
##'L6': 0.85
}

# NEURON simulation parameters
tstop = 8.0 # ms
dt = 0.0125 # ms
v_init = -70.0432010302 # mV

# range of stimulation intensities to be used in the binary search.
# In SimnNIBS, the stimulation intensity is in units of A/μs (with
# default value of 1*10^6 A/s = 1.0 A/μs).
thresholds_search_range = [0.0, 1500.0] # A/μs

# parameters of the RLC circuit for simulation of the TMS waveform (I)
# and electric field time course (dI/dt)
R = 0.5
L = 16
C = 200
wf_delay = 1.0
wf_dur = 1.0
wf_amp = 1.0

# indices of neurons to be simulated, if it is desired to do the
# simulations only for specific neurons in the population.
# set it to None if you want to simulate all N neurons in the population.
NEURONS_TO_SIMULATE = range(0, 101) 

# default electric field (E-field) vector to *override* the E-field
# vectors loaded from SimNIBS simulation files.
# set it to None if you want to use the E-field vectors calculated by SimNIBS.
DEFAULT_E_VECTOR = None #np.array([0.0, -1.0, 1.0])

# if you want to assume that E-field is uniform around each single neuron,
# leave this as False.
INTERPOLATE_AT_SEGMENTS = False

# ----- parameters
config_filename = 'simulation_config.json'
orig_folder = os.getcwd()


# which model to use: "pure" Blue Brain Project neuron (directive: 'bbp'),
# or my StimCell class with Aberra et al. 2018
# modifications (directive: 'StimCell.py').
model_template_directive = 'StimCell.py'

populations_names = ['pop_L23',
                     'pop_L23_NMDA',
                     'pop_L23_GABAa',
                     'pop_L23_GABAb',
                     'pop_L5']

### change this if you would like to set the path to the neurons positions
### file manually.
### if None, the default is used.
##CUSTOMIZED_DATA_FOLDER = '/home/master/Desktop/msh_inputs_old/data/'
##MAIN_DATA_FOLDER = CUSTOMIZED_DATA_FOLDER

# config files for synapse models
SYNAPSE_FILES = ['NMDA_ExcToExc.json',
                 'GABAb_InhToExc.json',
                 'InputSynapses_ExcToExc.json',
                 'AMPA_ExcToExc_modified.json',
                 'GABA_InhToExc_modified.json']

# -------------

# number of neurons in each population
# it's written like this to facilitate the automatic editing done by
# the parameter search script.
NUM_SCALE = 5
##N_pop_L23 = 19*NUM_SCALE
##N_pop_L23_GABAa = 4*NUM_SCALE
##N_pop_L23_NMDA = 5*NUM_SCALE
##N_pop_L23_GABAb = 2*NUM_SCALE
N_pop_L23 = 19*NUM_SCALE
N_pop_L23_GABAa = 4*NUM_SCALE
N_pop_L23_NMDA = 5*NUM_SCALE
N_pop_L23_GABAb = 2*NUM_SCALE

network_neurons_count = {
'pop_L23':          N_pop_L23,
'pop_L23_GABAa':    N_pop_L23_GABAa,
'pop_L23_NMDA':     N_pop_L23_NMDA,
'pop_L23_GABAb':    N_pop_L23_GABAb,
'pop_L5': 1
}

pop_input_N = 2500

# estimative of id of L5 node
L5_node_id = sum(network_neurons_count.values())-1

# synapse template to be used from L2/3 to L5 connections (e.g., Exp2Syn).
# ASyN_STD: from 'alpha_STD.mod'.
synapse_template = 'ASyN_STD'

# angle of coil orientation relative to central sulcus
# (0°: perpendicular to central sulcus)
global COIL_ANGLE
COIL_ANGLE = 270

# -------------

# NEURON simulation parameters
sim_tstop = 20.0
sim_dt = 0.025

# neurons positions
USE_POP_POSITIONS = 1
USE_DUMPED_POSITIONS = 0 # load positions from 'network_neurons_positions.dat'

# background activity
USE_PSG = 0

# TMS parameters
USE_TMS_STIM = 1
TMS_STIM_INTENSITY  = 400.0
TMS_PARAMETERS = {'stim_intensity': TMS_STIM_INTENSITY}
##wf_params = {'dt':sim_dt, 'tstop':sim_tstop, # use with psg
##             'Rc':0.5, 'Lc':16, 'Cc':200,
##             'DEL':55.0, 'DUR':1.0, 'AMP':1.0}
# dur and amp are to be kept 1.0; changing them has no effect
wf_params = {'dt':sim_dt, 'tstop':sim_tstop,
             'Rc':0.5, 'Lc':16, 'Cc':200,
             'DEL':5.0, 'DUR':1.0, 'AMP':1.0}

# synapses parameters
# synaptic weights following Rusu et al. 2014:
#   inh: log normal, (0.2, 0.5).
#   exc: log normal, (0.1, 0.5).
base_weight  = -1.3
weight_0  = base_weight
weight_1  = 2.0 * base_weight
weight_2  = base_weight
weight_3  = 2.0 * base_weight
weight_sigma  = 0.5

# synapses proportion (exc:inh):
# glutamatergic: 3 to 12 synapses; GABAergic: 5 to 30 synapses.
# Microcircuitry of the Neocortex, Rawasmamy et al.
nsyn_scale  = 1
inh_syn_num = 2 # equal number of synapses between exc and inh.
##exc_syn_num = 2
exc_syn_num = 4
nsyn_scale_0   = exc_syn_num*nsyn_scale
nsyn_scale_1   = inh_syn_num*nsyn_scale
nsyn_scale_2   = exc_syn_num*nsyn_scale
nsyn_scale_3   = inh_syn_num*nsyn_scale

# multiplicative factor for delay in synapses
# 1e-5 is too small of a delay and causes 'mindelay' NEURON errors.
delay_factor = 5.5e-05

# input network params
N_thalamus = 70
thalamus_net_name = 'mthalamus'

# name of each thalamus network
thalamus_net_names = [thalamus_net_name+'_0', thalamus_net_name+'_1']

# name of populations to which the thalamus networks will project
thalamus_target_pops = {thalamus_net_names[0]: \
                            ['pop_L23', 'pop_L23_NMDA'],
                        thalamus_net_names[1]: \
                            ['pop_L23_GABAa', 'pop_L23_GABAb']}

# input network
# A constant firing rate
##times = (0.00, sim_tstop) #sim_tstop-(sim_tstop*0.8))
##firing_rates = {thalamus_net_names[0]: 6.0,
##                thalamus_net_names[1]: 12.0}

base_rate = 6.0 # Hz
max_rate = 2*base_rate
phi = 0.0
mu_freq = 10.0 # Hz
const_freq = 1.0

# Rusu et al. 2014 (via Schaworonkow et al. 2018)
base_rates = {thalamus_net_names[0]: 6.0, 
                thalamus_net_names[1]: 12.0}

## Uncomment to model the input firing rates on a sine wave
# times = np.linspace(0.0, 3.0, 1000)
# firing_rate = 10.0*np.sin(times) + 10.0
N_times = sim_tstop/sim_dt
times = np.linspace(0.0, sim_tstop, N_times)

# from formula: A*sin(2*pi*freq*t-phi)
firing_rates = {thalamus_net_names[0]: \
                 base_rates[thalamus_net_names[0]] \
                 * (np.sin(2.0*np.pi*mu_freq*times-phi) + const_freq),
                thalamus_net_names[1]: \
                 base_rates[thalamus_net_names[1]] \
                 * (np.sin(2.0*np.pi*mu_freq*times-phi) + const_freq)}

#debug
nsynmax_thalamus = 10
synw_thalamus = 1.0
maxrang_thalamus = 150.0*1e3

# input network
input_network_model = {
    'input_network': {
        'N': pop_input_N,
        'ei': 'e',
        'pop_name': 'input_network',
        'model_type': 'virtual'
    }
}

# ----- parameters log
rank = MPI.COMM_WORLD.Get_rank()
if rank == 0:
    print(
        ['network_neurons_count: ', network_neurons_count,
        'USE_PSG: ', USE_PSG,
        'tstop: ', sim_tstop, '\t', 'dt: ', sim_dt,
        #'times: ', times,
        #'firing_rates: ', firing_rates,
        'USE_DUMPED_POSITIONS: ', USE_DUMPED_POSITIONS,
        'USE_TMS_STIM: ', USE_TMS_STIM,
        'TMS_PARAMETERS: ', TMS_PARAMETERS,
        'wf_params: ', wf_params,
        'COIL_ANGLE: ', COIL_ANGLE,
        'L5_node_id: ', L5_node_id,
         ]
        )

if not os.path.isdir(simnibs_bin_path):
    print 'simnibs path not found: ', simnibs_bin_path
    exit()
