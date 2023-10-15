# files
import os
import sys
import shutil
import pickle
import time
import ast

# math
import numpy as np
import random

# neuron
import neuron

# tms_networks
from tms_networks.parameters import *
from .core import StimCell, NeuronPopulation, TMSSimulation
from .core.codes import linear_alg_routines as linar

from tms_networks.debugging import *

# this is declared global to be accessible inside the
# 'create_stimcell' function, which is called internally
# by BMTK to create a cell of the 'StimCell' class.
global bcl
bcl = []

try:
    os.mkdir(ROOT_DIR + '/results/')
    print("created folder: " + ROOT_DIR + '/results/')
except OSError:
    pass

# ------------------------------------------
# population

def create_population(simnibs_session_name,
                      population_cell_name,
                      waveform_params={'dt':0.025,
                                       'tstop':20.0,
                                       'Rc':0.5,
                                       'Lc':16,
                                       'Cc':200,
                                       'DEL':1.0,
                                       'DUR':1.0,
                                       'AMP':1.0},
                      plot_folder=None,
                      interp_E_path=None):

    session_data_folder = ROOT_DIR + "/core/data/" \
                          + simnibs_session_name + "/"
    population_cell_type = population_cell_name.split('_')[0]

    # load files
    layer_coords_dict = pickle.load(open(session_data_folder + \
                            "population_parameters/" + population_cell_type \
                            + "_soma_centers.dat", "r"))
    soma_centers = layer_coords_dict['layer_soma_centers']
    layer_neurons_distance = layer_coords_dict['layer_neuron_distances_factor']
    gm_normals = pickle.load(open(session_data_folder + 
        'population_parameters/' + population_cell_type + '_normals.dat', 'r'))
    gm_centers = pickle.load(open(session_data_folder + 
        "population_parameters/" + population_cell_type + "_gm_centers.dat", "r"))
    neurons_polygon_indices = pickle.load(open(session_data_folder 
        + "population_parameters/" + population_cell_type
        + "_neurons_indices.dat", "rb"))

    print("Cells positions and orientations loaded.")

    # tmssimulation object: configure magnetic stimulation of the population
    tms_sim = TMSSimulation.TMSSimulation(simnibs_session_name)

    # create RLC waveform
    tms_sim.generate_rlc_waveform(**waveform_params)

    # create population
    population = NeuronPopulation.Population(
            cell_name=population_cell_name,
            soma_centers=soma_centers,
            gm_normals=gm_normals,
            gm_centers=gm_centers,
            neurons_polygon_indices=neurons_polygon_indices,
            tms_sim=tms_sim,
            v_init=-70.0,
            tstop=waveform_params['tstop'],
            dt=waveform_params['dt'],
            simnibs_session_name=simnibs_session_name,
            results_folder=ROOT_DIR+'/results/',
            plot_folder=plot_folder
    )

    return population

# create Population objects for L23 and L5 neurons.
# this represents a full group of thousands of neurons;
# not all of them will necessarily be used in the network simulations.
simnibs_session_name = 'simnibs_simulation_coilAngle' + \
                       str(float(COIL_ANGLE)) \
                       if COIL_ANGLE!=0.0 else 'simnibs_simulation'
write_log('parameters', ['simnibs_session_name: '+simnibs_session_name])


def create_populations(pop_folders=None):
        
        if pop_folders is None:
                pop_folders = {}
                pop_folders['pop_L23'] = None
                pop_folders['pop_L5'] = None
        
        global pop_L23
        pop_L23 = create_population(simnibs_session_name,
                            'L23_PC_cADpyr229_1',
                            waveform_params=wf_params,
                            plot_folder=pop_folders['pop_L23']
                            )
        time.sleep(1.0) #debug
        global pop_L5
##        pop_L5 = create_population(simnibs_session_name,
##                            'L5_TTPC2_cADpyr232_1',
##                            waveform_params=wf_params,
##                            plot_folder=pop_folders['pop_L5']
##                            )
        pop_L5 = create_population(simnibs_session_name,
                            'L5_simplified_cell',
                            waveform_params=wf_params,
                            plot_folder=pop_folders['pop_L5']
                            )
        
        return [pop_L23, pop_L5]

# obs: there are 3 meanings of 'population' here.
# one is the Population object representing a single layer (L23 or L5);
# other is the sub-group of the Population object, representing neurons used
# in the bmtk network simulations;
# and other is the sub-population used in the bmtk network, like 'pop_L23_AMPA'
# and 'pop_L23_GABAa'.

# ------------------------------------------
# cell models functions

def create_bluebrain_cell(bionode, cell_name, load_mechanisms=False):
    """
    FOR DEBUGGING: Instead of creating the cell modified by Aberra et al. 2018,
    creates only the original Blue Brain Project cell.
    """
    
##    if load_mechanisms:
##        neuron.load_mechanisms(BBP_MECHS_FOLDER)
    
    os.chdir(ROOT_DIR + '/core/AberraEtAl2018_edit/cells')
    
    sys.path.append(os.getcwd())
    exec('import ' + cell_name + '.run as cell_model')
    
    os.chdir(cell_name)
    
    # load cell definitions
    try:
        neuron.h.xopen("morphology.hoc")
        neuron.h.xopen("biophysics.hoc")
        neuron.h.xopen("template.hoc")
    except:
        print 'morphology/biophysics/template not redefined.'
    
    cell = cell_model.create_cell(False)
    os.chdir('..')
    
    # set position
    x,y,z = bionode.position

    # assumes the default soma position is the origin: (0,0,0)
    curr_x, curr_y, curr_z = [0.0, 0.0, 0.0] 
    for sec in cell.all:
        for i in range(sec.n3d()):
            neuron.h.pt3dchange(i,
                    x - curr_x + sec.x3d(i),
                    y - curr_y + sec.y3d(i),
                    z - curr_z + sec.z3d(i),
                    sec.diam3d(i), sec=sec)
    
    return cell

def create_stimcell(bionode, template_name, dynamics_params=None):
    '''
    Function called internally by BMTK to create a cell of the StimCell class.
    '''
    
    # the cell dynamics parameters are already defined in the
    # modified Aberra el al. 2018 code.
    # for now, we don't need to override them.
    if dynamics_params is not None:
            print("StimCell: dynamics_params ( ", \
                  dynamics_params, " ) were passed; ignoring.")
    
    if USE_POP_POSITIONS: #debug
        # debug
        pop = None
        if 'L23' in template_name:
                pop = pop_L23
        elif 'L5' in template_name:
                pop = pop_L5
        else:
                raise "template_name not understood; it must be 'L23' or 'L5'."
                quit()
        
        pop_name = bionode._node.node_type_properties['pop_name']
        chosen_neuron_idx = bionode._node._group_props['neurons_indices']
        
        assert(np.all(np.array(bionode.position) == \
                      np.array(pop.soma_centers[chosen_neuron_idx]))) #debug
        
        s0 = StimCell.StimCell(cell_name=template_name,
                               population=pop,
                               neuron_idx=chosen_neuron_idx)
        
        s0.position_neuron_at_relative_depth()
        s0.align_to_cortical_surface()
        s0.random_azimuthal_rotation()

        # test: this replaces the lines marked with '#debug#test'.
        # undo the comments in those lines if you remove this line.
        s0._update_pt3d()

        try:
            assert(map(int, bionode.position) == map(int, s0.somapos)) #debug
        except AssertionError:
            s0.set_pos(x=bionode.position[0],
            y=bionode.position[1],
            z=bionode.position[2])

            s0._update_pt3d()
                
        if USE_TMS_STIM:
                # set TMS stimulation
                v_seg_values = s0.set_stimulation(stim_intensity=TMS_PARAMETERS['stim_intensity'])

    else:
        s0 = StimCell.StimCell(cell_name=template_name, population=None, neuron_idx=0)
        
        # set soma position
        s0.set_pos(x=bionode.position[0],
                y=bionode.position[1],
                z=bionode.position[2])
        
        s0._update_pt3d()
        
        chosen_neuron_idx = None
    
    print('StimCell created: ', template_name,
          ' | index: ', chosen_neuron_idx)
    
    # appends the cell to a global list so that all cells
    # can be accessed after the end of the simulation.
    bcl.append(s0)
    
    return s0.cell

# my function, written modifying 'Exp2Syn' from bmtk code.
def ASyN_STD(syn_params, x, sec):
    """
    Create a list of ASyN_STD synapses
    
    :param syn_params: parameters of a synapse
    :param xs: list of normalized distances along the section
    :param secs: target sections
    :return: list of NEURON synapse objects
    """    
    syn = neuron.h.ASyN_STD(x, sec=sec)
    
    # synapse parameters	
    syn.e = syn_params['erev']
    syn.tau1 = syn_params['tau1']
    syn.tau2 = syn_params['tau2']
    syn.gpeak = syn_params['gpeak']

    # short-term depression parameters
    syn.d1 = syn_params['d1']
    syn.tau_D1 = syn_params['tau_D1']
    
    write_log('synapses', ['ASyN_STD created at sec: ', sec.name()])
    write_log('synapses', ['syn_params: ', syn_params])

    return syn

def lognormal_synweights(edge_props, src_props, trg_props):

    mean = edge_props['syn_weight']

    syn_weight_params = ast.literal_eval(edge_props['syn_weight_params'])
    sigma = syn_weight_params['sigma']

    weight = np.random.lognormal(mean=mean, sigma=sigma)

    write_log('synapses', ['lognormal : synapse weight: ', weight])
    
    return weight

# ------------------------------------------

# debug
# creates random positions inside a sphere
sphere_positions = lambda x, center, radius: \
    np.array([linar.random_point_inside_sphere(center=center, radius=radius) \
              for i in range(x)])

# mine
def choose_positions(positions, preserve_original=False):
    """ Randomly choose lines from a matrix. Returns the chosen
        line and the matrix with the chosen line removed
        (if preserve_original is False). """
    
    chosen_pos_index = random.randint(0, positions.shape[0])
    pos = positions[chosen_pos_index, :]
    
    if not preserve_original:
        positions = np.delete(positions, chosen_pos_index, axis=0)
    
    return pos, positions

# mine
def delay_proportional_to_distance(source, target, delay_factor=1e-3):
    """
    given two bionodes (objects from bmtk library), calculates distance
    between them.
    """
    return delay_factor*linar.distance(source.position, target.position)

# from bio_450cells bmtk example.
def random_connections(source, target, p=0.1):
    sid = source['node_id']  # Get source id
    tid = target['node_id']  # Get target id
    
    # Avoid self-connections.
    if sid == tid:
        return None
    
    return np.random.binomial(1, p)  # nsyns

# from bio_450cells bmtk example.
def n_connections(src, trg, prob=0.1, min_syns=1, max_syns=5):
    """
    Referenced by add_edges() and called by build() for every source/target
    pair. For every given target/source pair will connect the two with
    a probability prob (excludes self-connections).
    """
    if src.node_id == trg.node_id:
        return 0
    
    return 0 if np.random.uniform() > prob else np.random.randint(min_syns, max_syns)

def load_nodes_positions(L23_num_neurons,
                         L5_num_neurons,
                         synapse_types,
                         pop_params_path=ROOT_DIR + \
                             '/core/data/simnibs_simulation/population_parameters/'):
    
    POP_PARAMS_FOLDER = os.path.abspath(pop_params_path)
    soma_centers = {}
    L23_pop_positions = {}
    for layer_name in ('L23','L5'):
        soma_centers_file = POP_PARAMS_FOLDER + '/' + layer_name + '_soma_centers.dat'
        soma_centers[layer_name] = pickle.load(open(soma_centers_file, 'rb'))['layer_soma_centers']
    
    # L23: 4 populations
    for syn_type in synapse_types:
        L23_pop_positions[syn_type] = []
        for i in range(L23_num_neurons[syn_type]):
            pos, soma_centers['L23'] = choose_positions(soma_centers['L23'])
            L23_pop_positions[syn_type].append(pos)
        L23_pop_positions[syn_type] = np.array(L23_pop_positions[syn_type])
    
    # L5: 1 population
    L5_pop_positions = []
    for i in range(L5_num_neurons):
        pos, soma_centers['L5'] = choose_positions(soma_centers['L5'])
        L5_pop_positions.append(pos)
    L5_pop_positions = np.array(L5_pop_positions)
    
    return [L23_pop_positions, L5_pop_positions]

def get_time_label():
    from socket import gethostname
    from datetime import datetime
    
    label = gethostname() + '_' \
        + str(datetime.now())[:-7].replace('-','').\
        replace(' ','').replace(':','')
    return label

def create_dir_file(network_folder_name, filename='default_base_dir.txt'):
    # save path of default base dir of the bionet
    # (so it can be used by the 'run bionet' script),
    # since its name changes on different runs.
    base_dir_file = open(filename, "w+")
    base_dir_file.write(network_folder_name)
    base_dir_file.close()

def clean_network_folders(folders_to_remove=['components',
                                             'network',
                                             'output'],
##                         files_to_remove=[],
                         recreate=False,
                         ask_confirmation=False):

    write_log('files', ['current folder: ', os.getcwd()])

    if ask_confirmation:
            raw_input('cleaning current folder. press ctrl+c to stop; \
                      enter to continue. ')

    for folder in folders_to_remove:
        try:
            shutil.rmtree(folder)
            if recreate:
                os.mkdir(folder)
        except:
            write_log('files', ['Folder not removed: ', folder])

    write_log('files', ['removed folders: ', folders_to_remove])

##    for f in files_to_remove:
##        try:
##            os.remove(f)
##        except:
##            write_log('synapses', ['File not removed: ', f])
##    write_log('files', ['removed files: ', files_to_remove, '\n'])


def copy_mechanisms(current_modfiles_folder,
                    mechs_folder,
                    compile_mechs=True):
    
    orig_folder = os.getcwd()

    os.chdir(current_modfiles_folder)
    for fn in os.listdir('.'):
        shutil.copy(fn, "..")
        write_log('files', ['Copied mechanism: ', fn])
    
    os.chdir('..')
    for fn in os.listdir(mechs_folder):
        if fn[-4:] == ".mod":
            shutil.copy(mechs_folder+'/'+fn, '.')
            write_log('files', ['Copied mechanism: ', fn])
    
    if compile_mechs:
        # compiles mechanisms which are
        # inside the 'components/mechanisms' folder
        os.system('nrnivmodl')
    
    os.chdir(orig_folder)

def copy_synapse_files(synapse_files):
    for synapse_file in synapse_files:
        shutil.copy(ROOT_DIR + '/synapse_models/' + synapse_file,
                    BMTK_DIR + '/components/synaptic_models')

def load_mechanisms():
    ''' Copy cell model (Aberra et al 2018) mechanisms
        to BioNet folder. (mod files with same name are
        overwritten, since the cell model depends on them.
        besides, mod files are copied to the parent folder,
        'mechanisms', and compiled there.) '''

    modfiles_folder = BMTK_DIR + "/components/mechanisms/modfiles/"
    copy_mechanisms(current_modfiles_folder=modfiles_folder,
                    mechs_folder=BBP_MECHS_FOLDER,
                    compile_mechs=COMPILE_MECHS)
    neuron.load_mechanisms(modfiles_folder)
