# this file of the bmtk library was modified to allow custom
# netcons (whose source section is not the soma, but the axon terminals):
# /home/master/Desktop/bmtk-0.0.8-py2.7.egg/bmtk/simulator/bionet/biocell.py

# files
import sys
import os
import shutil
import pickle
import time
import datetime

# math
import math
import random
import numpy as np

# parallel
from mpi4py import MPI

# neuron
import neuron

# bmtk
from bmtk.builder.networks import NetworkBuilder
from bmtk.simulator import bionet
from bmtk.utils.sim_setup import build_env_bionet
from bmtk.utils.reports.spike_trains import plotting
from bmtk.analyzer.spike_trains import load_spikes_file
from bmtk.simulator.bionet.pyfunction_cache import *
from bmtk.builder.auxi.edge_connectors import distance_connector
from bmtk.builder.auxi.edge_connectors import connect_random
from bmtk.utils.reports.spike_trains import PoissonSpikeGenerator
from bmtk.utils.sim_setup import build_env_bionet

# tms_networks
from tms_networks.functions import *
from tms_networks.parameters import *

from tms_networks.debugging import *


# ------------------------ additional models
# add cell model directives

add_cell_model(func=create_stimcell,
               directive='StimCell.py',
               model_type='biophysical',
               overwrite=True)
add_synapse_model(func=ASyN_STD,
                  name='ASyN_STD',
                  overwrite=True)
add_weight_function(lognormal_synweights,
                    'lognormal_synweights',
                    overwrite=False)


# ----------------------------------------------------------------------------
def create_nodes(network, populations_names, network_neurons_count,
                 network_neurons_positions, network_neurons_indices):
    '''
    Creates excitatory (AMPA/NMDA-projecting) and inhibitory
    (GABA-A/GABA-B-projecting) neurons for TMS-activated network.
    '''
    
    L23_morphology_file = BBP_MODELS_FOLDER \
        + '/L23_PC_cADpyr229_1/morphology/' \
        + 'dend-C170897A-P3_axon-C260897C-P4_-_Clone_4.asc'
    L5_morphology_file = ROOT_DIR \
        + '/core/L5_simplified_cell/morphology/' \
        + 'dend-C060114A7_axon-C060116A3_-_Clone_2.asc'
    
    # L2/3 neurons
    for pop_name in populations_names:
        if 'L23' in pop_name:
            
            if 'GABA' in pop_name:
                ei = 'i'
            else:
                ei = 'e'
            
            network.add_nodes(N=network_neurons_count[pop_name],
                              pop_name=pop_name,
                              ei=ei,
                              positions=network_neurons_positions[pop_name],
                              neurons_indices=network_neurons_indices[pop_name], # custom parameter
                              model_type='biophysical',
                              model_template='StimCell.py:L23_PC_cADpyr229_1',
                              morphology=L23_morphology_file)
    
    # L5 neurons
    network.add_nodes(N=network_neurons_count['pop_L5'],
                     pop_name='pop_L5',
                     positions=network_neurons_positions['pop_L5'],
                     # custom parameter
                     neurons_indices=network_neurons_indices['pop_L5'],
                     model_type='biophysical',
                     model_template='StimCell.py:L5_simplified_cell',
                     morphology=L5_morphology_file)


def create_edges(network, nsyn_scale, weights,
                 weight_sigma, synapse_template):
    
    # todo: put all add_edges into a single loop.
    network.add_edges(source={'pop_name': 'pop_L23'},
                      target={'pop_name': 'pop_L5'},
                     connection_rule=nsyn_scale[0],
##                     connection_rule=distance_connector,
##                     connection_params={'d_weight_min': 0.0,
##                                        'd_weight_max': weight_0,
##                                        'd_max': 50.0*1e3, #test
##                                        'nsyn_min': 0,
##                                        'nsyn_max': int(0.6*nsyn_scale)}, #test
                     
                     syn_weight=weights[0],  # used as mean by lognormal_synweights
                     syn_weight_params={'sigma': weight_sigma},
                     weight_function='lognormal_synweights',
                     
##                     distance_range=[30.0, 150.0],
                     distance_range=[0.0, 1e20],
##                     target_sections=['basal', 'apical', 'soma'],
                     target_sections=['dend', 'apic'],

                     # delay proportional to distance:
##                     delay=delay_proportional_to_distance,
##                     delay_params={'delay_factor': delay_factor},
                     
                     # delay constant:
                     delay=1.0,
                     
                     dynamics_params='AMPA_ExcToExc_modified.json',
                     model_template=synapse_template) 
    
    network.add_edges(source={'pop_name': 'pop_L23_GABAa'},
                     target={'pop_name': 'pop_L5'},
                     connection_rule=nsyn_scale[1],
##                     connection_rule=distance_connector,
##                     connection_params={'d_weight_min': 0.0,
##                                        'd_weight_max': weight_1,
##                                        'd_max': 50.0*1e3, #test
##                                        'nsyn_min': 0,
##                                        'nsyn_max': int(0.2*nsyn_scale)},
                     syn_weight=weights[1],  # used as mean by lognormal_synweights
                     syn_weight_params={'sigma':weight_sigma},
                     weight_function='lognormal_synweights',
                     
##                     distance_range=[30.0, 150.0],
                     distance_range=[0.0, 1e20],
##                     target_sections=['basal', 'apical', 'soma'],
                     target_sections=['dend', 'apic'],
                     
##                     delay=delay_proportional_to_distance,
##                     delay_params={'delay_factor':delay_factor},
                     delay=1.0,
                     
                     dynamics_params='GABA_InhToExc_modified.json',
                     model_template=synapse_template)
    
    network.add_edges(source={'pop_name': 'pop_L23_NMDA'},
                     target={'pop_name': 'pop_L5'},
                     connection_rule=nsyn_scale[2],
##                     connection_rule=distance_connector,
##                     connection_params={'d_weight_min': 0.0,
##                                        'd_weight_max': weight_2,
##                                        'd_max': 50.0*1e3, #test
##                                        'nsyn_min': 0,
##                                        'nsyn_max': int(0.1*nsyn_scale)},
                     
                     syn_weight=weights[2], # used as mean by lognormal_synweights
                     syn_weight_params={'sigma':weight_sigma},
                     weight_function='lognormal_synweights',
                     
##                     distance_range=[30.0, 150.0],
                     distance_range=[0.0, 1e20],
##                     target_sections=['basal', 'apical', 'soma'],
                     target_sections=['dend', 'apic'],

##                     delay=delay_proportional_to_distance,
##                     delay_params={'delay_factor':delay_factor},
                     delay=1.0,
                     
                     dynamics_params='NMDA_ExcToExc.json',
                     model_template=synapse_template)
    
    network.add_edges(source={'pop_name': 'pop_L23_GABAb'},
                     target={'pop_name': 'pop_L5'},
                     connection_rule=nsyn_scale[3],
##                     connection_rule=distance_connector,
##                     connection_params={'d_weight_min': 0.0,
##                                        'd_weight_max': weight_3,
##                                        'd_max': 50.0*1e3, #test
##                                        'nsyn_min': 0,
##                                        'nsyn_max': int(0.05*nsyn_scale)},
                     
                     syn_weight=weights[3], # used as mean by lognormal_synweights
                     syn_weight_params={'sigma': weight_sigma},
                     weight_function='lognormal_synweights',
                     
##                     distance_range=[30.0, 150.0],
                     distance_range=[0.0, 1e20],
##                     target_sections=['basal', 'apical', 'soma'],
                     target_sections=['dend', 'apic'],
                    
##                     delay=delay_proportional_to_distance,
##                     delay_params={'delay_factor':delay_factor},
                     delay=1.0,
                     
                     dynamics_params='GABAb_InhToExc.json',
                     model_template=synapse_template)


def create_neuron_positions(populations_names,
                            layer_populations
                            ):
    ''' populations_names : list of strings naming all populations, grouped
        by synapse types: AMPA, GABA-A, etc.
        layer_populations : dictionary of population objects, grouped
        by layer: L23, L5.
    '''
    
    nnpos_path = ROOT_DIR + '/simulation_temp/network_neurons_positions.dat'
    nnind_path = ROOT_DIR + '/simulation_temp/network_neurons_indices.dat'
    
    # first case: loads the already-generated positions of the neurons
    if USE_DUMPED_POSITIONS:
        with open(nnpos_path, 'rb') as nnpos_file:
            network_neurons_positions = pickle.load(nnpos_file)
        with open(nnind_path, 'rb') as nnind_file:
            network_neurons_indices = pickle.load(nnind_file)
    else:
    # second case: generates the positions of the neurons
        network_neurons_positions = {}
        
        # indices of neurons at the Population object
        network_neurons_indices = {}
        
        # randomly selects neurons from the Population object
        # and gets their positions.
        for pop_name in populations_names:
                N = network_neurons_count[pop_name]
                network_neurons_positions[pop_name] = []
                network_neurons_indices[pop_name] = []
                
                if 'L23' in pop_name:
                        pop = layer_populations['L23']
                elif 'L5' in pop_name:
                        pop = layer_populations['L5']
                
                for i in range(N):
                        # problem: there is a low chance of two neurons being
                        # of the same index.
                        chosen_neuron_idx = random.randint(0, pop._N-1)
                        network_neurons_indices[pop_name].\
                            append(chosen_neuron_idx)
                        network_neurons_positions[pop_name].\
                            append(pop.soma_centers[chosen_neuron_idx])
                
                network_neurons_positions[pop_name] = np.array(
                    network_neurons_positions[pop_name])
                network_neurons_indices[pop_name] = np.array(
                    network_neurons_indices[pop_name])
        
        write_log('positions', ['generated positions: ' \
            + str(network_neurons_positions)])
        
        with open(nnpos_path, 'w+') as nnpos_file:
            pickle.dump(network_neurons_positions, nnpos_file)
        with open(nnind_path, 'w+') as nnind_file:
            pickle.dump(network_neurons_indices, nnind_file)
        print('written files in simulation_temp: ' \
              + nnpos_path + ' & ' + nnind_path)
        
    return {'positions': network_neurons_positions,
            'indices': network_neurons_indices}


def create_network_populations():
    
    # if an interpolated E-field file already exists,
    # then it is loaded and no interpolation is remade.
    if os.path.isfile(POP_FOLDERS_PATH):
        pop_folders_detected = True
        with open(POP_FOLDERS_PATH, 'rb') as pop_folders_file:
            pop_folders = pickle.load(pop_folders_file)
    else:
        pop_folders_detected = False
        pop_folders = None
    
    # E-field interpolation is done here if it was not loaded.
    pop_L23, pop_L5 = create_populations(pop_folders=pop_folders)
    pop_folders = {'pop_L23': pop_L23.plots_folder_path,
                   'pop_L5': pop_L5.plots_folder_path}
    
    # saves locations of folders of the populations
    # if it does not already exist.
    if pop_folders_detected is False:
        with open(POP_FOLDERS_PATH, 'w+') as pop_folders_file:
            pickle.dump(pop_folders, pop_folders_file)
            print('file written: ' + POP_FOLDERS_PATH)
    
    return {'L23': pop_L23, 'L5': pop_L5}


def create_psgs(thalamus_net_names, N_thalamus, thalamus_target_pops,
                connect_random, nsynmax_thalamus, synw_thalamus,
                maxrang_thalamus):
    '''
    Creates Poisson spikes generator which represents
    background activity arriving to the main network.
    '''
    
    thalamus_networks = {}
    psgs = {}
    
    for thalamus_net_name in thalamus_net_names:
        thalamus_networks[thalamus_net_name] = \
                                             NetworkBuilder(thalamus_net_name)
        thalamus = thalamus_networks[thalamus_net_name]
        
        thalamus.add_nodes(N=N_thalamus,
                           pop_name=thalamus_net_name,
                           potential='exc',
                           model_type='virtual')
        
        for target_pop in thalamus_target_pops[thalamus_net_name]:
            thalamus.add_edges(source=thalamus.nodes(),
                               target=cortex.nodes(pop_name=target_pop),
                               connection_rule=connect_random,
                               connection_params={'nsyn_min': 0,
                                   'nsyn_max': nsynmax_thalamus},
                               syn_weight=synw_thalamus, # 1.0
                               distance_range=[0.0, maxrang_thalamus], #[0.0, 1.5*1e5], #test
    ##                           target_sections=['basal', 'apical'],
                               target_sections=['dend', 'apic'],
                               delay=0.0,
                               dynamics_params='InputSynapses_ExcToExc.json',
                               model_template='exp2syn')
        
        write_log('connections', ['connections were set: '
                                  + thalamus_net_name \
                                  + " -> " + target_pop])
        
        thalamus.build()
        thalamus.save_nodes(output_dir=BMTK_DIR+'/network')
        thalamus.save_edges(output_dir=BMTK_DIR+'/network')
        
        # ------ psg's
        psgs[thalamus_net_name] = PoissonSpikeGenerator(\
            population=thalamus_net_name)
        psg = psgs[thalamus_net_name]
        psg.add(node_ids=range(N_thalamus),
                firing_rate=firing_rates[thalamus_net_name],
                times=times)
        psg.to_sonata(BMTK_DIR+'/inputs/'+thalamus_net_name+'_spikes.h5')
        
    return {'networks': thalamus_networks, 'psgs': psgs}

# ----------------------------------------------------------------------------
# functions for running a single network

def build():
    clean_network_folders(folders_to_remove=[BMTK_DIR])
    clean_network_folders(folders_to_remove=[
        ROOT_DIR+'/simulation_temp/connections',
        ROOT_DIR+'/simulation_temp/logs',
        ROOT_DIR+'/simulation_temp/morphologies'],
        recreate=True)

    layer_populations = create_network_populations()

    nndict = create_neuron_positions(populations_names,
                                     layer_populations
                                     )
    network_neurons_positions = nndict['positions']
    network_neurons_indices = nndict['indices']
    
    # --- configures cortex network
    nsyn_scale = [nsyn_scale_0, nsyn_scale_1, nsyn_scale_2, nsyn_scale_3]
    weights = [weight_0, weight_1, weight_2, weight_3]
    
    # todo: put all 'add_nodes' in the same loop.
    cortex = NetworkBuilder('mcortex')
    create_nodes(cortex, populations_names, network_neurons_count,
                 network_neurons_positions, network_neurons_indices)
    create_edges(cortex, nsyn_scale, weights,
                 weight_sigma, synapse_template)
    cortex.build()
    cortex.save_nodes(output_dir=BMTK_DIR+'/network')
    cortex.save_edges(output_dir=BMTK_DIR+'/network')
    
    print_network_info(cortex)
    
    # --- background activity
    if USE_PSG:
        create_psgs()
        spikes_inputs = []
        for thalamus_net_name in thalamus_net_names:
            spikes_inputs.append((thalamus_net_name,
                                  BMTK_DIR+'/inputs/'
                                    +thalamus_net_name+'_spikes.h5'))
    else:
        spikes_inputs = None
    
    # --- builds networks
    build_env_bionet(base_dir=BMTK_DIR,      
                     network_dir=BMTK_DIR+'/network',
                     tstop=sim_tstop,
                     dt=sim_dt,
                     #report_vars=['v'], # also 'cai'
                     spikes_inputs=spikes_inputs,
                     include_examples=True, # Copies components files
                     compile_mechanisms=False
                    )
    
    copy_synapse_files(SYNAPSE_FILES)
    load_mechanisms()


def run():
    '''
    Run a single simulation of the network.
    This function must be called with the following command, where you must
    replace 'NUM_CORES' with the number of desired parallel nodes:
        mpirun -n NUM_CORES nrniv -mpi -python net.py run
    
    '''
    
    __name__ = 'tms_networks.net'
    
    with open(POP_FOLDERS_PATH, 'rb') as pop_folders_file:
        pop_folders = pickle.load(pop_folders_file)
    
    t_start = time.time()
    
    # it is easier to recreate the population objects than deal with its
    # unpickable attributes. the interpolated E-field is loaded here from
    # the file created inside 'build()'.
    pop_L23, pop_L5 = create_populations(pop_folders=pop_folders)
    
    conf = bionet.Config.from_json(BMTK_DIR + '/simulation_config.json')
    conf.build_env()
    net = bionet.BioNetwork.from_config(conf)
    
    # 'create_stimcell' is called inside this 'from_config'.
    sim = bionet.BioSimulator.from_config(conf, network=net)
    set_recorders(bcl)
    sim.run()

    # dump results
    write_log('time', ['Total time', t_start, time.time()])
    dump_connections(net)
    dump_cells_vmem(bcl, pop_L23.tms_sim.stim_time_course)
    dump_morphologies(bcl)

    MPI.COMM_WORLD.barrier()
    bionet.nrn.quit_execution()

# ----------------------------------------------------------------------------
# functions for running a pool of networks


def run_pool(N_simulations=12):
    
    write_log('pool_log',
              'running pool with {} networks. '.format(N_simulations))
    
##    if os.path.isdir(POOL_DIR):
##        shutil.rmtree(POOL_DIR)
    os.mkdir(POOL_DIR)
    write_log('pool_log',
              'created folder: {}'.format(POOL_DIR))
    
    pool_start_time = time.time() #debug
    
    # run pool
    for sim_id in range(N_simulations):
        # create folder for the current simulation. 
        sim_folder = os.path.abspath(POOL_DIR + '/' + str(sim_id))
        os.mkdir(os.path.abspath(sim_folder))

        # build bmtk network
        os.chdir(ROOT_DIR+'/..')
        os.system('./run_tms_network.sh')

        # save raster plot
        config_file = BMTK_DIR + '/simulation_config.json'
        spikes_file = None

        try:
            spike_trains = load_spikes_file(config_file=config_file,
                                            spikes_file=spikes_file)
            st = spike_trains

            # calculate mean square sum for the resulting raster plot.
            st_times = list(spike_trains.to_dataframe()['timestamps'])

            # raster plot filename
            raster_plot_fn = sim_folder + '/raster_' + str(sim_id) + '.pdf'
            plotting.plot_raster(spike_trains,
                                 show_plot=False,
                                 save_as=raster_plot_fn)

            shutil.copy(ROOT_DIR+'/vmem_results/L5_vmem.dat', sim_folder)
            shutil.copy(ROOT_DIR+'/parameters.py',
                        sim_folder)

            # generates vmem file which will be used for simulating
            # the epidural response.
            vmem = pickle.load(open(ROOT_DIR+'/vmem_results/L5_vmem.dat',
                                    'r'))
            vmem_axon = vmem[L5_AXON_SEG_INDEX, :]
            pickle.dump(vmem_axon, open(
                POOL_DIR+'/'+str(sim_id)+'/L5_vmem_axon.dat', 'w+'))
        except:
            write_log('pool_log',
            'could not load spikes in simulation: {}'.format(sim_id))

        os.chdir(TMS_NETWORK_DIR)

    pool_end_time = time.time() #debug    

    write_log('pool_log', 'run pool: total time (h:min:s): {}'.\
        format(str(datetime.timedelta(seconds=pool_end_time-pool_start_time))))
    print('finished; ran {} simulations. '.format(N_simulations))

# ----------------------------------------------------------------------------
# main


if __name__ == '__main__':
    if sys.argv[-1] == 'build':
        print('building network: ')
        build()
    elif sys.argv[-1] == 'run':
        print('running network: ')
        run()
    elif sys.argv[-1] == 'pool':
        print('running pool: ')
        run_pool(N_simulations=12)
    else:
        print('building and running network: ')
        build()
        run()
