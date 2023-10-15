# -*- coding: utf-8 -*-
import os
import shutil
import itertools
import numpy as np
import pickle
import time
import sys

# import simulation parameters
from ..parameters import *
from ..preprocessing import *

import TMSSimulation
import NeuronPopulation

#debug
# for customizing path of results
##CURRENT_RESULTS_FOLDER = '/home/master/Desktop/TMSEffects/tms_isolated_neurons'
CURRENT_RESULTS_FOLDER = os.getcwd() + '/results' + '/'

# ---------------------------------------------------------------------

def runAllPopulations(msh_folders, simnibs_mat_fn, number_of_coil_orientations, simnibs_command):
    """ Run simulations for synaptically isolated neurons of multiple populations, which can vary in
        coil orientation (and therefore applied electric field), cortical layer, or cell model.
        The default option if for running all L23 and L5 populations for all 8 coil orientations. """
    
    # which steps will be run
    pipeline = {'generate_files': 0,
                'run_simnibs':0,
                'arranging':0,
                'preprocessing':0,
                'simulations':1}
    
    print pipeline
    
    # if you only want to run some of the sessions, then change this.
    msh_folders_mask = np.array([0,1,1,1, 1,1,1,1])
##    msh_folders_mask = np.ones(8)
    
    # this selects the desired sessions (set in msh_folders_mask).
    msh_folders = list(itertools.compress(msh_folders, msh_folders_mask))
    
    # -------------------
    if pipeline['generate_files']:
        
        print "===================== Generating .mat files: "
        angles = generate_mat_files(simnibs_mat_fn, number_of_coil_orientations, show_coil_basis=False)
    
    # -------------------
    if pipeline['run_simnibs']:
        print "===================== Running SimNIBS: "
        run_simnibs(simnibs_command, simnibs_mat_fn, angles)
    
    # -------------------
    if pipeline['arranging']:
        print "===================== Arranging: "
        
        
        if len(msh_folders) == 0:
            raise OSError("No folders found inside: " + os.getcwd() + MSH_INPUTS_FOLDER + " \
                            Please put one or more folder resulting from a SimNIBS simulation inside the aforementioned folder so they \
                          can be processed and used for the neuronal simulations.")
        
        # for each folder inside MSH_INPUTS_FOLDER:
        for folder in msh_folders:
            
            # get all .msh files inside this folder
            msh_fn = filter(lambda string: (string[-4:] == '.msh'), os.listdir(MSH_INPUTS_FOLDER + '/' + folder))
            
            if len(msh_fn) == 0:
                raise FileNotFoundError("No .msh files found inside the folder: " + folder)
            elif len(msh_fn) != 1:
                print "WARNING: there should be only one .msh file inside the folder inside " + MSH_INPUTS_FOLDER + " : " + folder
                print "Using the first .msh found instead."
            
            msh_fn = msh_fn[0][:-4]
            
            # msh_fn will be the name of the corresponding data folder
            session_data_folder = MAIN_DATA_FOLDER + '/' + folder
            if not os.path.isdir(session_data_folder):
                os.mkdir(session_data_folder)
                print "created folder: ", session_data_folder

            # create folders for data generated during preprocessing
            data_subfolders = [session_data_folder + '/fields',
                               session_data_folder + '/meshes',
                               session_data_folder + '/population_parameters']
            for dsfolder in data_subfolders:
                if not os.path.isdir(dsfolder):
                    os.mkdir(dsfolder)
                    print "created subfolder: ", dsfolder
            
            msh_path = MSH_INPUTS_FOLDER+'/'+folder+'/'+msh_fn+'.msh'
    
    # -------------------
    if pipeline['preprocessing']:
        print "===================== Preprocessing:"
        # run 'preprocessing.py'
        import preprocessing
        
        for simnibs_session_name in msh_folders:
            # this function will read the .msh file inside the 'input' folder and
            # write, inside data_folder, files to be used in the simulations.
            preprocessing.preprocessing(simnibs_session_name, relative_depths)
    
    # -------------------
    if pipeline['simulations']:
        print "===================== Simulations will start."
        
        print "NEURONS_TO_SIMULATE: ", NEURONS_TO_SIMULATE
        neurons_to_simulate = NEURONS_TO_SIMULATE
        
        for simnibs_session_name in msh_folders:
            for population_cell_name in ['L23_PC_cADpyr229_1']:#, 'L23_PC_cADpyr229_2', 'L23_PC_cADpyr229_3']: #['L23_PC_cADpyr229_1', 'L5_TTPC2_cADpyr232_1']:
                runPopulation(simnibs_session_name, population_cell_name, neurons_to_simulate)
                # TODO: make NeuronPopulation save string indicating data_folder
    
    #...
    # save file indicating coil orientation [and coil position matrix] inside the results folder.
    # (each population simulation will have one folder inside the 'results' folder.)
    #...

# ---------------------------------------------------------------------

def runPopulation(simnibs_session_name, population_cell_name, neurons_to_simulate=None):
    """ Run simulations for synaptically isolated neurons of a single population (i.e., of a single cortical layer).

    neurons_to_simulate : list of indices of neurons of the population to be simulated.
                          If None, all neurons of the population will be simulated. """
    
    # ---- load population parameters
    start_time = time.time()
    
    session_data_folder = ROOT_FOLDER + "/core/data/" + simnibs_session_name + "/"
    population_cell_type = population_cell_name.split('_')[0]
    
    # load files
    with open(session_data_folder + "population_parameters/" + population_cell_type + "_soma_centers.dat", "r") as fileobj:
        layer_coords_dict = pickle.load(fileobj)
    layer_center_coords = layer_coords_dict['layer_soma_centers']
    layer_neurons_distance = layer_coords_dict['layer_neuron_distances_factor']
    with open(session_data_folder + "population_parameters/" + population_cell_type + "_normals.dat", "r") as fileobj:
        gm_normals = pickle.load(fileobj)
    with open(session_data_folder + "population_parameters/" + population_cell_type + "_gm_centers.dat", "r") as fileobj:
        gm_centers = pickle.load(fileobj)
    with open(session_data_folder + "population_parameters/" + population_cell_type + "_neurons_indices.dat", "rb") as fileobj:
        neurons_polygon_indices = pickle.load(fileobj)
    
    print "Cells positions and orientations loaded."
    
    # --- tmssimulation object
    # configure the magnetic stimulation of the population
    tms_sim = TMSSimulation.TMSSimulation(simnibs_session_name)
    
    # --- waveform
    # create RLC waveform
    tms_sim.generate_rlc_waveform(dt=dt, tstop=tstop,
                                  DEL=wf_delay, DUR=wf_dur, AMP=wf_amp,
                                  Rc=R, Lc=L, Cc=C)
    
    #debug
    PLOT_ALL = False
    if PLOT_ALL:
        show_waveform(tms_sim, tstop, dt)
    
    # ------------------------------
    # population
    print "=========== Starting population simulation..."
    print "- Running population in layer: ", population_cell_type
    print "- With neurons at this relative cortical depth: ", layer_neurons_distance
    
    # create population, and run the simulations (with neurons under stimulation)
    population = NeuronPopulation.Population(cell_name=population_cell_name,
                                                    soma_centers=layer_center_coords,
                                                    gm_normals=gm_normals,
                                                    gm_centers=gm_centers,
						    neurons_polygon_indices=neurons_polygon_indices,
                                                    tms_sim=tms_sim,
                                                    v_init=v_init,
                                                    tstop=tstop,
                                                    dt=dt,
                                                    simnibs_session_name=simnibs_session_name,
                                                    results_folder=CURRENT_RESULTS_FOLDER) #debug
    
    # in case something went wrong
    if population is None:
        return
    
    # save information for this population simulation
    with open(population.plots_folder_path + "/simulation_info.dat", "w+") as dump_file:
        pickle.dump({'simnibs_session_name':simnibs_session_name,
                     'population_cell_name':population_cell_name,
                     'layer_neurons_distance':layer_neurons_distance,
                     'thresholds_search_range':thresholds_search_range,
                     'waveform':tms_sim.stim_time_course,
                     'tstop':tstop,
                     'dt':dt,
                     'DEFAULT_E_VECTOR':DEFAULT_E_VECTOR,
                     'neurons_to_simulate':neurons_to_simulate},
                    dump_file)
    
    # run simulations
    layer_thresholds = population.simulate_cells(thresholds_search_range, run_parallel, neurons_to_simulate)
    
    # ------------------------------ 
    end_time = time.time()
    print 'total simulation time (seconds): ', end_time - start_time

    if len(layer_thresholds) > 0:
        print "\nPopulation simulation finished.\n"
    else:
        print "\nWARNING: simulate_cells did not finish as expected: the length of 'thresholds' is zero. \n"

    print "====== " + population_cell_name + "_thresholds: "
    print layer_thresholds
    
    # ------------------------------
    # save the final results (activation thresholds) for this population
    os.chdir(population.plots_folder_path)
    with open(population_cell_type + "_thresholds.dat", "w+") as dump_file:
        pickle.dump(layer_thresholds, dump_file)
    print "\n" + population_cell_name + "_thresholds saved at " + os.getcwd() + " ."
    os.chdir('../..')

# ---------------------------------------------------------------------

if __name__ == '__main__':
    if sys.argv[-1] == 'all':
        print "Running all populations. "
        
        # these arguments are defined inside 'parameters.py'.
        runAllPopulations(msh_folders, simnibs_mat_fn, number_of_coil_orientations, simnibs_command)
        
    else:
        simnibs_session_name = 'simnibs_simulation_coilAngle45.0' #'simnibs_simulation_coilAngle270.0' #os.listdir(ROOT_FOLDER+'/core/data/')[3]
        population_cell_name = 'L23_PC_cADpyr229_1' #'L5_TTPC2_cADpyr232_1' #'L23_PC_cADpyr229_2' #'L1_NGC-DA_bNAC219_1'
        neurons_to_simulate = NEURONS_TO_SIMULATE
        
        print "· Session: ", simnibs_session_name
        print "· Cell type: ", population_cell_name
        print "\n"
        
        runPopulation(simnibs_session_name, population_cell_name, neurons_to_simulate)
