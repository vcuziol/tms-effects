# -*- coding: utf-8 -*-

import numpy as np
import sys
import pickle
import os
import traceback
import time
from socket import gethostname
from datetime import datetime

from neuron import h

from ..core.StimCell import StimCell
from ..core.codes.list_routines import *
from ..core.codes import save_data as save_data
from ..core.codes.parallel_functions import *
from ..parameters import *

class Population:
    """ A collection of synaptically isolated neurons, positioned
        at a distance from the centers of the polygons
        (i.e., mesh elements) of a polygon mesh representing a
        gray matter region.
        
       soma_centers : list of x,y,z coordinates which represent
       the position where each soma center will be.
    """
    
    def __init__(self, cell_name,
                 soma_centers, gm_normals, gm_centers,
                 neurons_polygon_indices,
                 tms_sim, v_init, tstop, dt,
                 simnibs_session_name, results_folder,
                 plot_folder=None):
        """
        :param cell_name: Name of the folder of the cell model. string.
        """
        self.cell_name = cell_name
        self.simnibs_session_name = simnibs_session_name
        self.data_folder = ROOT_DIR + '/core/data/' + simnibs_session_name
        self.msh_fn = self.find_msh_file(MSH_INPUTS_FOLDER + '/' + \
            simnibs_session_name)
        
        # list of 3D coordinates of the centers of the neurons' somas.
        self.soma_centers = np.array(soma_centers)
        
        try:
            assert(len(gm_normals) == len(soma_centers) == len(gm_centers))
        except:
            print("POPULATION SIM ABORTED: %s | %s" % (cell_name, \
                traceback.format_exc()))
            return None
        
        # list of normal vectors and centers for each corresponding mesh
        # polygon/triangle of the gray matter motor cortex mesh.
        self.gm_normals = np.array(gm_normals)
        self.gm_centers = np.array(gm_centers)
        
        # indices of GM mesh polygons corresponding to each neuron:
        # since the polygons of the 3D mesh of the gray matter were used
        #  for positioning and rotating the neurons, and some of them were
        # maybe eliminated because of the algorithm for positioning according
        # to cortical depth,  we need to store the ids of the polygons of
        # the gray matter mesh that correspond to the current population
        # neurons' to be able to find the electric field vector for each neuron.
        self.neurons_polygon_indices = np.array(neurons_polygon_indices)
        
        self.v_init = v_init
        self.tstop = tstop
        self.dt = dt
        
        self.tms_sim = tms_sim
        
        N = len(soma_centers)
        self._N = N

        # create folders for saving plots
        if not os.path.isdir(results_folder):
            os.mkdir(results_folder)
        
        if plot_folder is None:
            self.plots_folder_path = results_folder + \
                gethostname() + '_' + str(hex(id(self))) + \
                '_' + str(datetime.now())[:-7].replace('-','').\
                replace(' ','').replace(':','')
            os.mkdir(self.plots_folder_path)
        
        else:
            # if 'plot_folder' was passed, subfolders
            # are assumed to already exist.
            self.plots_folder_path = plot_folder
        
        self.neuron_3d_data_path = self.plots_folder_path + \
                                   '/neurons_3d_data/'
        self.neuron_Efield_vecs_path = self.plots_folder_path + \
                                       '/neurons_Efield_vecs/'
        
        if plot_folder is None:
            os.mkdir(self.neuron_3d_data_path)
            print "created folder: ", self.neuron_3d_data_path
            os.mkdir(self.neuron_Efield_vecs_path)
            print "created folder: ", self.neuron_Efield_vecs_path

        csv_file_path = self.neuron_Efield_vecs_path + 'soma_centers_mm.csv'
        interpolated_Efield_file_path = csv_file_path[:-4] + '_E.csv'
        
        # if plot_folder is not None, load interpolated E field from it.
        # if it is, then generate the E field.
        if plot_folder is not None:
            self.interpolated_E_field = np.loadtxt( \
                interpolated_Efield_file_path, delimiter=',')
            print("loaded interpolated_E_field: " + \
                  interpolated_Efield_file_path)
        else:
            # interpolate electric field vectors at center of each soma.
            # if INTERPOLATE_AT_SEGMENTS is True, then interpolation is
            # made inside 'stimcell.set_E_field'.
            if (not INTERPOLATE_AT_SEGMENTS):
                print 'Starting interpolation: ', self.cell_name, ' \t ', \
                      self.simnibs_session_name
                print('''E-field will be interpolated at each soma center, 
                    and then one single E-field vector will be used for all 
                    segments of each neuron.''')
                
                # '1e-3' converts from micrometers to millimeters.
                np.savetxt(csv_file_path, self.soma_centers*1e-3, delimiter=',')
                
                msh_file_path = MSH_INPUTS_FOLDER + '/' + self.simnibs_session_name + '/' + self.msh_fn
                command_string = get_fields_command + ' -s ' + csv_file_path + ' -m ' + msh_file_path
                print 'running get_fields_at_coordinates as: ', command_string
                os.system(command_string)
                
                # load new interpolated E-field file.
                self.interpolated_E_field = np.loadtxt(interpolated_Efield_file_path, delimiter=',')
                
                if VERBOSE:
                    print('loaded interpolated_E_field: ' + interpolated_Efield_file_path)
                    import core.codes.linear_alg_routines as linar
                    print 'max normE across all neurons, at base TMS intensity: ', np.max(linar.magnitudes(self.interpolated_E_field))
                
                #debug
    ##            try:
    ##                assert(len(soma_centers) == len(self.interpolated_E_field)) #debug
    ##                # max normE on motor cortex at base TMS intensity (1.0 A/micros.) is expected not to be too low.
    ##                assert(np.max(linar.magnitudes(self.interpolated_E_field)) > 1.7) #debug
    ##            except AssertionError:
    ##                import pdb; pdb.set_trace()  

    # ----------------------------------------------------------------------------
    def find_msh_file(self, msh_folder):
        
        # get all .msh files inside this folder
        msh_fn = filter(lambda string: (string[-4:] == '.msh'), os.listdir(msh_folder))
        
        if len(msh_fn) == 0:
            raise OSError("No .msh files found inside the folder: " + msh_folder)
        elif len(msh_fn) != 1:
            print "WARNING: there should be only one .msh file inside the folder inside: " + msh_folder
            print "Using the first .msh found instead."
        
        return msh_fn[0]
    
    # ------------------------------------------------------
    def print_neuron_header(self, neuron_idx, thresholds_search_range):
        print "\n"
        print "███████████████████████████████████████████████████████████████████████████"
        print "█████ ----> Search started for neuron ", str(neuron_idx), " --- range: ", thresholds_search_range, "██████"
        print "███████████████████████████████████████████████████████████████████████████"
    
    # ------------------------------------------------------
    def neuron_threshold_search_wrapped(self, stimcell, thresholds_search_range, n_iterations=10):
        try:
            result = self.neuron_threshold_search(stimcell, thresholds_search_range, n_iterations=10)
            return result
        except:
            print("SEARCH ABORTED: %s | %s" % (stimcell.neuron_idx, traceback.format_exc()))
            return None
    
    # ------------------------------------------------------
    def neuron_threshold_search(self, stimcell, thresholds_search_range, n_iterations=10):
        
        self.print_neuron_header(stimcell.neuron_idx, thresholds_search_range)
        
        low_end, high_end = thresholds_search_range
        
        for iter_idx in range(n_iterations):
            print "\n:::::::: iteration ", iter_idx, " started."
            
            middle = (low_end + high_end)/2.0
            
            # the stimulation intensity unit is assumed to be the default unit in SimNIBS, which is 1 A/μs = 10^6 A/s.
            # the default intensity in SimNIBS is 1.0 * 10^6 A/s.
            stim_intensity = middle
            
            # configuring stimulation:
            # "v_segments" are the quasipotentials at each segment.
            # With no arguments, "set_stimulation" uses the time course and E-field vector
            # stored in the Population object for this specific neuron.
            v_segments = stimcell.set_stimulation(stim_intensity=stim_intensity)
            
            # --------------------------------------------
            # simulate
            was_activated = stimcell.simulate_neuron()
            
            if was_activated:
                # go left
                high_end = middle
            else:
                # go right
                low_end = middle
            
            print "__θ__ current threshold: ", middle, " A/μs"
            
            print "stimcell vmem nanpercentage: ", nanpercentage(stimcell.vmem)
            if np.isnan(np.sum(stimcell.vmem)):
                print "***** stimcell.vmem has NaNs."
                raise ValueError("SEARCH ABORTED: stimcell.vmem has NaNs.")
            
            # --------------------------------------------
            # get information about the beginning site of the action potential for this neuron
            AP_beginning_dict = stimcell.detect_AP_beginning_site(stimcell.vmem, h.dt)
            AP_beginning_dict['stim_intensity'] = stim_intensity
            
            print 'AP_beginning_seg_idx: ', AP_beginning_dict['AP_beginning_seg_idx']
            
            if (AP_beginning_dict['AP_beginning_seg_idx'] is not None) and (len(AP_beginning_dict['AP_beginning_seg_idx']) != 0):
                AP_beginning_dict['AP_beginning_seg_name'] = map(stimcell.get_idx_name, AP_beginning_dict['AP_beginning_seg_idx'])
                
                print("++++ AP beginning segment: ")
                for seg_name, seg_i in zip(AP_beginning_dict['AP_beginning_seg_name'], range(len(AP_beginning_dict['AP_beginning_seg_name']))):
                    print "  ", seg_i, ") ", seg_name
                
    	        # save terminals data only if there was an action potential.
                if DUMP_AP_SITE_MEASURES:
                    save_data.save_terminals_info(stimcell, AP_beginning_dict)
            else:
                print("---- AP beginning segment: no APs in this iteration.")
            
            # -------------
            # saves files for only the first iteration.
            if iter_idx == 0:
                save_data.save_neuron_morphology(stimcell, DUMP_NEURON_3D_DATA)
                save_data.save_neuron_secs_info(stimcell)
                
                # save vmem matrix, with Vm time course of each segment.
                if DUMP_VMEM:
                    with open(self.plots_folder_path + '/' + "stimcell.vmem_" + str(stimcell.neuron_idx) + ".dat", "w+") as dump_file:
                        pickle.dump(stimcell.vmem, dump_file)
                    
            # --------------------------------------------
            # saves ALL vmem of ALL simulations and iterations.
            # be careful: this will take around 250MB for each of the 2000 neurons! (for a total of ~500GB.)
            # please reduce the total number of neurons to just a few before turning this on.
    ##        pickle.dump(stimcell.vmem, open(self.plots_folder_path + '/' + "stimcell.vmem_" + str(i) + "-" + str(stim_intensity) + ".dat", "w+"))
    
            #debug
            print 'np.max(stimcell.vmem):', np.max(stimcell.vmem)
	    
        from codes.linear_alg_routines import magnitudes
        print "..... Search finished for neuron", str(stimcell.neuron_idx), "."
        
        base_E_mean_magnitude = np.mean(magnitudes(flattenLL(stimcell.E_vectors)))
        print "base E-field, mean magnitude: ", base_E_mean_magnitude, " V/m"
        print "==> θ threshold found: ", stim_intensity, " A/μs ; \t", stim_intensity*base_E_mean_magnitude, " V/m"
        
        result = {'neuron_idx':stimcell.neuron_idx, 'threshold_found':stim_intensity}
        pickle.dump(result, open(self.plots_folder_path + "/threshold_neuron_" + str(stimcell.neuron_idx) + ".dat", "w+"))

        # if the 'stimcell.E_vectors' is a single vector but repeated (i.e., the E-field is considered uniform along the neuron),
        # then 'stim_intensity*base_E_mean_magnitude' is the threshold in V/m. otherwise, if the E-field is not uniform, then it is just an
        # approximation based on the mean of the E-field vectors, since we can't speak of a single 'threshold E-field value'.
        result_threshold_Efield = {'neuron_idx':stimcell.neuron_idx, 'threshold_found_(V/m)':stim_intensity*base_E_mean_magnitude}
        pickle.dump(result_threshold_Efield, open(self.plots_folder_path + "/threshold_Efield_neuron_" + str(stimcell.neuron_idx) + ".dat", "w+"))
        
        return result
        
    # ------------------------------------------------------
    def run_cell_simulation(self, search_parameters):
        i, thresholds_search_range = search_parameters
        
        print "\n---- cell number: ", i
        cell_t_st = time.clock()
	
        # --------------------
        # neuron imports
        from neuron import h#, gui
        import neuron
        h.load_file('stdrun.hoc')
        h.load_file('import3d.hoc')
        
        # neuron vectors (to be used in the network version)
        t_vec = h.Vector()   # Spike time of all cells
        id_vec = h.Vector()  # Ids of spike times
        
        # --------------------
        # neuron creation
        stimcell = StimCell(self.cell_name,
                            neuron_idx=i,
                            population=self,
                            load_mechanisms=True)
        
        cell_t_en = time.clock(); verboseprint('create stimcell time:', cell_t_en-cell_t_st)
        
        # must be run after creating StimCell.
        h.tstop = self.tstop
        h.dt = self.dt
        
        # --------------------
        # transformations
        try:
            stimcell.position_neuron_at_relative_depth()
            stimcell.align_to_cortical_surface()
            stimcell.random_azimuthal_rotation()

            stimcell._update_pt3d()
        except:
            print("WARNING: NEURON SIM ABORTED: %s | %s" % (i, traceback.format_exc()))
            return None
        
        # --------------------
        # run iterations for threshold binary search.
        if RUN_SIMULATION:
            print "[!] uniques values of cm: ", np.unique([sec.cm for sec in h.allsec()]) #debug
            
            result = self.neuron_threshold_search_wrapped(stimcell, thresholds_search_range)
        else:
            print "*** WARNING: simulations will be SKIPPED."
            result = None
            
      # --------------------
        return result
    	
    # ------------------------------------------------------
    def simulate_cells(self,
                       thresholds_search_range,
                       run_parallel=True,
                       neurons_to_simulate=None):
        """ Simulate cells/neurons in the population.
        
        neurons_to_simulate : list of indices of neurons for which simulations
        will be done. If None, all N neurons of the population are simulated.
        """
        
        self.cells = []
        N = self._N
        
        # for cell creation
        sys.path.append(ROOT_DIR + "/core/AberraEtAl2018_edit/cells")
        os.chdir('../core/AberraEtAl2018_edit/cells')
        
        neurons_thresholds = []
        
        if neurons_to_simulate is None:
            neurons_to_simulate = range(N)
            print "All neurons of the population will be simulated. "
        else:
            print "Ids of neurons to be simulated: ", neurons_to_simulate

        print "Total number of neurons to simulate: ", len(neurons_to_simulate)
        
        # run simulations
        if run_parallel:
            print "\nRunning parallel."
            print "n_procs: ", n_procs
            
            paral_arguments = []
            for i in neurons_to_simulate:
                paral_arguments.append([i, thresholds_search_range])
            
            # TODO: as it is, the parallelized function ('f') must receive a single list as argument.
            # this could be fixed to have many arguments, in the proper way.
            result = parallel_map(f=self.run_cell_simulation,
                                  X=paral_arguments,
                                  nprocs=n_procs)
            neurons_thresholds.append(result)
        else:
            print "\nRunning serial."
            
            for i in neurons_to_simulate:
                result = self.run_cell_simulation([i, thresholds_search_range])
                neurons_thresholds.append(result)
        
        return neurons_thresholds
