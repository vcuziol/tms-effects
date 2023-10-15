# -*- coding: utf-8 -*-
# files
import os
import time
import datetime
import sys

# math
import numpy as np
import copy

# NEURON
import neuron
from neuron import h

# tms_networks
from ..parameters import *
import codes.linear_alg_routines as linar
from codes.list_routines import *
from codes.neuron_routines import *
from codes.LFPyCell import Cell

from ..debugging import write_log


class StimCell(Cell):
    ''' Representation of a neuron under external stimulation
        by a uniform electric field. '''
    
    def __init__(self,
                 cell_name, population=None, neuron_idx=0, 
                 v_init=-70.0, tstart=0, tstop=20.0, dt=0.025,
                 nsegs_method='lambda100', lambda_f=100, d_lambda=0.1,
         max_nsegs_length=None, pt3d=True, load_mechanisms=False,
                 verbose=False, make_checks=False, **kwargs):
        
        t_start = time.time(); #debug
        
        # original cell/neuron model
        self.cell = self.create_cell(cell_name, load_mechanisms)
        
        write_log('time', ['create cell', t_start, time.time()]) #debug
        
        self.cell_name = cell_name
        
        t_start = time.time() # debug
        # creates a list of the sections of the cell in the order that they
        # are traversed in the topology tree.
        # this list is relevant to the calculations of the quasipotentials.
        self.section_list = traverse_topology(self.cell.soma[0])
        write_log('time', ['traverse topology', t_start, time.time()]) #debug
        
        self.allseclist = self.section_list # for compatibility with LFPy code.
        
        Cell.__init__(self,
                      cell_seclist=self.section_list,
                      v_init=v_init, tstart=tstart, tstop=tstop, dt=dt,
                      nsegs_method=nsegs_method, lambda_f=lambda_f,
                      d_lambda=d_lambda, max_nsegs_length=max_nsegs_length,
                      pt3d=pt3d, verbose=verbose, **kwargs)
        
        # index of neuron inside the population
        self.neuron_idx = neuron_idx
        
        self.E_vectors = None
        
        # population:
        # the population attributes are given priority over the attributes of
        #  the single cell.
        self.population = population
        if population is not None:
            self.v_init = population.v_init
            self.dt = population.dt
            self.tstop = population.tstop
        else:
            self.v_init = v_init
            self.dt = dt
            self.tstop = tstop

##        h.finitialize(self.v_init)
        
        write_log('stimcell', ['stimcell __init__ is finished.']) #debug
        
        # disabling of main axon terminal done by Aberra et al. 2019
        if 'L5' in cell_name:
            assert(len(h.secrefs[int(h.min_sec_ind)].sec.children()) == 0)
            h.secrefs[int(h.min_sec_ind)].sec(1.0).diam = 1000.0
            write_log('stimcell', ['main axon terminal disabled: ',
                h.secrefs[int(h.min_sec_ind)].sec.name()])
    
        if make_checks:
            # check if there are isolated sections (with no parent 
            # and no children).
            check_topology()
            
            if len(list(self.section_list)) != len(list(self.cell.all)):
                write_log('stimcell', ['''WARNING: self.section_list has
                    different number of sections from cell.all: ''',
                    len(list(self.section_list)), ' != ',
                    len(list(h.allsec()))])
    
    def create_cell(self, cell_name, load_mechanisms=False):
        ''' Creates the cell object in NEURON of the given cell type.
            The neuron models were adapted by Aberra et al. (2018)
            (ModelDB accession number: 241165) from BlueBrain project models.
            Obs.: To allow for correct instantiation of these cells,
            change NSTACK to 100000 and NFRAME to 20000 in 'nrn.defaults'
            (usually located in '/usr/local/nrn/share/nrn/lib/nrn.defaults').
        '''
        
        t_start = time.time(); #debug
        
        h.load_file('stdgui.hoc')
        h.load_file('import3d.hoc')
        
        if cell_name == 'L5_simplified_cell':

            from L5_simplified_cell import create_L5_simplified_cell
            
            os.chdir(ROOT_DIR + '/core/L5_simplified_cell')
            cell = create_L5_simplified_cell.create_cell(mechanisms_path=\
                ROOT_DIR + '/core/L5_simplified_cell/mechanisms')
            h('objref cell')
            h.cell = cell
            
            # prepare to load Aberra et al. 2018 modifications
            os.chdir(ROOT_DIR + '/core/AberraEtAl2018_edit')
        else:
            
            os.chdir(ROOT_DIR + '/core/AberraEtAl2018_edit')
            
            if load_mechanisms:
                neuron.load_mechanisms('mechanisms')
            
            os.chdir('cells/' + cell_name)
        
            try:
                # Load morphology
                neuron.h.xopen("morphology.hoc")
                # Load biophysics
                neuron.h.xopen("biophysics.hoc")
                # Load main cell template
                neuron.h.xopen("template.hoc")
            except:
                write_log('stimcell',
                          ['Template not redefined; this is expected.'])
            
            os.chdir('../..')
        
        # --
        t_start = time.time(); #debug
        
        # load 'init.hoc' from Aberra et al. 2018 code.
        h('objref cell')
        h.load_file("interpCoordinates.hoc")
        h.load_file("setPointers.hoc")
##        h.load_file("calcVe.hoc")
##        h.load_file("stimWaveform.hoc")
        h.load_file("cellChooser.hoc")
        h.load_file("setParams.hoc")
        
        # if we are using 'L5_simplified_cell', then we do not want the myelin
        # biophysics to be changed; we want the original parameters of Rusu
        # et al. 2014 to be maintained. 
        # 'axonEdit1_L5simplifiedCell.hoc' was modified for this purpose.
        if cell_name == 'L5_simplified_cell':
            h.load_file('axonEdit1_L5simplifiedCell.hoc')
        else:
            h.load_file('axonEdit1.hoc')
        
        h.celsius = 37
        h.setParamsAdultHuman()
        
        if cell_name != 'L5_simplified_cell':
            # TODO: the index of the chosen cell must correspond to the order
            # of the cell names given inside 'cellChooser.hoc';
            # the code of 'cellChooser.hoc' must be changed to use
            # automatically the list of the folders' names.
            all_cell_names = [stringobj.s for stringobj in list(h.cell_names)]
            chosen_cell_index = all_cell_names.index(cell_name) + 1
            
            write_log('time', ['load_aberra_files', t_start, time.time()])
            # ----------------
            t_start = time.time(); #debug 
            # create neuron and make modifications (such as geometry
            # and myelination).
            h.create_chosen_cell(chosen_cell_index)
            write_log('time', ['create cell', t_start, time.time()])
        else:
            h.cell = cell
            
            # snippet from 'create_chosen_cell'
            h('totalnseg = 0')
            h('forsec cell.all { totalnseg += nseg }')
            write_log('stimcell', ['totalnseg: ', h.totalnseg])
            h('forsec cell.all { \n insert xtra \n insert extracellular \n }')
            
        t_start = time.time(); #debug
        
        # make cell morphology modifications (from Aberra et al. 2018)
        h.modify_cell_geometry()
            
        h('iseg_secList = new SectionList()')
        h('myelin_secList = new SectionList()')
        h('node_secList = new SectionList()')
        h('unmyelin_secList = new SectionList()')
            
        # for storing the new sections of the axon (i.e., myelin, node and
        # unmyelin)
        h('axonal = new SectionList()')

        # increments numMyelin, numNode, and numUnmyelin
        h.count_myelin_sections(h.cell.axonal) 
        
        # temporary fix to L5 cell having remaining isolated sections
        # (i.e., no parent and no children) in memory.
        if 'L5' in cell_name: #debug
            h.numMyelin = h.numMyelin - 1
            h.numNode = h.numNode - 1
        
        write_log('stimcell', ['numMyelin, numNode, numUnmyelin: ',
            h.numMyelin, h.numNode, h.numUnmyelin])

        # axon_secList must exist after count_myelin_sections
##        write_log('stimcell', ["axon_secList: \n", list(h.axon_secList)])
##        write_log('stimcell', ["cell.axonal: \n", list(h.cell.axonal)])
    
        # --
    
        nMyelin, nNode, nUnmyelin = h.numMyelin, h.numNode, h.numUnmyelin
        h.cell.redefine_myelin(nMyelin, nNode, nUnmyelin)

        h('totalnsec=0'); h('forall { totalnsec=totalnsec+1 }')
        
        write_log('stimcell', ['cell.all:', len(list(h.cell.all)),
                               ' \t--\t totalnsec: ', h.totalnsec])
        
        h('for i=0,numMyelin-1 { cell.myelin[i] { axonal.append() }}')
        h('for i=0,numNode-1 { cell.node[i] { axonal.append() }}')
        h('for i=0,numUnmyelin-1 { cell.unmyelin[i] { axonal.append() }}')
        
        h('chdir("../..")')
        
        os.chdir(ROOT_DIR + '/core/AberraEtAl2018_edit')
        h.load_file("editMorphology.hoc")
        
        write_log('stimcell', 
            ['myelin_secList: ', len(list(h.myelin_secList)), '\t', 
            'node_secList: ', len(list(h.node_secList)), '\t',
            'unmyelin_secList: ', len(list(h.unmyelin_secList))])
        
        h('for i=0,numMyelin-1 { cell.myelin[i] {  \
        myelin_secList.append() \
        axonal.append() } }')
        h('for i=0,numNode-1 { cell.node[i] {  \
            node_secList.append() \
            axonal.append() } }')
        h('for i=0,numUnmyelin-1 { cell.unmyelin[i] { \
            unmyelin_secList.append() \
            axonal.append() } }')
        
        write_log('stimcell', ['len cell.all: ', len(list(h.cell.all))])
        
        h.add_cell_myelin()
        
        if cell_name == 'L5_simplified_cell':
            create_L5_simplified_cell.load_axon_biophysics(cell)
        
        # After any change to cell geometry or nseg, be sure to invoke setpointers().
        h.setpointers()

        write_log('stimcell', ['== create_cell (on Stimcell) has finished. =='])
        
        os.chdir('..')
        
        return h.cell
    
    #debug
    def create_bluebrain_cell(self, cell_name, load_mechanisms=True):
        ''' (FOR DEBUGGING) Instead of creating the cell modified by Aberra 
            et al. 2018, creates the original Blue Brain Project cell. '''
        
        import neuron
        
        if load_mechanisms:
            neuron.load_mechanisms(ROOT_DIR \
                + '/core/AberraEtAl2018_edit/mechanisms')
        
        os.chdir(ROOT_DIR + '/core/AberraEtAl2018_edit/cells')
        h.load_file('nrngui.hoc')
        h.load_file('import3d.hoc')
        sys.path.append(os.getcwd())
        exec('import ' + cell_name + '.run as cell_model')

        os.chdir(cell_name)
        try:
            # Load morphology
            h.xopen("morphology.hoc")
        except:
            write_log('stimcell', ['morphology.hoc not redefined.'])
        
        try:
            # Load biophysics
            h.xopen("biophysics.hoc")
        except:
            write_log('stimcell', ['biophysics not redefined.'])
        
        try:
            # Load main cell template
            h.xopen("template.hoc")
        except:
            write_log('stimcell', 
            ["Template not redefined; this is expected."])
        
        cell = cell_model.create_cell(False)
        os.chdir('..')
        
        return cell
    
    
    def set_rotation_vectorial(self, initial_vector, final_vector): #v
        '''
        Rotate morphology of cell according to the direction given by
        initial_vector and final_vector
        (given that the cell's somatodendritic axis is oriented
        according to 'initial_vector'.)
        '''
        
        theta = linar.ccw_angle_vectors(initial_vector, final_vector)
        M_rot = np.matrix(
            linar.rotation_matrix(-np.cross(initial_vector, final_vector),
                                theta))
        
        rel_start, rel_end = self._rel_positions()
        
        rel_start = rel_start * M_rot
        rel_end = rel_end * M_rot
        
        self._real_positions(rel_start, rel_end)

        # rotate the pt3d geometry accordingly
        if self.pt3d and hasattr(self, 'x3d'):
            self._set_pt3d_rotation_vectorial(initial_vector, final_vector)


    def _set_pt3d_rotation_vectorial(self, initial_vector, final_vector):
        ''' '''

        theta = linar.ccw_angle_vectors(initial_vector, final_vector)
        M_rot = np.matrix(
            linar.rotation_matrix(-np.cross(initial_vector, final_vector),
                                theta))
        
        for i in range(len(self.x3d)):
            rel_pos = self._rel_pt3d_positions(self.x3d[i],
                                               self.y3d[i],
                                               self.z3d[i])
            
            rel_pos = rel_pos * M_rot
            
            self.x3d[i], self.y3d[i], self.z3d[i] = \
                                        self._real_pt3d_positions(rel_pos)

        self._update_pt3d()


    def somatodend_axis(self):
        return linar.centroid(
            get_section_list_points(self.cell.apic)) - self.somapos
        

    def calculate_segments_centers(self, with_seg_names=False, flatten=False):
        # This is an adaptation of the "grindaway" function from
        # "extracellular_stim_and_rec".
        
        segments_centers = []
        seg_count = 0
        
        for sec in self.section_list:

            section_seg_centers = []
            
            # get data for the section
            # nn = sec_num_pts
            sec_num_pts = int(sec.n3d()) # number of 3d points for the current section
            xx = []
            yy = []
            zz = []
            length = []
            
            # for each point, get x,y,z coordinates and arc length position.
            for ii in range(0, sec_num_pts-1): 
                xx.append(sec.x3d(ii))
                yy.append(sec.y3d(ii))
                zz.append(sec.z3d(ii))
                length.append(sec.arc3d(ii))
            
            length = np.array(length)
            if int(length[-1]) != 0:
                length = length/(length[-1])
            
            # initialize the x-coordinates at which to evaluate
            # the interpolated values.
            rangev = []
            rangev_step = 1.0/sec.nseg
            rangev_length = sec.nseg+2
            rangev = np.linspace(0, 0+(rangev_step*rangev_length),
                                 rangev_length, endpoint=False)
            rangev = rangev - 1.0/(2.0*sec.nseg)
            rangev[0] = 0.0
            rangev[-1] = 1.0
            
            # numpy interp function: y = interp(x, xp, fp), where
            # y are the interpolated values.
            xint = np.interp(rangev, length, xx)
            yint = np.interp(rangev, length, yy)
            zint = np.interp(rangev, length, zz)
            
            # stores the segment centers separately, by section.
            for ii in range(0, sec.nseg):
                seg_count = seg_count + 1
                xr = rangev[ii]
                section_seg_centers.append([xint[ii], yint[ii], zint[ii]])

            segments_centers.append(np.array(section_seg_centers))
        
        # puts segments names in a list
        if with_seg_names:
            seg_names_grouped_by_sec = []
            
            for sec in self.section_list:
                seg_names = []
                for seg in sec:
                    seg_names.append(sec.name() + "(" + str(seg.x) + ")")
                    
                seg_names_grouped_by_sec.append(seg_names)
            
            segments_centers = [segments_centers, seg_names_grouped_by_sec]
        else:
            segments_centers = np.array(segments_centers)
        
        # returns centers, either flattened
        # (i.e., each element is a segment center) or grouped by section.
        if flatten:
            return np.array(flattenLL(segments_centers))
        else:
            return segments_centers
        

    def position_neuron_at_relative_depth(self):
        '''
        Translate the neuron to be at a distance from the cortical
        surface at a fraction of the cortical depth at that point.
        A population must be defined for this neuron so this method
        can to be used.
        '''

        # check soma position
        try:
            a = self.population.soma_centers[self.neuron_idx]
            b = self.population.gm_centers[self.neuron_idx]
            c = self.population.gm_centers[self.neuron_idx] \
                + self.population.gm_normals[self.neuron_idx]
            assert(linar.angle_vectors(b-a, c-a) < 0.1)
        except AssertionError as e:
            import pdb; pdb.set_trace() #v 

        soma_center = self.population.soma_centers[self.neuron_idx]

        # translate the entire neuron so that its soma center is
        # positioned in the given coordinates.
        self.set_pos(soma_center[0],
                     soma_center[1],
                     soma_center[2])

        # check
        try:
            #debug: if 3d pts are being updated after each transform
            # and not manually, uncomment this.
            assert(sum(self.somapos
                       -self.population.soma_centers[self.neuron_idx]) < 0.01)
        except AssertionError:
            import pdb; pdb.set_trace() #v


    def rotate_axis(self, final_axis):
        '''
        Set somatodendritic axis of the cell to be oriented in
        the same direction as the given vector.
        When cell/neuron is created, the soma location is near (0,0,0).
        The following code rotates the cell so that its somatodendritic axis
        is aligned with the normal of the corresponding cortex mesh polygon.
        Since the soma position is necessarily in the span of the vector
        normal to the polygon (because it was positioned using the normal
        vector), then the rotation of the cell to be aligned with the
        vector "final_sal" - which connects soma and center of polygon -
        will make the neuron's axis to be normal to the cortical surface
        (as it happens in the layers of the human brain cortex).
        '''

        # calculate the somatodendritic axis, i.e., the vector
        # that connects the soma center to the centroid
        # of the apical dendrites.
        old_axis = self.somatodend_axis()

        # rotate cell morphology, changing the somatodendritic axis
        self.set_rotation_vectorial(old_axis, final_axis)

        new_axis = self.somatodend_axis()

        # check
        epsilon = 0.1
        try:
            assert(linar.ccw_angle_vectors(final_axis, new_axis) < epsilon)
        except AssertionError as e:
            import traceback; traceback.print_exc()
            import pdb; pdb.set_trace() #v

    def align_to_cortical_surface(self):
        '''
        Rotate the neuron so its somatodendritic axis becomes orthogonal
        to the cortical surface. A population must be defined for this
        neuron for this method to be used.
        '''

        # the normals of the mesh are pointing 'inwards', so we invert them.
        final_axis = -self.population.gm_normals[self.neuron_idx]

        # check
        try:
            assert(abs(linar.angle_vectors(final_axis,
                self.population.gm_normals[self.neuron_idx]) - np.pi) < 0.1)
        except AssertionError:
            import pdb; pdb.set_trace() #v
        
        self.rotate_axis(final_axis)

        return 0


    def azimuthal_rotation(self, rotation_angle):
        '''
        Rotate counter-clockwise the neuron around its own somatodendritic
        axis ("somatodend_axis").
        '''

        old_axis = self.somatodend_axis()

        rot_initial_vector = linar.random_orthogonal_vector(old_axis)
        rot_final_vector = np.dot(linar.rotation_matrix(old_axis,
                                                        rotation_angle),
                                        rot_initial_vector)

        # apply azimuthal rotation
        self.set_rotation_vectorial(rot_initial_vector, rot_final_vector)

        new_axis = self.somatodend_axis()

        # check
        epsilon = 0.1
        assert(linar.angle_vectors(old_axis, new_axis) < epsilon)

        return 0


    def random_azimuthal_rotation(self):
        '''
        Rotate the neuron around its own somatodendritic axis 
        ("somatodend_axis") by a random angle so that, considering the 
        entire population of neurons,the axon collaterals are not 
        preferentially aligned to any specific direction. (If they were, 
        it could bias the simulation results, since they depend on
        the applied electric field direction.)
        '''

        random_azimrot_angle = 2.0*np.pi*np.random.random()
        self.azimuthal_rotation(random_azimrot_angle)
        
        return 0
        

    def calculate_cell_quasipotentials(self, E_field, segments_centers=None):
        '''
        Calculate quasipotentials by numerical integration of a given 
        eletric field's values,
        following the order of segments given by 'self.section_list'.

        E_field : E-field 3D vectors given as a list of lists, where each 
        list contains the vectors for the segments of a given section.
        '''
        
        segment_index = 0
        quasipotentials = []  # list of lists
        sec_names = []
        
        # 'centers': centers of the segments of the cell given as a list of 
        # lists, where each list contains the segment centers of 
        # a given section.
        if segments_centers is None:
            centers = self.calculate_segments_centers()
        else:
            centers = segments_centers
        
        for sec in self.section_list:
            
            segment_in_section_index = 0
            section_quasipotentials = []
            
            for seg in sec:
                # if the segment is the root segment (first segment of soma):
                if segment_index == 0:
                    section_quasipotentials.append(0.0) # phi_soma = 0
                    segment_in_section_index += 1
                    segment_index += 1
                    continue
                # if the segment is the first of the section:
                elif segment_in_section_index == 0:
                    # get previous section's id, for use as
                    # index with E_field and centers
                    previous_sec_name = sec.parentseg().sec.name()
                    previous_sec_id = sec_names.index(previous_sec_name)
                    
                    # displacement vector
                    s_pc = centers[len(sec_names)][segment_in_section_index] \
                        - centers[previous_sec_id][-1]
                    E_p = E_field[previous_sec_id][-1]
                    
                    # get the quasipotential of the
                    # previous section's last segment.
                    phi_p = quasipotentials[sec_names.index(previous_sec_name)][-1]
                # if the segment is other than the first of the section:
                else: 
                    s_pc = centers[len(sec_names)][segment_in_section_index] \
                        - centers[len(sec_names)][segment_in_section_index-1]
                    E_p = E_field[len(sec_names)][segment_in_section_index-1]
                    
                    phi_p = section_quasipotentials[-1]
                    
                E_c = E_field[len(sec_names)][segment_in_section_index]
                
                # converts from micrometers to milimeters, so that E-field
                # of unit mV/mm (which is equivalent to V/m) can be used.
                s_pc = s_pc * 1e-3
                
                # calculate quasipotential of current segment
                phi_c = phi_p - 0.5 * np.dot((E_c + E_p), s_pc)

                section_quasipotentials.append(phi_c)
                
                segment_in_section_index += 1
                segment_index += 1

            assert(len(list(sec)) == len(section_quasipotentials)) #debug

            sec_names.append(sec.name())
            quasipotentials.append(section_quasipotentials)

        return quasipotentials


    def set_E_field(self, E_vectors=None):
        '''
        Sets electric field vectors defining the stimulus over this cell,
        and calculates quasipotentials (i.e., electric potential under
        the quasistatic assumption).
        '''
        
        segments_centers = self.calculate_segments_centers()

        # this will become something other than None only if a uniform e-field
        # is meant to be used.
        E_vector = None
        
        # ------------ E-field
        if E_vectors is None:
            if self.population is not None:
                # 'neuron_idx' refers to the index of the neuron in the order of the current simulations.
                # 'neurons_polygon_indices' is the index of the polygon (of the GM mesh) that corresponds to the neuron.
        ##                neuron_E_vector = self.population.tms_sim.E_field_vectors[self.population.neurons_polygon_indices[self.neuron_idx]]

                if INTERPOLATE_AT_SEGMENTS:
                    # flattens segments centers; and converts from micrometers to milimeters.
                    fl_segcenters_mm = np.array(flattenLL(segments_centers)) * 1e-3
                    
                    # save, in the simulation results folder, a csv file containing the segments centers of this neuron.                
                    csv_file_path = self.population.neuron_Efield_vecs_path + 'segcenters_mm_neuron_' + str(self.neuron_idx) + '.csv'
                    np.savetxt(csv_file_path, fl_segcenters_mm, delimiter=',')

                    # files for the interpolated fields will be created inside the folder of the original csv file.
                    # obs: this command requires that simnibs is installed, since it uses the 'get_fields_at_coordinates' function.
                    msh_file_path = MSH_INPUTS_FOLDER + '/' + self.population.simnibs_session_name + '/' + self.population.msh_fn
                    command_string = get_fields_command + ' -s ' + csv_file_path + ' -m ' + msh_file_path
                    write_log('stimcell', 
            ['running get_fields_at_coordinates as: ', 
            command_string])
                    os.system(command_string)
                    
                    # load interpolated E-field file.
                    interpolated_Efield_file_path = csv_file_path[:-4] + '_E.csv'
                    E_vectors_flattened = np.loadtxt(interpolated_Efield_file_path, delimiter=',')
                    
                    # this 'unflattens' the flattened E-field (i.e, not grouped by sections), loaded from the csv file.
                    # each element of 'segments_centers' is a list of all segments centers of a single section.
                    E_vectors = []
                    all_seg_index = 0
                    for sec_centers in segments_centers:
                        E_vectors.append(E_vectors_flattened[all_seg_index:all_seg_index+len(sec_centers)])
                        all_seg_index += len(sec_centers)

                else:
                    # E-field was interpolated only at center of soma and was stored in 'interpolated_E_field' during init of NeuronPopulation.
                    E_vector = self.population.interpolated_E_field[self.neuron_idx]

            else:
                if DEFAULT_E_VECTOR is not None:
                    # default E-vector, defined in parameters.py, is used.
                    E_vector = DEFAULT_E_VECTOR
                    write_log('stimcell', 
            ["WARNING: DEFAULT_E_VECTOR is being used: ", 
                        DEFAULT_E_VECTOR])
                else:
                    raise TypeError('''Electric field vectors for StimCell are 
                        not defined anywhere (neuron_E_vectors,
                        DEFAULT_E_VECTOR, or self.population.). ''')
        
        if E_vector is not None:
            # assume E-field to be the same in all segments
            # (i.e., uniform E-field).
            # this replicates the E-field vector to create 'E_field', which
            #  is the E-field at each segment; the segments are grouped by
            # section and ordered according to 'self.section_list'.
            E_vectors = []
            for sec_sc in segments_centers:
                sec_ief = np.array([E_vector]*len(sec_sc))
                E_vectors.append(sec_ief)
            
        # check
        for i in range(len(E_vectors)):
            assert(len(E_vectors[i]) == len(segments_centers[i]))

        self.E_vectors = E_vectors
        
        # quasipotentials are calculated using the E-field at each segment's center.
        self.quasipotentials = self.calculate_cell_quasipotentials(E_field=\
            self.E_vectors, segments_centers=segments_centers)

        # values of electric potential at segments centers, but flattened, 
        # i.e., not grouped by sections.
        # (in 'self.quasipotentials', they are grouped by section.)
        self.v_segments = np.array(flattenLL(self.quasipotentials))
    

    def set_stimulation(self, stim_intensity=1.0, stim_time_course=None):
        ''' Sets stimulus intensity and time course, which modulate the 
            electric field vectors pre-stored in 'self.E_vectors'.
            If stim_time_course is None, the default time course given by the 
            Population object (self.tms_sim.stim_time_course) is used.
        
            There are two different methods - 'set_E_field' and 
            'set_stimulation' - because it is desirable to change only stimulus 
            intensity, which multiplies predefined E-field vectors.
        
        stim_intensity : intensity of external stimulus, to be applied by 
                         scaling the electric field time course. In the case 
                         of magnetic stimulation and any electric field 
                         calculated by SimNIBS (under quasi-static assumption), 
                         this intensity is given in units of A/μs.
             (ampere per microsecond).
        neuron_E_vector : vector of electric field for this neuron, assuming
                          the electric field is uniform across all of the 
                          neuron's morphology. '''
        
        import TMSSimulation
        tms_sim = TMSSimulation.TMSSimulation()
        
        h.dt = self.dt
        h.tstop = self.tstop
        
        if self.E_vectors is None:
            self.set_E_field()
        
        # change stimulation intensity by scaling the pulse waveform.
        if stim_time_course is None:
            if self.population is not None:
                stim_time_course = stim_intensity \
                    * self.population.tms_sim.stim_time_course
            else:
                # generate waveform with parameters given in 'parameters.py'
                from parameters import R,L,C
                stim_time_course = stim_intensity * tms_sim.\
                    generate_rlc_waveform(h.dt, h.tstop, 1.0, 1.0, 1.0, 
                                          Cc=C, Rc=R, Lc=L)
##                stim_time_course = tms_sim.generate_rlc_waveform(h.dt, h.tstop, 1.0, 1.0, 1.0, Cc=C, Rc=R, Lc=L) #debug
        else:
            stim_time_course = stim_intensity * stim_time_course
        self.stim_time_course = stim_time_course
        self.stim_intensity = stim_intensity

        # ------------ v_ext
        write_log('stimcell', ['ief (1.0 A/μs) magnitude: ', 
        linar.magnitudes([flattenLL(self.E_vectors)[0]]), ' V/m'])
        write_log('stimcell', ['ief (w/ current intensity) magnitude: ', 
            np.max(stim_time_course) * linar.magnitudes([flattenLL(\
            self.E_vectors)[0]]), ' V/m'])
        write_log('stimcell', ['max v_segments (at current intensity): ', 
            np.max(self.v_segments)*stim_intensity, ' mV'])
        
        # generate v_ext matrix
        v_ext = tms_sim.build_v_ext(self.v_segments, stim_time_course)
        t_ext = np.arange(h.tstop / h.dt + 1) * h.dt
        
        assert(np.isnan(np.sum(np.array(v_ext.ravel()))) == False)
        write_log('stimcell', ['max v_ext: ', np.max(v_ext)]) #debug
        
        # insert extracellular potentials
        self.insert_v_ext(v_ext, t_ext)
        
        # ------------
        return self.v_segments
    

    def simulate_neuron(self):
        ''' Run the simulation for this neuron, with the predefined h.dt 
            and h.tstop. '''
        
        # debug:
        # defined from 'setPointers.hoc', in the Aberra et al. models folder
#        h.setpointers()
        
        # run simulation
        write_log('stimcell', ['simulation started. ', '- tstop: ', h.tstop, 
            ' / ', 'dt: ', h.dt])
        sim_t_st = time.clock()
        self.simulate(rec_imem=True, rec_vmem=True, rec_ipas=True, 
                      rec_icap=True,
                      #rec_variables=['cai'],
                      tstart=0, tstop=h.tstop, dt=h.dt)
        sim_t_en = time.clock()
        write_log('stimcell', ['simulation finished; time: ', 
            sim_t_en - sim_t_st])
        
        # detects if an action potential occured in this neuron or not.
        was_activated = self.detect_neuron_activation(self.vmem)
        
        return was_activated
    

    def detect_neuron_activation(self, vmem):
        ''' After simulation, the 'vmem' matrix is checked to detect
            activation (i.e., occurrence of action potential).
            If more than 3 segments/compartments reached - at least once - 
            a membrane potential greater than 0.0 in a single neuron, the 
            neuron is considered as having activated. '''
        
        n_activated_segments = 0
        
        for i in range(vmem.shape[0]):
            if sum(vmem[i,:]>0.0) > 0:
                n_activated_segments += 1

        if n_activated_segments > 3:
            return 1
        else:
            return 0
        

    def detect_AP_beginning_site(self, vmem, dt):
        
        n_rows = vmem.shape[0]
        
        # stores the first time that Vm crosses 0 for each row (i.e., for 
        # each segment, since each row represents the time 
        # course of Vm for a single segment).
        rows_idx_first_crossing = []
        
        # for each segment, find the first time that Vm crosses 0 (if any).
        for i in range(n_rows):
            idces_with_v_over_zero = np.where(vmem[i]>0.0)[0]
            
            if len(idces_with_v_over_zero) > 0:
                idx_of_first_zero_crossing = idces_with_v_over_zero[0]
                rows_idx_first_crossing.append(idx_of_first_zero_crossing)
            else:
                # numpy.argmin (used later) ignores empty lists
                rows_idx_first_crossing.append([])
        
        assert(len(rows_idx_first_crossing) == n_rows)
        
        rows_idx_first_crossing = np.array(rows_idx_first_crossing, 
            dtype='object')
        
        # case when the membrane voltage never goes over 0 in any segment.
        if len(rows_idx_first_crossing) <= 0 or \
           np.array(rows_idx_first_crossing).size == 0:
            return {'AP_beginning_seg_idx':None, 'AP_start_t':None}
        
        # earliest time that Vm crosses 0 at any segment, given in number of 
        # time steps since time 0.
        AP_time_step = np.min(rows_idx_first_crossing)
        
        # get all segments where AP begins at the earliest found time.
        AP_beginning_segs = np.where(\
            rows_idx_first_crossing == AP_time_step)[0]
        
        # time instant (in ms) when AP first crosses zero in any segment.
        AP_beginning_t = np.array(AP_time_step) * dt
        
        # AP beginning segment indices must be used in the list of segments of
        # the cell to find the name of the sections and their exact position.
        return {'AP_beginning_seg_idx':AP_beginning_segs,
                'AP_start_t':AP_beginning_t}
        

    # aliases
    stim = set_stimulation
    run = simulate_neuron
