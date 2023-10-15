# -*- coding: utf-8 -*-

import os
import pickle
import numpy as np

from neuron import h

from . import linear_alg_routines as linar
from .list_routines import *
from .neuron_routines import *


def save_terminals_info(stimcell, AP_beginning_dict):
    """ Save info about the directions of the sections
    containing each terminal. """
    
    soma_apic_line = linar.centroid(get_section_list_points(stimcell.cell.apic)) - stimcell.somapos
    stim_intensity = AP_beginning_dict['stim_intensity']
    E_vectors = stimcell.E_vectors # vectors of the uniform eletric field.
    flattened_E_vectors = np.array(flattenLL(stimcell.E_vectors))
    
    axon_section_types = ['axon','Myelin','Node','Unmyelin']
    labs = get_last_axon_branches(stimcell.cell, axon_section_types) # gets all labs of the current neuron.
    
    base_E_mean_magnitude = np.mean(linar.magnitudes(np.array(flattenLL(E_vectors))))
    print "base mean magnitude E (for 1.0 A/μs) : ", base_E_mean_magnitude
    print "stim. intensity: ", stim_intensity, " A/μs"
    print "-> peak E-field mean magnitude: ", stim_intensity*base_E_mean_magnitude, " V/m"
    
    # ------- get info about AP lab (i.e., lab where AP initiated)
    all_APlab_vectors = []
    all_angle_E_with_APlab = []
    all_angle_E_with_sda = []
    all_APlab_lengths = []

    E_vector = np.mean(flattened_E_vectors, axis=0)
    
    # there are possibly many segments we can say 'at which the AP initiated', since the Vm can cross 0 for the first time at many 
    # segments and all at the same time. so each one of them is accounted for here.
    for AP_seg_name in AP_beginning_dict['AP_beginning_seg_name']:
        # AP_seg_name is actually not a string; it is a list composed by segment index, section name (string) and seg.x 
        # (a float between 0 and 1 indicating the position of the segment along the arc length of the section).
        
        exec("AP_seg = h." + AP_seg_name[1] + '(' + str(AP_seg_name[2]) + ')')
        
        # find in which LAB APsec appears (if any).
        APsecname = AP_seg.sec.name()
        APlab_idx = None
        for i in range(len(labs)):
            for sec in labs[i]:
                if APsecname in sec.name():
                    APlab_idx = i
        
        if APlab_idx is not None:
            # info about the LAB where APsec appears.
            APlab = labs[APlab_idx]
            APlab_ends = get_seclist_ends([APlab])[0]
            APlab_vector = APlab_ends[1] - APlab_ends[0]
            
            # angle of mean E-field vector with last axon branch where the AP initiated.
            angle_E_with_APlab = linar.angle_vectors(APlab_vector, E_vector)
            
            APlab_length = np.linalg.norm(APlab_vector)
        else:
            # for the case that the AP did not occur in a LAB.
            APlab = None
            APlab_vector = None
            angle_E_with_APlab = None
            APlab_length = None
        
        # angle of mean E-field vector with somatodendritic axis ('sda'), which points from the soma center to centroid of apical dendrites.
        angle_E_with_sda = linar.angle_vectors(soma_apic_line, E_vector)
        
        all_APlab_vectors.append(APlab_vector)
        all_angle_E_with_APlab.append(angle_E_with_APlab)
        all_angle_E_with_sda.append(angle_E_with_sda)
        all_APlab_lengths.append(APlab_length)
    
    # ------- get info about all LABs of this neuron
    
    # create list where each element is a list containing the names of the sections of the corresponding LAB.
    all_lab_sec_names = []
    for lab in labs:
        lab_sec_names = []
        for sec in lab:
            lab_sec_names.append(sec.name())
        all_lab_sec_names.append(lab_sec_names)
    
    all_lab_ends = get_seclist_ends(labs)
    
##    # another definition of 'section vector'. we could define a vector for checking if
##    # a terminal-cointaining-section (i.e. the last branch of an axon collateral) is aligned with the electric field by:
##    # 1) getting the first and final points of the section (above);
##    # 2) getting the second-last and last points of the section (below).
##    alternative_AP_term_ends = np.array(get_section_points(AP_seg.sec))[[-2,-1]]
##    alternative_AP_section_vector = alternative_AP_term_ends[1] - alternative_AP_term_ends[0]
    
    # debug
    ##    # get labs pts
    ##    labs_pts = np.zeros((1,3))
    ##    for lab in labs:
    ##      for sec in lab:
    ##        labs_pts = np.concatenate((labs_pts,getsec3d(sec)),axis=0)
    ##    labs_pts = np.delete(labs_pts, (0), axis=0)
    ##
    ##    import my_helper_code.mayavi_functions as myv; from mayavi import mlab ; myv.draw_points(getcell3d(h.cell), scl=20.0, color=(0,0,0)); myv.draw_points(labs_pts, scl=50.0, color=(1,0,0)) ; mlab.show()
    ##
    
    # get angle between vectors of all LABs and vector of electric field.
    all_lab_vectors = []
    all_lab_angE = []
    for lab_ends in all_lab_ends:
        # vector connecting first and last points of the LAB (last axon branch).
        lab_vector = lab_ends[1] - lab_ends[0]
        all_lab_vectors.append(lab_vector)
        
        # store the angle of the E-field with the LAB vector.
        all_lab_angE.append(linar.angle_vectors(E_vector, lab_vector))

    # dictionary storing info about AP and LABs of this neuron.
    # if the AP did not occur in a LAB, then APlab_vector and the 3 subsequent fields will be 'None'.
    # however, the fields with the prefix 'labs_' will always be defined and will store data related to
    # the LABs of the current neuron.
    all_labs_dict = {'E_vector':E_vector,
                     
                     'APlab_vector':all_APlab_vectors,
                     'angle_E_with_APlab':all_angle_E_with_APlab,
                     'angle_E_with_sda':all_angle_E_with_sda,
                     'APlab_length':all_APlab_lengths,
                     
                     'labs_ends':all_lab_ends,
                     'labs_vectors':all_lab_vectors,
                     'labs_angE':all_lab_angE,
                     'labs_sec_names':all_lab_sec_names,
                     
                     'neuron_index':stimcell.neuron_idx,
                     'stim_intensity':AP_beginning_dict['stim_intensity'],
                     'AP_beginning_seg_name':AP_beginning_dict['AP_beginning_seg_name']
                    }
    
    # --------------------------------
    with open(stimcell.population.plots_folder_path + '/' + 'AP_beginning_site_' + str(stimcell.neuron_idx) + '.dat', 'w+') as dump_file:
        pickle.dump(AP_beginning_dict, dump_file)
    
    with open(stimcell.population.plots_folder_path + '/' + 'labs_dict_' + str(stimcell.neuron_idx) + '.dat', 'w+') as dump_file:
        pickle.dump(all_labs_dict, dump_file)
    
    return 0


def save_neuron_morphology(stimcell, filetype, save_path):
    '''
    Saves the neuron morphology at the given path.
    '''

    # check
    if not os.path.isdir(save_path):
        os.mkdir(save_path)

    # save morphology
    cell = stimcell.cell
    filename = save_path \
                + '/neuron_morph_' + str(stimcell.neuron_idx)
    if filetype != '':
        if filetype == 'obj':
            save_morphology_obj(cell, filename+'.obj')
        elif filetype == 'pts':
            with open(filename+'.dat', 'w+') as dump_file:
                pickle.dump(get_cell_points(cell), dump_file)

##    soma_apic_line = linar.centroid(\
##        get_section_list_points(stimcell.cell.apic)) - stimcell.somapos
##    with open(save_path +
##              'soma_apic_line_' +
##              str(stimcell.neuron_idx) + '.dat', 'w+') as dump_file:
##        pickle.dump({'somapos':stimcell.somapos,
##                     'soma_apic_line_vector':(stimcell.somapos + \
##                                              soma_apic_line),
##                     'polygon_center':pop.gm_centers[stimcell.neuron_idx],
##                     'polygon_normal':pop.gm_normals[stimcell.neuron_idx]},
##                    dump_file)

def save_neuron_secs_info(stimcell):
    ''' Saves information about the neuron's sections. '''
    
    sections_info = []
    cell_seg_centers = stimcell.calculate_segments_centers() # centers of all segments of the cell
    
    for idx, sec in enumerate(stimcell.section_list):
        parentseg = sec.parentseg()
        parentsec_name = parentseg.sec.name() if parentseg is not None else None
        
        sec_info = {'name': sec.name(),
                    'index_in_section_list': idx,
                    'nseg': sec.nseg,
                    'n3d': sec.n3d(),
                    'parent_sec_name': parentsec_name,
                    'children_sec_name': [childrensec.name() for childrensec in sec.children()],
                    '3d_pts': [[sec.x3d(i), sec.y3d(i), sec.z3d(i)] for i in range(sec.n3d())],
                    'seg_centers': cell_seg_centers[idx],
                    'seg_x': [seg.x for seg in sec]
                    }
        sections_info.append(sec_info)
    
    fn = stimcell.population.neuron_3d_data_path + 'sections_info_of_neuron_'+str(stimcell.neuron_idx)+'.dat'
    
    with open(fn, 'w+') as dump_file:
        pickle.dump(sections_info, dump_file)
    
    print 'secs info saved at: ', fn


point_index = 0
points = []
lines = []

def save_morphology_obj(cell, filename):
    """ Writes the 3D points of a given NEURON cell to an '.obj' file. """
    
    global point_index
    global morph_file
    
    point_index = 0
    
    morph_file = open(filename, "w")
    
    root_section = cell.soma[0]
    write_section(root_section)
    
    morph_file.close()


def write_section(section):
    '''
    Writes the 3D points of a given section and of all its
    children sections (i.e., recursively) to a file opened as 'morph_file'.
    '''

    global point_index

    # writes all points of current section
    for i in range(0, section.n3d()):
        # write point
        morph_file.write("v %f %f %f\n" %
                         tuple([section.x3d(i),
                                section.y3d(i),
                                section.z3d(i)]))

        if i != 0:
            # when i == 0, there is no need to write a line.

            # write line: current point (point_index) and
            # previous point (point_index-1)
            # '+ 1' to adjust to obj file format's convention
            # to start at 1, instead of starting at 0 like Python
            # and NEURON indices.
            morph_file.write("l %d %d\n" % tuple([point_index-1 + 1,
                                                  point_index + 1]))

        point_index += 1

    if section.children() != None:
        current_section_last_point_index = point_index - 1

    # writes children of current section
    for child_section in section.children():
        if child_section != None:
            child_section_first_point_index = point_index

            # '+ 1' to adjust to obj file format's convention to start at 1,
            # instead of starting at 0 like Python and NEURON indices.
            morph_file.write("l %d %d\n" %
                             tuple([current_section_last_point_index + 1,
                                    child_section_first_point_index + 1]))

            write_section(child_section)


def save_obj(vertices, lines, obj_file):
    save_obj_vertices(vertices, obj_file)
    save_obj_lines(lines, obj_file)


def save_obj_vertices(vertices, obj_file):
    for i in range(0, len(vertices)):
        # write vertex
        obj_file.write("v %f %f %f\n" % tuple(vertices[i]))


def save_obj_lines(lines, obj_file):
    for i in range(0, len(lines)):
        obj_file.write("l %d %d\n" % tuple(lines[i]))
