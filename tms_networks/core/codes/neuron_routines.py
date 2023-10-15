import numpy as np
from list_routines import *
from neuron import h


def get_section_points(sec):
    L = []
    for i in range(int(sec.n3d())):
        L.append([sec.x3d(i), sec.y3d(i), sec.z3d(i)])
    return np.array(L)


def get_section_list_points(seclist):
    '''
    (after calling "pts = get_section_list_points(cell.dend)",
    use pts[0] to access the points.)
    '''
    L = []
    for sec in list(seclist):
        L.append(get_section_points(sec))
    return np.array(flattenLL(L))


def get_cell_points(cell):
    '''
    Returns the 3D points that make up the morphology of
    the given neuron.
    Each element of the returned array is an array representing
    the points of a single section.
    '''
    
    L = []
    for sec in cell.allsec():
        L.append(get_section_points(sec))
        
    return np.array(L)


def check_topology():
    isolated_sections = []
    for sec in h.allsec():
        if (sec.parentseg() is None) and (len(sec.children()) == 0):
            print "WARNING: section " + sec.name() \
                  + " has neither parent nor children!"
            isolated_sections.append(sec)
    return isolated_sections


def get_terminals(cell=None, section_types=[]):
    terminals = []

    if isinstance(section_types, str):
        section_types = [section_types]
    if not isinstance(section_types, list):
        raise TypeError("'section_types' must be a list"
                        + "of strings, or a single string.")
    
    if cell is None:
        allsections = h.allsec()
    else:
        allsections = cell.allsec()
    
    for sec in allsections:
        # for a segment to be a terminal,
        # 1) its section must have no children, and
        # 2) it must be the last segment of the section.
        if len(sec.children()) == 0:
            if len(section_types) != 0:
                for section_type in section_types:
                    if section_type in sec.__str__():
                        terminals.append(list(sec)[-1])
            else:
                terminals.append(list(sec)[-1])
    return terminals


def traverse_topology(section, section_list=None):
    ''' Returns a list of the sections in the tree starting at the given root section. '''

    # prevents 'section_list' from accumulating from one base call to another
    if section_list is None:
        section_list = []
    
    section_list.append(section)
    
    # traverse children of current section
    for child_section in section.children():
        if child_section != None:
            section_list = traverse_topology(child_section, section_list)
    
    return section_list


def get_last_axon_branches(cell, section_types=[]):
    ''' Returns a list of the sections that make up the last branches of the morphology's tree,
        i.e., branches that do not bifurcate and contain a terminal. '''
    
    LAB_list = [] # last axon branches list
    sec_temp_list = []

    sec_ordered_by_tree = traverse_topology(cell.soma[0])

    for sec in sec_ordered_by_tree:

        if (sec.name().split('[')[0] in section_types) or (len(section_types) == 0):

            sec_temp_list.append(sec)

            # if sec has 2 or more children (i.e., is branching):
            if len(sec.children()) >= 2:
                sec_temp_list = []

            # if sec has 0 children (i.e., is a terminal):
            if len(sec.children()) == 0:
                LAB_list.append(sec_temp_list)
                sec_temp_list = []

    return LAB_list


def get_seclist_ends(seclists):
    ''' Get first and last points for each of the section lists.
        The sections are assumed to be connected in the order
        given by the list.'''
    
    all_seclist_ends = []
    for seclist in seclists:
        init_pt = None; final_pt = None
      
        for sec in seclist:
            if len(sec.children()) == 0:
              init_pt = get_section_points(sec)[0]
            if sec.parentseg().sec.__str__() != sec.name():
              final_pt = get_section_points(sec)[-1]
              
        if init_pt is None or final_pt is None:
            raise ValueError("no endpoints were found for this section. ")

        seclist_ends = [init_pt, final_pt]
        all_seclist_ends.append(seclist_ends)
        
    return all_seclist_ends


def load_obj_mesh(filename):
    """ Loads only vertices coordinates and face ids from obj file.
        Materials, vertex normals, texture coordinates, etc are ignored.

        Note that if you want to use the returned ids in 'lines' and 'faces' as
        indices in numpy arrays or python lists, you must subtract 1 from all indices,
        because the indices in OBJ file format start at 1, and not at 0.
        So do this to the returned lines (and faces):
            nm0['lines'] = np.array(nm0['lines'])-1
    """
    
    vertices = []
    edges = []
    faces = []
    
    with open(filename, 'rb') as f:
        for line in f:
            line = line.replace('\n','').replace('\r','').split(' ')

            if line[0] == 'v':
                line.remove(line[0])
                vertex_coords = [float(coord) for coord in line]
                vertices.append(vertex_coords)

            elif line[0] == 'l':
                line.remove(line[0])
                vertices_ids = [int(vertex_id) for vertex_id in line]
                edges.append(vertices_ids)

            elif line[0] == 'f':
                line.remove(line[0])

                face_vertices = []

                # remove empty strings
                line = [x for x in line if x!='']
                
                # extract only vertices ids
                for elem in line:
                    vertex_id = int(elem.split('/')[0])
                    face_vertices.append(vertex_id)

                faces.append(face_vertices)

            else:
                pass # ignore vt, vn, etc.

    return {'vertices': vertices, 'edges': edges, 'faces': faces}


