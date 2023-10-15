import pickle
import os
import time
import datetime
import sys
sys.path.append('..') # ROOT_DIR

from mpi4py import MPI

import math
import numpy as np
import matplotlib.pyplot as plt

# tms_networks
from .parameters import *
from .core.codes.save_data import *
from .core.codes.neuron_routines import *
import tms_networks.core.codes.linear_alg_routines as linar

# ------------------------------------------
# debug functions

def flatten_dict_list(list0):
    final_dict = {}
    for dict0 in list0:
        for key, value in dict0.items():
            final_dict[key] = value
    return final_dict

flatten_list = lambda l: [item for sublist in l for item in sublist]

def write_log(log_name, params):

    file_path = ROOT_DIR + '/simulation_temp/logs/' + log_name + '.txt'
    file = open(file_path, 'a+')
    
    if log_name == 'time':
        # params: made up by name of code section, t_start, t_end
        file.write(params[0] + ' - time (h:min:s): {}'.\
            format(str(datetime.timedelta(seconds=params[1]-params[2]))))
    else:
        log_string = ''
        for item in params:
            
            log_string += str(params) + ' '
        
        file.write(log_string)
    
    file.write('\n')
    file.close()

def write_file(string, filename):
    f = open(filename, 'a')
    f.write(string)
    f.close()

def dump_cells_vmem(bcl, tms_wf):
    '''
    Save membrane potential matrix for each cell.
    This function will be called by each of the nodes.
    '''

    rank = MPI.COMM_WORLD.Get_rank()

    # save current waveform
    if rank == 0:
        dat_path = ROOT_DIR + '/vmem_results/pulse_waveform.dat'
        with open(dat_path, 'w+') as wf_file:
            pickle.dump(tms_wf, wf_file)

        fig_path = ROOT_DIR + '/vmem_results/pulse_waveform.pdf'
        plt.plot(np.arange(len(tms_wf)), tms_wf)
        plt.title('Current waveform')
        plt.savefig(fig_path, dpi=600)
        print('saved waveform: ', fig_path)

    # save vmem of all cells (each vmem: ~50MB)
    try:
        for i in range(len(bcl)):
            bcl[i]._collect_vmem()
            cell_layer = bcl[i].cell_name.split('_')[0]
            
            if cell_layer == 'L5':
                filename = 'L5_vmem.dat'
            else:
                filename = cell_layer + '_vmem_node' + str(rank) \
                    + '_cell' + str(i) + '.dat'
            
            with open(ROOT_DIR + '/vmem_results/' + filename, 'w+') as vmem_file:
                pickle.dump(bcl[i].vmem, vmem_file)
                print('vmem was dumped: ', bcl[i].cell_name, filename)
            
            # save vmem figure
            plt.clf()
            plt.imshow(bcl[i].vmem)
            #ax = fig.get_axes()
            plt.ylabel('Segment index')
            plt.xlabel('Time steps')
            plt.clim(-175., 60.)
            plt.colorbar()
            #ax.set_xticks([0,2,4,6])
            plt.title('Membrane potential at segments')
            fig_path = ROOT_DIR + '/vmem_results/' + filename[:-4] + '.pdf'
            plt.savefig(fig_path)
            print('saved fig: ', fig_path)
    except Exception as e:
        print('could not save vmem: ', str(e))

def dump_connections(net):
    '''
    write netcons to file
    '''
    
    write_log('connections', ['edge populations: '])
    for edge_pop in net._edge_populations:
        write_log('connections', [edge_pop.name + '\t' +
                  'number of edges:' + \
                  str(len(list(edge_pop.get_edges()))) + '\n'])

    rank = MPI.COMM_WORLD.Get_rank()
    bcs = net.get_local_cells()
##    string = ''

    # for creating obj representing connections between neurons of the network
    conn_vertices = []
    conn_lines = []
    pt_idx = 1 # obj indices start at 1

    for cell_idx in bcs.keys():
##            string += '\n' + '===== cell_idx: ' + str(cell_idx) + '; ' + \
##                str(bcs[cell_idx].morphology_file) + '\n'           

        for nc in bcs[cell_idx].netcons:

            if nc.preseg() is None or nc.postseg() is None:
                pass
            else:
                approx_preseg_pt = get_section_points(
                    nc.preseg().sec)[-1]
                approx_postseg_pt = linar.centroid(
                    get_section_points(nc.postseg().sec))

                conn_vertices.append(approx_preseg_pt)
                conn_vertices.append(approx_postseg_pt)
                conn_lines.append([pt_idx, pt_idx+1])
                pt_idx += 2

##                string += 'precell' + str(nc.precell()) + ' ; preseg ' + \
##                    str(nc.preseg()) + '\t->\t' + 'postseg: ' + \
##                    str(nc.postseg()) + '\n'

    conn_vertices = np.array(conn_vertices)
    conn_lines = np.array(conn_lines)

    save_obj(conn_vertices, conn_lines, open(ROOT_DIR
        +'/simulation_temp/connections/connections'
        +str(rank)+'.obj', 'w+'))

##        netcons_text_file = os.path.abspath(ROOT_DIR + '/simulation_temp'+ \
##                                            '/netcons_node0.txt')
##        write_file(string=string, filename=netcons_text_file)
##        write_log('connections',
##                  ['netcons_node0.txt dumped: ', netcons_text_file])

def dump_morphologies(cell_list):
    # for storing the longitudinal axes of cells (i.e., the axis
    # pointing from the soma to the centroid of the 3D points of
    # the apical dendrites).
    cells_axes = []

    soma_centers = []

    # dump obj files representing morphologies
    for cell in cell_list:
        save_neuron_morphology(cell, 'obj',
            save_path=ROOT_DIR+'/simulation_temp/morphologies')
        cells_axes.append(cell.somatodend_axis())
        soma_centers.append(cell.somapos)

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    cells_axes = comm.gather(cells_axes, root=0)
    soma_centers = comm.gather(soma_centers, root=0)

    # dump cell axes
    if rank == 0:
        cells_axes_f, soma_centers_f = map(flatten_list,
                                           [cells_axes, soma_centers])

        pickle.dump(cells_axes_f, open(ROOT_DIR+
            '/simulation_temp/morphologies/cells_axes.dat', 'w+'))
        pickle.dump(soma_centers_f, open(ROOT_DIR+
            '/simulation_temp/morphologies/soma_centers.dat', 'w+'))


def set_recorders(bcl):    
    for i in range(len(bcl)):
        bcl[i]._set_voltage_recorders()
        write_log('vmem',
                  ['voltage recorders were set for ' \
                   + str(bcl[i].cell_name) + ': bcl[' \
                   + str(i) + ']'])


def check_ncs(ncs):
        search_str = raw_input('Type the search string for precell (e.g., "[22]"): ')
        if len(search_str) == 0:
                search_str = '[22]'
        print "search_str: ", search_str
        
        for i in range(len(ncs)):
                if search_str in str(ncs[i].precell()):
                        print 'found : ', i

def check_bcl():
    ''' for checking blue brain cells, not StimCell. '''
    
    for i in range(len(bcl)):
        print 
        print " ------ bcl ", i
        print len(list(bcl[i].all))

        for sec in bcl[i].all:
            print sec.n3d(),

    for i in range(len(bcl)):
        print "bcl ", i, ": \t", bcl[i].soma[0].x3d(0), bcl[i].soma[0].y3d(0), bcl[i].soma[0].z3d(0)

def check_scells(cl):
    """ for checking StimCells. """
    from .core.neuron_routines import *

    print " \n number of sections for each cell:"
    for i in range(len(cl)):
        print ""
        print " ------ cl ", i, "\t", "n sections: ", len(list(cl[i].cell.all))
    
    print ' \n n3d for each cell: '
    for i in range(len(cl)):
        print 'cell ', i, ': '
        for sec in cl[i].cell.all:
            print sec.n3d(),

    raw_input('press enter to continue.')
    
    print ' \n all somapos: '
    for c0 in bcl:
        print c0.somapos
    
    print ' \n first point of each cell: '
    for c0 in cl:
        print get_cell_points(c0.cell)[0][0]

def check_cm(cl):
	i=0
	for s0 in cl:
	    i=i+1
	    print ""
	    print " cell ", i, ':'

	    nsec = 0
	    for sec in s0.cell.all:
	        nsec=nsec+1
	        print sec.cm,

chkcm = lambda cell: [sec.cm for sec in cell.all]

def check_nodes(graph, node_pop_name):
        print "nodes list: "
        nodes = list(graph.get_node_population(node_pop_name).get_nodes())
        for node in nodes: print node.node_id, '\t', node.model_template, '\t', node.position

        print "edge populations: "
        for edge_pop in graph._edge_populations: print edge_pop.name, '\t', 'number of edges:', len(list(edge_pop.get_edges()))
        print "\n"

        print "edges list: "
        edges = list(graph._edge_populations[0].get_edges())
        for edge, i in zip(edges,range(len(edges))): print i, '\t', edge.source_population, edge.source_node_id, '\t\t -> \t', edge.target_population, edge.target_node_id, edge.target_sections

def print_syn_table(net):
    ''' net: NetworkBuilder instance. '''
        
    syn_table = list(net._DenseNetwork__edges_tables)[0]['syn_table']
    nsyns = np.sum(syn_table.nsyn_table)
    print '\n nsyn table: '
    print syn_table.nsyn_table
    print 'nsyns: ', nsyns
    print '\n'

def plot_cells(cl, show_cortex=False):
    """ cl: cell list; list of StimCell objects. """
    
    import my_helper_code.mayavi_functions as myv; from mayavi import mlab
    import pymesh as pm
    from .core.list_routines import *
    
    l23_pts = []; l5_pts = []
    for c0 in cl:
        centers = c0.calculate_segments_centers(flatten=True)
        if 'L23' in c0.cell_name:
            l23_pts.append(centers)
        elif 'L5' in c0.cell_name:
            l5_pts.append(centers)
        else:
            print 'cell pt3ds not appended.'
    
    if show_cortex:
            #debug
            resize_scale = 1e3
            gm_mesh = pm.load_mesh('/home/master/Desktop/TMSEffects/core/data/simnibs_simulation/meshes/motor_cortex_meshes/gm_motor.obj')
            rescaled_vertices = resize_scale * gm_mesh.vertices
            gm_mesh = pm.form_mesh(rescaled_vertices, gm_mesh.faces)
            wm_mesh = pm.load_mesh('/home/master/Desktop/TMSEffects/core/data/simnibs_simulation/meshes/motor_cortex_meshes/wm_motor.obj')
            rescaled_vertices = resize_scale * wm_mesh.vertices
            wm_mesh = pm.form_mesh(rescaled_vertices, wm_mesh.faces)
    
    myv.draw_points(flattenLL(l23_pts), scl=100.0, color=(1,0,0))
    myv.draw_points(flattenLL(l5_pts), scl=100.0, color=(0,0,1))

    if show_cortex:
        myv.draw_tri_mesh(gm_mesh, opacity=0.1)
        myv.draw_tri_mesh(wm_mesh, opacity=0.1)
    
    mlab.show()

# debugging

def print_network_info(network):
    print "network nnodes:", network.nnodes, '\t', 'built: ', network.nodes_built
    print "network nedges:", network.nedges, '\t', 'built: ', network.edges_built
    
    print "network edges table: "
    print network.edges_table()
    
    print "network.edges_table()[0]  nsyn table: "
    print network.edges_table()[0]['syn_table'].nsyn_table
    
    print "network 1st node: ", list(network.nodes())[0]
##    print "network 1st edge: ", network.edges()[0]
    
