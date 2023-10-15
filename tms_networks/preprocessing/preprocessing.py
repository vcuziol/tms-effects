import os
import sys

import numpy as np
import scipy.io
import pickle

try:
    import pymesh as pm
except ImportError:
    pass

import linear_alg_routines as linar

# simnibs codes
sys.path.append(os.path.abspath('../core/libraries'))
import simnibs.msh.gmsh_numpy as gmsh_numpy
import simnibs.msh.surface as surface

sys.path.append(os.path.abspath('../core'))
from parameters import *

# ----------------- FUNCTIONS
def print_magnitudes():
    magnitudes_E = []
    for i in range(len(gm_motor_E)):
        magnitudes_E.append(np.linalg.norm( gm_motor_E[i] ))
    magnitudes_E = np.array(magnitudes_E)

    print "Magnitudes of the normal vectors:"
    print "np.mean(magnitudes): ", np.mean(magnitudes_E)
    print "np.std(magnitudes): ", np.std(magnitudes_E)

def create_mesh_geom(fname, resize_scale):
    original_mesh = pm.load_mesh(fname)
    rescaled_vertices = resize_scale * original_mesh.vertices
    mesh = pm.form_mesh(rescaled_vertices, original_mesh.faces)

    mesh.add_attribute("face_centroid")
    face_centers = mesh.get_attribute("face_centroid").reshape(mesh.num_faces, 3)
    
    mesh.add_attribute("face_normal")
    # by default, the normals point inside the mesh.
    face_normals = mesh.get_attribute("face_normal").reshape(mesh.num_faces, 3)
    
    mesh_geom = {'nodes': mesh.vertices,
                 'trnodes': mesh.faces,
                 'centers': face_centers,
                 'normals': face_normals}

    return mesh_geom

def calculate_faces_centers(pm_mesh):
    pm_mesh.add_attribute("face_centroid")
    face_centers = pm_mesh.get_attribute("face_centroid").reshape(pm_mesh.num_faces, 3)
    return face_centers

# ----------------------------------------------------------------------------------------------------------------
def get_pathfem(mat_file_path):
    return scipy.io.loadmat(mat_file_path)['pathfem']

def read_coil_matrix(simnibs_mat_fn, show=False):
    ''' Reads coil position matrix from a .mat file generated by SimNIBS. '''

    # ------- reading from file

    mdict = scipy.io.loadmat(simnibs_mat_fn)

    coil_position_matrix = mdict['poslist'][0][0][0][0][4][0][0][6]

    # x,y,z axes in homogeneous coordinates
    origin_and_axes = np.concatenate((np.zeros((3,1)), np.identity(3)), axis=1)
    origin_and_axes = np.concatenate((origin_and_axes, np.ones((1,4))), axis=0)

    # gets the coil center and axes, as visualized in simnibs.
    # the coil 'origin' (or position) is the first column of 'coil_points'.
    coil_points = np.dot(coil_position_matrix, origin_and_axes)

    # ------- plot
    if show:
        fig = plt.figure(); ax = Axes3D(fig)
        cpts=coil_points

        #v1 = [cpts[:,0][:-1], cpts[:,1][:-1]]; v2 = [cpts[:,0][:-1], cpts[:,2][:-1]]; v3 = [cpts[:,0][:-1], cpts[:,3][:-1]]
        v1 = cpts[:,1][:-1] - cpts[:,0][:-1]
        v2 = cpts[:,2][:-1] - cpts[:,0][:-1]
        v3 = cpts[:,3][:-1] - cpts[:,0][:-1]

        cpts = coil_points[:-1].T

        draw3dscatter(ax, cpts)

        for v, c in zip([v1,v2,v3], ['r','g','b']):
            draw_3d_one_quiver(ax, cpts[0], v, color=c, length=1.0)

        plt.title("coil orientation as given by the matrix")
        plt.show()

    return {'coil_points':coil_points, 'coil_position_matrix':coil_position_matrix}

def write_coil_matrix(new_coil_matrix, old_simnibs_mat_fn, new_simnibs_mat_fn=None, new_pathfem_suffix=None):
    ''' Given a SimNIBS .mat file, creates a copy of the file (or overwrites it) but
        with the coil position matrix overwritten with the matrix given as parameter.
        If only 'old_simnibs_mat_fn' is given, the file is overwritten.
        If both 'old_simnibs_mat_fn' and 'new_simnibs_mat_fn' are given, makes a copy of the old file.

        'new_pathfem_suffix' is a suffix to be appended to the path of the simulation results folder, so that
        each .mat file, after run in SimNIBS, generates the results in a different folder. '''

    if new_simnibs_mat_fn is None:
        new_simnibs_mat_fn = old_simnibs_mat_fn
    
    assert(new_coil_matrix.shape == (4,4))
    
    import scipy.io
    
    mdict = scipy.io.loadmat(old_simnibs_mat_fn)
    
    mdict['poslist'][0][0][0][0][4][0][0][6] = new_coil_matrix
    
    if new_pathfem_suffix is not None:
        mdict['pathfem'] = np.array([mdict['pathfem'][0] + new_pathfem_suffix])
    
    scipy.io.savemat(new_simnibs_mat_fn, mdict)
    
    return 0
    
# for debugging
def check_coil_axes(folder=os.getcwd()):
    ''' Given a folder with SimNIBS .mat files, prints info about the coil axes of each .mat file.
        The angles printed should be always close to 90 degrees and the vector printed (which is the z-axis of the coil) should not change. '''
    
    matfiles=filter(lambda x: x[-3:]=='mat', os.listdir(folder))
    mdic=read_coil_matrix(mat_fn); coilM = mdic['coil_position_matrix']; v1=coilM[:,0][:-1]; v2=coilM[:,1][:-1]; v3=coilM[:,2][:-1]; print np.rad2deg(ang(v1,v2)),np.rad2deg(ang(v2,v3)), np.rad2deg(ang(v1,v3)),'\n'; print v3, '\n\n'
    
# ----------------------------------------------------------------------------------------------------------------
def extract_simnibs_msh(simnibs_session_name, SURFACE_TYPE):
    # --- reading files
    # wm: white matter; gm: gray matter
    surface_code_dict = {'wm':1001, 'gm':1002}
    volume_code_dict = {'wm':1, 'gm':2}
    
    surface_code = surface_code_dict[SURFACE_TYPE]
    volume_code = volume_code_dict[SURFACE_TYPE]
    
    # selects only the first .msh file
    msh_fn = filter(lambda string: (string[-4:] == '.msh'), os.listdir(MSH_INPUTS_FOLDER + '/' + simnibs_session_name))[0]
    msh_fn = MSH_INPUTS_FOLDER + '/' + simnibs_session_name + '/' + msh_fn
    
    # this line gives an error on Windows 10 ("TypeError: ord() expected a character, but string of length 0 found"). run on linux instead.
    mesh_struct = gmsh_numpy.read_msh(msh_fn)
    
    # read surface
    m_surf = surface.Surface(mesh_struct, [volume_code, surface_code])
    
    print "mesh_struct.elmdata: ", type(mesh_struct.elmdata)
    
    m_elmdata_normE = []; m_elmdata_E = []; m_elmdata_v = []
    # get element data ('normE', 'E' and 'v') for the chosen surface
    for i in range(len(mesh_struct.elmdata)):
        if mesh_struct.elmdata[i].field_name == 'normE':
                m_elmdata_normE = mesh_struct.elmdata[i].value[mesh_struct.elm.tag1 == surface_code]
        elif mesh_struct.elmdata[i].field_name == 'E':
                m_elmdata_E = mesh_struct.elmdata[i].value[mesh_struct.elm.tag1 == surface_code]
        elif mesh_struct.elmdata[i].field_name == 'v':
                m_elmdata_v = mesh_struct.elmdata[i].value[mesh_struct.elm.tag1 == surface_code]
    
    # --- saving the data using pickle
    session_data_folder = ROOT_FOLDER + '/core/data/' + simnibs_session_name + '/'
    meshes_folder_name = session_data_folder + 'meshes/' + SURFACE_TYPE + '_files/'
    fields_folder_name = session_data_folder + 'fields/' + SURFACE_TYPE + '_files/'
    
    for folder in [meshes_folder_name, fields_folder_name]:
        if not os.path.isdir(folder):
            os.mkdir(folder)
            print "created folder: ", folder
    
    print "Writing files..."
    
    pickle.dump(m_elmdata_v, open(fields_folder_name + SURFACE_TYPE + "_elmdata_v.dat","w+"))
    pickle.dump(m_elmdata_E, open(fields_folder_name + SURFACE_TYPE + "_elmdata_E.dat","w+"))
    pickle.dump(m_elmdata_normE, open(fields_folder_name + SURFACE_TYPE + "_elmdata_normE.dat","w+"))
    
    pickle.dump(m_surf.tr_nodes, open(meshes_folder_name + SURFACE_TYPE + "_surf-tr_nodes.dat","w+"))
    pickle.dump(m_surf.nodes, open(meshes_folder_name + SURFACE_TYPE + "_surf-nodes.dat","w+"))
    
    m_mesh = pm.form_mesh(m_surf.nodes, m_surf.tr_nodes)
    pm.save_mesh(meshes_folder_name + SURFACE_TYPE + ".obj", m_mesh)

    print "Fields files were written. "
    
    return 0
    
# ----------------------------------------------------------------------------------------------------------------
def generate_neurons_positions_with_relative_depth(wm_pm_mesh, gm_centers, gm_normals, relative_depth, max_distance=5000.0, remove_nans=True):
    """ Given a white matter mesh (PyMesh Mesh object), the centers and normals of the polygons of a gray matter mesh,
       and a parameter denoting relative depth of a cortical layer (i.e., the percentage of the cortical thickness),
       creates points indicating centers of the neurons' somas.
        
       'relative_depth': depth of the neuron relative to the total cortical thickness in the cortex point.
       'wm_pm_mesh': White matter PyMesh Mesh object.
        
       If 'remove_nans' is True, the polygons of the GM for which no intersection is found are simply removed from the final positions list,
       instead of being replaced with a new position calculated using the median of the cortical depths. It is useful activating this because
       using the median could introduce artifacts, like having a few neuron positions above the cortex surface.
       """
    
    from vtk import *
    import vtk_functions as vtkf
    
    wm_centers = calculate_faces_centers(wm_pm_mesh) # wm polygon centers
    mesh_verts_range = np.round(np.abs(np.max(wm_pm_mesh.vertices)) + np.abs(np.min(wm_pm_mesh.vertices)))
    
    pd = vtkf.form_polydata(wm_pm_mesh.vertices, wm_pm_mesh.faces)
    
    # vtk object for finding points of intersection between a mesh and a line.
    obbTree = vtkOBBTree()
    obbTree.SetDataSet(pd)
    obbTree.BuildLocator()
    
    neurons_positions = []
    cortical_depths = []
    neurons_indices = []
    
    for i in range(gm_centers.shape[0]):
        
        ray_origin = gm_centers[i]
        ray_direction = gm_normals[i]
        ray_endpoint = ray_origin + mesh_verts_range * ray_direction
        
        # find intersection of the prolonged normal vector with the white matter mesh,
        # to calculate cortical depth (a.k.a. cortical thickness).
        isecpts = vtkPoints()
        obbTree.IntersectWithLine(ray_origin, ray_endpoint, isecpts, None)
        intersection_points = np.array(vtkf.vtkpoints_to_numpy(isecpts))
        
        if len(intersection_points) == 0:
            pos = None
        elif intersection_points[0].size == 0:
            pos = None
        else:
            # gets first intersection point
            intersection_point = np.array(intersection_points[0])
            
            # distance between 'gm_centers[i]' and 'intersection_point'.
            cortical_depth = np.linalg.norm(gm_centers[i] - intersection_point)
            
            # since the relative depth is a percentage (between 0 and 1) and the cortical depth is in micrometers,
            # the result will be the distance from the cortex surface, in micrometers.
            distance_to_surface = relative_depth * cortical_depth
            
            # point plus vector means the point is displaced to the 'tip' of the vector.
            pos = gm_centers[i] + distance_to_surface * linar.normalize_vector(gm_normals[i])

            assert(linar.collinear(pos, gm_centers[i], gm_centers[i]+gm_normals[i]))
            
        if pos is not None:
            neurons_indices.append(i)
            neurons_positions.append(pos)
            cortical_depths.append(cortical_depth)
        else:
            if not remove_nans:
                neurons_positions.append(np.array([np.nan,np.nan,np.nan]))
                cortical_depths.append(np.nan)
    
    if max_distance is None:
        max_distance = np.inf
    
    # --- replace outliers
    neurons_positions = np.array(neurons_positions)
    cortical_depths = np.array(cortical_depths)
    
    # get indices of values to be replaced (either nans or values over the given max)
    
    if remove_nans:
        outliers_mask = cortical_depths>max_distance
    else:
        overmax_mask = cortical_depths>max_distance
        nan_mask = np.isnan(cortical_depths)
        outliers_mask = np.logical_or(overmax_mask, nan_mask)
    outliers_indices = np.arange(len(cortical_depths))[outliers_mask]
    
    cortical_depths_without_outliers = cortical_depths[~outliers_mask]
    cortical_depth_replace_value = np.median(cortical_depths_without_outliers)
    
    neurons_positions[outliers_indices] = gm_centers[outliers_indices] + (relative_depth * cortical_depth_replace_value) * gm_normals[outliers_indices]
    cortical_depths[outliers_indices] = cortical_depth_replace_value
    
    assert(sum(cortical_depths>max_distance) == 0)
    
    return {'neurons_positions':neurons_positions,
            'cortical_depths':cortical_depths,
            'neurons_indices':neurons_indices
            }

def generate_neurons_positions_with_constant_depth(pm_mesh, scale_factor, centers, normals):
    """ Given a mesh (PyMesh Mesh object), the centers and normals of its polygons,
     and a parameter denoting depth from the surface, creates points indicating centers of neurons.
     The "scale factor" adjusts the distance from each neuron to the cortical surface.
     (each neuron will have the same distance to the cortical surface.) """
    
    neurons_positions = []
    for i in range(pm_mesh.num_faces):
        # point plus vector: the point is displaced to the 'tip' of the vector.
        pos = centers[i] + float(scale_factor)*normals[i]
        neurons_positions.append(pos)

    return neurons_positions

# ========================================================================================================================
def generate_mat_files(simnibs_mat_fn, number_of_coil_orientations, show_coil_basis):
    ''' Code for generating .mat files (to be used as input to SimNIBS simulations) for a given number of coil orientations.
        The necessary files are generated from a given base .mat SimNIBS file. '''

    coil_dict = read_coil_matrix(simnibs_mat_fn, show=show_coil_basis)
    
    # 'coil_points' refers to 4 points: the coil position ('coil_origin') and the tips of 3 vectors representing the coil axes ('coil_basis').
    # in 'coil_points', they are stored in world space; i.e., the tips of the coil axes are considered as points measured from the global origin.
    # in 'coil_position_matrix', they are stored as usual vectors; i.e., as unitary vectors meant to be measured from the 'coil_origin' point.
    # (that is, if we add 'coil_origin' to the vectors in 'coil_position_matrix', we get the vectors in 'coil_points'.)
    coil_origin = coil_dict['coil_points'][:,0][:-1] 
    coil_basis = coil_dict['coil_position_matrix'][:-1,:-1] # each column is an axis.
    
    # -------------- creates simnibs files with different coil orientations

    angles = np.linspace(0, 2.0*np.pi, number_of_coil_orientations+2)
    angles = angles[1:-1] # removes angles 0 and 360 degrees

    # list of y-axis of all calculated coil bases, because the coil's handle is pointed in the direction of the y-axis
    # (with the handle-coil direction along the +y axis.)
    y_axis_list = []
    y_axis_list.append(coil_basis[:,1])
    
    for angle in angles:
        new_coil_basis = linar.rotate_ccw(coil_basis, angle, coil_basis[:,2])
        
        y_axis_list.append(new_coil_basis[:,1])
        
        new_coil_position_matrix = linar.add_matrix_homog(np.concatenate((new_coil_basis, coil_origin.reshape((3,1))), axis=1))

        #debug
##        gm_mesh = "../core/data/simnibs_simulation/meshes/gm_files/gm.obj"; import linear_alg_routines as linar; ang=linar.angle_vectors; from helper import *
##        v1=new_coil_basis[:,0]; v2=new_coil_basis[:,1]; v3=new_coil_basis[:,2];
##        print ang(v1,v2), ang(v1,v3), ang(v2,v3)
##        print new_coil_position_matrix[:,2][:-1], '\n', coil_basis[:,2]
        
        new_fn = simnibs_mat_fn[:-4] + '_ang' + str(np.round(np.rad2deg(angle))) + '.mat'
        result = write_coil_matrix(new_coil_position_matrix,
                                   old_simnibs_mat_fn = simnibs_mat_fn,
                                   new_simnibs_mat_fn = new_fn,
                                   new_pathfem_suffix = '_coilAngle' + str(np.round(np.rad2deg(angle),2)))
        
        if result == 0:
            print "Successfully wrote new coil matrix at: ", new_fn

    return angles

# ----------------------------------------------------------------------------------------------------------------
def run_simnibs(simnibs_command, simnibs_mat_fn, angles):
    # runs simulations with all created files
    from subprocess import call
    
    # for each generated .mat file, run SimNIBS simulation (serially)
    for angle in angles:
        mat_fn = simnibs_mat_fn[:-4] + '_ang' + str(np.round(np.rad2deg(angle))) + '.mat'
        print "\nRunning SimNIBS with parameter: ", mat_fn
        print "..."
        call([simnibs_command, mat_fn])
        print "(OK) Simulation finished. Results saved at: ", get_pathfem(mat_fn)

# ----------------------------------------------------------------------------------------------------------------
def preprocessing(simnibs_session_name, relative_depths):
    root_folder = os.getcwd()
    os.chdir(root_folder)

    print "----- preprocessing files for: ", simnibs_session_name
    print "Working at", root_folder, "..."

    # mesh processing parameters:
    # the original gray matter mesh is in milimeters (mm, the same unit as the original nifti file),
    # and we want it to be micrometers; so rescale of 1000.
    CORTICAL_MESH_RESCALE = 1e3

    # --- reads simnibs files and writes meshes and fields
    
    # read simnibs msh and create files
    extract_simnibs_msh(simnibs_session_name, 'gm')
    extract_simnibs_msh(simnibs_session_name, 'wm')

    print 'GM/WM meshes and fields written.'
    
    # --- loads the submesh of the motor cortex (region of interest)
    # TODO: create MC (motor cortex) indices for WM
    session_data_folder = ROOT_FOLDER + '/core/data/' + simnibs_session_name + '/'
    mc_base_folder = ROOT_FOLDER + '/preprocessing/motor_cortex_indices/'
    meshes_base_folder = session_data_folder + 'meshes/'
    fields_folder = session_data_folder + 'fields/'
    
    gm_fn = 'gm.obj'
    wm_fn = 'wm.obj'
    
    fn_dic = {'gm': gm_fn, 'wm': wm_fn}
    
    # extract data (E-field and potential v) relative to motor cortex submesh (for gm and wm)
    for tissue_type in fn_dic:
        tissue_field_folder = fields_folder + tissue_type + "_files/" + tissue_type
        
        # path to GM or WM obj file
        obj_mesh_filename = meshes_base_folder + tissue_type + '_files/' + fn_dic[tissue_type]
        
        # load whole GM/WM mesh and create submesh of GM/WM motor cortex
        m_mesh = pm.load_mesh(obj_mesh_filename)
        
        # reads file containing the polygon indices of the motor cortex (m.c.) submesh
        motor_indices = open(mc_base_folder + tissue_type + '_motor_indices.txt', 'rb').readline().split(',')
        motor_indices = map(int, motor_indices) # if python3, change to list(map(int, motor_indices))
        
        motor_mesh = pm.submesh(m_mesh, motor_indices, 0)

        if tissue_type == 'gm':
            if not os.path.isdir(meshes_base_folder + 'motor_cortex_meshes/'):
                os.mkdir(meshes_base_folder + 'motor_cortex_meshes/')
                print "created folder: ", meshes_base_folder + 'motor_cortex_meshes/'
        
##        if tissue_type == 'gm':
##            # generates ids for "lesser" GM motor cortex (MC) mesh, to be used as basis for population ("motor_indices");
##            # and ids for "greater" MC mesh, to be used for interpolation (since it must contain an area greater than the "lesser" MC's for better interpolation).
##            motor_mesh_for_interp = motor_mesh
##            motor_mesh, faces_kept_bool = helper_meshproc.remove_boundary(motor_mesh_for_interp) # removes vertices from the mesh boundary.
##            
##            motor_indices_for_interp = motor_indices
##
##            print "motor_mesh_for_interp.faces shape: ", motor_mesh_for_interp.faces.shape
##            
##            if not os.path.isdir(meshes_base_folder + 'motor_cortex_meshes/'):
##                    os.mkdir(meshes_base_folder + 'motor_cortex_meshes/')
##
##            pm.save_mesh(meshes_base_folder + 'motor_cortex_meshes/' + tissue_type + '_motor_for_interp.obj', motor_mesh_for_interp)
##            motor_E_for_interp = pickle.load(open(tissue_field_folder + '_elmdata_E.dat', 'rb'))[motor_indices_for_interp]
##            pickle.dump(motor_E_for_interp, open(tissue_field_folder + '_motor_E_for_interp.dat','wb'))
##            
##            motor_indices = np.array(motor_indices)
##            motor_indices = motor_indices[faces_kept_bool] # TODO: check if faces ids are conserved
##
##            print "motor_mesh.faces shape: ", motor_mesh.faces.shape
        
        # motor mesh is saved with vertices in millimeters, not micrometers.
        pm.save_mesh(meshes_base_folder + 'motor_cortex_meshes/' + tissue_type + '_motor.obj', motor_mesh)

        
        # select values of E and v for the motor cortex
        motor_E = np.array(pickle.load(open(tissue_field_folder + "_elmdata_E.dat","rb")))[motor_indices]
    

        elmdata_v = pickle.load(open(tissue_field_folder + "_elmdata_v.dat","rb"))
        if len(elmdata_v) > 0:
            motor_v = elmdata_v[motor_indices]
        else:
            motor_v = []
        
        #
        pickle.dump(motor_E, open(tissue_field_folder + "_motor_E.dat",'wb'))
        pickle.dump(motor_v, open(tissue_field_folder + "_motor_v.dat",'wb'))

    # load both motor cortex meshes
    mc_meshes_folder = meshes_base_folder + 'motor_cortex_meshes/'
    motor_meshes = filter(lambda string: (string[-4:] == '.obj'), os.listdir(mc_meshes_folder))

    # generate and save geometry for motor cortex meshes 
    for mesh_fn in motor_meshes:
        # gmmc: gray matter motor cortex
        mesh_geom = create_mesh_geom(mc_meshes_folder + mesh_fn, CORTICAL_MESH_RESCALE)

        geom_fn = session_data_folder + 'meshes/motor_cortex_meshes/' + mesh_fn[:-4] + '_geom.dat'
        pickle.dump(mesh_geom, open(geom_fn, 'wb'))

        print "Saved", mesh_fn, " geometry at: ", geom_fn
        
        # TODO: for wm, swap normals direction?

    # -----------------------------------------------------
    # this code loads a mesh of the gray matter motor cortex (in obj)
    # and calculates neurons positions and orientations according to it.
    # the results are saved as 'populations_parameters'.

    for population_layer_number in relative_depths.keys():

        # -- load GMMC geometry
        gmmc_geom = pickle.load(open(mc_meshes_folder + "gm_motor_geom.dat", "rb"))
        wmmc_geom = pickle.load(open(mc_meshes_folder + "wm_motor_geom.dat", "rb"))
        
        #debug
    ##        pm.save_mesh("gmmc_rescaled.obj", pm.form_mesh(gmmc_geom['nodes'], gmmc_geom['trnodes']))
    ##        pm.save_mesh("wmmc_rescaled.obj", pm.form_mesh(wmmc_geom['nodes'], wmmc_geom['trnodes']))

        # -- create neurons positions
        gmmc = pm.form_mesh(gmmc_geom['nodes'], gmmc_geom['trnodes'])
        wmmc = pm.form_mesh(wmmc_geom['nodes'], wmmc_geom['trnodes'])
        
        neuron_pos_dic = generate_neurons_positions_with_relative_depth(wm_pm_mesh=wmmc,
                                                                        gm_centers=gmmc_geom['centers'],
                                                                        gm_normals=gmmc_geom['normals'],
                                                                        relative_depth=relative_depths[population_layer_number],
                                                                        max_distance=6000.0)

        layer_soma_centers = neuron_pos_dic['neurons_positions']
        layer_neurons_indices = neuron_pos_dic['neurons_indices']
        
        # -- save population parameters
        pickle.dump({"layer_soma_centers":layer_soma_centers, "layer_neuron_distances_factor":relative_depths[population_layer_number]},
                                open(session_data_folder + "population_parameters/" + population_layer_number + "_soma_centers.dat", "w+"))
        pickle.dump(gmmc_geom['normals'][layer_neurons_indices], open(session_data_folder + "population_parameters/" + population_layer_number + "_normals.dat", "w+"))
        pickle.dump(gmmc_geom['centers'][layer_neurons_indices], open(session_data_folder + "population_parameters/" + population_layer_number + "_gm_centers.dat", "w+"))
        pickle.dump(layer_neurons_indices, open(session_data_folder + "population_parameters/" + population_layer_number + "_neurons_indices.dat", "w+"))
        
        print population_layer_number, "population parameters files written."

    print "preprocessing finished. "