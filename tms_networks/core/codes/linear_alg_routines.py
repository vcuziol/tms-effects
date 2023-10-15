import numpy as np
import math

def magnitudes(vectors):
    return np.array([np.linalg.norm(v) for v in vectors])

def distance(a, b):
    a = np.array(a)
    b = np.array(b)
    return np.linalg.norm(b-a)

def centroid(points):
    ''' Calculates the barycenter (a.k.a. 'center of mass' or
        'centroid') of the 3D points given by in an array. '''
    barycenter = [0.0, 0.0, 0.0]
    n = len(points)
    for i in range(n):
        point = points[i]
        barycenter[0] += point[0]
        barycenter[1] += point[1]
        barycenter[2] += point[2]
    barycenter[0] /= n
    barycenter[1] /= n
    barycenter[2] /= n
    
    return np.array(barycenter)

def collinear(p0, p1, p2, epsilon=1e-6):
    # https://stackoverflow.com/questions/9608148/python-script-to-determine-if-x-y-coordinates-are-colinear-getting-some-e
    x1, y1 = p1[0] - p0[0], p1[1] - p0[1]
    x2, y2 = p2[0] - p0[0], p2[1] - p0[1]
    return abs(x1 * y2 - x2 * y1) < epsilon

def angle_vectors(v1, v2):
    """ Returns the angle (in radians) between vectors 'v1' and 'v2'    """
    cosang = np.dot(v1, v2)
    sinang = np.linalg.norm(np.cross(v1, v2))
    return np.arctan2(sinang, cosang)

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.

    https://stackoverflow.com/questions/6802577/rotation-of-3d-vector
    """
    if np.all(axis == np.array([0.0, 0.0, 0.0])):
        raise ValueError("'axis' cannot be a null vector.")
    
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

def ccw_angle_vectors(v1, v2):
    """ Returns the angle between vectors v1 and v2 as
    measured in the counterclockwise direction
    (therefore the return value can range from 0 to [not-inclusive] 2*pi)."""

    # if the two vectors are the same, the angle is 0.
    if np.all(v1 == v2):
        return 0.0
    
    dot = np.dot(v1, v2)
    v_norm = np.linalg.norm(v1) * np.linalg.norm(v2)
    angle_t = np.arccos(dot/v_norm)

    n = np.cross(v1,v2)
    n = n/np.linalg.norm(n)
    tri_prod = np.dot(v1,np.cross(v2,n).T).sum(0)
    
    if tri_prod > 0:
        angle = angle_t
    elif tri_prod < 0:
        angle = 2.0*np.pi - angle_t
    else:
        raise ValueError("tri_prod has an invalid value.")
    
    return angle

def ccw_angle_vectors2(v1,v2):
    """ Instead of returning the angle as measured from v1
    to v2 counterclockwise (and thus always positive),
    returns a positive angle if ccw and negative if cw."""
    dot = np.dot(v1, v2)
    v_norm = np.linalg.norm(v1) * np.linalg.norm(v2)
    angle_t = np.arccos(dot/v_norm)

    n = np.cross(v1,v2)
    n = n/np.linalg.norm(n)
    tri_prod = np.dot(v1,np.cross(v2,n).T).sum(0)
    if tri_prod > 0:
        angle = angle_t
    elif tri_prod < 0:
        angle = -angle_t
    
    return angle

def normalize_vector(v):
    norm = np.linalg.norm(v)
    if norm == 0:
       return v
    return v / norm

def vector_mag(v):
    return np.linalg.norm(v)

def rotate_ccw(points, theta, axis):
        """
        Rotate points around an axis, for a given angle theta.
        The points must be represented by a matrix in which each column is a 3D point.
        """
        
        M_rot = np.array(rotation_matrix(axis, theta))

        return np.dot(M_rot, points)

def random_orthogonal_vector(v1):
    ''' Finds a vector that makes a right angle with the given one. '''
    v2 = v1 + random_point_inside_sphere()
    v3 = np.cross(v1, v2)

    try:
        assert(np.rad2deg(angle_vectors(v1,v3)) - 90.0 < 1e-5)
    except AssertionError:
        print 'assertion failed: resulting angle is not approximately a right angle.'
        import pdb; pdb.set_trace()

    return v3

def random_points_inside_cube(number_of_random_points):
    # https://stackoverflow.com/questions/26895097/generate-a-random-cluster-of-points-around-a-given-point-python
    point = [0.0, 0.0, 0.0]
    deviation = 10.0

    random_points = []
    for i in range(number_of_random_points):
        random_points.append([point[i] + random.random() * deviation for i in range(3)])

    random_points = np.array(random_points)

    return random_points

def random_point_inside_sphere(center=[0, 0, 0], radius=1.0):
    # https://stackoverflow.com/questions/51031077/generate-a-point-within-a-sphere-of-a-certain-radius-thats-around-a-specific-gi
    r = radius * np.random.random()**(1.0/3.0)

    phi = np.random.uniform(0, 2*np.pi) 
    costheta = np.random.uniform(-1.0, 1.0)
    theta = np.arccos(costheta) 

    x = np.sin(theta) * np.cos(phi) 
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)

    point = np.array([x,y,z]) * r

    new_point = center + point
    
    return new_point

def rotate_vectorial(pts, initial_vector, final_vector):
    ''' Rotate geometry of cell object according to the direction
        and angle given by initial_vector and final_vector
        (given that the cell is oriented according to 'initial_vector'.)
    '''

    # no rotation necessary
    if np.all(initial_vector == final_vector):
        return np.array(pts)
    
    theta = ccw_angle_vectors(initial_vector, final_vector)
    M_rot = np.matrix(
            rotation_matrix(-np.cross(initial_vector, final_vector),
                            theta))

    new_pts = pts * M_rot

    return np.array(new_pts)

def add_pts_homog(M):
    ''' Appends a row of ones to a matrix. This is to convert point coordinates (arranged with one point in each column) to homogeneous coordinates. '''
    return np.concatenate((M, np.ones((1,M.shape[1]))), axis=0)

def add_matrix_homog(M):
    ''' Appends a [0,0,0,1] row to a (3,4)-shaped matrix. This is to convert a linear transformation matrix to be used with homogeneous coordinates. '''
    return np.concatenate((M, np.array([0.0,0.0,0.0,1.0]).reshape(1,4)), axis=0)

##def _set_pt3d_rotation_vectorial(initial_vector, final_vector):
##    """ """
##
##    theta = lin.ccw_angle_vectors(initial_vector, final_vector)
##    M_rot = np.matrix(
##        lin.rotation_matrix(np.cross(initial_vector, final_vector),
##                            theta))
##    
##    for i in range(len(self.x3d)):
##        rel_pos = self._rel_pt3d_positions(self.x3d[i],
##                                           self.y3d[i],
##                                           self.z3d[i])
##        
##        rel_pos = rel_pos * M_rot
##        
##        self.x3d[i], self.y3d[i], self.z3d[i] = \
##                                    self._real_pt3d_positions(rel_pos)
##
##    self._update_pt3d()

##def updatecell3d(cell, newpts):
##    point_counter = 0
##
##    for sec in cell.allsec():
##        for i in range(sec.n3d()):
##            sec.
##            
