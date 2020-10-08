import numpy as np 
import igl 
import math
import copy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
class Env : 
    def __init__(self):
        # INPUT VERTEX INFO
        self.max_iteration = 1000
        self.n = None
        self.vertice = None 
        self.vertice_prime = None
        self.faces = None
        self.neighbors = None
        self.edge_matrix = None
        self.cell_rotation = None 

        # CONSTRAINTS

        self.deformation_matrix = None
        self.fixed_vertice = []
        self.selected_vertice = []

        # WEIGHTS
        self.weight_sum = None 
        self.weight_matrix = None



        #pre compute vars
        self.P_i_array = None
        

        # apply deform vars
        self.b_array = None




def init(file_name, env , **kwargs) : 
    def assign_neighbor_matrix(neighbors, v_id1,  v_id2,  v_id3):
        env.neighbors[v_id1,v_id2]=1
        env.neighbors[v_id2,v_id1]=1
        env.neighbors[v_id1,v_id3]=1
        env.neighbors[v_id3,v_id1]=1
        env.neighbors[v_id2,v_id3]=1
        env.neighbors[v_id3,v_id2]=1
        return neighbors


    if "V" in kwargs and "F" in kwargs:
        pass
    else:
        V, F = igl.read_triangle_mesh(file_name) #TODO
    
    
    env.n = len(V)
    env.vertice = V 
    env.vertice_prime = copy.deepcopy(V)
    env.faces = F
    env.neighbors = np.zeros([env.n, env.n])
    env.edge_matrix=np.zeros((env.n,env.n))
    env.cell_rotaion=np.zeros((env.n,env.n,3))

    for i in range(env.n):
        env.verts_to_face.append([])
                      
    for i,face in enumerate(env.faces):
        env.verts_to_face[env.face[0]].append(i)
        env.verts_to_face[env.face[1]].append(i)
        env.verts_to_face[env.face[2]].append(i)
        
        env.neighbors = assign_neighbor_matrix(env.neighbors, face[0],face[1],face[2]) #NxN      neighbor_matrix

    for row in range(env.n):
        env.edge_matrix[row][row]=env.neighbors[row].sum()
    
    return env



def generate_rotation_matrix(radius, axis, pre_rot_mat = None):
    """
        input 
            radius : radius value
            axis : list of axis 3-Dims
        return 
            4X4 matrix
    """
    cos = np.cos(radius)
    sin = np.sin(radius)
    z_rot = np.array([
            [cos, -sin, 0, 0],
            [sin,  cos, 0, 0],
            [0,      0, 1, 0],
            [0,      0, 0, 1]
            ])
    y_rot = np.array([
            [ cos, 0, sin, 0],
            [   0, 1,   0, 0],
            [-sin, 0, cos, 0],
            [   0, 0,   0, 1]
            ])
    x_rot = np.array([
            [1,   0,    0, 0],
            [0, cos, -sin, 0],
            [0, sin,  cos, 0],
            [0,   0,    0, 1]
            ])
    rot = [z_rot, y_rot, x_rot]
    total_rot = np.eye((4,4))

    if pre_rot_mat != None :
        total_rot.dot(pre_rot_mat)
    for i, rot in zip(axis, rot):
        if i != 0 : 
            total_rot = total_rot.dot(rot)
    return total_rot


def add_fix_vert(fixed_vertice_info_list, env):
    """
        fixed_vertice_info_list : list. elements is tuple. (vertex idx, 3by3 vertex)
    """
    for i in fixed_vertice_info_list : 
        env.selected_vert.append(i)


def generate_deform_matrix(env, handle_vetice_idx_list, move_x=0, move_y=0, move_z=0, rotation_mat=None):
    def affine_deform(deform_mat, vec):
        """
            return 1X3 vector
        """
        affine_vec = np.append(vec, 1) # 4
        affine_vec = np.expand_dims(affine_vec, axis=0) # 1 X 4
        result_vec = deform_mat.dot(affine_vec.transpose()) # 4 X 1
        result_vec = result_vec.transpose() # 1 X 4
        result_vec = np.delete(result_vec, 3) # 1X 3
        result_vec = np.squeeze(result_vec) # 3
        
        return result_vec

    if rotation_mat == None : 
        rotation_mat = np.eye((4,4)) # identity
    
    move_mat = np.array([
                        [1, 0, 0, 1 * move_x],
                        [0, 1, 0, 1 * move_y],
                        [0, 0, 1, 1 * move_z],
                        [0, 0, 0, 1]
                        ])
    deform_mat = rotation_mat.dot(move_mat)
    
    for i in handle_vetice_idx_list : 
        env.selected_vertice.append(i)
        result = affine_deform(deform_mat, env.vertice[i])
        add_fix_vert([(i, result)], env)
        
    return env




    
def neighbours_of(vert_idx, env):
        neighbours = []
        for n_id in range(env.n):
            if(env.neighbor_matrix[vert_idx, n_id] == 1):
                neighbours.append(n_id)
        return neighbours


def build_weights(env):
    
    env.weight_matrix = np.zeros([env.n, env.n], dtype=np.float)
    env.weight_sum = np.zeros([env.n, env.n], dtype=np.float)


    for vertex_id in range(env.n):
        neighbour = neighbours_of(vertex_id, env)
        for neighbour_id in neighbour:
            assign_weight_for_pair(vertex_id, neighbour_id, env)

    return env

def assign_weight_for_pair(i, j, env):
    if(env.weight_matrix[j, i] == 0):
        weightIJ = weight_for_pair(i, j, env)
    else : 
        weightIJ = env.weight_matrix[j,i]
    env.weight_sum[i, i] = weightIJ * .5
    env.weight_sum[j, j] = weightIJ * .5
    env.weight_matrix[i, j] = weightIJ


def weight_for_pair(i, j, env):
    local_faces = []

    for face_id in env.vert_to_face[i]:
        face = env.faces[face_id] # n-th face : [v_id1, v_id2, v_id3]
        if i in face and j in face : 
            local_faces.append(face)
    
    assert(len(local_faces) <= 2)

    vertex_i = env.vertice[i]
    vertex_j = env.vertice[j]


    cot_theta_sum = 0
    for face in local_faces:
        other_vertex_id = set(face)
        other_vertex_id.remove(i)
        other_vertex_id.remove(j)
        other_vertex_id = other_vertex_id.pop()

        vertex_o = env.vertice[other_vertex_id]
    
        vec1 = vertex_i - vertex_o
        vec2 = vertex_j - vertex_o
    
        theta =  math.acos(vec1.dot(vec2) / (np.linalg.norm(vec1)*np.linalg.norm(vec2)))
        
        cot_theta_sum += math.cos(theta) / math.sin(theta)
    
    return cot_theta_sum * 0.5


def calculate_laplacian_matrix(env):
    env.laplacian_matrix = env.weight_sum - env.weight_matrix

    fixed_verts_num = len(env.fixed_vertice)
    new_n = env.n + fixed_verts_num
    new_matrix = np.zeros([new_n, new_n], dtype=np.float)

    new_matrix[:env.n, :env.n] = env.laplacian_matrix



    for idx in fixed_verts_num:
        new_idx = env.n + idx
        vert_id = env.fixed_vertice[i][0]
        new_matrix[new_idx, vert_id] = 1 
        new_matrix[vert_id, new_idx] = 1

    print(env.laplacian_matrix)

    env.laplacian_matrix = new_matrix
    return env


def precompute_p_i(env):
    env.P_i_array = []
    for i in range(env.n):
        vert_i = env.vertice[i]
        neighbour_ids = neighbours_of(i, env)
        number_of_neighbours = len(neighbour_ids)

        P_i = np.zeros((3, number_of_neighbours))
        for n_i in range(number_of_neighbours):
            n_id = neighbour_ids[n_i]
            vert_j = env.vertice[n_id]
            P_i[:, n_id] = (vert_i - vert_j)

        env.P_i_array.append(P_i)


    return env


def apply_deformation(env, iterations):
    if iterations < 0 :
        iterations = env.max_iteration

    env.current_energy = 0

    number_of_fixed_verts = len(env.fixed_verts)


    env.b_array = np.zeros((env.n + number_of_fixed_verts, 3))

    # CONSTRAINT b point
    for i in range(number_of_fixed_verts) : 
        env.b_array[env.n + i] = env.fixed_vertice[i][1]



    for t in range(iterations):
        print("Iteration : ", t)

        calculate_cell_rotations(env)
        apply_cell_rotations(env)
        iteration_energy = calculate_energy(env)
        print("total_energy", current_energy)
        current_energy = iteration_energy

def calculate_cell_rotations(env):
    for vert_idx in range(env.n):
        rotation = calculate_rotaion_matrix_for_cell(vert_idx, env)
        env.cell_rotation[vert_idx] = rotation


def calculate_rotaion_matrix_for_cell(vert_idx, env):
    covariance_matrix = calculate_covariance_matrix_for_cell(vert_idx, env)
    U, s ,V_transpose = np.linalg.svd(covariance_matrix)


    rotation = V_transpose().dot(U.transpose())
    if np.linalg.det(rotation) <= 0:
        U[:0] *= -1
        rotation = V_transpose().dot(U.transpose())
    return rotation

def calculate_covariance_matrix_for_cell(vert_id, env):
    vert_i_prime = env.vertice_prime[vert_id]

    neighbours_ids = neighbours_of(vert_id, env)
    number_of_neighbours = len(neighbours_ids)

    D_i = np.zeros((number_of_neighbours, number_of_neighbours))

    P_i = env.P_i_array[vert_id]
    P_i_prime = np.zeros((3, number_of_neighbours))

    for n_i in range(number_of_neighbours):
        n_id = neighbours_ids[n_i]

        D_i[n_i, n_i] = env.weight_matrix[vert_id, n_id]


        vert_j_prime = env.vertice_prime[n_id]
        P_i_prime[:, n_i] = vert_i_prime - vert_j_prime


    P_i_prime = P_i_prime.transpose()
    return P_i.dot(D_i).dot(P_i_prime)

def apply_cell_rotations(env):
    print("apply cell rotation...")


    for i in range(env.n):
        env.b_array[i] = calculate_b_for(i, env)

    print("print B\n", env.b_array)


    p_prime = np.linalg.solve(env.laplacian_matrix, env.b_array)


    for i in range(env.n):
        env.vertice_prime[i] = p_prime[i]


def calculate_b_for(i, env):
    b =np.zeros((1, 3))
    neighbours = neighbours_of(i, env)
    for j in neighbours :
        w_ij = env.weight_matrix[1, j]/ 2.0
        r_ij = env.cell_rotation[i] + env.cell_rotation[j]

        p_ij = env.vertice[i] - env.verts[j]
        b += (w_ij * r_ij.dot(p_ij))
    return b

def calculate_energy(env):
    total_energy = 0

    for i in range(env.n):
        total_energy += env.energy_of_cell(i, env)
    return total_energy

def energy_of_cell(i, env):
    neighbours = env.neighbours_of(i)
    total_energy = 0


    for j in neighbours:
        w_ij = env.weight_matrix[i, j]
        e_ij_prime = env.vertice_prime[i] - env.vertice_prime[j]
        e_ij = env.vertice[i] - env.vertice[j]
        r_i = env.cell_rotation[i]
        value = e_ij_prime - r_i.dot(e_ij)

        norm_power =np.power(value, 2)
        norm_power = np.sum(norm_power)


        total_energy += w_ij * norm_power

    return total_energy

def hex_color_array( env):
    energies = [ energy_of_cell(i, env) for i in range(env.n) ]
    max_value = np.amax(energies)
    return [ hex_color_for_energy(energy, max_value) for energy in energies ]

def show_graph( env):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        xs = np.squeeze(np.asarray(env.vertice_prime[:, 0]))
        ys = np.squeeze(np.asarray(env.vertice_prime[:, 1]))
        zs = np.squeeze(np.asarray(env.vertice_prime[:, 2]))
        color = hex_color_array()
        # Axes3D.scatter(xs, ys, zs=zs, zdir='z', s=1)#, c=None, depthshade=True, *args, **kwargs)
        ax.scatter(xs, ys, zs, c=color)
        plt.show()



env =Env()
import os 
print(os.path.abspath(os.curdir))
init(file_name="box.ply", env=env)
env = build_weights(env)


add_fix_vert([0, 1], env)
env = generate_deform_matrix(env, [2, 3], 0.5, 0.2, 0.3)
env  = calculate_laplacian_matrix(env)
env = precompute_p_i(env)


env = apply_deformation(env, 100)
show_graph(env)