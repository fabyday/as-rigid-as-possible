import numpy as np 
import igl 
import math
import copy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class ArapData:
    def __init__(self, V, F):
        assert len(V.shape) == 2, "shape error"
        
        self.v = V
        self.v_prime = copy.deepcopy(V)
        self.v_size = self.v.shape[0]

        self.f = F


        self.fixed_verts = []
        

    def add_fixed_vertex(self, vert_num):
        pass

    def add_selected_vertex(self, vert_num):
        pass


    def compile(self, iteration=1000):
        self.iter_num = 1000
        self.__init_precompute()
        self._build_weight()

    def __init_precompute(self):

        def _vert_idx_neigbor(face_matrix):
            neighbor_of_vert_idx = [[] for i in range(self.v_size)]
            vert_in_face = [[] for i in range(self.v_size)]
            assert(len(neighbor_of_vert_idx)) == self.v_size, "test"
            
            for f_id, faces in enumerate(face_matrix):
                for i in range(len(faces)):
                    neighbor_of_vert_idx[faces[i]].extend([faces[(i+1)%3], faces[(i+2)%3]])
                    vert_in_facefaces[i].append(f_id)
            
                    

            return neighbor_of_vert_idx, vert_in_face
            
        def _make_adj_mat(adj_list):
            adj_mat = np.zero((self.v_size, self.v_size))
            for row_idx in range(len(adj_list)):
                for col_idx in adj_list[row_idx]:
                    adj_mat[row_idx][col_idx] = 1.0
            return adj_mat
                
        def _make_degree_mat(adj_list):
            degree_mat = np.zero((self.v_size, self.v_size))
            for row_idx in range(len(adj_list)):
                degree = len(adj_list[row_idx])
                degree_mat[row_idx][row_idx] = degree
            return degree_mat

        def _make_laplacian(degree_mat, adj_mat):
            return (degree_mat - adj_mat)

     
                    
        self.neighbor_of_vert_idx, self.vert_in_face = _vert_idx_neigbor(self.f)
        self.adj_mat = _make_adj_mat(self.neighbor_of_vert_idx)
        self.degree_mat = _make_degree_mat(self.neighbor_of_vert_idx)
        self.laplacian_mat = _make_laplacian(self.degree_mat, self.adj_mat)


    def _build_weight(self):
        def _assign_weight_for_pair(i, j):
            if(self.weight_matirx[j, i] == 0 ):
                weightIJ = _weight_for_pair(i, j)
            else : 
                weightIJ = self.weight_matirx[j, i]
            self.weight_sum[i, i] = weightIJ * .5
            self.weight_sum[j, j] = weightIJ * .5
            self.weight_matirx[i, j] = weightIJ

        def _weight_for_pair(i, j):
            local_face = []

            for face_id in self.vert_in_face[i]:
                face = self.f[face_id]
                if i in face and j in face : 
                    local_face.append(face)
            
            vertex_i = self.v[i]
            vertex_j = self.v[j]

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




        self.weight_sum = np.zeros([self.v_size, self.v_size], dtype=np.float)
        self.weight_matirx = np.zeros([self.v_size, self.v_size], dtype=np.float)


        for v_idx in range(self.v_size):
            neighbors = self.neighbor_of_vert_idx[v_idx]
            for neighbor_id in neighbors:
                _assign_weight_for_pair(v_idx, neighbor_id)        
    
    




    def apply_deformation(self):
        #CALC ROTATIONS
        def calculate_cell_rotation():
            for v_idx in range(self.v_size):
                rotation = caculate_rotation_matrix_for_cell(v_idx):
                self.cell_rotation[v_idx] = rotation
        
        def caculate_rotation_matrix_for_cell(v_idx):
            covariance_matrix = calculate_covariance_matrix_for_cell(v_idx)
            U, s , V_transpose = np.linalg.svd(covariance_matrix)


            rotation = V_transpose.dot(U.transpose())
            if np.linalg.det(rotation) <= 0 :
                U[:0] *= -1
                rotaiton = V_transpose.dot(U.transpose())
            return rotation

        
        def calculate_covariance_matrix_for_cell(vert_id):
            vert_i_prime = self.v_prime[vert_id]

            neighbor = self.neighbor_of_vert_idx[vert_id]
            number_of_neighbor = len(neighbor)

            D_i = np.zeros((number_of_neighbor, number_of_neighbor))


            P_i = P_i_array[vert_id]
            P_i_prime = np.zeros((3, number_of_neighbor))


            for n_i in range(number_of_neighbor):
                n_id = self.neighbor_of_vert_idx[n_i]
                D_i[n_i, n_i] = self.weight_matirx[vert_id, n_id]
                vert_j_prime = self.v_prime[n_id]
                P_i_prime[:, n_i] = vert_i_prime - vert_j_prime

            P_i_prime = P_i_prime.transpose()
            return P_i.dot(D_i).dot(P_i_prime)
        # APPLY
        def apply_cell_rotation():
            print("apply cell rotation...")
            for i in range(self.v_size):
                self.b_array[i] = calculate_b_for(i)
            p_prime = np.linalg.solve(self.laplacian_mat, self.b_array)
            for i in range(self.v_size):
                self.v_prime[i] = p_prime[i]
    
        self.current_energy = 0.0

        number_of_fixed_verts = len(self.fixed_verts)

        self.b_array = np.zeros((self.v_size + number_of_fixed_verts, 3))


        for i in range(number_of_fixed_verts):
            self.b_array[self.fixed_verts + i] = self.v[self.fixed_verts[i]]


        for t in range(self.iter_num):
            print("iter num : ", t)
            self.calculate_cell_rotation()
            self.apply_cell_rotation()
            iteration_energy = self.calculate_energy()
            self.current_energy = iteration_energy

    
    
