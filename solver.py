import numpy as np 
import copy
import mesh as mobj
import logging 

logging.basicConfig(level=logging.DEBUG)
slvr_logger = logging.getLogger("Solver")




class Sovler:
    """
        How Much Cell exists?
        cell mean 

            |
          --O--
            |
        each Vertex will become Cell.
        So Number of Cell, mean Number of Vertex.
    """
    def __init__(self, mesh, constraint, iteration = 100):
        # INPUTS DATA (PROBLEM)
        self.mesh = mesh # Main REFERNCE 
        self.v = mesh.get_vertex()
        self.f = mesh.get_face()
        self.laplacian = mesh.get_laplacian_matrix()


        
        
        # SOLVER vals
        self.constraint = constraint
        self.iteration = iteration
        
    def _build_weight_per_cell(self, vert_face_match_list, neighbor_list ):
        """
            #TODO Need To Refactoring.
            Cell Weights is 1.0. 

            vert_face_match_list : Face list that involved vertex 
            neighbor_list : neighbor list
            laplacian : L := D - A 
        """
        cell_weights = np.eye(self.f.shape[0], self.f.shape[0]) # M X M 
        edge_weights = np.zeros((self.v.shape[0], self.v.shape[0])) # N X N 


        
        # Set Usable Variables for Local Functions  
        v = self.v #vertex values
        vert_f_list =  vert_face_match_list  #vertex face match list
        f = self.f  # this is Face

        def _build_weight_per_edge(i, j):
            """
                i : starting point of vertex index.
                j : 1-ring neighbor vertex index of vertex i.
            """
            v_i = v[i]
            v_j = v[j]            
            vec_ij = v_j - v_i 

            # list of face associated with vertex
            # associated face, more than 2
            face_list_associated_with_i = vert_f_list[i]
            face_list_associated_with_j = vert_f_list[j]
            
            # Find only Common List
            # But I Expect it JUST 2
            face_list_associated_with_common = face_list_associated_with_i and face_list_associated_with_j  # [ 1, 2, 3 ] and [ 2, 3 ] := [ 2, 3 ]

            other_ne_vector_list = []
            for face_idx in face_list_associated_with_common:
                v_list = set(f[face_idx])
                result = v_list.difference( set([i , j]) ) 
                logging.debug("result set size ({}) and, result elems {} ".format(len(result), result) )
                
                
                result_vec = v_i - v[list(result)[0]] # TODO it's looks complicated.
                other_ne_vector_list.append(result_vec)
                # I Expect just 1 Argument.
                assert len(result) == 1, "we expect 1 Element Per each Triangle. but this result is {}".format(len(result))


                

            alpha_cos = vec_ij.dot(other_ne_vector_list[1])/ (np.linalg.norm(other_ne_vector_list[0]) * np.linalg.norm(vec_ij))
            beta_cos  = vec_ij.dot(other_ne_vector_list[1]) / (np.linalg.norm(other_ne_vector_list[0]) * np.linalg.norm(vec_ij))

            # Arccos, It is iverse Function of Cos. A.K.A AK-47~
            alpha_theta = np.arccos(alpha_cos)
            beta_theta = np.arccos(beta_cos)

            weight = 0.5 *(np.tan(beta_theta) + np.tan(alpha_theta))
            return weight

        # self.cell_weights
        # self.edge_weights
        # Num{n} := number of -1. 
        # that mean Num1 := 3
        # [Num1 -1      -1      -1  ]
        # [0    Num2    -1      -1  ]
        # [-1   0       Num3    -1  ]
        # [-1   0       0       Num4]
        for vert_idx, vert_ne_list in enumerate(neighbor_list) : 
            weight_sum = 0
            for neighbor_idx in vert_ne_list : 
                weightIJ = _build_weight_per_edge(vert_idx, neighbor_idx) 
                edge_weights[vert_idx, neighbor_idx] = -weightIJ
                weight_sum += weightIJ
            edge_weights[vert_idx, vert_idx] = weight_sum
            


        return cell_weights, edge_weights

    def _build_weight(self):
        """
            for every triangle <weights>
                for each Edge
            per Cell : number of Face?
            per Edge : NoF X 


        """

        #Sum(cell Energy Weight Sum(edge Energy Weight) )
        self.neighbor_list = self.mesh.get_vertex_neighbor_list()
        self.neighbor_matrix = self.mesh.get_neighbor_matrix()
        self.vert_f_list = self.mesh.get_face_vertex_mathcing_list()

        # TODO Robustly, doesn't Using Sparse Matrix. Only Dense Type Matrix is used.
        cell_weights, edge_weights = self._build_weight_per_cell(self.vert_f_list, self.neighbor_list)
        return cell_weights, edge_weights

    def _cell_rotation_init(self):
        """
            cell roation := eye(3,3)
        """
        reval = np.zeros( (len(self.v[0]),3,3) )
        for idx in len(self.v.shape[0]):
            reval[idx, ...] = np.eye(3,3)
        
        return reval


    def _concat_constraint(self, mat, cons):
        #[]     []
        #[] +   []
        #
        #   []
        #   []
        #   []
        #   []
        
        assert mat.shape[-1] == cons.shape[-1], "2 matrix shape is different"


        return np.concatenate( [mat, cons], axis = 0 )
    
        



    ############### PUBLIC METHOD #########################
    def precompute(self):
        #build Weight
        self.cell_weights, self.edge_weights = self._build_weight()
        self.b = np.zeros_like(self.v)
        # self.cell_ratations = self._cell_rotation_init()

        #Add Constraint
        Cons_A, Cons_b = self.constraint.to_dense_mat(len(self.edge_weights.shape[0]), len(self.v.shape[-1]))
        self.edge_weights   = self._concat_constraint(self.edge_weights, Cons_A)
        self.b              = self._concat_constraint(self.b, Cons_b)








        # TODO It need but not now. it's disable to having intuition.
        # self.v = self.v.transpose() # N X 3 -> 3 X N
        # self.v = self.v.reshape(-1, self.v.shape[0]) # [ X X X ... Y Y Y ... Z Z Z ] := 3 X N -> 3*N 
        # self.v_prime = copy.deepcopy(self.v) 


    def solve(self):
        """
            solve data S' and S
        """

        #initialize vertex, bias ... 
        self._cell_rotation_init()
        self.vertex_prime = self.v
        

        crruent_energy = 0
        for t in range(self.iteration):
            self.calc_cell_rotation()
            self.apply_cell_rotation()

            iteration_energey = self.calc_total_energy()
            print("total energy : ",  crruent_energy)
            self.b = self.apply_bias()


            self.vertex_prime = self.ls_solve()

        return self.vertex_prime, self.f





    def calc_cell_rotation(self):
        
        for idx in range(len(self.cell_ratations)) : # cell_roations := ( N, 3, 3 )
            rots = self.apply_cell_rotation(idx)
            self.cell_rotations[idx] = rots # update it!

        
        
    def apply_cell_rotation(self, idx):
        """
            Cell Rotation Calculation.
            Use SVD decomposition.
            U and V_T's multiplcation mean Rotations.
        """
        
        U, Sig, V_T = np.linalg.svd(self.cell_rotations[idx])
        new_rotation = V_T.T.dot(U.T)
        det = np.linalg.det(new_rotation)
        # S, _, D = np.linalg.svd()
        if det < 0. : #Orthogonal that is not to bocome 0.
            new_rotation *= -1

        return new_rotation


    def calc_total_energy(self):
        # Tr(L * P') = total Energy
        np.trace(self.edge_weights.dot(self.vertex_prime))
        pass
    
    def apply_bias(self):
        """
            b := (N, 3)
            b := 0.5 *  w_ij * (R_i + R_j) * (p_i - p_j) := 0.5 * (Ri + Rj) L(w_ij included in L) * x
            i
            iterative solution.
            one time solution...
        """
        N = self.neighbor_matrix.shape[-1] # Matrix N X N
        rots = self.cell_ratations.reshape(N, 9) # N X N :-> N X 9
        total_rots_sums = self.neighbor_matrix.dot(rots)  # := R_i + R_j, N X N * N X 9 :-> N X 9 
        total_rots_sums = total_rots_sums.reshape(N, 3, 3) # N X 3 X 3 
        
        W_L_x = self.edge_weights.dot(self.v) # N X 3 
        W_L_x_transpose = W_L_x.transpose()
        # np.tensordot(total_rots_sums, W_L_x, [[1, 2], [1]])
        for i in range(N): # TODO parell
            W_L_x_transpose[..., i, None] = total_rots_sums[i].dot(W_L_x_transpose[..., i, None])


        return W_L_x_transpose.transpose() # return N X 3
        

        def ls_solve(self):
            vertex_prime = np.linalg.solve(self.edge_weights, self.b)
            return vertex_prime




    




