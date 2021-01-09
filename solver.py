import numpy as np 
import copy
import mesh as mobj




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
        cell_weights = np.eye((self.f.shape[0], self.f.shape[0])) # M X M 
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
                result = v_i - v[result] # TODO it's looks complicated.
                other_ne_vector_list.append(result)
                # I Expect just 1 Argument.
                assert len(result) == 1, "we expect 1 Element Per each Triangle."

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
        neighbor_list = self.mesh.get_vertex_neighbor_list()
        neighbor_matrix = self.mesh.get_neighbor_matrix()
        vert_f_list = self.mesh.get_face_vertex_mathcing_list()

        # TODO Robustly, doesn't Using Sparse Matrix. Only Dense Type Matrix is used.
        cell_weights, edge_weights = self._build_weight_per_cell(vert_f_list, neighbor_list)
        return cell_weights, edge_weights

    def _cell_rotation_init(self):
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
        self.cell_ratations = self._cell_rotation_init()

        #Add Constraint
        Cons_A, Cons_b = self.constraint.to_dense_mat(len(self.edge_weights.shape[0]), len(self.v.shape[-1]))
        self.cell_weights   = self._concat_constraint(self.cell_weights, Cons_A)
        self.b              = self._concat_constraint(self.b, Cons_b)








        # TODO It need but not now. it's disable to having intuition.
        # self.v = self.v.transpose() # N X 3 -> 3 X N
        # self.v = self.v.reshape(-1, self.v.shape[0]) # [ X X X ... Y Y Y ... Z Z Z ] := 3 X N -> 3*N 
        # self.v_prime = copy.deepcopy(self.v) 


    def solve(self):
        pass





