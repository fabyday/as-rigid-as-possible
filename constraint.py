class Constraint : 
    """
        For ... SPARSE LIN SOLVER
    """

    class Triplet :
        def __init__(self, vertex_idx, coord):
            self.vertex_idx = vertex_idx
            self.value = coord
        
        def __call__(self):
            return vertex_idx, value

    def __init__(self):
        self.contraint_list = []


    def reset_constraint(self):
        self.contraint_list = []

    def add_constraint(self, vertex_idx, coord):
        """
            
            vertex_idx : index of vertex
            coord : constraint coordinate. [X, Y, Z]
        """
        self.contraint_list.append(Constraint.Triplet(vertex_idx, coord))
    
    def update_constraint(self, vertex_idx, update_coord):
        #TODO
        pass



    

    def to_sparse_mat(self):
        pass
    

    def to_dense_mat(self, v_rows, v_cols):
        # Constraint_list_size, v_rows(N)
        # Constraint_list_size, v_cols(3)
        Cons_A = np.zeros(len(self.constraint_list), v_rows).astype("np.float32")
        Cons_b = np.zeros(len(self.constraint_list), v_cols).astype("np.float32")
        

        for idx, data in enumerate(constraint_list) : 
            row_idx, coord = data()
            Cons_A[idx, row_idx] = 1
            Cons_b[idx] = coord

        return Cons_A, Cons_b






