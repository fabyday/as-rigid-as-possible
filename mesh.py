
import numpy as np 
class Mesh:
    """
        Class Mesh 
            mesh를 원하는 형태의 행렬로 변환해준다.
        TODO precomute가 필요하면 따로 만들것.
    """
    # i = vertex index.
    # j = neighbor vertex index.
    # p'_i - p'_j  = R(p_i - p_j)
    

    def __init__(self, V, F):
        """
            V : ndArray Object. (N, 3)
            F : ndArray Object. (M, 3)
        """
        self.V = V
        self.F = F

    def get_vertex(self):
        return self.V

    def get_face(self):
        return self.F
    

    def get_neighbor_matrix(self):
        """
            리턴 값. 두개
            
        """
        
        vertex_num = self.V.shape[0]
        reval = np.zeros((vertex_num, vertex_num)) # N X N Matrix.
        dims = 3
        # reval_2 = [[] for _ in len(self.F) ] 임시 보류  TODO

        for row_idx, face_mem in enumerate(self.F) : # row_idx, face_mem := idx, (v_id_1 , v_id_2, v_id_3)
            
            for i in range(len(face_mem)) : 
                reval [ face_mem[i] ] [ face_mem[ (i + 1)%dims ] ] = 1
                reval [ face_mem[i] ] [ face_mem[ (i + 2)%dims ] ] = 1

            
                
                

        return reval



    def get_face_vertex_mathcing_list(self):
        """ 
            v_idx1 : f_idx1 f_idx2
            v_idx2 : f_idxN f_idx22 .... f_idxN
            .
            .
            .

            return Just Python Built-in List(TODO change it numpy array if it was needed)
        """ 

        V = self.V
        F = self.F 
        N, _ = V.shape
        M, _ = F.shape

        reval = [[] for _ in range(N)] # list size is same as vertex row size.
        for face_row_idx, faces in enumerate(F) :
            for v_id in faces : 
                reval[v_id].append(face_row_idx)
        return reval
        


    def get_vertex_neighbor_list(self):
        """
            different neighbor Matrix. that's Adj Matrix.
            neighbor List : 
                v1 : v2 v3 ... vn
                v2 : v4 v5 ... n? 
                .
                .
                .
        """    
        vertex = self.V
        face = self.F
        assert len(vertex.shape) == 2, "vertex Shape is not 2-dims"
        assert len(face.shape) == 2, "face Shape is not 2-dims"

        neighbour = [ [] for _ in range(len(vertex)) ]
        
        for f_data_list in face:
            for i, vertex_num in enumerate(f_data_list):
                    if (f_data_list[(i+1)%3]+1) not in neighbour[vertex_num]:
                        neighbour[vertex_num].append(f_data_list[(i+1)%3] + 1)
                    if (f_data_list[(i+2)%3]+1) not in neighbour[vertex_num]:
                        neighbour[vertex_num].append(f_data_list[(i+2)%3] + 1)
                    
        pointnum = len(vertex)
        maxdegree = 0 
        degree = np.zeros(pointnum, dtype=np.int32)
        
        for i in range(pointnum):
            degree[i] = len(neighbour[i])
            if degree[i] > maxdegree:
                maxdegree = degree[i]
        after_neighbour = np.zeros((pointnum, maxdegree), dtype=np.int32)
        
        for i in range(pointnum):
            zero_adding_size = maxdegree - degree[i]
            after_neighbour[i] = neighbour[i]+[0]*zero_adding_size

        return after_neighbour        


    def get_laplacian_matrix(self):
        """
            return Laplacian Matrix.
            [[]]
        """
        
        neighbor = self.get_neighbor_matrix() # N X N 
        sums = np.sum(neighbor, axis = -1 ) #[N, 1]
        reval = np.eye(neighbor.shape[0]) * sums # [ N X N ]
        reval -= neighbor # Degree Mat - Adj Mat

        return reval
        