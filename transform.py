import numpy as np 




class Tramsform:

    """
        SYNTHE ROT
    """
    def __init__(self):
        pass 



    
    


    def get_rotation_mat(angular, axis):
        """

            angluar : angular, real number
            axis : x, y, z Axis

            return SO(3) Matrix
            a d g 
            b e h
            c f i
        """
        # np.pi 
        
        cos_rot = np.cos(angular)
        sin_rot = np.sin(angular)
        
        reval = np.eye(3,3)
        if axis in ["x", "X"]: 
            reval  = np.array([[1, 0, 0],
                                [ 0, cos_rot, -sin_rot]
                                [0, sin_rot, cos_rot]
                                ])
        elif axis in ["y", "Y"]: 
            reval  = np.array([[cos_rot, 0, sin_rot],
                                [ 0, 0, 0]
                                [-sin_rot, 0, cos_rot ]
                                ])
        elif axis in ["z", "Z"]: 
            reval  = np.array([[cos_rot, -sin_rot, 0],
                                [ sin_rot, cos_rot, 0]
                                [0, 0, 1]
                                ])
        
        


        return reval


