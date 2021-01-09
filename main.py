import igl

import solver 
import constraint
import transform
import open3d as o3d


cons = constraint.Constraint()
cons.add_constraint()
slvr = solver.Sovler()


if  __name__ == "__main__":
   slvr.precompute()
   V, F = slvr.solve()

