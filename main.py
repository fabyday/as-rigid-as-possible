import igl

import solver 
import constraint
import transform
import open3d as o3d
import mesh as m


cons = constraint.Constraint()
cons.add_constraint(0, [0,0,0])
# slvr = solver.Sovler()


if  __name__ == "__main__":
   V, F = igl.read_triangle_mesh("cube.obj")
   
   mesh = m.Mesh(V, F)

   slvr = solver.Sovler(mesh, cons)
   slvr.precompute()
   V, F = slvr.solve()

   igl.write_obj("cube_remake", V, F)

