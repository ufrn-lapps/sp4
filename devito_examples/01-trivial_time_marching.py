from devito import Eq, Grid, TimeFunction, Operator
from devito.ir.iet import  find_affine_trees, retrieve_iteration_tree

grid_2d = Grid(shape=(4, 4))
v = TimeFunction(name='v', grid=grid_2d, time_order=2, save=2)
equation = Eq(v.forward, v+1)
operator = Operator(equation)