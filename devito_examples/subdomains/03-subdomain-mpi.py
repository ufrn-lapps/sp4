## ------------------------------------------------------------------------
## Example to test the implementation of Functions Subdomain over MPI.
## This example should be executed with 4 process.
## ------------------------------------------------------------------------

import numpy as np
from mpi4py import MPI

from devito import (Grid, Function, TimeFunction, Eq, solve, Operator, SubDomain,
                    SubDomainSet, Dimension)
from devito.data import LEFT, RIGHT

class SubDomain_Middle(SubDomain):
    name = 'sd_middle'

    def define(self, dimensions):
        x, y = dimensions
        # Create a 4x4 subdomain in the middle of the grid.
        return {x: ('middle', 3, 3), y: ('middle', 3, 3)}

sd_middle = SubDomain_Middle()

grid = Grid((10, 10), subdomains=(sd_middle))
x, y = grid.dimensions

u = Function(name='u', grid=grid, subdomain=sd_middle)
v = Function(name='v', grid=grid)

# Initialize data with its rank
# u.data[:] = MPI.COMM_WORLD.rank
u.data[:] = 421.0

eqn = Eq(v, u + 1, subdomain=grid.subdomains['sd_middle'])

op = Operator(eqn)
op.apply()
print(op)
print(v.data[:])
print(u.data[:])

#if LEFT in glb_pos_map[x] and LEFT in glb_pos_map[y]:
#    expected = np.array([[0., 0., 0., 0., 0.,],
#                         [0., 0., 0., 0., 0.,],
#                         [0., 0., 0., 0., 0.,],
#                         [0., 0., 0., 2., 2.,],
#                         [0., 0., 0., 2., 2.,]], dtype=np.float)
#    assert ((np.array(v.data[:]) == expected).all())
#
#elif LEFT in glb_pos_map[x] and RIGHT in glb_pos_map[y]:
#    expected = np.array([[0., 0., 0., 0., 0.,],
#                         [0., 0., 0., 0., 0.,],
#                         [0., 0., 0., 0., 0.,],
#                         [2., 2., 0., 0., 0.,],
#                         [2., 2., 0., 0., 0.,]], dtype=np.float)
#    assert ((np.array(v.data[:]) == expected).all())
#
#elif RIGHT in glb_pos_map[x] and LEFT in glb_pos_map[y]:
#    expected = np.array([[0., 0., 0., 2., 2.,],
#                         [0., 0., 0., 2., 2.,],
#                         [0., 0., 0., 0., 0.,],
#                         [0., 0., 0., 0., 0.,],
#                         [0., 0., 0., 0., 0.,]], dtype=np.float)
#    assert ((np.array(v.data[:]) == expected).all())
#
#elif RIGHT in glb_pos_map[x] and RIGHT in glb_pos_map[y]:
#    expected = np.array([[2., 2., 0., 0., 0.,],
#                         [2., 2., 0., 0., 0.,],
#                         [0., 0., 0., 0., 0.,],
#                         [0., 0., 0., 0., 0.,],
#                         [0., 0., 0., 0., 0.,]], dtype=np.float)
#    assert ((np.array(v.data[:]) == expected).all())
