import numpy as np
from devito import (Grid, TimeFunction, Eq, VectorTimeFunction, TensorTimeFunction,
                    div, grad, curl, diag, Operator, SubDomain, NODE)
from examples.seismic import ModelElastic, plot_velocity, TimeAxis, RickerSource, plot_image

extent = (200., 100., 100.) # 200 x 100 x 100 m domain
h = 1.0 # Desired grid spacing
shape = (int(extent[0]/h+1), int(extent[1]/h+1), int(extent[2]/h+1))

# Model physical parameters:
vp = np.zeros(shape)
qp = np.zeros(shape)
vs = np.zeros(shape)
qs = np.zeros(shape)
rho = np.zeros(shape)

# Set up three horizontally separated layers:
vp[:,:,:int(0.5*shape[2])+1] = 1.52
vs[:,:,:int(0.5*shape[2])+1] = 0.
rho[:,:,:int(0.5*shape[2])+1] = 1.05

vp[:,:,int(0.5*shape[2])+1:int(0.5*shape[2])+1+int(4/h)] = 1.6
vs[:,:,int(0.5*shape[2])+1:int(0.5*shape[2])+1+int(4/h)] = 0.4
rho[:,:,int(0.5*shape[2])+1:int(0.5*shape[2])+1+int(4/h)] = 1.3

vp[:,:,int(0.5*shape[2])+1+int(4/h):] = 2.2
vs[:,:,int(0.5*shape[2])+1+int(4/h):] = 1.2
rho[:,:,int(0.5*shape[2])+1+int(4/h):] = 2.

origin = (0, 0, 0)
spacing = (h, h, h)
so = 4 # FD space order (Note that the time order is by default 1).
nbl = 20 # Number of absorbing boundary layers cells

# Define our 'upper' and 'lower' SubDomains:
class Upper(SubDomain):
    name = 'upper'
    def define(self, dimensions):
        x, y, z = dimensions
        return {x: x, y: y, z: ('left', int(0.5*shape[2])+1+nbl)}

class Lower(SubDomain):
    name = 'lower'
    def define(self, dimensions):
        x, y, z = dimensions
        return {x: x, y: y, z: ('right', shape[2]+nbl-(int(0.5*shape[2])+1))}

# Create these subdomains:
ur = Upper()
lr = Lower()

model = ModelElastic(space_order=so, vp=vp, vs=vs, rho=rho, origin=origin, shape=shape,
                     spacing=spacing, nbl=nbl, subdomains=(ur,lr))
s = model.grid.stepping_dim.spacing

# Source freq. in MHz (note that the source is defined below):
f0 = 0.12

# Thorbecke's parameter notation
l = model.lam
mu = model.mu
ro = model.irho

t0, tn = 0., 30.
dt = 0.9*model.critical_dt
time_range = TimeAxis(start=t0, stop=tn, step=dt)

# PDE fn's:
x, y, z = model.grid.dimensions
damp = model.damp

v = VectorTimeFunction(name='v', grid=model.grid, space_order=so, time_order=1)
tau = TensorTimeFunction(name='t', grid=model.grid, space_order=so, time_order=1)

# Source
src = RickerSource(name='src', grid=model.grid, f0=f0, time_range=time_range)
src.coordinates.data[:] = np.array([100., 50., 35.])

# The source injection term
src_xx = src.inject(field=tau[0, 0].forward, expr=src*s)
src_yy = src.inject(field=tau[1, 1].forward, expr=src*s)
src_zz = src.inject(field=tau[2, 2].forward, expr=src*s)

u_v_u = Eq(v.forward, model.damp * (v + s*ro*div(tau)), subdomain = model.grid.subdomains['upper'])
u_t_u = Eq(tau.forward, model.damp * (tau + s * l * diag(div(v.forward))),
           subdomain = model.grid.subdomains['upper'])

u_v_l = Eq(v.forward, model.damp * (v + s*ro*div(tau)), subdomain = model.grid.subdomains['lower'])
u_t_l = Eq(tau.forward, model.damp * (tau + s * l * diag(div(v.forward))
                                      + s * mu * (grad(v.forward) + grad(v.forward).T)),
           subdomain = model.grid.subdomains['lower'])

op = Operator([u_v_u, u_v_l, u_t_u, u_t_l] + src_xx + src_yy + src_zz, subs=model.spacing_map)

op(dt=dt)

## Plots

# Mid-points:
mid_x = int(0.5*(v[0].data.shape[1]-1))+1
mid_y = int(0.5*(v[0].data.shape[2]-1))+1

# Plot some selected results:

plot_image(v[0].data[1, :, mid_y, :], cmap="seismic")
plot_image(v[0].data[1, mid_x, :, :], cmap="seismic")

plot_image(tau[2, 2].data[1, :, mid_y, :], cmap="seismic")
plot_image(tau[2, 2].data[1, mid_x, :, :], cmap="seismic")

from IPython import embed; embed()

