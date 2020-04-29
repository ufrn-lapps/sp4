##########################################################################################
# This code is based on devito user api tutorial 03_subdomains.
# 
# Information:
# Grid:             200 x 100 x 100
# Grid spacing:     1.0
# Stencil Order:    4
# Number NBL cells: 20
#
# Upper Domain 
# Water layer
# vp:               1500 m.s⁻¹
# vs:               0 m.s⁻¹
# rho:              1020 kg.m⁻³
#
# Lower Domain
# Solid layer
# vp:               3400 m.s⁻¹
# vs:               1963 m.s⁻¹
# rho:              2599 m.s⁻¹
##########################################################################################



import numpy as np
from sympy import Symbol
from devito import (Grid, TimeFunction, Eq, VectorTimeFunction, TensorTimeFunction,
                    div, grad, curl, diag, solve, Operator, SubDomain, NODE)
from examples.seismic import Model, ModelElastic, RickerSource, TimeAxis, plot_image, plot_velocity


def experiment_subdomain(water_level=0.5):

    # Configuration
    extent = (200., 100., 100.)  # 200 x 100 x 100 m domain
    h = 1.0                      # Desired grid spacing
    spacing = (h, h, h)
    shape = (int(extent[0]/h+1),
             int(extent[1]/h+1),
             int(extent[2]/h+1))
    origin = (0, 0, 0)
    # FD space order (Note that the time order is by default 1).
    so = 4
    nbl = 20                     # Number of absorbing boundary layers cells

    # Model physical parameters:
    vp = np.zeros(shape)
    # qp = np.zeros(shape)
    vs = np.zeros(shape)
    # qs = np.zeros(shape)
    rho = np.zeros(shape)

    # Set up two horizontally separated layers:
    # Water layer
    vp[:, :, :int(water_level*shape[2])+1] = 1.50   # Vp 1500 m.s⁻¹
    vs[:, :, :int(water_level*shape[2])+1] = 0.     # Vs 0 m.s⁻¹
    rho[:, :, :int(water_level*shape[2])+1] = 1.02  # Water density 1020 kg.m⁻³

    # Solid layer
    vp[:, :, int(water_level*shape[2])+1:] = 3.4   # Vp 3400 m.s⁻¹
    vs[:, :, int(water_level*shape[2])+1:] = 1.93  # Vs 1963 m.s⁻¹
    rho[:, :, int(water_level*shape[2])+1:] = 2.5  # Solid density 2500 kg.m⁻³

    # Define our 'upper' and 'lower' SubDomains:
    class Upper(SubDomain):
        name = 'upper'

        def define(self, dimensions):
            x, y, z = dimensions
            return {x: x, y: y, z: ('left', int(water_level*shape[2])+1+nbl)}

    class Lower(SubDomain):
        name = 'lower'

        def define(self, dimensions):
            x, y, z = dimensions
            return {x: x, y: y, z: ('right', shape[2]+nbl-(int(water_level*shape[2])+1))}

    class Interface(SubDomain):
        """
        A small interface between fluid and solid domains.
        """
        name = 'interface'

        def define(self, dimensions):
            x, y, z = dimensions
            return {x: x, y: y, z: ('middle', int(water_level*shape[2])+1+nbl - 2, int(water_level*shape[2])+1+nbl + 2)}

    # Create these subdomains:
    ur = Upper()
    lr = Lower()
    interface = Interface()

    model_elastic = ModelElastic(space_order=so, vp=vp, vs=vs, rho=rho, origin=origin, shape=shape,
                                 spacing=spacing, nbl=nbl, subdomains=(ur, lr, interface))
    f0 = 0.12

    # Thorbecke's parameter notation
    l = model_elastic.lam                        # (vp² - 2vs²)*rho
    mu = model_elastic.mu                        # vs²*rho
    irho = model_elastic.irho                    # 1/rho
    s = model_elastic.grid.stepping_dim.spacing

    t0, tn = np.float32(0.), np.float32(30.)
    dt = np.float32(0.9*model_elastic.critical_dt)
    # dt = model_elastic.critical_dt
    time_range = TimeAxis(start=t0, stop=tn, step=dt)

    # PDE fn's:
    # x, y, z = model_elastic.grid.dimensions
    # damp = model_elastic.damp

    v = VectorTimeFunction(
        name='v', grid=model_elastic.grid, space_order=so, time_order=1)
    tau = TensorTimeFunction(
        name='t', grid=model_elastic.grid, space_order=so, time_order=1)

    # Source
    src = RickerSource(name='src', grid=model_elastic.grid,
                       f0=f0, time_range=time_range)
    src.coordinates.data[:] = np.array([100., 50., 35.])

    # The source injection term
    src_xx = src.inject(field=tau[0, 0].forward, expr=src*s)
    src_yy = src.inject(field=tau[1, 1].forward, expr=src*s)
    src_zz = src.inject(field=tau[2, 2].forward, expr=src*s)

    # print(solve(v.forward, v + s*ro*div(tau)))

    # Upper domain equations
    # Update velocity in upper domain
    u_v_u = Eq(v.forward, model_elastic.damp * (v + s*irho*div(tau)),
               subdomain=model_elastic.grid.subdomains['upper'])
    # Update tau in upper domain
    u_t_u = Eq(tau.forward, model_elastic.damp * (tau +
                                                  s * l * diag(div(v.forward))),
               subdomain=model_elastic.grid.subdomains['upper'])

    # Lower domain equations
    # Update velocity in lower domain
    u_v_l = Eq(v.forward, model_elastic.damp * (v + s*irho*div(tau)),
               subdomain=model_elastic.grid.subdomains['lower'])
    # Update tau in lower domain
    u_t_l = Eq(tau.forward, model_elastic.damp * (tau +
                                                  s * l * diag(div(v.forward)) +
                                                  s * mu * (grad(v.forward) + grad(v.forward).T)),
               subdomain=model_elastic.grid.subdomains['lower'])

    op = Operator([u_v_u, u_v_l, u_t_u, u_t_l] + src_xx + src_yy + src_zz, subs=model_elastic.spacing_map)
    op.apply(dt=dt)

    # Plots
    # Mid-points:
    #mid_x = int(0.5*(v[0].data.shape[1]-1))+1
    #mid_y = int(0.5*(v[0].data.shape[2]-1))+1

    # Plot some selected results:
    #plot_image(v[0].data[1, :, mid_y, :], cmap="seismic")
    #plot_image(v[0].data[1, mid_x, :, :], cmap="seismic")

    #plot_image(tau[2, 2].data[1, :, mid_y, :], cmap="seismic")
    #plot_image(tau[2, 2].data[1, mid_x, :, :], cmap="seismic")


for i in range(5):
    experiment_subdomain(0.5)
