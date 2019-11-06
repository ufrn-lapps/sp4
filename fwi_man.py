import numpy as np
import scipy
from sympy import Symbol
from devito import (configuration, Function, clear_cache, TimeFunction, Eq,
                    solve, Operator, Inc)
from examples.seismic import (AcquisitionGeometry, RickerSource, Receiver, TimeAxis,
                              plot_shotrecord, plot_velocity, plot_image)
from examples.seismic.acoustic import AcousticWaveSolver
from model import demo_model

configuration['log-level'] = 'WARNING'

# space
shape = (101, 101)
spacing = (10., 10.)
origin = (0., 0.)
model = demo_model('circle-isotropic', vp=3.0, vp_background=2.5,
                   origin=origin, shape=shape, spacing=spacing, nbpml=40)
model0 = demo_model('circle-isotropic', vp=2.5, vp_background=2.5,
                    origin=origin, shape=shape, spacing=spacing, nbpml=40,
                    grid=model.grid)

# time
t0 = 0.
tn = 1000.
f0 = 0.010
time_axis = TimeAxis(start=t0, stop=tn, step=model.critical_dt)
nt = time_axis.num

# source
nshots = 9
src_coordinates = np.empty((1, 2))
source_locations = np.empty((nshots, 2), dtype=np.float32)
source_locations[:, 0] = 30.
source_locations[:, 1] = np.linspace(0., 1000, num=nshots)

# receiver
nreceivers = 101
rec_coordinates = np.empty((nreceivers, 2))
rec_coordinates[:, 1] = np.linspace(0, model.domain_size[0], num=nreceivers)
rec_coordinates[:, 0] = 980.

# for comparison with 03_fwi.ipiynb tutorial:
geometry = AcquisitionGeometry(
    model, rec_coordinates, src_coordinates, t0, tn, f0=f0, src_type='Ricker')
solver = AcousticWaveSolver(model, geometry, space_order=4)

# stencil
u = TimeFunction(name='u', grid=model.grid,
                 time_order=2, space_order=4, save=nt)
H = Symbol('H')
eq = model.m * u.dt2 + model.damp * u.dt - H
eq_time = solve(eq, u.forward).subs({H: u.laplace})
eqn = [Eq(u.forward, eq_time)]

# stencil-back
v = TimeFunction(name='v', grid=model.grid,
                 time_order=2, space_order=4, save=nt)
eq_back = model.m * v.dt2 - H - model.damp * v.dt
eq_time_back = solve(eq_back, v.backward).subs({H: v.laplace})
eqn_back = [Eq(v.backward, eq_time_back)]

# True Data
u0_man = TimeFunction(name='u', grid=model.grid, time_order=2, space_order=4,
                      save=nt)
receivers = []
for i in range(9):
    src_true = RickerSource(name='src', grid=model.grid, time_range=time_axis,
                            coordinates=source_locations[i, :], npoint=1, f0=f0)
    src_term = src_true.inject(
        field=u.forward, expr=src_true * model.grid.stepping_dim.spacing**2 / model.m)

    # receiver
    receiver = Receiver(name='rec', grid=model.grid, time_range=time_axis,
                        coordinates=rec_coordinates, npoint=nreceivers)
    rec_term = receiver.interpolate(expr=u)
    receivers.append(receiver)

    # operator
    op_fwd = Operator(eqn + src_term + rec_term,
                      subs=model.spacing_map, name='Forward')

    op_fwd.apply(src=src_true, rec=receivers[i], u=u0_man,
                 vp=model.vp, dt=model.critical_dt)


def fwi_gradient(vp_in):
    # MANUAL
    grad_manual = Function(name="grad", grid=model.grid)
    residual_man = Receiver(name='rec', grid=model.grid,
                            time_range=time_axis,
                            coordinates=rec_coordinates)
    objective_manual = 0.

    for i in range(9):
        # MANUAL
        clear_cache()

        # source
        src_true = RickerSource(name='src', grid=model.grid, time_range=time_axis,
                                coordinates=source_locations[i, :], npoint=1, f0=f0)
        src_term = src_true.inject(
            field=u.forward, expr=src_true * model.grid.stepping_dim.spacing**2 / model.m)

        # receiver
        rec_term = receivers[i].interpolate(expr=u)

        # operator
        op_fwd = Operator(eqn + src_term + rec_term,
                          subs=model.spacing_map, name='Forward')

        u0_man.data.fill(0.)
        rec_smooth = Receiver(name='rec', grid=model.grid, time_range=time_axis,
                              coordinates=rec_coordinates, npoint=nreceivers)
        op_fwd.apply(src=src_true, rec=rec_smooth,
                     u=u0_man, vp=vp_in, dt=model.critical_dt)

        # back-receiver
        rec_back = Receiver(name='rec', grid=model.grid, time_range=time_axis,
                            coordinates=rec_coordinates, npoint=nreceivers)
        rec_back_term = rec_back.inject(field=v.backward,
                                        expr=rec_back * model.grid.stepping_dim.spacing**2 / model.m)

        # gradient
        gradient_update = Inc(grad_manual, - u.dt2 * v)
        op_grad = Operator(eqn_back + rec_back_term + [gradient_update], subs=model.spacing_map,
                           name='Gradient')
        residual_man.data[:] = rec_smooth.data[:] - receivers[i].data[:]
        objective_manual += .5*np.linalg.norm(residual_man.data.flatten())**2
        op_grad.apply(rec=residual_man, u=u0_man, vp=vp_in,
                      dt=model.critical_dt, grad=grad_manual)

        # sanity-check -> expect for 0!
        # plot_shotrecord(true_d.data[:] - rec_true.data[:], model, t0, tn)
        # plot_shotrecord(smooth_d.data[:] - rec_smooth.data[:], model, t0, tn)
    return objective_manual, -grad_manual.data

# ff, update,ff_man,update_man = x(model0.m)
# print('Objective AUTO is %f ' % ff)
# print('Objective MANUAL is %f ' % ff_man)
# plot_image(update, vmin=-1e4, vmax=1e4, cmap="jet")
# plot_image(update_man, vmin=-1e4, vmax=1e4, cmap="jet")


def apply_box_constraint(vp):
    return np.clip(vp, 2.0, 3.5)


fwi_iterations = 5
history_man = np.zeros((fwi_iterations, 1))
print('MANUAL')
for i in range(0, fwi_iterations):
    phi, direction = fwi_gradient(model0.vp)
    history_man[i] = phi
    alpha = .05 / np.abs(direction).max()
    model0.vp = apply_box_constraint(model0.vp.data - alpha * direction)
    print('Objective value is %f at iteration %d' % (phi, i+1))
