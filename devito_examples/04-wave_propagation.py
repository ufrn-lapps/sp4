import numpy as np
import sys

from devito.tools import pprint
from devito import Eq, Operator, TimeFunction, solve
from examples.seismic import Model, Receiver, RickerSource, TimeAxis, plot_image, plot_shotrecord


# DEFINE MODEL

# Define a physical size
shape = (101, 101, 101)    # Number of grid point (nx, nz)
spacing = (1., 1., 1.)  # Grid spacing in m. The domain size is now 1km by 1km

origin = (0., 0., 0.)
# the absolute location of the source and receivers

# Define a velocity profile. The velocity is in km/s
v = np.empty(shape, dtype=np.float32)
v[:, :, :] = 1.5

# With the velocity and model size defined, we can create the seismic model that
# encapsulates this properties. We also define the size of the absorbing layer as 10 grid points
model = Model(vp=v, origin=origin, shape=shape, spacing=spacing,
              space_order=2, nbl=40, bcs="damp")

# plot_velocity(model)

# TIME

t0 = 0.                 # Simulation starts a t=0
tn = 1000.              # Simulation last 1 second (1000 ms)
dt = model.critical_dt  # Time step from model grid spacing

time_range = TimeAxis(start=t0, stop=tn, step=dt)
print(dt)

# SOURCE

f0 = 0.010  # Source peak frequency is 10Hz (0.010 kHz)
src = RickerSource(name='src', grid=model.grid, f0=f0,
                   npoint=1, time_range=time_range)


# Src in the middle of the cube
src.coordinates.data[0, :] = np.array(model.domain_size) * .5

# We can plot the time signature to see the wavelet
# src.show()

# RECEIVERS


# Create symbol for 101 receivers
rec = Receiver(name='rec', grid=model.grid, npoint=101, time_range=time_range)

# Prescribe even spacing for receivers along the x-axis
rec.coordinates.data[:, 0] = np.linspace(0, model.domain_size[0], num=101)
rec.coordinates.data[:, 1:] = 20.  # Depth is 20m


# Define the wavefield with the size of the model and the time dimension
u = TimeFunction(name="u", grid=model.grid, time_order=2, space_order=2)

# We can now write the PDE
pde = model.m * u.dt2 - u.laplace + model.damp * u.dt

stencil = Eq(u.forward, solve(pde, u.forward))

# Finally we define the source injection and receiver read function to generate the corresponding code
src_term = src.inject(field=u.forward, expr=src * dt**2 / model.m)

# Create interpolation expression for receivers
rec_term = rec.interpolate(expr=u.forward)


op = Operator([stencil] + src_term + rec_term, subs=model.spacing_map, dle='noop')

print("Argumentos passados para o código C")
for k, v in op.arguments().items():
    print(k, v)


op.apply(time=time_range.num - 1, dt=model.critical_dt)
# import pdb; pdb.set_trace()
with open('out-time-1000.txt', 'w') as f:
    print(u.data[0][0][0][:], file=f)

# plot_shotrecord(rec.data, model, t0, tn)
print(u.data.shape)
# print(u.data[0])
