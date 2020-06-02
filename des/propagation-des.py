import logging
import matplotlib.pyplot as plt
import numpy as np
from rich.logging import RichHandler
from rich.progress import track
import sys

from deswaves import Ricker, WaveObserver, Wave2DOrder2, DESScheduler

# Increase recursion limit. Default: 10000
sys.setrecursionlimit(50000)

if __name__ == "__main__":

    # Create logger
    logging.basicConfig(level=logging.INFO, format="%(message)s", datefmt="[%X]", handlers=[RichHandler()]) 
    logger = logging.getLogger("DES")


    # --------------------------------- PARAMETERS ----------------------------------
    D = 0.8             # Element size (Cell size ?)
    nx = 500            # Grid elements in the X dimension
    ny = 500            # Grid elements in the Y dimension

    order = 2           # space order
    c = 3000.           # propagation speed [m/sec]
    f_max = 80          # Source max frequency [hz]
    # -------------------------------------------------------------------------------


    padding = max(1, int(order/2))  # Padding in each dimension/direction 
    lambda_min = c/(2*np.pi*f_max)  # Some source parameter...
    L = min((nx-1)*D, (ny-1)*D)     # Distance until the closest border 
    t0 = 0.04                       # Initial time step - same as reference
    Tend = L/2/c+t0                 # Time to reach closest boundary ??

    logger.info("Propagation speed = %g m/s" % (c))
    logger.info("Max. frequency in source = %g Hz" % (f_max))
    logger.info("Min. wave length = %g m" % (lambda_min))
    logger.info("Element size = %g m" % (D))
    logger.info("Grid = [ %d  x %d ]" % (nx, ny))
    logger.info("Padding = %d " % (padding))
    logger.info("Dimension = [ %g m  x %g m ]" % ((nx-1)*D, (ny-1)*D))
    logger.info("Time to reach boundary = %g s" % (Tend))
    logger.info("dt CFL = %g s" % (D/c))

    assert lambda_min > D

    src = Ricker(f_max/2.7, sampleRate=1./f_max*.1, center=t0)
    logger.info("source %s built. value(t=0 sec) = %g." % (src, src.getValue(0)))
    assert abs(src.getValue(0)) < 1e-2


    # Plot source data
    if True:
        plt.clf()
        tt = np.arange(0, Tend, src.getSampleRate(0.))
        plt.plot(tt, [src.getValue(t) for t in tt])
        plt.xlabel('q')
        plt.ylabel('t [sec]')
        plt.savefig("source.png")
        logger.info("source.png created.")


    obs0 = WaveObserver(name="p@0")
    obs1 = WaveObserver(name="p@5")
    obs2 = WaveObserver(name="p@10")

    # Create a grid of events
    eventArray = []
    for ix in range(0, nx+2*padding):
        eventArray.append([])
        for iy in range(0, ny+2*padding):
            interior = padding <= ix and ix < padding + \
                nx and padding <= iy and iy < padding+ny
            eventArray[-1].append(Wave2DOrder2(ix, iy, c=c,
                                               d0=D, interior=interior))

    # Turn them into a numpy array:
    eventArray = np.array(eventArray)

    # Add a source at the center of the domain:
    eventArray[int(padding+nx/2), int(padding+ny/2)].addSource(src)
    eventArray[int(padding+nx/2), int(padding+ny/2)].addObserver(obs0)
    eventArray[int(padding+nx/2-5), int(padding+ny/2)].addObserver(obs1)
    eventArray[int(padding+nx/2-10), int(padding+ny/2)].addObserver(obs2)

    # open the DES schedule
    des = DESScheduler(eventArray, numStages=1, thresholds=np.array([1e-6]))

    des.run(T_end=Tend)
    logger.info("T_end=%e reached. All done." % Tend)

    if True:
        plt.clf()
        # ,  linewidths=0.5, linestyles='solid',  colors='black')
        plt.plot(obs0.t, obs0.p, label=obs0.name)
        # ,  linewidths=0.5, linestyles='solid',  colors='black')
        plt.plot(obs1.t, obs1.p, label=obs1.name)
        # ,  linewidths=0.5, linestyles='solid',  colors='black')
        plt.plot(obs2.t, obs2.p, label=obs2.name)
        plt.xlabel('p')
        plt.ylabel('t [sec]')
        plt.legend()
        plt.savefig("plot.png")
        logger.info("plot.png created.")
        plt.clf()
        plt.imshow(des.U[:, :, 0])
        plt.colorbar()
        plt.savefig("image.png")
        logger.info("image.png created.")
