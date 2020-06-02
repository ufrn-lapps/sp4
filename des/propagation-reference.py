import logging
import matplotlib.pyplot as plt
import numpy as np
from rich.logging import RichHandler
from rich.progress import track
import sys

from deswaves import Ricker

sys.setrecursionlimit(50000)


if __name__ == "__main__":

    # Create logger
    logging.basicConfig(level=logging.INFO, format="%(message)s", datefmt="[%X]", handlers=[RichHandler()]) 
    logger = logging.getLogger("Propagation Reference")

    # --------------------------------- PARAMETERS ----------------------------------
    D = 1       # Element Size
    NX = 500    # Grid X dimension size
    NY = 500    # Grid Y dimension size
    Order = 2   # this is all we have at the moment - Space Order ? 
    c = 3000.0  # Wave velocity in m/sec
    f_max = 80  # Source max frequency in Hz
    # -------------------------------------------------------------------------------

    Padding = max(1, Order/2)
    lambda_min = c/(2*np.pi*f_max)
    L = min((NX-1)*D, (NY-1)*D)
    t0 = 0.04
    Tend = L/2/c+t0

    logger.info("propagation speed = %g m/s" % (c))
    logger.info("max. frequency in source = %g Hz" % (f_max))
    logger.info("min. wave length = %g m" % (lambda_min))
    logger.info("element size = %g m" % (D))
    logger.info("grid = [ %d  x %d ]" % (NX, NY))
    logger.info("padding = %d " % (Padding))
    logger.info("dimension = [ %g m  x %g m ]" % ((NX-1)*D, (NY-1)*D))
    logger.info("time to reach boundary = %g s" % (Tend))
    logger.info("dt CFL = %g s" % (D/c))

    assert lambda_min > D

    src = Ricker(f_max/2.7, sampleRate=1./f_max*.05, center=t0)
    logger.info("source %s built. value(t=0 sec) = %g." % (src, src.getValue(0)))
    assert abs(src.getValue(0)) < 1e-2


    dt = 0.2*min(D/c, 1./f_max*.05)
    # Source position: Middle of the grid
    iS = int(Padding + NX/2) 
    jS = int(Padding + NY/2)

    T = 0
    C = (dt*c/D)**2
    TT = []
    UU0 = []
    UU1 = []
    UU2 = []
    U = np.zeros((NX+2*Padding, NY+2*Padding))
    U1 = np.zeros((NX+2*Padding, NY+2*Padding))

    FROG = True
    if FROG:
        V = np.zeros((NX+2*Padding, NY+2*Padding))
        A = np.zeros((NX+2*Padding, NY+2*Padding))
        V1 = np.zeros((NX+2*Padding, NY+2*Padding))
        A1 = np.zeros((NX+2*Padding, NY+2*Padding))

        for T in track(np.arange(0, Tend, dt), description="Propagating..."):

            # Predictor step
            for j in range(1, NY):
                for i in range(1, NX):
                    U[i, j] = U1[i, j]+dt*V1[i, j]+dt**2/2*A1[i, j]
                    A[i, j] = -c**2/D**2 * (4*U[i, j]-U[i-1, j]-U[i+1, j]-U[i, j-1]-U[i, j+1])


            q = src.getValue(T+dt)
            A[iS, jS] += c**2*q

            for j in range(1, NY):
                for i in range(1, NX):
                    V[i, j] = V1[i, j] + dt/2*(A1[i, j]+A[i, j])

            TT.append(T)
            UU0.append(U[iS, jS])
            UU1.append(U[int(Padding+NX/2-5), int(Padding+NY/2)])
            UU2.append(U[int(Padding+NX/2-NX/10), int(Padding+NY/2)])

            logging.debug(T, q, U[iS, jS], U[int(Padding+NX/2-5), int(Padding+NY/2)], U[int(Padding+NX/2-NX/10), int(Padding+NY/2)])

            # Save values for next iteration
            U1, V1, A1 = U, V, A
        pass
    else:
        U1 = np.zeros((NX+2*Padding, NY+2*Padding))
        U2 = np.zeros((NX+2*Padding, NY+2*Padding))
        while T < Tend:

            for j in range(1, NY):
                for i in range(1, NX):
                    R = C*(4*U1[i, j]-U1[i-1, j]-U1[i+1, j]-U1[i, j-1]-U1[i, j+1])
                    U[i, j] = 2*U1[i, j]-U2[i, j]-R

            T += dt
            q = src.getValue(T)
            U[iS, jS] += (dt*c)**2*q

            TT.append(T)
            UU0.append(U[iS, jS])
            UU1.append(U[int(Padding+NX/2-5), int(Padding+NY/2)])
            UU2.append(U[int(Padding+NX/2-NX/10), int(Padding+NY/2)])

            print(T, q, U[iS, jS], U[int(Padding+NX/2-5), int(Padding+NY/2)],
                  U[int(Padding+NX/2-NX/10), int(Padding+NY/2)])

            U1, U2 = U, U1


    if True:
        import matplotlib.pyplot as plt
        plt.clf()
        # ,  linewidths=0.5, linestyles='solid',  colors='black')
        plt.plot(TT, UU0, label='pos1')
        # ,  linewidths=0.5, linestyles='solid',  colors='black')
        plt.plot(TT, UU1, label='pos2')
        # ,  linewidths=0.5, linestyles='solid',  colors='black')
        plt.plot(TT, UU2, label='pos3')
        plt.xlabel('p')
        plt.ylabel('t [sec]')
        plt.legend()
        plt.savefig("plot.png")
        logger.info("plot.png created.")
        plt.clf()
        plt.imshow(U[:, :])
        plt.colorbar()
        plt.savefig("image.png")
        logger.info("image.png created.")
