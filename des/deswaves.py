import numpy as np
import sys

import logging

Dlogger = logging.getLogger('DES')
Wlogger = logging.getLogger('DES.Waves')


class DESSourceTemplate(object):
    """
    template for source term whch can be added to DESEvent
    """

    def __init__(self, **kwargs):
        """
        """
        pass

    def __repr__(self):
        """
        a nice representation of the DESEvent (can be overwritten)
        """
        return "DESSource"

    def getSampleRate(self, t):
        """
        Returns the rate at which the source term is sampled. (needs to be overwritten)
        """
        raise NotImplemented

    def getValue(self, t):
        """
        get value of wavelet at time `t` (needs to be overwritten)
        """
        raise NotImplemented

    def getIntegral(self, t0, t1):
        """
        get value of wavelet at time `t` (needs to be overwritten)
        """
        raise NotImplemented


class Ricker(DESSourceTemplate):
    """
    The Ricker Wavelet w=f(t)
    """

    def __init__(self, dominateF=40, center=None, sampleRate=None):
        """
        Sets up a Ricker wavelet with dominant frequence `dominateF` [Hz] and
        center at time `center` [sec]. If `center` is not given an estimate
        for suitable `center` is calculated so f(0)~0.

        :note: maximum frequence is about 2 x the dominant frequence.
        """
        drop = 2.5
        self.maxF = np.sqrt(7.)*dominateF
        self.dominateF = dominateF

        if sampleRate:
            self.sampleRate = sampleRate
        else:
            self.sampleRate = 1./self.maxF

        self.__s = np.pi*self.dominateF
        if center == None:
            center = drop**2/self.__s
        self.center = center
        self.width = 2*np.sqrt(6)/np.pi/self.dominateF

    def __repr__(self):
        """
        a nice representation of the DESEvent (can be overwritten)
        """
        return "Ricker(%g hz@%g sec, %g ms)" % (self.dominateF, self.center, self.sampleRate*1000.)

    def getSampleRate(self, t):
        """
        Returns the rate at which the source term is sampled.
        """
        if t > self.center+self.width/2:
            return sys.float_info.max
        else:
            return self.sampleRate

    def getValue(self, t):
        """
        get value of wavelet at time `t`
        """
        e2 = (self.__s*(t-self.center))**2
        return (1-2*e2)*np.exp(-e2)

    def getIntegral(self, t0, t1):
        """
        get integral of wavelet from time `t0` to `t1`
        """
        I1 = np.exp((self.__s*(t1-self.center))**2)*(t1-self.center)
        I0 = np.exp((self.__s*(t0-self.center))**2)*(t0-self.center)

        return I1-I0


class DESObserverTemplate(object):
    """
    class template for  point value observer
    """

    def __init__(self, name="<none>"):
        self.name = name

    def addValue(self, t, v):
        raise NotImplemented


class DESEventTemplate(object):
    """
    this is a generic DES Event. On observer can be attached to the event to track the 
    value of the event over time. The solution is updated by a predictor-corrector type scheme.
    The predictor is given as 

       U^p_n = U^0  = U_{n-1} + dt * r_0 

    with 

          r_0= R(t_{n-1}, U_{n-1})

    followed by multi-stage corrector step:

        for i=1,...,numStages:
            r_i= R(t_{n-1}+c_i*dt,  U_{i-1})
            U^{i} = U_{n-1}+dt * ( a_{i,0} r_0 + a_{i,1} r_1 + ...+ a_{i,i-1} r_{i-1})           

        with U_n = U^{numStages}

    Notice that it is assumed that in the predictor step the rate r_0 does not require any knowledge of dt.
    Also, as the rate calculation involves a spatial stencil in some form the calculation of U^{i} in stage i need to be completed for all
    grid points before the next stage can be started (when threading is used a synchronization of threads need to be applied.)

    """

    def __init__(self, i0, i1, i2=None, T0=0.):
        """
        creates a DES Event connected with a grid point (i0,i1 [,i2]) and initial time T0.
        """
        self.has_valid_scheduled_time = False
        self.event_in_eternity = False          # Event that last forever? 
        self.is_synchronized = False
        self.inPEP = False
        self.i0 = i0                            # Grid point location
        self.i1 = i1
        self.i2 = i2
        self.dt_CFL = sys.float_info.max        # The local save time step size
        self.t_current = T0                     # Time of event
        self.dU = 0.                            # Accumulated solution increment
        self.observer = None                    # An observer (if available)
        self.source = None

    def addObserver(self, observer=None):
        """
        adds an observer to the event. The `observer` needs to be None or an instance of `DESObserverTemplate`
        """
        if observer:
            assert isinstance(observer, DESObserverTemplate)
        self.observer = observer

    def __repr__(self):
        """
        a nice representation of the DESEvent (can be overwritten)
        """
        return "DESEvent@(%s,%s)" % (self.i0, self.i1)

    def addSource(self, source=None):
        """
        adds a source `source` to the event. If `source` is None the source is removed.
        `source` needs to be an instance of `DESSourceTemplate`.
        """
        if source:
            assert isinstance(source, DESSourceTemplate)
            Wlogger.info("Source %s added to %s." % (source, self))
        self.source = source
        return self

    def synchronize(self, U, V, A, pepList):
        """
        this updates the predicted solution for all naighbours. if the predicted change 
        is greater than the target change dU_target synchronization for this naigbour is also triggered.  
        """
        self.has_valid_scheduled_time = False
        self.is_synchronized = True
        self.dU = 0.
        for e in self.neigbours:
            if not e.inPEP:
                e.makeSolutionPrediction(self.t_current, U, V, A, pepList)
                # print self," : ", e, " is updated. dU, dU_target=",abs(e.dU), e.dU_target, abs(e.dU) >e.dU_target
                if not e.is_synchronized and (abs(e.dU) > e.dU_target).any():
                    # print self," : synconize ", e,
                    e.synchronize(U, V, A, pepList)

    def makeSolutionPrediction(self, T_clock, U, V, A,  pepList):
        """
        applies a prediction step to update the solution U <- U+dt*R

        this may be corrected through a multi-stage corrector step.
        """
        self.dt = T_clock-self.t_current
        DU = V[self.i0, self.i1]*self.dt+A[self.i0, self.i1]*self.dt**2/2
        U[self.i0, self.i1] += DU
        self.dU += DU
        # print "SS: ",self,U[self.i0,self.i1], self.dU, self.has_valid_scheduled_time
        self.t_current = T_clock
        self.inPEP = True
        pepList.append(self)

        # print "added to PEP ",self, U[self.i0,self.i1]

    def setScheduledTime(self, T_end, V, A, thresholds, OMEGA_CFL=0.5):
        """
        sets the new scheduled time for the event and marks the scheduled time as valid. 
        If the scheduled time is bigger than the T_end the event is maked as inactive

        (was called Schedule)
        """
        dt = OMEGA_CFL*self.dt_CFL
        self.dU_target = dt*abs(V[self.i0, self.i1]+dt/2*A[self.i0, self.i1])
        if self.source:
            dtmax = self.source.getSampleRate(self.t_current)
        else:
            dtmax = sys.float_info.max
        # print self, self.dU_target, R[self.i0,self.i1]
        if (self.dU_target < thresholds).all():
            self.dt_target = dtmax
            self.dU_target = thresholds
        else:
            self.dt_target = min(dtmax, dt)
        self.t_scheduled = self.t_current+self.dt_target
        self.has_valid_scheduled_time = True
        if self.t_scheduled > T_end:
            self.event_in_eternity = True
        else:
            self.event_in_eternity = False
        return self

    def computeNewRate(self, U, V, A, first=False):
        """
        by default we assume no change but this should calculate a new rate of change and a corresponding CFL consistemt time step size
        this function needs to be overwritten
        Note: would like to run this function as a numpy operation
        """
        self.V[self.i0, self.i1] = 0
        self.A[self.i0, self.i1] = 0
        self.dt_CFL = sys.float_info.max

    def setNeigbours(self, eventArray):
        """
        select the naigbour  from a grid. This method is called only once.
        this does nothing and needs to be overwritten.
        """
        self.neigbours = []

    def updateSolution(self, stage, U, V, A):
        """
        this applies the stage `stage` update to the solution U. 
        at input U is solution 
        """
        pass


class WaveObserver(DESObserverTemplate):
    """
    class template for  point value observer
    """

    def __init__(self, dt=None, t0=0., name="WaveObserver"):
        self.name = name
        self.t = []
        self.p = []
        self.dt = dt
        if self.dt:
            self.tnext = t0+dt

    def addValue(self, t, u):
        if self.dt:
            if t >= self.tnext:
                self.tnext += self.dt
                self.t.append(t)
                self.p.append(u[0])
        else:
            self.t.append(t)
            self.p.append(u[0])


class Wave2DOrder2(DESEventTemplate):
    """
    this is a generic desevent connected with a grid point (i0,i1)
    """

    def __init__(self, i0, i1, c=1., d0=1., interior=True, T0=0, alpha=1):
        """
        as grid point '(i0,i1)' with cell size `d0` and initial time `t0`.
        `c` is the wave propagation speed at grid point `c`.
        if `interior` is set the grid point is updated. 
        `alpha` controles the time integration scheme with 
        alpha=1 -> Heun's method 
        alpha=0.5 -> midpoint method (not supported yet, needs additional R)
        alpha=2./3. -> two-stage second-order Runge-Kutta (not supported yet, needs additional R) 
        """

        DESEventTemplate.__init__(self, i0, i1, i2=None, T0=T0)
        self.c = c
        self.d0 = d0
        self._interior = interior
        assert alpha == 1
        self.alpha = alpha
        self.R = np.zeros((2,))
        self.dt_CFL_0 = 0.2*self.d0/self.c

    def __repr__(self):
        """
        a nice representation of the Event
        """
        if self._interior:
            if self.source:
                return "Wave2DOrder2@(%s, %s, s=%s)" % (self.i0, self.i1, self.source)
            else:
                return "Wave2DOrder2@(%s, %s)" % (self.i0, self.i1)

        else:
            return "Wave2DOrder2@Padding(%s,%s)" % (self.i0, self.i1)

    def setNeigbours(self, eventArray):
        """
        select the neigbours from a grid:
        """
        if self._interior:
            self.neigbours = [e for e in [eventArray[self.i0-1, self.i1], eventArray[self.i0+1, self.i1],
                                          eventArray[self.i0, self.i1-1], eventArray[self.i0, self.i1+1]] if e._interior]
        else:
            self.neigbours = []

    def computeNewRate(self, U, V, A, first=False):
        if self._interior:  # Don't update in the padding zone!

            a = -self.c**2*(4*U[self.i0, self.i1, 0]-U[self.i0-1, self.i1, 0]-U[self.i0+1,
                                                                                self.i1, 0]-U[self.i0, self.i1-1, 0]-U[self.i0, self.i1+1, 0])/self.d0**2
            if self.source:
                q = self.source.getValue(self.t_current)
                a += q*self.c**2
            if not first:
                dt = self.dt
                V[self.i0, self.i1, 0] += dt/2*(A[self.i0, self.i1, 0]+a)
            else:
                V[self.i0, self.i1, 0] = 0
            A[self.i0, self.i1] = a
            self.dt_CFL = self.dt_CFL_0
            if self.source:
                self.dt_CFL = min(
                    self.dt_CFL_0, self.source.getSampleRate(self.t_current))
                Wlogger.info("%s: source value at t= %s sec : %g (dt_CFL = %g, u=%g)" % (
                    self, self.t_current, q, self.dt_CFL, U[self.i0, self.i1, 0]))

    def updateSolution(self, stage, U, V, A):
        """
        NOT USED!
        this applies the stage `stage` update to the solution U. 
        at input U and R are solution from stage `stage-1` .
        At the inital all (stage=1) `R` is the rate at self.t_current-1 and U is the predicted solution

        At the output 
        """
        if stage == 1:
            self.R[0] = self.c**2*U[self.i0, self.i1, 1]
            self.R[1] = -(4 * U[self.i0, self.i1, 0] - U[self.i0-1, self.i1, 0] - U[self.i0+1, self.i1, 0] - U[self.i0, self.i1-1, 0] - U[self.i0, self.i1+1, 0]) / self.d0**2
            if self.source:
                q = self.source.getValue(t=self.t_current)
                self.R[1] += q
        elif stage == 2:
            dt = self.dt
            U[self.i0, self.i1] = U[self.i0, self.i1] + dt/2 * (self.R - R[self.i0, self.i1]) # DONT EXIST R, only self.R


class DESBoudaryEvent(DESEventTemplate):
    """
    this a special DESEvent to deal with the boundary to be added later
    """
    pass


class DESScheduler(object):
    def __init__(self, eventArray, thresholds, numStages=1, T_clock=0., OMEGA_CFL=0.9, OMEGA_PEP=0.9):
        self.T_clock = T_clock
        self.OMEGA_CFL = OMEGA_CFL
        self.OMEGA_PEP = OMEGA_PEP
        numSol = len(thresholds)
        Dlogger.info("Number of solution components is %d." % numSol)
        Dlogger.info("thresholds are %s." % str(thresholds))
        Dlogger.info("number of stages %d." % numStages)

        self.eventArray = eventArray
        for E in np.nditer(self.eventArray, flags=["refs_ok"]):
            E.tolist().setNeigbours(self.eventArray)

        self.thresholds = thresholds
        self.V = np.zeros(eventArray.shape+(numSol,))
        self.A = np.zeros(eventArray.shape+(numSol,))
        self.U = np.zeros(eventArray.shape+(numSol,))
        # update model
        # U{n} = U{n-1} + dt * V{n-1} + dt**2 * A_{n-1}
        # A{n},V{n} = f(t_n, U{n}, V{n-1}, A_{n-1})
        # e.g. V{n} = V_{n-1} + dt/2* (A_{n-1}+A_{n})

        self.numStages = numStages

    def run(self, T_end):
        """
        this runs the DESScheduler until `t_current` for all events have passed `T_end` 
        """
        # Filled by self.processPEP
        self.eventQueue = []  

        # Initialize all events by adding them to the pepList:
        self.pepList = []
        for E in np.nditer(self.eventArray, flags=["refs_ok"]):
            e = E.tolist()
            e.has_valid_scheduled_time = False
            e.inPEP = True
            e.event_in_eternity = False
            self.pepList.append(e)

        Dlogger.info("Initial PEP list with %s events is processed." % len(self.pepList))
        self.processPEP(T_end, first=True)  # fills self.eventQueue
        iterCount = 0
        Dlogger.info("Initial eventQueue has %s elements." % len(self.eventQueue))

        # now we start to progress in time until T_end is reached.
        # all events are processed and pepList is empty - nothing to do anymore
        while len(self.eventQueue) > 0:

            # Sort the eventQueue by dt_target set by e.setScheduledTime
            self.eventQueue.sort(key=lambda e: e.dt_target)

            # Pick the event at the top of the list
            self.T_clock = self.eventQueue[0].t_scheduled
            dt_prep = self.OMEGA_PEP*self.eventQueue[0].dt_target
            Dlogger.info("Event %s, T_clock = %g sec" % (iterCount, self.T_clock))
            niter = 0
            while (len(self.eventQueue) > 0 and
                   self.eventQueue[0].t_scheduled <= self.T_clock + dt_prep):

                e_top = self.eventQueue[0]
                Dlogger.debug("%s selected from event queue, scheduled time = %g sec" % (e_top, e_top.t_scheduled))

                e_top.makeSolutionPrediction(self.T_clock, self.U, self.V, self.A, self.pepList)
                e_top.synchronize(self.U, self.V, self.A, self.pepList)
                dt_prep = min(dt_prep, self.OMEGA_PEP * e_top.dt_target)

                # Make sure that we are looking at events that are still valid
                self.eventQueue = [e for e in self.eventQueue if e.has_valid_scheduled_time and not e.inPEP]
                niter += 1

            Dlogger.info("PEP list created in  %s steps" % (niter))
            Dlogger.info("Start to process PEP list with %s elements" % (len(self.pepList)))
            self.processPEP(T_end)
            iterCount += 1
            if iterCount > 30:
                print((abs(self.U) < (self.thresholds[0])).sum(), 
                      (abs(self.U) > (self.thresholds[0])).sum())

    def processPEP(self, T_end, first=False):
        """
        process the events in the pepList.
        the events are removed form the pepList (empty on exit)
        """
        # first update the solution through various stages:
        # at this point R is the still the rate R from the previous times step
        # and U is holding the predictor for T_clock
        if not first:
            for i in range(1, self.numStages):
                # omp parallel for
                for e in self.pepList:
                    e.updateSolution(i, self.U, self.V, self.R)
                # with treading add a barrier synchronization HERE!

        # At this point U is the solution at T_clock and a new rate R(T_clock, U)
        # for the next prediction step is calculated:
        for e in self.pepList:
            # This compute the rates R for a predictor step
            if e.observer:
                e.observer.addValue(e.t_current, self.U[e.i0, e.i1])

            e.computeNewRate(self.U, self.V, self.A, first)

            if not e.has_valid_scheduled_time:
                e.setScheduledTime(T_end, self.V, self.A, self.thresholds, self.OMEGA_CFL)

                if not e.event_in_eternity:
                    self.eventQueue.append(e)

            e.inPEP = False
            e.is_synchronized = False

        # With treading add a barrier synchronization HERE!
        self.pepList = []
