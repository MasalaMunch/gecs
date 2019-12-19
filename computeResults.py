import numpy as np
import random
import time
import json
import cvxopt
cvxopt.solvers.options['show_progress'] = False

timestep = 1/30
inverseTimestep = 1/timestep

class Gec:

    def __init__(self, mass=None, radius=None, color=None, 
                 position=None, velocity=(0, 0), acceleration=(0, 0)):

        assert mass is not None
        self.inverseMass = 1 / mass

        assert radius is not None
        self.radius = radius

        assert color is not None
        self.color = color

        assert position is not None
        self.position = np.array(position).astype(float)

        self.velocity = np.array(velocity).astype(float)

        self.acceleration = np.array(acceleration).astype(float)

    def advance(self):

        self.position += self.velocity * timestep
        self.velocity += self.acceleration * timestep

    def applyImpulse(self, impulseVector):

        self.velocity += impulseVector * self.inverseMass

class GecPair:

    def __init__(self, aGec, bGec):

        self.aGec = aGec
        self.bGec = bGec
        self.inverseMassScalar = self.aGec.inverseMass + self.bGec.inverseMass
        self._massScalar = 1 / self.inverseMassScalar
        self._radiiSum = self.aGec.radius + self.bGec.radius
        self.advance()

    def advance(self):

        centerDifference = self.bGec.position - self.aGec.position
        centerDistance = np.linalg.norm(centerDifference)
        self._direction = centerDifference / centerDistance
        self._distanceComponentOfVelocityCushion = (
            (centerDistance - self._radiiSum) * inverseTimestep
            )
        self._appliedImpulse = 0

    def getVelocityCushion(self):
        # if negative: how much you need to fix it in velocity units
        # if positive: how much extra velocity

        return (
            # distance component: vel needed to colide ina time stop

            self._distanceComponentOfVelocityCushion + 
            np.dot(self.bGec.velocity - self.aGec.velocity, self._direction)
            #^ negative => moving closer to each other, 
            #  positive => moving farther away from each other
            )

    def applyImpulse(self, impulseScalar):

        impulseVector = self._direction * impulseScalar
        self.aGec.applyImpulse(impulseVector)
        self.bGec.applyImpulse(-impulseVector)
        self._appliedImpulse += impulseScalar

    def fixVelocity(self):
    #^ used to compute approximate results

        self.applyImpulse(min(
            self._massScalar * self.getVelocityCushion() + self._appliedImpulse, 
            0.0,
            ) - self._appliedImpulse)

    def getImpulseToVelocityScalar(self, other, sharedGec):
    #^ used to compute optimal results

        return (
            (1.0 if ((sharedGec == self.aGec and sharedGec == other.aGec) 
            or (sharedGec == self.bGec and sharedGec == other.bGec)) 
            else -1.0) * 
            sharedGec.inverseMass * 
            np.dot(self._direction, other._direction)
            )

def getGecPairs(gecs):
    return tuple(
        GecPair(gecs[i], gecs[j]) 
        for i in range(len(gecs)) 
        for j in range(i+1, len(gecs))
        if gecs[i].inverseMass != 0 or gecs[j].inverseMass != 0
        )

resolution = (1280, 720)

def getTestGecs():
    white = (1.0, 1.0, 1.0, 1.0)
    purple = (0.9, 0.0, 0.9, 1.0)
    green = (0.0, 0.9, 0.0, 1.0)
    bigRadius = resolution[0]/4
    bigOffset = resolution[0]/8
    gecs = [
        Gec(mass=float("Infinity"), radius=bigRadius-bigOffset, color=white,
            position=(resolution[0]/2, 0)),
        ]
    radius = 15
    for i in range(int((resolution[1]-(bigRadius-bigOffset))/(2*radius))):
        gecs.append(Gec(mass=1, radius=radius, color=purple, 
            position=(resolution[0]/2, (bigRadius-bigOffset)+(i+1)*(2*radius+75)),
            acceleration=(0, -250),
            ))
    return gecs

def storeResults(gecs, stepGecPositions, testName, computeTime):
    with open("results/" + testName + ".json", "w") as resultFile:
        json.dump({
            "computeTime": computeTime,
            "gecCount": len(gecs),
            "gecRadii": [g.radius for g in gecs],
            "gecColors": [g.color for g in gecs],
            "timestep": timestep,
            "stepCount": simulationSteps,
            "resolution": resolution,
            "stepGecPositions": stepGecPositions,
            }, resultFile)

simulationTime = 20
simulationSteps = int(simulationTime / timestep)

testName = "coarseApproximation"
startTime = time.time()
gecs = getTestGecs()
gecPairs = getGecPairs(gecs)
stepGecPositions = []
for i in range(simulationSteps):
    print(str(round(i/simulationSteps*100, 2))+"%", "through", testName)
    stepGecPositions.append([])
    for g in gecs:
        stepGecPositions[-1].append(tuple(g.position))
    for j in range(1):
    #^ for the coarsest approximation, only do 1 iteration
        for p in gecPairs:
            p.fixVelocity()
    for g in gecs:
        g.advance()
    for p in gecPairs:
        p.advance()
storeResults(gecs, stepGecPositions, testName, time.time()-startTime)

testName = "fineApproximation"
startTime = time.time()
gecs = getTestGecs()
gecPairs = getGecPairs(gecs)
stepGecPositions = []
for i in range(simulationSteps):
    print(str(round(i/simulationSteps*100, 2))+"%", "through", testName)
    stepGecPositions.append([])
    for g in gecs:
        stepGecPositions[-1].append(tuple(g.position))
    for j in range(4):
    #^ for the finer approximation, do 4 iterations
        for p in gecPairs:
            p.fixVelocity()
    for g in gecs:
        g.advance()
    for p in gecPairs:
        p.advance()
storeResults(gecs, stepGecPositions, testName, time.time()-startTime)

testName = "optimal"
# the part done in gams
startTime = time.time()
gecs = getTestGecs()
gecPairs = getGecPairs(gecs)
gecPairsWithGec = {
    g: tuple((i, p) for i, p in enumerate(gecPairs) if p.aGec == g or p.bGec == g)
    for g in gecs
    }

# 7 is arbitrary
q = cvxopt.matrix(7.0, (len(gecPairs), 1))

# 0 is not
h = cvxopt.matrix(0.0, (2*len(gecPairs), 1))
P = cvxopt.spmatrix([], [], [], (len(gecPairs), len(gecPairs)))
G = cvxopt.spmatrix([], [], [], (len(gecPairs)*2, len(gecPairs)))
for i, p in enumerate(gecPairs):
    P[i, i] = 2*p.inverseMassScalar

    # div by mass converts impuse into velocity change
    G[i, i] = -p.inverseMassScalar
    G[i+len(gecPairs), i] = -1.0
stepGecPositions = []
for i in range(simulationSteps):
    print(str(round(i/simulationSteps*100, 2))+"%", "through", testName)
    stepGecPositions.append([])
    for g in gecs:
        stepGecPositions[-1].append(tuple(g.position))
    for j, p in enumerate(gecPairs):
        cushion = p.getVelocityCushion()
        q[j] = cushion
        h[j] = cushion
        for k, otherP in gecPairsWithGec[p.aGec]:
            if j != k:
                scalar = p.getImpulseToVelocityScalar(otherP, p.aGec)
                doubleScalar = 2*scalar
                P[j, k] = doubleScalar
                P[k, j] = -doubleScalar
                G[j, k] = -scalar
                G[k, j] = scalar
        for k, otherP in gecPairsWithGec[p.bGec]:
            if j != k:
                scalar = p.getImpulseToVelocityScalar(otherP, p.bGec)
                doubleScalar = 2*scalar
                P[j, k] = doubleScalar
                P[k, j] = -doubleScalar
                G[j, k] = -scalar
                G[k, j] = scalar
    for p, impulse in zip(gecPairs, cvxopt.solvers.qp(P, q, G=G, h=h)["x"]):
        p.applyImpulse(-impulse)
    for g in gecs:
        g.advance()
    for p in gecPairs:
        p.advance()
storeResults(gecs, stepGecPositions, testName, time.time()-startTime)