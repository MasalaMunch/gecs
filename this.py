import numpy as np
import pyglet
import pygletDraw
import random

timestep = 1/60
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

        self.position += self.velocity*timestep
        self.velocity += self.acceleration*timestep

    def applyImpulse(self, impulseVector):

        self.velocity += impulseVector * self.inverseMass

    def render(self, position=None):

        assert position is not None
        pygletDraw.Circle(x=position[0], y=position[1], 
            width=self.radius*2, color=self.color).render()

class GecPair:

    def __init__(self, aGec, bGec):

        self.aGec = aGec
        self.bGec = bGec

    def advance(self):

        centerDifference = self.bGec.position - self.aGec.position
        centerDistance = np.linalg.norm(centerDifference)
        self._direction = centerDifference / centerDistance
        self._distanceComponentOfVelocityCushion = (
            (centerDistance - (self.aGec.radius + self.bGec.radius))
            * inverseTimestep
            )
        self._massScalar = 1 / (self.aGec.inverseMass + self.bGec.inverseMass)
        self._appliedImpulse = 0

    def getVelocityCushion(self):

        return (
            self._distanceComponentOfVelocityCushion + 
            np.dot(self.bGec.velocity - self.aGec.velocity, self._direction)
            #^ (negative => moving closer to each other, 
            #   positive => moving farther away from each other)
            )

    def applyImpulse(self, impulseScalar):

        impulseVector = self._direction * impulseScalar
        self.aGec.applyImpulse(impulseVector)
        self.bGec.applyImpulse(-impulseVector)
        self._appliedImpulse += impulseScalar

    def fixVelocity(self):

        impulse = min(
            self._massScalar*self.getVelocityCushion() + self._appliedImpulse, 
            0.0,
            ) - self._appliedImpulse
        if impulse != 0.0:
            self.applyImpulse(impulse)

    def getImpulseToVelocityScalar(self, otherGecPair=None, sharedGec=None):

        assert otherGecPair is not None
        assert sharedGec is not None
        return (
            (1.0 if ((sharedGec == self.aGec and sharedGec == otherGecPair.aGec) 
            or (sharedGec == self.bGec and sharedGec == otherGecPair.bGec)) 
            else -1.0) * 
            sharedGec.inverseMass * 
            np.dot(self._direction, otherGecPair._direction)
            )

resolution = (1280, 720)

white = (1.0, 1.0, 1.0, 1.0)
purple = (0.9, 0.0, 0.9, 1.0)

def get100gecs():
    gecs = [
        Gec(mass=float("Infinity"), radius=resolution[0]/2, color=white, 
            position=(0, 0)),
        Gec(mass=float("Infinity"), radius=resolution[0]/2, color=white,
            position=(resolution[0], 0)),    
        ]
    minX = resolution[0]/2*np.sqrt(2)/2
    maxX = resolution[0]-minX
    minY = minX
    maxY = resolution[1] + 1000
    random.seed(a=2)
    for i in range(100-len(gecs)):
        gecs.append(Gec(mass=1, radius=5, color=purple, 
            position=(random.uniform(minX, maxX), random.uniform(minY, maxY)),
            acceleration=(0,-200)))
    return tuple(gecs)

def getGecPairs(gecs):
    return tuple(
        GecPair(gecs[i], gecs[j]) 
        for i in range(len(gecs)) 
        for j in range(i+1, len(gecs))
        if gecs[i].inverseMass != 0 or gecs[j].inverseMass != 0
        )

simulationTime = 10
simulationSteps = int(simulationTime / timestep)

approximateGecs = get100gecs()
approximateGecPairs = getGecPairs(approximateGecs)
approximateResults = []
for i in range(simulationSteps):
    print(i/simulationSteps)
    for p in approximateGecPairs:
        p.advance()
    for j in range(5):
        for p in approximateGecPairs:
            p.fixVelocity()
    approximateResults.append([])
    for g in approximateGecs:
        approximateResults[-1].append(np.array(g.position))
    for g in approximateGecs:
        g.advance()

#TODO separate into result computation and storage file and result display file 

if __name__ == "__main__":
    
    stepIndex = 0

    stage = pyglet.window.Window(resolution[0], resolution[1], vsync=True)

    @stage.event
    def on_draw():
        global stepIndex
        stage.clear()
        for i, g in enumerate(approximateGecs):
            g.render(approximateResults[stepIndex][i])

    def advanceStep(dt):
        global stepIndex
        if stepIndex == simulationSteps-1:
            pyglet.app.exit()
        else:
            stepIndex += 1

    pyglet.clock.schedule_interval(advanceStep, timestep)

    pyglet.app.run()