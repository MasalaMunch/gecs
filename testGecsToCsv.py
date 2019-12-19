import numpy as np
import random
import time
import json
import csv

resolution = (1280, 720)

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

        return (
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

gecs = getTestGecs()

with open('testGecsData.csv', mode='w') as csv_file:
    fieldnames = ['', 'inverseMass', 'radius', 'color', 'velocityx', 'velocityy', 'accelerationx', 'accelerationy', 'positionx', 'positiony']
    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
    writer.writeheader()

    for g in gecs:    
        writer.writerow({'':'g1', 'inverseMass': g.inverseMass, 'radius': g.radius, 'color': g.color[0], 
        'velocityx': g.velocity[0], 'velocityy': g.velocity[1], 
        'accelerationx': g.acceleration[0], 'accelerationy': g.acceleration[1], 
        'positionx': g.position[0], 'positiony': g.position[1]})

    # dummy row since gams deletes columns with all-zero values, causing errors
    writer.writerow({'':'gbad', 'inverseMass': 100000000000000000000, 'radius': 0, 'color': 0.9, 
    'velocityx': 0.000001, 'velocityy': 0.000001, 
    'accelerationx': 0.000001, 'accelerationy': 0.000001, 
    'positionx': 0, 'positiony': 0})
