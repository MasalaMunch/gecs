import sys
import json
import pyglet
import pygletDraw

with open(sys.argv[1]) as resultFile:
    result = json.load(resultFile)

stage = pyglet.window.Window(*result["resolution"], vsync=True)

stepIndex = 0

@stage.event
def on_draw():

    stage.clear()

    for i in range(result["gecCount"]):

        x, y = result["stepGecPositions"][stepIndex][i]

        pygletDraw.Circle(x=x, y=y, width=result["gecRadii"][i]*2, 
            color=result["gecColors"][i]).render()

def advanceStep(_):

    global stepIndex

    if stepIndex == result["stepCount"]-1:
        pyglet.app.exit()
    else:
        stepIndex += 1

pyglet.clock.schedule_interval(advanceStep, result["timestep"])

pyglet.app.run()