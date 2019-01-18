# This example demonstrates how to rotate a camera in SAMSON,
# capture screenshots, and create a movie

# See: https://documentation.samson-connect.net/scripting-guide/using-camera-and-producing-animations/

# Open any molecule

import os
import numpy as np
from math import pi, log10
import gif_using_imageio
#import gif_using_moviepy

def rotate(camera, velocity, center):
        '''
        rotate the camera with a given velocity and update the viewport
        '''
        camera.rotate(velocity, center)
        SAMSON.requestViewportUpdate()
        SAMSON.processEvents()


camera = SAMSON.getActiveCamera()                                           # get the active camera
camera.center()                                                             # centers the camera
#camera.frontView()                                                          # set the view of the camera
SAMSON.processEvents()                                                      # process events: updates viewport

indexer = SAMSON.getNodes('n.t a')                                          # get all atoms

centerOfMass = np.zeros(3)                                                  # compute center of mass
for a in indexer:
    centerOfMass[0] += a.getX().value;
    centerOfMass[1] += a.getY().value;
    centerOfMass[2] += a.getZ().value;
    
centerOfMass /= indexer.size

center = Type.length3(Quantity.length(centerOfMass[0]), Quantity.length(centerOfMass[1]), Quantity.length(centerOfMass[2]))  # center of the molecule

velocity = pi * 1./180.;                                                    # rotation velocity
velocity3 = Type.radianPerSecond3(Quantity.radPerS(0), Quantity.radPerS(velocity), Quantity.radPerS(0)) # rotation is only in y-direction

pathToFile = '/tmp/camera/'
try:                                                                        # try to create a tmp folder
    os.mkdir(pathToFile)
except:
    pass

ext = '.png'                                                                    # extension of files
numberOfRorations = 36                                                          # rotate 36 times

SAMSON.showProgressBar('Rotating and capturing', 0, numberOfRorations)          # show a progress bar to inform user about the progress

for i in range(numberOfRorations):
        rotate(camera, velocity3, center)                                       # rotate the camera around the center of mass with a given velocity
        i_str = str(i).zfill(int(log10(numberOfRorations)) + 1)
        filename = pathToFile + i_str + ext
        SAMSON.captureViewportToFile(filename, 800, 800)                        # save a capture of the viewport in a file with 800x800 resolution
        SAMSON.setProgressBarValue(i)                                           # update the progress bar
        if SAMSON.isProgressBarStopped(): break                                 # break if the rotation is finished

os.chdir(pathToFile)                                                            # change directory to the directory with images
img_files = sorted((fn for fn in os.listdir('.') if fn.endswith('.png')))       # get a list of created images

gif_using_imageio.create_gif(img_files, 0.1)                                    # create a gif and mp4 files using imageio
#gif_using_moviepy.create_gif(img_files, 12)                                    # create a gif and mp4 files using moviepy

SAMSON.hideProgressBar()
