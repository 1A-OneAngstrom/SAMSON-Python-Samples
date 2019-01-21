'''
This example demonstrates how to make an animation of a path in SAMSON
by capturing screenshots and creating a movie

See: https://documentation.samson-connect.net/scripting-guide/using-camera-and-producing-animations/

Open a molecule in SAMSON
'''

import os
import numpy as np
from math import pi, log10
import gif_using_imageio
#import gif_using_moviepy

class HaltException(Exception): pass

# Get path to the script folder
path_to_pyscript = os.path.realpath(__file__)
path_to_files = os.path.dirname(path_to_pyscript) + '/tmp/path/'
try: os.makedirs(path_to_files)                                                 # try to create a tmp folder
except: HaltException('Could not create a folder')

camera = SAMSON.getActiveCamera()                                               # get the active camera
#camera.frontView()                                                             # set the view of the camera
SAMSON.processEvents()                                                          # process events: updates viewport

paths = SAMSON.getNodes('n.t path')                                             # get the path
path = paths[0]

ext = '.png'                                                                    # extension of files
path.forwardFlag = True                                                         # set a change in positions in a forward direction
path.animationFlag = True                                                       # set an animation
path.animationType = sam.Modeling.StructuralModel.Path.AnimationType.Loop       # set Loop type of an animation

SAMSON.showProgressBar('steps', 0, path.numberOfSteps)                          # show a progress bar to inform user about the progress

for i in range(path.numberOfSteps):
        path.updateState()                                                      # update the state - change positions to the next step
        SAMSON.processEvents()
        i_str = str(i).zfill(int(log10(path.numberOfSteps)) + 1)
        filename = path_to_files + i_str + ext
        SAMSON.captureViewportToFile(filename, 800, 600)                        # save a capture of the viewport in a file with 800x600 resolution
        SAMSON.setProgressBarValue(i)                                           # update the progress bar
        if SAMSON.isProgressBarStopped(): break                                 # break if the rotation is finished

os.chdir(path_to_files)                                                         # change directory to the directory with images
img_files = sorted((fn for fn in os.listdir('.') if fn.endswith('.png')))       # get a list of created images

gif_using_imageio.create_gif(img_files, 0.1)                                    # create a gif and mp4 files using imageio
#gif_using_moviepy.create_gif(img_files, 12)                                    # create a gif and mp4 files using moviepy

SAMSON.hideProgressBar()
