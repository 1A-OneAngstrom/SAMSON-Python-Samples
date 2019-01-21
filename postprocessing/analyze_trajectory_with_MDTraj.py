'''
This example shows how to use MDTraj package to analyze MD trajectories, 
for example, compute RMSD and the radius of gyration

See: https://documentation.samson-connect.net/scripting-guide/using-mdtraj-in-samson/
'''

import os
import mdtraj

import export_trajectory

class HaltException(Exception): pass

def load_trajectory_in_MDTraj(trajectory_files):
    '''
    Loads trajectories in MDTraj
    '''
    traj = mdtraj.load(trajectory_files[0])                             # load the first trajectory
    for i in range(1, len(trajectory_files)):                           # loop over trajectories
        traj += mdtraj.load(trajectory_files[i])                        # append new trajectory to the previous one
        
    return traj
    
def compute_RMSD(trajectory):
    '''
    Compute the RMSD
    '''
    rmsd = mdtraj.rmsd(traj[1:], traj[0])                               # get RMSD: 1st arg - target frames, 2nd arg - reference frame
    
    for i in range(rmsd.size):
        print("RMSD [frame: {0}] = {1} nm".format(i, rmsd[i]))
        
    return rmsd
        
def compute_Rgyr(trajectory):
    '''
    Compute the radius of gyration
    '''
    rgyr = mdtraj.compute_rg(trajectory)                                # compute radius of gyration for each frame
    for i in range(rgyr.size):
        print("Rgyr [frame: {0}] = {1} nm".format(i, rgyr[i]))
        
    return rgyr


# Get path to the script folder
path_to_pyscript = os.path.realpath(__file__)
path_to_files = os.path.dirname(path_to_pyscript) + '/tmp'
try: os.mkdir(path_to_files)                                            # try to create a tmp folder
except: HaltException('Could not create a folder')

# Open a molecule with a path
SAMSON.importFromFile(os.path.dirname(path_to_pyscript) + '/../samples/1VPK-path.sam')
SAMSON.processEvents()

sbpaths = SAMSON.getNodes('n.t path')                                   # get all paths from the active document in SAMSON

if sbpaths.size == 0:
    raise HaltException('No paths found in the active document')
    
sbpath = sbpaths[0]                                                     # get the first Path node

trajectory_files = export_trajectory.export_to_PDB(sbpath, path_to_files + '/traj-')

traj = load_trajectory_in_MDTraj(trajectory_files)                      # read from a list of trajectories

rmsd = compute_RMSD(traj)
rgyr = compute_Rgyr(traj)

# plot using matplotlib
from IPython import get_ipython
ipython = get_ipython()
ipython.magic('matplotlib inline')      # make matplotlib plot inlined into Jupyter QtConsole            

import matplotlib.pyplot as plt

plt.plot(range(1, rmsd.size+1), rmsd)
plt.ylabel('RMSD, nm')
plt.xlabel('steps')
plt.grid()
plt.show()

plt.plot(rgyr)
plt.ylabel('Rgyr, nm')
plt.xlabel('steps')
plt.grid()
plt.show()
