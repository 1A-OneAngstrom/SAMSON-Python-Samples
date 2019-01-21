'''
This example shows how to use MDAnalysis package to analyze MD trajectories, 
for example, compute RMSD and the radius of gyration

See: https://documentation.samson-connect.net/scripting-guide/analyzing-trajectories-with-mdanalysis/
'''

import os
import MDAnalysis
import MDAnalysis.analysis.rms

import export_trajectory

class HaltException(Exception): pass

def compute_RMSD(u):
    '''
    Compute the RMSD
    '''
    rmsd_path = []
    bb = u.select_atoms('backbone')                                     # a selection (AtomGroup)
    u.trajectory[0]                                                     # choose the first frame
    frame_0 = bb.positions                                              # coordinates of first frame
    for i in range(1, u.trajectory.n_frames):                           # loop over trajectories
        u.trajectory[i]                                                 # forward to the next frame
        frame_i = bb.positions                                          # coordinates of the current frame
        rmsd = MDAnalysis.analysis.rms.rmsd(frame_0, frame_i)           # get RMSD
        rmsd_path.append(rmsd)
        print("RMSD [frame: {0}] = {1} A".format(i, rmsd))
        
    return rmsd_path

def compute_Rgyr(u):
    '''
    Compute the radius of gyration
    '''
    rgyr_path = []
    bb = u.select_atoms('protein and backbone')                         # a selection (AtomGroup)
    for ts in u.trajectory:                                             # loop over trajectories
        rgyr = bb.radius_of_gyration()                                  # get the radius of gyration
        rgyr_path.append(rgyr)
        print("Rgyr [frame: {0}] = {1} A".format(ts.frame, rgyr))
    
    return rgyr_path


# Get path to the script folder
path_to_pyscript = os.path.realpath(__file__)
path_to_files = os.path.dirname(path_to_pyscript) + '/tmp'
try: os.mkdir(path_to_files)                                            # try to create a tmp folder
except: HaltException('Could not create a folder')

# Open a molecule with a path
SAMSON.importFromFile(os.path.dirname(path_to_pyscript) + '/../samples/1VPK-path.sam')
SAMSON.processEvents()
topology_file = os.path.dirname(path_to_pyscript) + '/../samples/1VPK.pdb'

sbpaths = SAMSON.getNodes('n.t path')                                   # get all paths from the active document in SAMSON

if sbpaths.size == 0:
    raise HaltException('No paths found in the active document')
    
sbpath = sbpaths[0]                                                     # get the first Path node

trajectory_files = export_trajectory.export_to_XYZ(sbpath, path_to_files + '/traj-')

u = MDAnalysis.Universe(topology_file, trajectory_files)                # read from a list of trajectories

rmsd = compute_RMSD(u)
rgyr = compute_Rgyr(u)

# plot using matplotlib
from IPython import get_ipython
ipython = get_ipython()
ipython.magic('matplotlib inline')      # make matplotlib plot inlined into Jupyter QtConsole            

import matplotlib.pyplot as plt

plt.plot(range(1, len(rmsd)+1), rmsd)
plt.ylabel('RMSD, Angstrom')
plt.xlabel('steps')
plt.grid()
plt.show()

plt.plot(rgyr)
plt.ylabel('Rgyr, Angstrom')
plt.xlabel('steps')
plt.grid()
plt.show()
