'''
Let us consider a case, when we have two structures (e.g., two structural models in terms of SAMSON) 
in the active Document in SAMSON and we want to superimpose them in a way that minimizes the RMSD.
For that we can use MDTraj module.

Let the first molecule be the reference one, and the second one - the one which should be moved. 

See: https://documentation.samson-connect.net/scripting-guide/using-mdtraj-in-samson/
'''

import os
import mdtraj

class HaltException(Exception): pass

# Get path to the script folder
path_to_pyscript = os.path.realpath(__file__)
path_to_files = os.path.dirname(path_to_pyscript) + '/tmp'
try: os.mkdir(path_to_files)                                            # try to create a tmp folder
except: HaltException('Could not create a tmp folder')

# Open a molecule with
SAMSON.importFromFile(os.path.dirname(path_to_pyscript) + '/../samples/superimpose-sample.sam')
SAMSON.processEvents()

indexer = SAMSON.getNodes('n.t sm')                                     # get indexer of structural models    

# export molecules each in a different pdb-file
if indexer.size < 2:
	HaltException('The active document contains less than 2 structural models')
        
ref = indexer[0]                                                        # get the first structural model    
mob = indexer[1]                                                        # get the second structural model    
ref_indexer = ref.getNodes()                                            # get an indexer of nodes for the first structural model
mob_indexer = mob.getNodes()

filename_ref = path_to_files + '/ref.pdb'
filename_mob = path_to_files + '/mob.pdb'

SAMSON.exportToFile(ref_indexer, filename_ref, '')                      # export the reference structure
                                                                        # the third input parameter is for options used for importing: '' is for default
SAMSON.exportToFile(mob_indexer, filename_mob, '')                      # export the structure which should be rotated

t_ref = mdtraj.load(filename_ref)                                       # load the reference into MDTraj
t_mob = mdtraj.load(filename_mob)                                       # load the molecule which needs to be superposed

t_mob.superpose(t_ref)                                                  # superpose each conformation in the trajectory 'traj' upon a reference

filename_mob_superposed = path_to_files + '/mob_superposed_on_ref_mdtraj.pdb'
t_mob.save(filename_mob_superposed)

SAMSON.importFromFile(filename_mob_superposed, '')                      # import the superposed structure into SAMSON
SAMSON.processEvents()
