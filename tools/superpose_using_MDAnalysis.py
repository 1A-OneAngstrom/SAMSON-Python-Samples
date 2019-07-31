'''
Let us consider a case, when we have two structures (e.g., two structural models in terms of SAMSON) 
in the active Document in SAMSON and we want to superimpose
in a way that minimizes the RMSD. For that we can use functions in the MDAnalysis.analysis.align module.
Let the first molecule be the reference one, and the second one - one which should be moved. 

See: https://documentation.samson-connect.net/scripting-guide/analyzing-trajectories-with-mdanalysis/
Superposition of structure: https://www.mdanalysis.org/MDAnalysisTutorial/analysismodule.html

'''
import os
import MDAnalysis
from MDAnalysis.analysis import align, rms

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

if indexer.size < 2:
	HaltException('The active document contains less than 2 structural models')
    
# export molecules each in a different pdb-file    
ref = indexer[0]                                                        # get the first structural model    
mob = indexer[1]                                                        # get the second structural model    
ref_indexer = ref.getNodes()                                            # get an indexer of nodes for the first structural model
mob_indexer = mob.getNodes()

filename_ref = path_to_files + '/ref.pdb'
filename_mob = path_to_files + '/mob.pdb'

SAMSON.exportToFile(ref_indexer, filename_ref, '')                      # export the reference structure
                                                                        # the third input parameter is for options used for importing: '' is for default
SAMSON.exportToFile(mob_indexer, filename_mob, '')                      # export the structure which should be rotated

u_ref = MDAnalysis.Universe(filename_ref)
u_mob = MDAnalysis.Universe(filename_mob)

rms.rmsd(u_mob.atoms.CA.positions, u_ref.atoms.CA.positions)

# Note that in this example translations have not been removed.
# In order to look at the pure rotation one needs to superimpose the centres of mass (or geometry) first:
u_ref0 = u_ref.atoms.CA.positions - u_ref.atoms.CA.center_of_mass()     # get a structure shifted by its center of mass
u_mob0 = u_mob.atoms.CA.positions - u_mob.atoms.CA.center_of_mass()     # get a structure shifted by its center of mass

rms.rmsd(u_mob0, u_ref0)

# The rotation matrix that superimposes mob on ref while minimizing the CA-RMSD is obtained with the rotation_matrix() function
R, rmsd = align.rotation_matrix(u_mob0, u_ref0)

# Putting all this together one can superimpose all of mob onto ref:
if 0:
    # do translations and rotations in SAMSON on the mob itself
    u_ref_com = u_ref.atoms.CA.center_of_mass()
    u_mob_com = u_mob.atoms.CA.center_of_mass()
    ref_center_of_mass = Type.vector3(Quantity.angstrom(u_ref_com[0]), Quantity.angstrom(u_ref_com[1]), Quantity.angstrom(u_ref_com[2]))
    mob_center_of_mass = Type.vector3(Quantity.angstrom(u_mob_com[0]), Quantity.angstrom(u_mob_com[1]), Quantity.angstrom(u_mob_com[2]))
    rotation_matrix = Type.matrix33(R.tolist())                             # create a rotation matrix in SAMSON. R is an ndarray and should be converted into list

    mob_indexer_a = mob.getNodes('n.t a')                                   # get an indexer of all atoms in the mob

    for a in mob_indexer_a:
        a.setPosition( a.getPosition() - mob_center_of_mass)                # translate by mob's center of mass
        a.setPosition( rotation_matrix * a.getPosition() )                  # rotate
        a.setPosition( a.getPosition() + ref_center_of_mass)                # translate by ref's center of mass

else:
    # superimpose using MDAnalysis and import the resulting structure in SAMSON
    u_mob.atoms.translate(-u_mob.atoms.CA.center_of_mass())
    u_mob.atoms.rotate(R)
    u_mob.atoms.translate(u_ref.atoms.CA.center_of_mass())
    filename_mob_superposed = path_to_files + '/mob_superposed_on_ref_mdanalysis.pdb'
    u_mob.atoms.write(filename_mob_superposed)

    SAMSON.importFromFile(filename_mob_superposed, '')
    SAMSON.processEvents()
