'''
Export trajectory from a SAMSON path into files
'''

from samson.Facade import SAMSON

def export_trajectory(sbpath, filename, fileformat):
    '''
    Exports trajectories from a Path node 'sbpath' to files with a name starting with 'filename' and extension 'fileformat'
    '''
    sbpath.currentStep = 0                                              # set currentStep for the Path to 0
    trajectory_files = []
        
    indexer = SAMSON.getNodes('n.t sm')                                 # get a node indexer for all atoms
    
    if indexer.size > 0:
        for step in range(sbpath.numberOfSteps):                        # loop over steps in the Path
            sbpath.currentStep = step                                   # increment currentStep
            
            fn = filename + str(step) + '.' + fileformat                # a name of a file 
            trajectory_files.append(fn)                                 # append list of trajectory files
            
            SAMSON.exportToFile(indexer, fn, '')                        # export current trajectory into a file 'fn'
            
    return trajectory_files                                             # return list of trajectory files

def export_to_PDB(sbpath, filename):
    '''
    Export trajectory from path into pdb format
    '''
    return export_trajectory(sbpath, filename, 'pdb')
    
def export_to_XYZ(sbpath, filename):
    '''
    Export trajectory from path into xyz format
    '''
    return export_trajectory(sbpath, filename, 'xyz')
