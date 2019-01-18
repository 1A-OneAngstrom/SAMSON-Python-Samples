from math import sqrt

def compute_rmsd(mol1_indexer, mol2_indexer):
        '''
        Computes RMSD between two conformations of a molecule
        This is not an optimized implementation, it is only for the tutorial purpose
        '''
        if mol1_indexer.size != mol2_indexer.size:                # check if indexers are of the same size
                print('molecules are not the same')
                return -1
                
        rmsd = 0.0;
        
        for i in range(mol1_indexer.size):                        # go through all atoms in the indexer
                if mol1_indexer[i].elementType != mol2_indexer[i].elementType:
                        print('atoms are not the same')
                        return -1
                pos1 = mol1_indexer[i].getPosition()              # get position of an atom from an indexer
                pos2 = mol2_indexer[i].getPosition()
                distance = pos1 - pos2                            # compute distance between atoms
                rmsd += distance.norm2().angstrom                 # accumulate squared norm for distances [in angstrom]
                
        rmsd = sqrt(rmsd / mol1_indexer.size)                     # take a square root of the sum of squared distances divided by the number of atoms
        
        return rmsd

def compute_rmsd_2mol(mol1, mol2):
        '''
        Computes several RMSD between two conformations:
         1) between all atoms; 2) between backbone atoms; 3) between alpha-Carbon atoms
        '''
        mol1_all_indexer      = mol1.getNodes('n.t a')            # get all atoms
        mol2_all_indexer      = mol2.getNodes('n.t a')
        mol1_backbone_indexer = mol1.getNodes('n.t a in n.t bb')  # get all atoms in the backbone
        mol2_backbone_indexer = mol2.getNodes('n.t a in n.t bb')
        mol1_CA_indexer       = mol1.getNodes('"CA" in n.t bb')   # get all alpha Carbons in the backbone
        mol2_CA_indexer       = mol2.getNodes('"CA" in n.t bb')

        rmsd_all      = compute_rmsd(mol1_all_indexer, mol2_all_indexer)
        rmsd_backbone = compute_rmsd(mol1_backbone_indexer, mol2_backbone_indexer)
        rmsd_CA       = compute_rmsd(mol1_CA_indexer, mol2_CA_indexer)

        if (rmsd_all != -1):      print("RMSD (all atoms): %e Angstrom" % rmsd_all)
        if (rmsd_backbone != -1): print("RMSD (backbone) : %e Angstrom" % rmsd_backbone)
        if (rmsd_CA != -1):       print("RMSD (C_alpha)  : %e Angstrom" % rmsd_CA)
        
        return [rmsd_all, rmsd_backbone, rmsd_CA], ['all', 'backbone', 'Calpha']                 # return list of different RMSD
