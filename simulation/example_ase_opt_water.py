'''
In this example, we will show how to use ASE Python package together with SAMSON.
We will create a water molecule using ASE package, copy it to SAMSON, 
optimize it using ase (BFGS), and update atoms positions in SAMSON accordingly.

See: https://documentation.samson-connect.net/scripting-guide/usage-of-other-packages-ase/
'''

import numpy as np
from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.emt import EMT

class ExampleASEOptimizationH20:
        
        molecule = None
        indexer = None

        def __init__(self):
                '''
                Initialize ASE and SAMSON
                '''
                self.initialize_system_ase()                                    # initialize the system in ASE
                self.initialize_system_samson()                                 # initialize the system in SAMSON

        def compute(self):
                '''
                Do computations using ASE and update positions in SAMSON
                '''
                dyn = BFGS(self.molecule)                                       # set up computations in ASE
                nsteps = 20
                for i in range(nsteps):
                        dyn.run(fmax = 0.005, steps = 1)                        # run the computations in ASE
                        self.update_positions_samson()                          # update positions in SAMSON
        
        def initialize_system_ase(self):
                '''
                Create an H2O molecule using ASE package; set a calculator
                '''
                d = 0.9575                                                      # set parameters
                t = np.pi / 180.0 * 104.51
                # create a molecule in ASE (positions in ASE are in Angstrom); set calculator
                self.molecule = Atoms('H2O',
                                positions = [(d, 0, 0.1),
                                             (d * np.cos(t), d * np.sin(t), 0),
                                             (0, 0, 0)],
                                calculator = EMT())
        
        def initialize_system_samson(self):
                '''
                Adds atoms created by ASE into a new layer in SAMSON
                '''
                document = SAMSON.getActiveDocument()                           # get the active document
                layer = sam.DataModel.Document.Layer('ase H2O opt')             # a new layer
                layer.create()                                                  # create it in SAMSON
                document.addChild(layer)                                        # add this layer to the active document
                structural_model = sam.Modeling.StructuralModel.StructuralModel() # a new structural model
                structural_model.create()                                       # create this structural model in SAMSON
                layer.addChild(structural_model)                                # add the structural model into the layer
                structural_root = structural_model.getStructuralRoot()          # get the root of the structural model

                self.add_atoms_in_samson(structural_root, self.molecule)        # add atoms in SAMSON
                
                self.indexer = structural_root.getNodes('n.t a')                # get a list of all atoms in the structural_root
                
                SAMSON.getActiveCamera().topView()                              # set a camera view in SAMSON
        
        def add_atoms_in_samson(self, structural_root, molecule):
                '''
                Add atoms in SAMSON according to their type and positions
                '''
                for a in molecule:                                             # loop through a list of atoms created by ase
                        sa = sam.Modeling.StructuralModel.Atom(SAMSON.getElementTypeBySymbol(a.symbol),
                                sam.DataModel.Quantity.angstrom(a.position[0]),
                                sam.DataModel.Quantity.angstrom(a.position[1]),
                                sam.DataModel.Quantity.angstrom(a.position[2])) # construct an atom in SAMSON with a given element type and position
                        sa.create()                                             # create an atom in SAMSON
                        structural_root.addChild( sa )                          # add an atom to the root of the structural model

        def update_positions_samson(self):
                '''
                Update positions of atoms in SAMSON based on new positions of these atoms in ASE
                '''
                new_positions = self.molecule.get_positions()                   # get the new positions from ASE

                for i in range(self.indexer.size):
                        self.indexer[i].setX(sam.DataModel.Quantity.angstrom(new_positions[i, 0]))
                        self.indexer[i].setY(sam.DataModel.Quantity.angstrom(new_positions[i, 1]))
                        self.indexer[i].setZ(sam.DataModel.Quantity.angstrom(new_positions[i, 2]))
                
                SAMSON.processEvents()                                          # process all events in SAMSON (to update positions in the viewport)


# Create an object of the ExampleASEOptimizationH20 class, which will initialize ASE and SAMSON.
# After this step, you will see a water molecule in the SAMSON’s viewport in a newly created layer
ase_h20_opt = ExampleASEOptimizationH20()
SAMSON.processEvents()          # process events: updates viewport
ase_h20_opt.compute()           # run the computations
