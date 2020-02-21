'''
Usage of the SAMSON’s Interaction Model
In this example, we will be using an Interaction Model provided by SAMSON together with ASE for minimization of the system.
This Interaction Model will be used to compute forces in SAMSON and ASE will be used for minimization of the system by using FIRE.
For that, we implement our own calculator in ASE (see ASE calculators).

See: https://documentation.samson-connect.net/scripting-guide/usage-of-other-packages-ase/
Calculator for ASE: https://wiki.fysik.dtu.dk/ase/development/calculators.html
'''

import os
import numpy as np
from ase import Atoms, Atom
from ase.optimize import FIRE
from ase.calculators.calculator import Calculator, all_changes

class IMWrap(Calculator):
        '''
        Calculator for ASE, see: https://wiki.fysik.dtu.dk/ase/development/calculators.html
        '''
        
        implemented_properties = ['energy', 'forces']                   # properties that calculator can handle 
        
        im = None
        molecule = None
        indexer = None
        
        def __init__(self, **kwargs):
                Calculator.__init__(self, **kwargs)
                if kwargs is not None:
                        print("set IMPS")
                        self.im = kwargs["im"]
                        self.molecule = kwargs["molecule"]
                        self.indexer = kwargs["indexer"]
                else:
                        print("Error: input parameters are not provided")
                        pass
                if (self.im == None or self.molecule == None or self.indexer == None):
                        print("Error: improper input parameters")
                        pass
        
        def calculate(self, atoms = None, properties = ['energy', 'forces'], system_changes = all_changes):
                '''
                The main function for Calculator that does the calculations
                '''
                
                self.im.updateInteractions()                            # update interactions in SAMSON
                
                forces = np.zeros((len(self.molecule), 3))
                for i in range(self.indexer.size):                      # get forces from the SAMSON's InteractionModelParticleSystem
                        forces[i, 0] = self.im.getForce(i)[0].N
                        forces[i, 1] = self.im.getForce(i)[1].N
                        forces[i, 2] = self.im.getForce(i)[2].N
                
                force_convertion = 1. / 8.238722514e-8                  # 1 atomic unit force = Eh / a0 = 8.2387225(14) x 10^(-8) N
                # set calculated properties for ASE
                self.results['forces'] = forces * force_convertion      # convert forces from N to atomic unit force
                self.results['energy'] = self.im.energy.eV.value        # get energy from the SAMSON's InteractionModelParticleSystem
                
                self.update_positions_samson()                          # update positions of atoms in SAMSON
        
        def update_positions_samson(self):
                '''
                Update positions of atoms in SAMSON based on new positions of these atoms in ASE
                '''
                new_positions = self.molecule.get_positions()           # get the new positions from ASE
                
                for i in range(self.indexer.size):
                        self.indexer[i].setX(sam.DataModel.Quantity.angstrom(new_positions[i, 0]))
                        self.indexer[i].setY(sam.DataModel.Quantity.angstrom(new_positions[i, 1]))
                        self.indexer[i].setZ(sam.DataModel.Quantity.angstrom(new_positions[i, 2]))
                
                SAMSON.processEvents()                                  # process all events in SAMSON (to update positions in the viewport)

class ExampleASECalculator:
        '''
        An example of the usage of a SAMSON's InteractionModelParticleSystem together with ASE
        We assume, that the InteractionModelParticleSystem has been already created in SAMSON for the atoms in the SAMSON's data graph
        '''
        
        molecule = None
        
        def __init__(self):
        
                indexer = SAMSON.getNodes('n.t a')                      # get list of all atoms in SAMSON
                
                self.create_atoms_in_ase(indexer)                       # create atoms in ASE
                        
                imps_indexer = SAMSON.getNodes('n.t imps')              # get all interaction models from the active document in SAMSON
                imps = imps_indexer[0].castToInteractionModelParticleSystem()  # get the first interaction model
                
                self.molecule.set_calculator(IMWrap(im = imps, molecule = self.molecule, indexer = indexer))  # set an ASE calculator
                
        def create_atoms_in_ase(self, indexer):
                '''create atoms in ASE'''
                self.molecule = Atoms()                                 # create an empty molecule in ASE
                
                for a in indexer:                                       # add atoms in the molecule
                        self.molecule.append( Atom( a.elementSymbol, (a.getX().angstrom, a.getY().angstrom, a.getZ().angstrom) ) )
                
        def compute(self):
                dyn = FIRE(self.molecule)                               # set the computation
                dyn.run(fmax = 0.005, steps = 100)                      # run the computation


def samson_add_simulator(simulator):
        '''
        Initialize and add simulator in the active document in SAMSON 
        To add a simulator in SAMSON we need to create and add a dynamical model, an interaction model, and a state updater.
        '''
        dynamicalModel = simulator.getDynamicalModel()                  # get the dynamical model from the simulator
        dynamicalModel.create()                                         # create the dynamical model

        interactionModel = simulator.getInteractionModel()              # get the interaction model from the simulator
        interactionModel.create()                                       # create the interaction model
        SAMSON.showProperties(interactionModel)                         # show a window with interaction model properties
        interactionModel.initializeInteractions()                       # initialize interaction model

        stateUpdater = simulator.getStateUpdater()                      # get the state updater from the simulator
        stateUpdater.create()                                           # create the state updater
        SAMSON.showProperties(stateUpdater)                             # show a window with state updator properties

        simulator.create()                                              # create the simulator

        document = SAMSON.getActiveDocument()                           # get the active document in SAMSON
        SAMSON.beginHolding('Add simulator')                            # holding: take care of undo/redo in SAMSON
        SAMSON.hold(dynamicalModel)                                     # hold an object for undo/redo in SAMSON
        document.addChild(dynamicalModel)                               # add the dynamical model to the active document in SAMSON
        SAMSON.hold(interactionModel)
        document.addChild(interactionModel)                             # add the interaction model to the active document in SAMSON
        SAMSON.hold(simulator)
        document.addChild(simulator)                                    # add the simulator to the active document in SAMSON
        SAMSON.endHolding()


# Open a hydrocarbon molecule
path_to_pyscript = os.path.realpath(__file__)
SAMSON.importFromFile(os.path.dirname(path_to_pyscript) + '/../samples/Cyclopropane.pdb', [])
SAMSON.processEvents()

# Apply a simulator to this molecule
# with Brenner interaction model and Interactive modeling state updater
indexer = SAMSON.getNodes('n.t sm')                                             # get an indexer of all structural models in the active document in SAMSON
IMuuid = sam.Core.Container.UUID('AD608CB6-6971-7CD4-6FCC-34531998E743')        # UUID of the Brenner SAMSON Element
SUuuid = sam.Core.Container.UUID('F912F119-7CBB-B5BD-972A-0A02DFCF683D')        # UUID of the State Updater pack SAMSON Element
simulator = SAMSON.makeSimulator(indexer, 'SMMBrennerInteractionModel', IMuuid, 'SESInteractiveModelingUpdater', SUuuid)        # make an instance of a simulator

samson_add_simulator(simulator)         # add the simulator in the active document in SAMSON

# Run ASE calculations
ase_calc = ExampleASECalculator()       # Create an object of the ExampleASECalculator class, which will initialize system in ASE and in SAMSON
ase_calc.compute()                      # Run computations: force computation using SAMSON and minimization using ASE
