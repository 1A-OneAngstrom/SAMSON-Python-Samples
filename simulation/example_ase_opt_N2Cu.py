'''
In this example, we will show how to use ASE Python package together with SAMSON.
We will create in ASE a system consisten of Nitrogen molecule and a slab of Copper atoms using ASE package,
copy it to SAMSON, compute adsorption energy, simulate it in ASE using velocity Verlet,
and update atoms positions in SAMSON accordingly.

See: https://documentation.samson-connect.net/scripting-guide/usage-of-other-packages-ase/
ASE tutorial: https://wiki.fysik.dtu.dk/ase/tutorials/surface.html
'''

from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import fcc111, add_adsorbate
from ase.md.verlet import VelocityVerlet
from ase import units

class ExampleASEOptimizationN2Cu:
        
        slab = None
        molecule = None
        indexer = None

        def __init__(self):
                '''
                Initialize ASE and SAMSON
                '''
                self.initialize_system_ase()                                    # initialize the system in ASE
                self.initialize_system_samson()                                 # initialize the system in SAMSON
            
        def compute_adsorption_energy(self):
                '''
                Computations using ASE; update positions in SAMSON
                Optimize the structure of the N2 molecule adsorbed on the Cu surface
                '''
                # set up computations in ASE
                h = 1.85
                add_adsorbate(self.slab, self.molecule, h, 'ontop')
                constraint = FixAtoms(mask = [a.symbol != 'N' for a in self.slab])
                self.slab.set_constraint(constraint)
                dyn = QuasiNewton(self.slab) #, trajectory = '/tmp/N2Cu.traj'
                dyn.run(fmax=0.05)

                self.update_positions_samson()                                  # update positions in SAMSON
                        
                print('Adsorption energy:', self.e_slab + self.e_N2 - self.slab.get_potential_energy())

        def simulate(self):
                '''
                MD computations using ASE; update positions in SAMSON
                '''
                # set up computations in ASE
                dyn = VelocityVerlet(self.molecule, dt = 1.0 * units.fs)
                for i in range(20):
                        pot = self.molecule.get_potential_energy()
                        kin = self.molecule.get_kinetic_energy()
                        print('%2d: %.5f eV, %.5f eV, %.5f eV' % (i, pot + kin, pot, kin))
                        dyn.run(steps=20)                                       # run the computations in ASE
                        self.update_positions_samson()                          # update positions in SAMSON
        
        def initialize_system_ase(self):
                '''
                Create an H2O molecule using ASE package; set a calculator
                '''
                d = 1.10                                                        # set parameters
                # create atoms using ASE (positions in ASE are in Angstrom); set calculator
                self.slab = fcc111('Cu', size = (4, 4, 2), vacuum = 10.0)       # build Cu crystal
                self.slab.set_calculator(EMT())
                self.e_slab = self.slab.get_potential_energy()
                
                self.molecule = Atoms('2N', positions = [(0., 0., 0.), (0., 0., d)])
                self.molecule.set_calculator(EMT())
                self.e_N2 = self.molecule.get_potential_energy()
                
        def initialize_system_samson(self):
                '''
                Adds atoms created by ASE into a new layer in SAMSON
                '''
                document = SAMSON.getActiveDocument()                           # get the active document
                layer = sam.DataModel.Document.Layer('ase slab opt')            # a new layer
                layer.create()                                                  # create it in SAMSON
                document.addChild(layer)                                        # add this layer to the active document
                structural_model = sam.Modeling.StructuralModel.StructuralModel() # a new structural model
                structural_model.create()                                       # create this structural model in SAMSON
                layer.addChild(structural_model)                                # add the structural model into the layer
                structural_root = structural_model.getStructuralRoot()          # get the root of the structural model

                self.add_atoms_in_samson(structural_root, self.slab)            # add slab atoms in SAMSON
                self.add_atoms_in_samson(structural_root, self.molecule)        # add molecule atoms in SAMSON
                
                self.indexer = structural_root.getNodes('n.t a and N')          # get a list of all Nitrogen atoms in the structural_root
                
                SAMSON.getActiveCamera().leftView()                             # set a camera view in SAMSON
        
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
# After this step, you will see Nitrogen and Copper molecules in the SAMSON’s viewport in a newly created layer
ase_N2Cu_opt = ExampleASEOptimizationN2Cu()
SAMSON.processEvents()                          # process events: updates viewport
ase_N2Cu_opt.compute_adsorption_energy()        # compute the adsorption energy
ase_N2Cu_opt.simulate()                         # perform MD simulation
