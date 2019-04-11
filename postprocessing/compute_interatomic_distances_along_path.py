'''
Compute minimum, maximum, and average bond length for each conformation in each path
and plot the result
'''

import numpy as np
import matplotlib.pyplot as plt

def compute_min_interatomic_distance(atom_indexer, neighbor_search, rcutoff = Quantity.angstrom(2)):
	'''
	Compute the minimal interatomic distances between atoms in atom_indexer
	'''
	neighbor_search.updateNeighborLists()
	
	min_distance = rcutoff														# initialize min distance to rcutoff
	for j in range(atom_indexer.size):											# go though all atoms
		a = atom_indexer[j]
		neighbor_indexer = neighbor_search.getNeighborIndexer(j)
		for n in neighbor_indexer:
			distance = (a.getPosition() - n.getPosition()).norm()
			if distance < min_distance : min_distance = distance
	
	return min_distance

def compute_min_interatomic_distance_per_path(path, atom_indexer, rcutoff = Quantity.angstrom(2.0)):
	'''
	Compute interatomic distances between atoms in atom_indexer along the path
	'''
	# Create a particle system from the atom indexer	
	particle_system = sam.Modeling.DynamicalModel.ParticleSystem(atom_indexer)
	# Create the neighbor search for this particle system
	neighbor_search = sam.Simulation.NeighborSearch.ParticleSystemGrid(particle_system, rcutoff)
	
	interatomic_length = []
	
	for i in range(path.numberOfSteps):											# loop over path steps
		path.currentStep = i													# set the current step of the path
		min_distance = compute_min_interatomic_distance(atom_indexer, neighbor_search, rcutoff)
		interatomic_length.append(min_distance.angstrom)
		# If you want to see in the viewport how the path changes, force SAMSON to process events
		# Comment it for bigger systems, since it may slow down the computations
		SAMSON.processEvents()
	
	return interatomic_length


path_indexer = SAMSON.getNodes('n.t path')										# get an indexer of all paths in the active document
atom_indexer = SAMSON.getNodes('n.t atom')										# get an indexer of all atoms in the active document
rcutoff = Quantity.angstrom(1.2)												# set the rcutoff to 1.2 Angstrom

# Compute min interatomic distance for each path in the document
interatomic_length_per_path = []
for path in path_indexer:														# loop over paths in the path indexer
	interatomic_length_per_path.append( compute_min_interatomic_distance_per_path(path, atom_indexer, rcutoff) )

# Plot results using matplotlib
legend = []
for i in range(0, path_indexer.size):
	plt.plot(interatomic_length_per_path[i])
	legend.append(path_indexer[i].name)

plt.legend(legend)
plt.ylabel('min interatomic distance, Å')
plt.xlabel('conformation')
plt.show()
