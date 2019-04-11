'''
Compute minimum interatomic distance for each conformation in each path
and plot the result
'''

import numpy as np
import matplotlib.pyplot as plt

def compute_bond_length_stats(bond_indexer):
	'''
	Compute min, max, and avg length of bonds
	'''
	min_bl = bond_indexer[0].length.angstrom;									# initialize min length
	max_bl = bond_indexer[0].length.angstrom;									# initialize max length
	avg_bl = 0.0;
	for b in bond_indexer:														# loop over bonds in the bond indexer
		bl = b.length.angstrom													# get the bond length in Angstroms	
		avg_bl = avg_bl + bl
		if bl < min_bl : min_bl = bl
		if bl > max_bl : max_bl = bl
	
	avg_bl = avg_bl / bond_indexer.size
	
	return [max_bl, avg_bl, min_bl]

def compute_bond_length_stats_per_path(path, bond_indexer):
	'''
	Compute min, max, and avg length of bonds along the path
	'''
	bond_length_along_path = []
	
	for i in range(path.numberOfSteps):											# loop over path steps
		path.currentStep = i													# set the current step of the path
		bond_length_along_path.append(compute_bond_length_stats(bond_indexer))	# compute bond length stats for the current step
		# If you want to see in the viewport how the path changes, force SAMSON to process events
		# Comment it for bigger systems, since it may slow down the computations
		SAMSON.processEvents()
	
	return bond_length_along_path


path_indexer = SAMSON.getNodes('n.t path')										# get an indexer of all paths in the active document
bond_indexer = SAMSON.getNodes('n.t bond')										# get an indexer of all bands in the active document

# Compute bond length stats for each path in the document
bond_length_along_paths = []
for path in path_indexer:														# loop over paths in the path indexer
	bond_length_along_paths.append( compute_bond_length_stats_per_path(path, bond_indexer) )

# Plot results using matplotlib
for i in range(0, path_indexer.size):
	plt.figure(i)
	plt.plot(bond_length_along_paths[i])
	plt.title(path_indexer[i].name)
	plt.legend(['max','mean','min'])
	plt.ylabel('bond length, Å')
	plt.xlabel('conformation')
	plt.show()
