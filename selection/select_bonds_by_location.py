'''
This example shows how to select bonds which center lies in between of two planes: z_min and z_max

If you want to select bonds which center lies on a specific plane, just set z_min=z_max
'''

# Get an indexer of all bonds in the active document. Here we use SAMSON Node Specification Language to get only bonds
allBondsIndexer = SAMSON.getNodes('node.type bond')

# two planes: 1.7A and 2.3A
z_min = Quantity.angstrom(1.7)
z_max = Quantity.angstrom(2.3)

selection = sam.DataModel.DataGraph.NodeIndexer()				# an indexer in which we will be adding the desired bonds

for bond in allBondsIndexer:							# a loop over an indexer
	bondCenterZ = (bond.leftAtom.getZ() + bond.rightAtom.getZ()) / 2.0	# compute the z-axis center of the bond 
	if bondCenterZ >= z_min and bondCenterZ <= z_max:			# check whether the bond lies in between of desired planes
		'''
		Add the bond to the selection indexer.
		We can set the selectionFlag to True for the bond here.
		But for the sake of the example we will add the selected bonds in an indexer which can be used later.
		'''
		selection.addNode(bond)						

allBondsIndexer.clear()								# clear an indexer

for bond in selection:								# loop over the desired bonds
    bond.selectionFlag = True							# select these bonds by setting the selectionFlag to True for them

selection.clear()								# clear an indexer
