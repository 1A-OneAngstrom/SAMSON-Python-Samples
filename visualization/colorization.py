'''
This example shows how to colorize atoms
based on their position (x, y, or z)

see: https://documentation.samson-connect.net/scripting-guide/colorization/

Open any molecule in SAMSON
'''

class HaltException(Exception): pass

indexer = SAMSON.getNodes('n.t atom')                   # get all atom nodes

if indexer.size == 0: HaltException('No atoms found in the active document')
 
direction = ['x','y','z']
ipos = 0;                                               # parameter that determines the coordinate according to which the colorization should be done
print('Colorize set of %d nodes in %s-direction' % (indexer.size, direction[ipos]))
 
SAMSON.setBusy(True)                                    # notify user that SAMSON is busy
 
minx = maxx = indexer[0].getPosition()[ipos].value      # set an initial value for minimal and maximal coordinates

for atom in indexer:                                    # compute minimal and maximal coordinates
    val = atom.getPosition()[ipos].value
    if val > maxx: maxx = val
    if val < minx: minx = val
 
for atom in indexer:                                    # colorize each atom according to its coordinate
                                                        # get RGB color for a node based on HSV color
    color = sam.DataModel.Color.Color.fromHSV(239.5 / 360.0 * (atom.getPosition()[ipos].value - minx) / (maxx - minx), 205.0/255.0, 1.0);
    atom.setColor(color);                               # set the color of the node
 
SAMSON.setBusy(False)                                   # SAMSON is not busy anymore
