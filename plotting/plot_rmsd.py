import pandas as pd
from compute_rmsd import compute_rmsd_2mol

# In this example, we assume that the active document in SAMSON
# contains a number of structural models representing a trajectory of a bio molecule

indexer = SAMSON.getNodes('n.t m')                              # get all molecules in the active document

df = pd.DataFrame()                                             # create pandas data frame

if indexer.size < 2:
        print("less than two molecules are present in the data graph")
else:                                                           # we assume here that we need to compute RMSD between the first two molecules
        goal = indexer[-1]                                      # get the goal molecule [assume the last molecule to be the goal one]
        indexer.removeNode(goal)                                # remove the target node from the indexer
        for t in range(indexer.size):
                mol = indexer[t]
                rmsd, rmsd_info = compute_rmsd_2mol(goal, mol)  # compute RMSD between two conformations
                # fill in the data frame
                for j in range(len(rmsd)):
                        df = df.append( pd.DataFrame({'name': [mol.name], 'info': rmsd_info[j], 'RMSD': rmsd[j], 'timestep': t }) )

df.head()                                                       # show the first several rows

#######################
# plot using matplotlib
import matplotlib.pyplot as plt

if df.size > 0 :
        df_all = df.groupby('info').get_group('all')            # group by info of RMSD and get RMSD computed based on all atoms
        df_bb  = df.groupby('info').get_group('backbone')       # group by info of RMSD and get RMSD computed based on backbone atoms
        df_CA  = df.groupby('info').get_group('Calpha')         # group by info of RMSD and get RMSD computed based on Calpha atoms
        line_RMSD_all, = plt.plot(df_all["timestep"], df_all["RMSD"], label="all")
        line_RMSD_bb,  = plt.plot(df_bb ["timestep"], df_bb ["RMSD"], label="backbone")
        line_RMSD_CA,  = plt.plot(df_CA ["timestep"], df_CA ["RMSD"], label="Calpha")
        plt.ylabel('RMSD, Angstrom')
        plt.xlabel('timestep')
        plt.legend(handles = [line_RMSD_all, line_RMSD_bb, line_RMSD_CA])
        plt.grid()
        plt.show()

#######################
# plot using bokeh 
# the plot will be openned in a default browser
from bokeh.plotting import figure, show

if df.size > 0 :
        p = figure(title="RMSD")
        p.line(x=df_all["timestep"], y=df_all["RMSD"], color='red',   legend='all')
        p.line(x=df_bb ["timestep"], y=df_bb ["RMSD"], color='blue',  legend='backbone')
        p.line(x=df_CA ["timestep"], y=df_CA ["RMSD"], color='green', legend='Calpha')
        show(p)
