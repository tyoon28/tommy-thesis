'''Unused testing'''
from graphlets import *

def main():
    for c in ['30','15']:
        for r in ['R1','R2','R3']:
            xtcs = []
            for file in os.listdir(f'{r}-{c}-closed'):
                if file.endswith('.xtc'):
                    xtcs.append(f'{r}-{c}-closed/'+file)
            xtcs.sort(key=lambda x: int(x.split('-')[1]))
            u = mda.Universe(f'{r}-{c}-closed/{r}-0-start-membrane-3JYC.pdb',*xtcs)
            protein = u.select_atoms('not resname CHOL and not resname POPC').residues
            for t in u.trajectory:
                
                res = protein.atoms.center_of_mass(compound='residues')
                mat = distances.contact_matrix(res, cutoff=6)
                
                # Graph from adjacency matrix
                G = nx.from_numpy_array(mat)
                deg = G.degree()
                for n in G.nodes:
                    if deg[n] == 3:
                        pass






    
    return