from network_generation import *
from network_statistics import *
import os




u = mda.Universe('R1-30-closed/R1-0-start-membrane-3JYC.pdb','R1-30-closed/R1-0-1000-3JYC.xtc')

#WRONG 
with open('30bonds.txt','w') as f:
    for i,atom in enumerate(u.atoms):
        if i != len(u.atoms)-1:
            nextatom = u.atoms[i+1]
        if atom.resname != 'POPC' and atom.resname != 'CHOL':
            if atom.resid % 337 != 0:
                print(f'{atom.resid-1} {atom.name} {nextatom.resid-1} {nextatom.name}',file=f)
        elif atom.resname == 'POPC':
            if atom.name != 'C4B':
                print(f'{atom.resid-1} {atom.name} {nextatom.resid-1} {nextatom.name}',file=f)
        elif atom.resname == 'CHOL':
            if atom.name != 'RC2':
                print(f'{atom.resid-1} {atom.name} {nextatom.resid-1} {nextatom.name}',file=f)

        

seen = []
for i,atom in enumerate(u.atoms):
    if atom.name in seen: break
    if atom.resname == 'CHOL':
        print(atom.name,atom.resid)
        seen.append(atom.name)
    