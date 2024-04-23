from network_generation import *
from network_statistics import *
import os




u = mda.Universe('R1-15-closed/R1-0-start-membrane-3JYC.pdb','R1-15-closed/R1-0-1000-3JYC.xtc')

#WRONG 
with open('15bonds.txt','w') as f:
    for i,atom in enumerate(u.atoms):
        if i != len(u.atoms)-1:
            nextatom = u.atoms[i+1]
        if atom.resname == 'POPC':
            if atom.name != 'C4B':
                print(f'{atom.resid} {atom.name} {nextatom.resid} {nextatom.name}',file=f)
        elif atom.resname == 'CHOL':
            if atom.name != 'RC2':
                print(f'{atom.resid} {atom.name} {nextatom.resid} {nextatom.name}',file=f)
    
    back = u.select_atoms('name BB')
    for i,bb in enumerate(back):
        if i != 1348-1:
            nextatom = back[i+1]
        if (i+1)%337 != 0:
            print(f'{bb.resid} {bb.name} {nextatom.resid} {nextatom.name}',file=f)
    
    for r in u.residues:
        if r.resname != 'CHOL' and r.resname != 'POPC':
            if len(r.atoms) < 4:
                for i,a in enumerate(r.atoms):
                    if i < len(r.atoms)-1:
                        nextatom = r.atoms[i+1]
                        print(f'{a.resid} {a.name} {nextatom.resid} {nextatom.name}',file=f)
            elif r.resname == 'PHE' or r.resname == 'TYR' or r.resname == 'GLY':
                bb = r.atoms[0]
                SC1 = r.atoms[1]
                SC2 = r.atoms[2]
                SC3 = r.atoms[3]
                print(f'{bb.resid} {bb.name} {SC1.resid} {SC1.name}',file=f)
                print(f'{SC1.resid} {SC1.name} {SC2.resid} {SC2.name}',file=f)
                print(f'{SC2.resid} {SC2.name} {SC3.resid} {SC3.name}',file=f)
                print(f'{SC3.resid} {SC3.name} {SC1.resid} {SC1.name}',file=f)
            elif r.resname =='TRP':
                bb = r.atoms[0]
                SC1 = r.atoms[1]
                SC2 = r.atoms[2]
                SC3 = r.atoms[3]
                SC4 = r.atoms[4]
                print(f'{bb.resid} {bb.name} {SC1.resid} {SC1.name}',file=f)
                print(f'{SC1.resid} {SC1.name} {SC2.resid} {SC2.name}',file=f)
                print(f'{SC2.resid} {SC2.name} {SC3.resid} {SC3.name}',file=f)
                print(f'{SC3.resid} {SC3.name} {SC1.resid} {SC1.name}',file=f)
                print(f'{SC3.resid} {SC3.name} {SC4.resid} {SC4.name}',file=f)
                print(f'{SC2.resid} {SC2.name} {SC4.resid} {SC4.name}',file=f)




        

seen = []
for i,atom in enumerate(u.atoms):
    if atom.name in seen: break
    if atom.resname == 'CHOL':
        print(atom.name,atom.resid)
        seen.append(atom.name)
    