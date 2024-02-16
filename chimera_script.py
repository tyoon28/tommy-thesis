'''
Makes a chimera script to color and show surfaces of all lipids/sterols.
Requires the .pdb file, but will be changed later to not need it
'''
import pandas as pd
import numpy as np
from matplotlib.pyplot import cm
from network_generation import *


pdb_file = '/Users/Tommy/Desktop/thesis/tommy-thesis/R1-30-closed/R1-0-start-membrane-3JYC.pdb'

n = ['type','atomid','name','res','resid','x','y','z','b','t']
x = pd.read_csv(pdb_file,delim_whitespace=True,header=4,names = n)


fn = pdb_file[:-4] + '_surface.cxc'


with open(fn,'w') as f:
    print('select ~:CHOL & ~:POPC','color sel plum','surface sel',sep='\n',file=f)
    print('select :CHOL','color sel orange',sep='\n',file=f)
    print('select :POPC','color sel white',sep='\n',file=f)
        
    for i in x.loc[x['res'] == 'CHOL']['resid'].unique():
        print(f'surface :{int(i)}',file=f)


    for i in x.loc[x['res'] == 'POPC']['resid'].unique():
        print(f'surface :{int(i)}',file=f)


################
        
        
u.select_atoms('resid 1646').atoms.center_of_mass(compound='residues')
u.select_atoms('resid 2018').atoms.center_of_mass(compound='residues')
u.select_atoms('resid 2095').atoms.center_of_mass(compound='residues')

with open('selecthide.cxc','w') as f:
    for r in u.residues:
        po = r.atoms.center_of_mass()[0]
        if po < 115.7066 and (r.resname == 'CHOL' or r.resname =='POPC'):
            print(f'hide :{r.resid} surfaces',file=f)
            print(f'hide :{r.resid}',file=f)




def color_by_centrality(d,fn):
    d_sort = dict(sorted(d.items(), key=lambda x: x[1]))

    a = np.array(list(d_sort.values()),dtype=np.longdouble)
    a = a/a.max()
    color = iter(cm.bwr(a))

    with open(f'{fn}.cxc','w') as f:
        for i in d_sort:
            c = next(color)
            c = c * 255
            s = MD_to_chimerax(i)
            print(f'color {s} rgba({round(c[0])},{round(c[1])},{round(c[2])},1)',file=f)
            
    return