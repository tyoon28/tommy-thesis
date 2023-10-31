'''
Makes a chimera script to color and show surfaces of all lipids/sterols.
Requires the .pdb file, but will be changed later to not need it
'''
import pandas as pd
import numpy as np


pdb_file = 'test.pdb'

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



