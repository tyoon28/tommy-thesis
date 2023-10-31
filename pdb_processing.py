'''
add chain IDs to .pdb files. each molecule has its own chain. does not work and is useless.
'''

import pandas as pd
import numpy as np

n = ['type','atomid','name','res','resid','x','y','z','b','t']
x = pd.read_csv('test.pdb',delim_whitespace=True,header=4,names = n)


x.iloc[int(3064/4) -1]
x = x.dropna()

l = 766 * ['A'] + 766 * ['B'] + 766 * ['C'] + 766 * ['D']

x['chain'] = x['res'] + x["resid"].astype(str)

x.loc[:3063,'chain'] = l

x['type'] = np.where(x['res'] == 'CHOL', 'HETATM', x['type'])
x['type'] = np.where(x['res'] == 'POPC', 'HETATM', x['type'])

x['atomid'] = x["atomid"].astype(int)
x['resid'] = x["resid"].astype(int)


cols = list(x.columns.values)
cols.insert(4, cols.pop(-1))

x = x[cols]

x.to_csv('testout.pdb',sep='\t',header=False,index=False)

