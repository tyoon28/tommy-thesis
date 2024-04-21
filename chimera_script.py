'''
Makes a chimera script to color and show surfaces of all lipids/sterols.
Requires the .pdb file, but will be changed later to not need it
'''
import pandas as pd
import numpy as np
from matplotlib.pyplot import cm
from network_generation import *
import csv
from matplotlib.colors import LinearSegmentedColormap


one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
'GLY':'G', 'PRO':'P', 'CYS':'C'}

# pdb_file = '/Users/Tommy/Desktop/thesis/tommy-thesis/R1-30-closed/R1-0-start-membrane-3JYC.pdb'

# n = ['type','atomid','name','res','resid','x','y','z','b','t']
# x = pd.read_csv(pdb_file,delim_whitespace=True,header=4,names = n)


# fn = pdb_file[:-4] + '_surface.cxc'


# with open(fn,'w') as f:
#     print('select ~:CHOL & ~:POPC','color sel plum','surface sel',sep='\n',file=f)
#     print('select :CHOL','color sel orange',sep='\n',file=f)
#     print('select :POPC','color sel white',sep='\n',file=f)
        
#     for i in x.loc[x['res'] == 'CHOL']['resid'].unique():
#         print(f'surface :{int(i)}',file=f)


#     for i in x.loc[x['res'] == 'POPC']['resid'].unique():
#         print(f'surface :{int(i)}',file=f)


# ################
        
        
# u.select_atoms('resid 1646').atoms.center_of_mass(compound='residues')
# u.select_atoms('resid 2018').atoms.center_of_mass(compound='residues')
# u.select_atoms('resid 2095').atoms.center_of_mass(compound='residues')

# with open('selecthide.cxc','w') as f:
#     for r in u.residues:
#         po = r.atoms.center_of_mass()[0]
#         if po < 115.7066 and (r.resname == 'CHOL' or r.resname =='POPC'):
#             print(f'hide :{r.resid} surfaces',file=f)
#             print(f'hide :{r.resid}',file=f)




def color_by_centrality(d,fn,normalize=True):
    d_sort = dict(sorted(d.items(), key=lambda x: x[1]))
    a = np.array(list(d_sort.values()),dtype=np.longdouble)

    if normalize:
        a = a/a.max()
    cmap = LinearSegmentedColormap.from_list("", ["darkgray","red"])
    color = iter(cmap(a))

    with open(f'{fn}.cxc','w') as f:
        for i in d_sort:
            c = next(color)
            c = c * 255
            s = MD_to_chimerax(int(i))
            print(f'color {s} rgba({round(c[0])},{round(c[1])},{round(c[2])},1)',file=f)
            
    return




def dd_to_csv(d,fn,u):
    d_sort = dict(sorted(d.items(), key=lambda x: x[1]))

    a = np.array(list(d_sort.values()),dtype=np.longdouble)
    a = a/a.max()
    ind = 0
    with open(f'{fn}.csv', 'w', newline='') as csvfile:
        
        wr = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        sh = wr.writerow(['subunit','residue ID','MD simulation ID','residue name','node PCA distance','PCA distance normalized to max'])
        for i in d_sort:
            resname = u.atoms.select_atoms(f'resid {int(i)}').residues.resnames[0]
            s = MD_to_chimerax(i)
            subunit = s[3]
            residue = s[5:]
            value = d_sort[i]
            value_normalized = a[ind]
            
            mdsimid = i
            ind+=1
            
            sh = wr.writerow([subunit,residue,mdsimid,resname, value,value_normalized])
    return


def MD_to_chimerax(s):
    if s%337 == 0:
        e = 337
    else: e = s%337
    c = (s-1)//337 + 1
    if c == 1: chain = 'A'
    elif c == 2: chain = 'B'
    elif c == 3: chain = 'C'
    else: chain = 'D'
    
    return f'#1/{chain}:{e + 35}'



def MD_to_resid(s,u=None):
    one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
    'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
    'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
    'GLY':'G', 'PRO':'P', 'CYS':'C'}

    if s%337 == 0:
        e = 337
    else: e = s%337
    c = (s-1)//337 + 1
    if u is not None:
        res = one_letter[u.select_atoms(f'resid {s}').residues.resnames[0]]
    else: res = ''

    
    
    return res+str(e + 36)



def u_to_top(u):
    pass

def symmetrize(fn):
    fn = '/Users/Tommy/Desktop/thesis/figures/out/all-dynnodepca.csv'
    do = pd.read_csv(fn)
    pa = do.groupby('residue ID')['PCA distance normalized to max'].mean()
    dd = pa.to_dict()

    ddd = {}
    for k in dd:
        for key in [(i * 337) + (k - 35) for i in range(0,4)]:
            ddd[key] = dd[k]


    print(dd[265])
    print(dd[265])
    print(dd[178])
    print(dd[179])
    print(dd[182])
    
    color_by_centrality(ddd,'all-dynnodepca-symmetric',normalize=False)