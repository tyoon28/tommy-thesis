# testing teneto package. Verdict: not very useful

from teneto import TemporalNetwork
import teneto
import numpy as np
import os
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PSF, DCD, GRO, XTC
from MDAnalysis.analysis import contacts, distances
from MDAnalysis.analysis.base import AnalysisBase
import networkx as nx



xtcs = []
for file in os.listdir('R1-30-closed'):
    if file.endswith('.xtc'):
        xtcs.append('R1-30-closed/'+file)
xtcs.sort(key=lambda x: int(x.split('-')[1]))
u = mda.Universe('R1-30-closed/R1-0-start-membrane-3JYC.pdb',*xtcs)


res = u.select_atoms('not resname CHOL and not resname POPC')
r = res.atoms.center_of_mass(compound='residues')
a = np.zeros((100,len(r),len(r)))

for i, ts in enumerate(u.trajectory[:100]):
    frame = u.trajectory.frame
    r = res.atoms.center_of_mass(compound='residues')
    mat = distances.contact_matrix(r, cutoff=6)
    np.fill_diagonal(mat, 0)
    a[i] = mat
a = np.transpose(a,(1,2,0))

tnet = TemporalNetwork(from_array=a)
tnet.calc_networkmeasure('volatility',distance_func='hamming')