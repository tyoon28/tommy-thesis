'''
Functions/AnalysisBase classes for calculating basic statistics from MD trajectory
'''
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PSF, DCD, GRO, XTC
from MDAnalysis.analysis import contacts, distances
from MDAnalysis.analysis.base import AnalysisBase
import numpy as np
from matplotlib import pyplot as plt
import networkx as nx
import itertools
import warnings
import time
import tqdm

class Cholesterol_contact(AnalysisBase):
    '''
    atomgroup: selection of cholesterol(s) and protein or all atoms
    results: 
        each item of results[id] is:
            contacts - list of frames that contact the channel
            contact_time - total contact time
            longest_contact - duration (frames) of longest contact with channel
            binding - T/F does it contact the channel for longer than 25000 unique residue contacts
    '''

    def __init__(self, atomgroup, **kwargs):
        super(Cholesterol_contact, self).__init__(atomgroup.universe.trajectory,
                                          **kwargs)
        self.ag = atomgroup

    def _prepare(self):
        # OPTIONAL
        # Called before iteration on the trajectory has begun.
        # Data structures can be set up at this time
        # get all cholesterol ids
        self.sterols = self.ag.select_atoms(f'(resname CHOL and not (name RC1 or name RC2))').residues
        self.protein = self.ag.select_atoms('not resname CHOL and not resname POPC').residues
        self.lenprotein = len(self.protein)
        self.lensterols = len(self.sterols)

        for i in self.sterols.resids:
            self.results[i] = {'contacts':[]}


        self.protein_sterols = self.sterols + self.protein
        self.contact_threshold = 6

    def _single_frame(self):
        # REQUIRED
        # Called after the trajectory is moved onto each new frame.
        # store an example_result of `some_function` for a single frame
        res = self.protein_sterols.atoms.center_of_mass(compound='residues')
        contact_mat = distances.contact_matrix(res, cutoff=self.contact_threshold)
        for i, sterol in enumerate(self.sterols):
            resid = sterol.resid
            mat_ind = self.lenprotein + i
            if np.any(contact_mat[mat_ind][self.lenprotein]):
                self.results[resid]['contacts'].append(self.ag.universe.trajectory.frame)

    def _conclude(self):
        # OPTIONAL
        for i in self.results:
            self.results[i]['contact_time'] = len(self.results[i]['contacts'])
            self.results[i]['longest_contact'] = self._longest(self.results[i]['contacts'])

        # Called once iteration on the trajectory is finished.
        # Apply normalisation and averaging to results here.
        pass

    def _longest(self, a):
        longest = 0
        last = -1
        s = 0
        for i in a:
            if i-1 == last:
                s += 1
            else:
                s = 0
            if s > longest: longest = s
            last = i

        return longest



def chol_contact_duration(u,sterol_id):
    '''returns longest continuous duration of cholesterol's contact with the channel'''
    contact_threshold = 6
    protein = u.select_atoms('not resname CHOL and not resname POPC')

    n_res = len(protein.atoms.center_of_mass(compound='residues'))

    longest_contact = 0
    contact = 0

    start = time.time()
    for ts in tqdm.tqdm(u.trajectory[:1000]):
        

        # excluding cholesterol tails from center of mass calculations
        sterol = (u.select_atoms(f'resname CHOL and resid {sterol_id} and not (name RC1 or name RC2)')
                  .atoms.center_of_mass(compound='residues'))


        res = u.select_atoms('not resname CHOL and not resname POPC').atoms.center_of_mass(compound='residues')

        # worried about efficiency using distance array but can't find built in contact function for two selections...
        #distance array between sterol and all other residues
        dist_arr = distances.distance_array(sterol, 
                                    res, 
                                    box=u.dimensions)
        
        # if cholesterol is contacting any residues
        if np.any(dist_arr < contact_threshold):
            contact += 1
        else:
            longest_contact = contact
            contact = 0
        
    print(time.time()-start)
    print(longest_contact)

    return longest_contact

    

def find_binding_sterols(u):
    sterols = set(u.select_atoms(f'resname CHOL').resids)
    touches = []

    for id in tqdm.tqdm(sterols):
        duration = chol_contact_duration(u,id)
        if duration != 0:
            touches.append(id)
        
    print(touches)
    print(len(touches))

def find_binding_sterols2(u):
    sterols = set(u.select_atoms(f'resname CHOL').resids)
    touches = []
    binders = {}
    for id in tqdm.tqdm(sterols):
        binders[id] = binder(u,id)
        
    return binders


def binder(u,sterol_id):
    '''does sterol contact channel for more than 25000 unique residue contacts'''
    contact_threshold = 6
    protein = u.select_atoms('not resname CHOL and not resname POPC')

    n_res = len(protein.atoms.center_of_mass(compound='residues'))

    longest_contact = 0
    contact = 0

    start = time.time()
    unique_residue_contacts = np.zeros((n_res,n_res+1))
    for ts in u.trajectory:
        
        # all residues and sterol
        res = (u.select_atoms(f'(not resname CHOL and not resname POPC) or \
                                  (resname CHOL and resid {sterol_id} and not (name RC1 or name RC2))')
                                  .atoms.center_of_mass(compound='residues'))
        
        contact_mat = distances.contact_matrix(res, cutoff=contact_threshold)

        # last row of contact matrix should be the sterol (?)
        if np.any(contact_mat[-1]):
            contact += distances.contact_matrix(res, cutoff=contact_threshold)
            unique_residue_contacts += contact_mat[:-1]
            if np.count_nonzero(unique_residue_contacts) > 25000:
                return True
            
        else:
            unique_residue_contacts = np.zeros((n_res,n_res))
        
    return False



def chol_contact_site(u,sterol_id):
    '''duration (or which frames) of cholesterol's contact with given site'''

    pass

# make function to find contact duration (percent contact) and maximal occupancy (longest contact)