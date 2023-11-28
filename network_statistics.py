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
import os
from scipy import stats
import seaborn as sns



# TODO: to find binding events I think will have to use a separate function?
class Cholesterol_contact(AnalysisBase):
    '''
    atomgroup: selection of cholesterol(s) and protein or all atoms
    results[id] is a dictionary containing:
            contacts - list of frames that contact the channel
            contact_time - total contact time
            longest_contact - duration (frames) of longest contact with channel
            binding_events - ranges of continuous interaction
            mostcontacts - most residues that were touched in one interaction
            residuecontacts - # residue contacts in each interaction
            need to add start of longest contact.
            
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


        self.sterol_residuecontacts = {i:0 for i in self.sterols.resids}
        self.currentinteractions = {i:[] for i in self.sterols.resids}

        for i in self.sterols.resids:
            self.results[i] = {'contacts':[],'binding_events':[],'mostcontacts':0,'residuecontacts':[]}



        self.protein_sterols =  self.protein + self.sterols
        self.contact_threshold = 6
        

    def _single_frame(self):
        # Called after the trajectory is moved onto each new frame.
        res = self.protein_sterols.atoms.center_of_mass(compound='residues')
        contact_mat = distances.contact_matrix(res, cutoff=self.contact_threshold)
        frame = self.ag.universe.trajectory.frame


        for i, sterol in enumerate(self.sterols):
            resid = sterol.resid

            # contact matrix index
            mat_ind = self.lenprotein + i

            # row of contact matrix corresponding to this sterol and the protein
            c = contact_mat[mat_ind][:self.lenprotein]
            nonzero = np.count_nonzero(c)

            if nonzero:
                # store that this frame is a contact
                self.results[resid]['contacts'].append(frame)

                # count how many residue contacts are happening
                self.sterol_residuecontacts[resid] += np.sum(c)
                self.currentinteractions[resid].append(frame)

            else:
                # end of interaction. store binding event and reset counters
                if self.currentinteractions[resid]:
                    startint = self.currentinteractions[resid][0]
                    endint = self.currentinteractions[resid][-1]
                    self.results[resid]['binding_events'].append((startint,endint))
                    self.results[resid]['residuecontacts'].append(self.sterol_residuecontacts[resid])
                
                
                if self.sterol_residuecontacts[resid] > self.results[resid]['mostcontacts']:
                    self.results[resid]['mostcontacts'] = self.sterol_residuecontacts[resid]

                self.sterol_residuecontacts[resid] = 0
                self.currentinteractions[resid] = []

    def _conclude(self):

        for sterol in self.sterols:
            resid = sterol.resid
            if self.currentinteractions[resid]:

                startint = self.currentinteractions[resid][0]
                endint = self.currentinteractions[resid][-1]
                self.results[resid]['binding_events'].append((startint,endint))
                self.results[resid]['residuecontacts'].append(self.sterol_residuecontacts[resid])
                
            if self.sterol_residuecontacts[resid] > self.results[resid]['mostcontacts']:
                    self.results[resid]['mostcontacts'] = self.sterol_residuecontacts[resid]

            self.currentinteractions[resid] = []
            self.sterol_residuecontacts[resid] = 0

        
        for i in self.results:
            self.results[i]['contact_time'] = len(self.results[i]['contacts'])
            self.results[i]['longest_contact'] = self._longest(self.results[i]['contacts'])

        # Called once iteration on the trajectory is finished.
        # Apply normalisation and averaging to results here.
        pass

    def _longest(self, a):
        # return duration of longest continuous contact given list of contact frames. assumes list is sorted.
        if not a: 
            return 0
        longest = 0
        last = a[0]-1
        s = 0
        for i in a:
            if i-1 == last:
                s += 1
            else:
                s = 1
            if s > longest: longest = s
            last = i

        return longest


def unique_residue_contacts(u,sterol_ids):
    #find number of unique residue contacts over each binding event (continuous frame stretch) for each sterol in sterol_ids.
    pass

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

    

def chol_contact_site(u,sterol_id):
    '''duration (or which frames) of cholesterol's contact with given site'''

    pass

def binding_sites(site,include_subunits = True):
    '''return lists of resids of binding sites'''
    if site == 'Ia':
        l =  [47, 51, 131, 134, 135]
    if site == 'Ib':
        l =  [48, 51, 52, 55, 127, 128, 131]
    if site == 'II':
        l =  [53, 56, 57, 132, 135, 136, 139]
    if site == 'III':
        l = [34, 47, 48, 50, 51, 54, 134, 139, 140, 144]
    if site == 'IV':
        l = [34, 35, 38, 46, 49, 53, 136, 139, 140]
    if site == 'all':
        l = [34, 35, 38, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 127, 128, 131, 132, 134, 135, 136, 139, 140, 144]
    
    if not l: return
    if include_subunits:
        x = l + [i+337 for i in l] + [i+674 for i in l] + [i+1011 for i in l]
    else: x = l
    return x
# make function to find contact duration (percent contact) and maximal occupancy (longest contact)


def cholesterol_occupancy(u,sites):
    '''
    return the number of cholesterols in contact with binding sites.
    could do this with the rog object...........
    '''

    contact_threshold = 6
    selectstring = 'resid ' +  ' or resid '.join(map(str,sites))
    l = len(sites)
    s = u.select_atoms(selectstring).residues
    chol = u.select_atoms(f'(resname CHOL and not (name RC1 or name RC2))').residues

    protein_sterols = s + chol

    r = protein_sterols.atoms.center_of_mass(compound='residues')
    contact_mat = distances.contact_matrix(r, cutoff=contact_threshold)

    # take only rows representing sites
    z = contact_mat[:len(sites),len(sites):]
    occupancy = np.sum(np.any(z,0))

    return occupancy

def cholesterol_occupancy2(rog):
    # for given timepoint return number of cholesterols bound. 
    # have to put a bound (touching binding site) result in the analysisbase object for this to work.
    pass


def mutual_exclusivity(u):
    '''
    Finds correlation between contact edges (residue-residue contacts only) over timeseries.
    '''
    contact_threshold = 6
    sterols = u.select_atoms(f'(resname CHOL and not (name RC1 or name RC2))').residues
    protein = u.select_atoms('not resname CHOL and not resname POPC').residues
    lenp = len(protein)

    #each row is a timestep. each column is an edge.
    numedges = int(((lenp**2)-lenp)/2)
    results = np.zeros((len(u.trajectory),numedges))

    for ts in tqdm.tqdm(u.trajectory):
        t = u.trajectory.frame
        res = protein.atoms.center_of_mass(compound='residues')
        contact_mat = distances.contact_matrix(res, cutoff= contact_threshold)
        #flatten contact matrix
        ind = np.triu_indices(lenp,k=1)
        flattened = contact_mat[ind]

        # mapping from residue pairs to edges:
        # edge between resiudes i, j is edge number int((n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1)
        # where n is length of protein.

        results[t] = flattened

    # get edge names, so we can label.
    edges = []
    for i in range(len(flattened)):
        edges.append((ind[0][i],ind[1][i]))
        

    # delete zero columns (edges that never exist)
    idx = np.argwhere(np.all(results[..., :] == 0, axis=0))
    results2 = np.delete(results, idx, axis=1)
    e = np.delete(edges,idx,axis=0)

    # get spearman correlation matrix and associated p values
    cormat, pvals = stats.spearmanr(results2)

    # get p values of all negative correlations
    l = pvals * (cormat<0)
    signeg = np.sum(cormat<0)


    # 57 and 63 are mutually exclusive?
    sub = cormat[57:100,57:100]
    ax = sns.heatmap(sub, linewidth=0)
    plt.xlabel('edge #')
    plt.ylabel('edge #')
    plt.title('correlation between edges')
    plt.show()



    
    return results
