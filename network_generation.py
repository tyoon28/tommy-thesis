'''
Functions for making networks from MD data.
'''


import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PSF, DCD, GRO, XTC
from MDAnalysis.analysis import contacts, distances
import numpy as np
from matplotlib import pyplot as plt
import networkx as nx
import itertools
import warnings
# suppress some MDAnalysis warnings about PSF files
warnings.filterwarnings('ignore')


def network_1_frame(u):
    '''protein contact network from 1 frame.'''

    #Contact is if centers of mass are within 6A'
    contact_threshold = 6
    # select all residues except cholesterol and POPC
    res = u.select_atoms('not resname CHOL and not resname POPC').atoms.center_of_mass(compound='residues')

    # contact matrix
    mat = distances.contact_matrix(res, cutoff=contact_threshold)

    # Graph from adjacency matrix
    G = nx.from_numpy_array(mat)

    return G

    # # plot residue distance
    # a = [(x,y) for x,y in itertools.product(range(n_res),range(n_res)) if sq_dist_res[x][y] < 6]
    # fig, ax = plt.subplots()
    
    # ax.set_aspect('equal')
    # # add figure labels and titles
    # plt.ylabel('Residue IDs')
    # plt.xlabel('Residue IDs')
    # ax.grid(True)

    # plt.scatter(*zip(*a),s=2,c='black')

    # plt.show()

def contact_percentage_network(u, threshold = 0.9):
    '''protein contact graph preserving only contacts present above a certain percentage of all frames'''

    # 6A is contact definition
    contact_threshold = 6

    # empty flattened triangular distance array
    n_res = len(u.select_atoms('not resname CHOL and not resname POPC').atoms.center_of_mass(compound='residues'))
    contact = np.zeros((n_res,n_res))

    
    for ts in u.trajectory:
        time = u.trajectory.time
        res = u.select_atoms('not resname CHOL and not resname POPC').atoms.center_of_mass(compound='residues')

        contact += distances.contact_matrix(res, cutoff=contact_threshold)

    contact /= len(u.trajectory)

    # graph from adjacency matrix, cutoff at threshold
    G = nx.from_numpy_array((contact > threshold))

    return G

def find_binding_sterols(u):
    '''
    Find cholesterols with binding events. 
    Barbera 2018: "for each cholesterol interaction, 
    the number of unique residue contacts at each time step during the interaction was counted. 
    Binding events were defined as those interactions in which the total number 
    of unique residue contacts over the duration of the cholesterol interaction 
    exceeded a threshold of 25,000."
    '''

    pass

 # each network should be an event. this should return a list of networks.
 # find stretches of contact with protein and then 
def chol_net(u,sterol_id):
    sel_sterol = f'resname CHOL and resid {sterol_id} and not (name RC1 or name RC2)'
    sel_protein = 'not resname CHOL and not resname POPC'

    # why do i need a reference group.
    ca = contacts.Contacts(u,
                       select=(sel_sterol, sel_protein),
                       refgroup=(acidic, basic),
                       method=fraction_contacts_between,
                       radius=5.0,
                       kwargs={'radius': 5.0,
                               'min_radius': 2.4}).run()
    
    pass


def cholesterol_network(u, sterol_id):
    '''
    Node size = relative duration of contact between the sterol of interest and the residue in question
    Edges connecting nodes indicate concurrent contacts with the sterol. 
    The edges are weighted according to the relevant correlation coefficient
    '''

    contact_threshold = 6
    protein = u.select_atoms('not resname CHOL and not resname POPC')

    n_res = len(protein.atoms.center_of_mass(compound='residues'))

    # each columns is a frame and each row is a residue. 1 means residue contacts sterol at given frame.
    res_mat = np.zeros((len(u.trajectory),n_res))

    for ts in u.trajectory:
        time = u.trajectory.time

        # excluding cholesterol tails from center of mass calculations
        sterol = (u.select_atoms(f'resname CHOL and resid {sterol_id} and not (name RC1 or name RC2)')
                  .atoms.center_of_mass(compound='residues'))


        res = u.select_atoms('not resname CHOL and not resname POPC').atoms.center_of_mass(compound='residues')

        # worried about efficiency using distance array but can't find built in contact function for two selections...
        #distance array between sterol and all other residues
        dist_arr = distances.distance_array(sterol, 
                                    res, 
                                    box=u.dimensions)
        
        # store contacts in res_mat
        res_mat[ts.frame] = dist_arr < contact_threshold

    # correlation matrix for contacts with sterol.
    x = np.corrcoef(res_mat,rowvar=False)

    
    # make graph from correlation matrix. put this in a different function
    G = nx.from_numpy_array(x)
    
    nx.set_node_attributes(G, node_duration, name=None)


    # can iterate over all positions here since we only do this 100s-1000s of times.

    pass