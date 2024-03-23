'''
Functions for making networks from MD data.
'''

from network_statistics import *
import dynetx as dn


# suppress some MDAnalysis warnings about PSF files
warnings.filterwarnings('ignore')
def output_consensus_graph(u,filename,s=0,d=None,threshold = 0.9,**kwargs):
    res = u.select_atoms('not resname CHOL and not resname POPC')
    lenr = len(res.residues)
    edgesmat = np.zeros(shape=(lenr,lenr))

    for ts in u.trajectory[s:d]:
        frame = u.trajectory.frame
        r = res.atoms.center_of_mass(compound='residues')
        mat = distances.contact_matrix(r, cutoff=6)
        np.fill_diagonal(mat, 0)
        edgesmat += mat

    
    t = len(u.trajectory[s:d]) * threshold
    G = nx.from_numpy_array((edgesmat >= t))

    #make folder if doesn't exist
    p = os.path.dirname(filename)
    Path(p).mkdir(parents=True, exist_ok=True)

    with open(filename, 'w') as f:
        f.write(f'{len(G.nodes)} {len(G.edges)}\n') # for orca input
        for line in nx.generate_edgelist(G,delimiter = ' ',data=False):
            linesp = line.split()
            f.write(f'{linesp[0]} {linesp[1]}\n')
    return


def network_1_frame(u):
    '''protein contact network from 1 frame.'''

    #Contact is if centers of mass are within 6A'
    contact_threshold = 6
    # select all residues except cholesterol and POPC
    res = u.select_atoms('not resname CHOL and not resname POPC').atoms.center_of_mass(compound='residues')

    # contact matrix
    mat = distances.contact_matrix(res, cutoff=6)

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



def contact_network_residues(u,res):
    '''
    static contact network including only select residue(s) and its neighbors (and edges between neighbors.)
    res: list of residue id's to focus on
    '''
    contact_threshold = 6
    chol = u.select_atoms(f'(resname CHOL and not (name RC1 or name RC2))').residues
    protein = u.select_atoms('not resname CHOL and not resname POPC').residues

    protein_sterols = protein + chol

    r = protein_sterols.atoms.center_of_mass(compound='residues')
    contact_mat = distances.contact_matrix(r, cutoff=contact_threshold)

    # take only selected residues and the residues they contact 
    z = contact_mat[res,:].any(0)
    adjacency = contact_mat[z][:,z]

    resids = [i for i in range(len(z)) if z[i] ]

    labels = {i:j for i,j in enumerate(resids)}
    G = nx.from_numpy_array(adjacency)
    nx.relabel_nodes(G, labels, copy=False)

    return G

def network_1_cholesterol(u,sterol_id):
    '''make a contact network centered on specified cholesterol'''
    contact_threshold = 6
    chol = u.select_atoms(f'resid {sterol_id}').residues
    protein = u.select_atoms('not resname CHOL and not resname POPC').residues

    protein_sterols = protein + chol

    r = protein_sterols.atoms.center_of_mass(compound='residues')
    contact_mat = distances.contact_matrix(r, cutoff=contact_threshold)

    # take only selected residues  
    z = contact_mat[len(protein),:]
    adjacency = contact_mat[z][:,z]

    resids = [i for i in range(len(z)) if z[i] ]

    labels = {i:j for i,j in enumerate(resids)}
    labels[len(labels)-1] = sterol_id
    
    G = nx.from_numpy_array(adjacency)
    nx.relabel_nodes(G, labels, copy=False)

    return G



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
        t = u.trajectory.time

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


def network_selected(u,residues):
    '''make a contact network including only the selected residues'''
    contact_threshold = 6
    residues.sort()
    selstring = 'resid ' + ' or resid '.join(map(str,residues))
    protein = u.select_atoms(selstring).residues

    r = protein.atoms.center_of_mass(compound='residues')
    contact_mat = distances.contact_matrix(r, cutoff=contact_threshold)

    

    labels = {i:j for i,j in enumerate(residues)}
    
    G = nx.from_numpy_array(contact_mat)
    nx.relabel_nodes(G, labels, copy=False)

    return G

def viz_consensus_graph(u,start=0,end=None,threshold = 0.9,**kwargs):
    # for gephi input
    res = u.select_atoms('not resname CHOL and not resname POPC')
    lenr = len(res.residues)
    edgesmat = np.zeros(shape=(lenr,lenr))

    for ts in tqdm.tqdm(u.trajectory[start:end]):
        frame = u.trajectory.frame
        r = res.atoms.center_of_mass(compound='residues')
        mat = distances.contact_matrix(r, cutoff=6)
        np.fill_diagonal(mat, 0)
        edgesmat += mat

    
    t = len(u.trajectory[start:end]) * threshold
    G = nx.from_numpy_array((edgesmat >= t))
    nx.write_gexf(G,'closed_15_90.gexf')

    # for i in nx.connected_components(G):
    #     if len(i) > 30:
    #         k = MD_list_to_chimerax(i)
    #         print(k)
    #         print()
    return

def MD_list_to_chimerax(s):

    l = sorted(s)
    model = ''
    i = -1
    out = []
    for n in l:
        m = f'#1.{n//337+1}'
        if m != model:
            model = m
            out.append(f'{model}:{n%337 + 35}')
            i +=1
        else:
            out[i] += f',{n%337 + 35}'

    return ''.join(out)





def hypergraph_1_frame(u):
    '''hypergraph inferred using '''
    return

def gephi_input(u):
    # for visualizing. not done
    res = u.select_atoms('not resname CHOL and not resname POPC')
    lenr = len(res.residues)
    edgesmat = np.zeros(shape=(lenr,lenr))

    for ts in u.trajectory:
        frame = u.trajectory.frame
        r = res.atoms.center_of_mass(compound='residues')
        mat = distances.contact_matrix(r, cutoff=6)
        np.fill_diagonal(mat, 0)
        edgesmat += mat
    
    t = len(u.trajectory) * 0.9
    G = nx.from_numpy_array((edgesmat >= t))
    two = nx.shortest_path(G,265)
    for i in binding_sites('all'):
        try:
            print(two[i])
        except KeyError:
            pass
    d = nx.to_pandas_edgelist(G,source='Source',target='Target')
    d['Type'] = 'undirected'
    d['Id'] = d.index
    d['Weight'] = 1
    d['binding'] = 1
    nx.write_gexf(G,'closed_15_90.gexf')