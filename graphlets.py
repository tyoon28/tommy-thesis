from network_statistics import *
import random
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoLocator
from sklearn.model_selection import train_test_split 
from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression
from subprocess import call
from matplotlib.pyplot import cm
from scipy import spatial
from chimera_script import *
from pathlib import Path
import statsmodels.api as sm 
import hypernetx as hnx
from matplotlib.pyplot import cm
import pickle 





one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
'GLY':'G', 'PRO':'P', 'CYS':'C'}


def output_sampled_graphs(u,numgraphs,basename,lengraph=1000):
    starts = random.sample(range(len(u.trajectory)-lengraph), numgraphs)
    for i in starts:
        fn = f'{basename}-{i}-{i+lengraph}.in'
        output_consensus_graph(u,fn,s=i,d=i+lengraph)
    return

def output_uniform_graphs(u,numgraphs,basename,lengraph=1000):
    starts = np.linspace(0,len(u.trajectory)-lengraph,numgraphs)
    starts = np.rint(starts).astype(np.int32) # round integer
    for i in tqdm.tqdm(starts):
        fn = f'{basename}-{i}-{i+lengraph}.in'
        output_consensus_graph(u,fn,s=i,d=i+lengraph)
    return



def graphlet_degree_distribution(filename):
    # get graphlet degree distribution of orca output
    gdd = np.zeros(73)
    with open(filename) as f:
        for line in f:
            l = np.fromstring(line,dtype = int,sep = ' ')
            gdd += l
    return gdd

def PCA_gdd(ldirs):
    # ldirs: list of directories containing separate classes of networks
    #ldirs = ['/Users/Tommy/Desktop/thesis/orca/output/R1-15-closed-uniform','/Users/Tommy/Desktop/thesis/orca/output/R1-30-closed-uniform']
    gdds = []
    for d in ldirs:
        for f in os.listdir(d):
            if f == '.DS_Store': continue
            gdd = graphlet_degree_distribution(os.path.join(d, f))

            if '30' in f: chol = 30
            else: chol = 15
            if 'closed' in f: state = 'closed'
            else: state = 'open'
            gdds.append([f,chol,state] + list(gdd))
    
    df = pd.DataFrame(gdds,columns=["name", "chol", "state"] + list(range(73)))
    # https://builtin.com/machine-learning/pca-in-python
    features = list(range(73))
    x = df.loc[:, features].values
    y = df.loc[:,['chol']].values
    x = StandardScaler().fit_transform(x)
    pca = PCA()
    principalComponents = pca.fit_transform(x)
    evr = pca.explained_variance_ratio_.cumsum()
    for i,j in enumerate(evr):
        if j > 0.99:
            nPCs = i + 1
    pca = PCA(n_components=nPCs)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data = principalComponents
             , columns = [f'PC{x}' for x in range(1,nPCs+1)])
    finalDf = pd.concat([principalDf, df[['chol','name']]], axis = 1)
    finalDf['start'] = finalDf['name'].str.split('-').str[3]
    finalDf['start'] = finalDf['start'].apply(int)

    return finalDf,pca


def plot_PCA_gdd(finalDf,out):
    # PLOTTING PCA. only for 15 and 30
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1) 
    # ax.set_xlabel('Principal Component 1', fontsize = 15)
    ax.set_xlabel('PC 1', fontsize = 15)

    ax.set_ylabel('PC 2', fontsize = 15)
    ax.set_title('PCA graphlet degree distribution', fontsize = 20)

    targets = sorted(finalDf['chol'].unique())
    color = iter(cm.rainbow(np.linspace(0, 1, len(targets))))


    # ax.scatter(finalDf['start']
    #             , finalDf['PC2']
    #             , c = 'r')
    for target in targets:
        indicesToKeep = finalDf['chol'] == target
        if len(indicesToKeep) > 10:
            c = next(color)
            ax.scatter(finalDf.loc[indicesToKeep, 'PC1']
                    , finalDf.loc[indicesToKeep, 'PC2'],
                    color=c)
            
    ax.legend(targets)
    ax.grid()
    ax.yaxis.set_major_locator(AutoLocator())
    ax.xaxis.set_major_locator(AutoLocator())
    plt.savefig(f'{out}.png')

def PCA_logistic_selection(finalDf,pca,nPCs):
    #TODO: do this including all replicates. maybe do it with just graphlets instead of PCs.
    X = finalDf[[f'PC{x}' for x in range(1,nPCs+1)]]
    y = finalDf['chol']
    X_train, X_test, y_train, y_test = train_test_split(X,y , 
                                    random_state=104,  
                                    train_size=0.8,  
                                    shuffle=True) 
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    model = LogisticRegression(solver='liblinear', random_state=0)
    model.fit(X_train, y_train)
    model.coef_
    y_pred = model.predict(X_test)
    print('Accuracy of logistic regression classifier on test set: {:.2f}'.format(model.score(X_test, y_test)))
    coefficients = model.coef_[0]

    feature_importance = pd.DataFrame({'Feature': X.columns, 'Importance': np.abs(coefficients)})
    feature_importance = feature_importance.sort_values('Importance', ascending=True)
    feature_importance.plot(x='Feature', y='Importance', kind='barh', figsize=(10, 6))
    plt.show()




def node_degree_signatures(filename):
    # get node degree signatures of orca output
    nodes = []
    with open(filename) as f:
        for line in f:
            l = np.fromstring(line,dtype = int,sep = ' ')
            nodes.append(l)
    return nodes

def do_30_15(replicate):
    '''make uniform snapshots for 30 and 15 conditions (length 500)'''
    for i in ['30','15']:
        xtcs = []
        for file in os.listdir(f'{replicate}-{i}-closed'):
            if file.endswith('.xtc'):
                xtcs.append(f'{replicate}-{i}-closed/'+file)
        xtcs.sort(key=lambda x: int(x.split('-')[1]))
        u = mda.Universe(f'{replicate}-{i}-closed/{replicate}-0-start-membrane-3JYC.pdb',*xtcs)
        output_uniform_graphs(u,400,f'../orca/input/{replicate}-{i}-closed/{replicate}-{i}-closed',lengraph=500)
    return 



def find_pca_loadings(pca,component):
    c = pca.components_[component-1]
    indices = np.argsort(c)
    loadings = [c[i] for i in indices]
    return

def nodes_in_graphlet(g,d):
    # which nodes participate in graphlet g?
    # which nodes are most likely to be in the graphlet together
    # read orca files
    # input: g = graphlet number, d = directory containing orca outputs

    # count appearances of nodes in graphlet
    nodes = {}

    for f in os.listdir(d):
        if f == '.DS_Store': continue
        filename = os.path.join(d,f)
        with open(filename) as fil:
            n = 0 # node number
            for line in fil:
                l = np.fromstring(line,dtype = int,sep = ' ')
                if l[g] != 0:
                    if n in nodes:
                        nodes[n] += l[g]
                    else:
                        nodes[n] = l[g]
                n += 1

    return



def PCA_node_signature(ldirs):
    # filenames - list of orca output files
    #ldirs = ['/Users/Tommy/Desktop/thesis/orca/output/R1-closed-15','/Users/Tommy/Desktop/thesis/orca/output/R1-closed-30']
    gdds = []
    c = 0
    for d in ldirs[:-1]:
        for f in os.listdir(d):
            if f == '.DS_Store': continue
            nds = node_degree_signatures(os.path.join(d, f))
            for i,j in enumerate(nds):
                if j[57]: 
                    c+=1
                    print(i)
    




def output_graphs_graphlets_cholesterol(replicate,thresh=0.4):
    '''make orca input for a replicate'''
    # which graphlets are associates with cholesterol occupancy (magnitude and boolean)?
    # find frames where cholesterol is bound
    #random sample of a bunch of frames
    winlen = 100
    for i in ['30','15']:
        print(f'making cholesterol graphs for {replicate}-{i}')
        xtcs = []
        for file in os.listdir(f'{replicate}-{i}-closed'):
            if file.endswith('.xtc'):
                xtcs.append(f'{replicate}-{i}-closed/'+file)
        xtcs.sort(key=lambda x: int(x.split('-')[1]))
        u = mda.Universe(f'{replicate}-{i}-closed/{replicate}-0-start-membrane-3JYC.pdb',*xtcs)
        
        starts = random.sample(range(len(u.trajectory)-winlen), 1000)
        cholesterol = {}
        sites = binding_sites('closed')
        basename = f'{replicate}-{i}-closed'

        # find how many cholesterols are binding in each sampled window
        print(f'calculating cholesterol binding for {replicate}-{i}')
        rog = Cholesterol_contact(u)
        rog.run(start=0,verbose=True)
        print(f'outputting graphs of length {winlen} for cholesterol binding for {replicate}-{i}')
        for t in tqdm.tqdm(starts):
            numbinding = 0
            numcontacts = 0
            for c in rog.results:
                a = np.array(rog.results[c]['binding_events_actual'])
                nb = np.count_nonzero(np.logical_and(a >= t, a <= t+winlen))

                b = np.array(rog.results[c]['contacts'])
                nc = np.array(rog.results[c]['contacts'])
        
                if nb >= winlen * thresh:
                    numbinding += 1
                if nc >= winlen * thresh:
                    numcontacts += 1
                output_consensus_graph(u,f'../orca/input/{basename}-contact/{basename}-contact-c{numcontacts}-cc{numbinding}-f{f}.in',s=f,d=f+winlen)
        print(f'done with cholesterol graphs for {replicate}-{i}')
    return

def graphlets_cholesterol_pca(r):
    '''analyze orca output. r is R1-R3. Orca output must be in ../orca/output/{r}-<x>-closed'''
    #PCA
    if r == 'all':
        ldirs = [f'../orca/output/{r}-{c}-closed-contact' for r in ['R1','R2','R3'] for c in ['15,30']]
    else:
        ldirs = [f'../orca/output/{r}-15-closed-contact',f'../orca/output/{r}-30-closed-contact']
    gdds = []
    for d in ldirs:
        for f in os.listdir(d):
            if f == '.DS_Store': continue
            gdd = graphlet_degree_distribution(os.path.join(d, f))

            chol = f.split('-')[4][1:]
            chol_bind = f.split('-')[5][2:]

            if 'closed' in f: state = 'closed'
            else: state = 'open'
            if '-15-' in f:
                chol_condition = 15
            else:
                chol_condition = 30

            gdds.append([f,chol,chol_bind,chol_condition,state] + list(gdd))
    
    
    df = pd.DataFrame(gdds,columns=["name", "chol",'chol_bind','chol_condition', "state"] + list(map(str,range(73))))
    df['chol'] = df['chol'].apply(int)
    df['chol_bind'] = df['chol_bind'].apply(int)


    # https://builtin.com/machine-learning/pca-in-python
    # see difference between cholesterol contacts and cholesterol "binding"
    for column in ['chol','chol_bind']:
        features = list(map(str,range(73)))
        x = df.loc[:, features].values
        y = df.loc[:,['chol']].values
        x = StandardScaler().fit_transform(x)

        # take PCA components with 99% variance explained
        pca = PCA()
        principalComponents = pca.fit_transform(x)
        evr = pca.explained_variance_ratio_.cumsum()
        for i,j in enumerate(evr):
            if j > 0.99:
                nPCs = i + 1
        pca = PCA(n_components=nPCs)
        principalComponents = pca.fit_transform(x)
        principalDf = pd.DataFrame(data = principalComponents
                , columns = [f'PC{x}' for x in range(1,nPCs+1)])
        
        finalDf = pd.concat([principalDf, df[['chol','name']]], axis = 1)
        finalDf['start'] = finalDf['name'].str.split('-').str[5].str[1:-4]
        finalDf['start'] = finalDf['start'].apply(int)

        plot_PCA_gdd(finalDf,f'{r}_graphlet_cholesterol_{column}')



    #bar chart to show prevalence of graphlet x in each cholesterol state
    # result_group_chol= df.groupby(['chol'])
    # total_by_chol = result_group_chol['14'].agg([np.mean])

    # pl = total_by_chol.plot(kind='bar' ,y='mean',rot=0)

    # plt.xlabel('# cholesterol "bound"')
    # plt.ylabel('average # graphlet 8')
    # plt.savefig('graphlet_cholesterol')
    # plt.clf()

    # bar chart to show mean euclidean distance from zero in each cholesterol state
    # do an ANOVA here...
    result_group_chol= df.groupby(['chol'])
    total_by_chol = result_group_chol[[str(i) for i in range(73)]].mean()
    tbc = np.linalg.norm(total_by_chol.values,axis=1)
    tbc = tbc-tbc[0]
    total_by_chol['graphlet vector L2 norm'] = tbc
    total_by_chol.plot(kind='bar' ,y='graphlet vector L2 norm',rot=0)
    plt.savefig(f'{r}_movement_by_chol')
    return

def output_single_frame_graph(u,frame,basename):
    res = u.select_atoms('not resname CHOL and not resname POPC')
    lenr = len(res.residues)

    frame = u.trajectory[frame]
    r = res.atoms.center_of_mass(compound='residues')
    mat = distances.contact_matrix(r, cutoff=6)
    np.fill_diagonal(mat, 0)

    G = nx.from_numpy_array(mat)
    
    with open(f'{basename}-{u.trajectory.frame}.in', 'w') as f:
        f.write(f'{len(G.nodes)} {len(G.edges)}\n') # for orca input
        for line in nx.generate_edgelist(G,delimiter = ' ',data=False):
            linesp = line.split()
            f.write(f'{linesp[0]} {linesp[1]}\n')
    return


def output_single_sampled_graphs(u,numgraphs,basename,exclude):
    x = [i for i in range(len(u.trajectory)-1) if i not in exclude]
    
    starts = random.sample(x, numgraphs)
    for i in starts:
        output_single_frame_graph(u,i,basename)

    return

def output_for_node_pca(r):
    # output 2 consensus graph for the entire thing
    for conc in ['15','30']:
        xtcs = []
        for file in os.listdir(f'{r}-{conc}-closed'):
            if file.endswith('.xtc'):
                xtcs.append(f'{r}-{conc}-closed/'+file)
        xtcs.sort(key=lambda x: int(x.split('-')[1]))
        u = mda.Universe(f'{r}-{conc}-closed/{r}-0-start-membrane-3JYC.pdb',*xtcs,continuous=True)

        filename = f'../orca/input/{r}-full/{r}-{conc}-closed-full.in'
        output_consensus_graph(u,filename)
    ## RUN ORCA
    return

def node_pca_analysis(r,output=False):
    # get graphlet composition of each node
    if r == 'all':
        outdirs = [f'../orca/output/{r}-full' for r in ['R1','R2','R3']]
    else: outdirs = [f'../orca/output/{r}-full']
    #outdirs = ['/Users/Tommy/Desktop/thesis/orca/output/R1-15-closed-uniform','/Users/Tommy/Desktop/thesis/orca/output/R1-30-closed-uniform']
    rows = []
    for d in outdirs:
        for fn in os.listdir(d):
            if fn == '.DS_Store': continue
            if '-30-' in fn: chol = 30
            else: chol = 15
            if 'closed' in fn: state = 'closed'
            else: state = 'open'

            with open(os.path.join(d, fn)) as f:
                for i,line in enumerate(f): 
                    l = line.split(' ')
                    rows.append([fn,i,chol,state] + l)
    df = pd.DataFrame(rows,columns=["name", 'node',"chol", "state"] + list(range(73)))

    # load mda universe. just for getting residue names
    #u = mda.Universe(f'{r}-30-closed/{r}-0-start-membrane-3JYC.pdb',f'{r}-30-closed/{r}-0-1000-3JYC.xtc')
    u = mda.Universe(f'R1-30-closed/R1-0-start-membrane-3JYC.pdb',f'R1-30-closed/R1-0-1000-3JYC.xtc')


    #PCA
    features = list(range(73))
    x = df.loc[:, features].values
    y = df.loc[:,['chol']].values
    x = StandardScaler().fit_transform(x)
    pca = PCA(n_components=37)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data = principalComponents
             , columns = [f'PC{x}' for x in range(1,38)])
    finalDf = pd.concat([principalDf, df[['chol','name','node']]], axis = 1)
    PCs = [f'PC{x}' for x in range(1,38)]
    
    # PCA see which nodes move the most
    df[features] = df[features].astype(int)
    result_group_node= finalDf.groupby('node')
    distance_by_node = result_group_node[PCs].apply(lambda x: spatial.distance.pdist(x.to_numpy(),metric='euclidean').mean())
    if output:
        ax = distance_by_node.plot(kind='bar' ,y='distance_by_node',rot=0,ec='blue')
        thresh = 25
        for p in ax.patches:
            resid = int(p.get_x() + 1.25)
            res = int(MD_to_chimerax(resid)[5:])
            resname = one_letter[u.residues[resid-1].resname]
            label = resname + str(res)
            if p.get_height() > thresh:
                ax.annotate(
                    label, xy=(p.get_x(), p.get_height() + 0.1), fontsize=4
                )
        plt.tight_layout()
        plt.savefig(f'{r}-nodemovement.png')


    # These results are zero-indexed. MD results are 1-indexed so change to 1
    distance_by_node.index += 1 
    d_s = distance_by_node.sort_values()

    dd=d_s.to_dict()
    if output: dd_to_csv(dd,f'{r}-nodepca',u)


    return dd

def replicates_node_graphlet():
    bigd = {}
    for r in ['R1','R2','R3']:
        d = node_pca_analysis(r)
        for i in d:
            if i not in bigd:
                bigd[i] = d[i]
            else:
                bigd[i] += d[i]
    dd_to_csv(bigd,f'R1-3combined',u)



def node_graphlets_cholesterol(u,basename):
    outdirs = ['/Users/Tommy/Desktop/thesis/orca/output/R1-15-closed-cholesterol-threshold40','/Users/Tommy/Desktop/thesis/orca/output/R1-30-closed-cholesterol-threshold40']
    rows = []
    for d in outdirs:
        for fn in os.listdir(d):
            if fn == '.DS_Store': continue
            if '-30-' in fn: chol = 30
            else: chol = 15
            if 'closed' in fn: state = 'closed'
            else: state = 'open'

            with open(os.path.join(d, fn)) as f:
                for i,line in enumerate(f): # count resid. does it start at 0...?
                    l = line.split(' ')
                    rows.append([fn,i,chol,state] + l)
    df = pd.DataFrame(rows,columns=["name", 'node',"chol", "state"] + list(range(73)))

    #PCA
    features = list(range(73))
    x = df.loc[:, features].values
    y = df.loc[:,['chol']].values
    x = StandardScaler().fit_transform(x)
    pca = PCA(n_components=50)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data = principalComponents
             , columns = [f'PC{x}' for x in range(1,51)])
    finalDf = pd.concat([principalDf, df[['chol','name','node']]], axis = 1)
    PCs = [f'PC{x}' for x in range(1,51)]

    # PCA see which nodes move the most
    df[features] = df[features].astype(int)
    result_group_node= finalDf.groupby('node')
    distance_by_node = result_group_node[PCs].apply(lambda x: spatial.distance.pdist(x.to_numpy(),metric = 'cosine').mean())
    distance_by_node = distance_by_node[~(distance_by_node >= 1)]
    distance_by_node.plot(kind='bar' ,y='distance_by_node',rot=0)

    # These results are zero-indexed. MD results are 1-indexed so change to 1
    distance_by_node.index += 1 
    d_s = distance_by_node.sort_values()

    dd=d_s.to_dict()
    color_by_centrality(dd,'pca_node_movement_cosine')
    MD_to_chimerax()
    return


def L264(r='R1'):
    # L 264 moves a lot in the graphlet analysis. so does I177.(Transmembrane)
    # (it's 229 in the MD resid)
    outdirs = [f'../orca/output/{r}-full']
    #outdirs = ['/Users/Tommy/Desktop/thesis/orca/output/R1-15-closed-uniform','/Users/Tommy/Desktop/thesis/orca/output/R1-30-closed-uniform']
    rows = []
    for d in outdirs:
        for fn in os.listdir(d):
            if fn == '.DS_Store': continue
            if '-30-' in fn: chol = 30
            else: chol = 15
            if 'closed' in fn: state = 'closed'
            else: state = 'open'

            with open(os.path.join(d, fn)) as f:
                for i,line in enumerate(f): 
                    l = line.split(' ')
                    rows.append([fn,i,chol,state] + l)
    df = pd.DataFrame(rows,columns=["name", 'node',"chol", "state"] + list(range(73)))

    u.select_atoms('resid 143').residues.resnames
    dff = df.loc[df['node'].isin([143+337*i for i in range(4)])] 
    dff.to_csv('G178_nodegraphlet.csv')
    # G4, 5,15!!,16,19,


    # load mda universe. just for getting residue names
    #u = mda.Universe(f'{r}-30-closed/{r}-0-start-membrane-3JYC.pdb',f'{r}-30-closed/{r}-0-1000-3JYC.xtc')
    u = mda.Universe(f'R1-30-closed/R1-0-start-membrane-3JYC.pdb',f'R1-30-closed/R1-0-1000-3JYC.xtc')

    dff = df.loc[df['node'].isin([228+337*i for i in range(4)])] 
    dff['node'] = dff['node'] + 1
    dff
    dff.to_csv('L264_nodegraphlet.csv')
    return



def residue_investigation(residues,r):
    #residues: list of residues to look at
    # r: replicate
    # what does it touch in 15 vs 30 condition?

    d = {} # key: residue; item: dict of concentrations with dict of contacts and frequencies
     
    

    for i in ['30','15']:
        d[i] = {}
        xtcs = []
        for file in os.listdir(f'{r}-{i}-closed'):
            if file.endswith('.xtc'):
                xtcs.append(f'{r}-{i}-closed/'+file)
        xtcs.sort(key=lambda x: int(x.split('-')[1]))
        u = mda.Universe(f'{r}-{i}-closed/{r}-0-start-membrane-3JYC.pdb',*xtcs,continuous=True)

        protein = u.select_atoms('not resname CHOL and not resname POPC')
        for res in residues:
            # record frequency of each unique set of contacts
            d[i][res] = {}
            print(res,i)
            for ts in tqdm.tqdm(u.trajectory):
                frame = u.trajectory.frame
                r_compound = protein.atoms.center_of_mass(compound='residues')
                mat = distances.contact_matrix(r_compound, cutoff=6)
                np.fill_diagonal(mat, 0) 
                # +/- 1 from resid for proper indexing
                nz = np.nonzero(mat[res-1])[0] +1
                fs = frozenset(nz)
                if fs in d[i][res].keys():
                    d[i][res][fs] += 1 
                else: d[i][res][fs] = 1 
    
        
    d_difference = {}
    for res in residues:
        d_difference[res] = {}

        s = set(d['15'][res].keys()).union(set(d['30'][res].keys()))

        for fs in s:
            if fs in d['30'][res].keys() and fs in d['15'][res].keys():
                d_difference[res][fs] = d['30'][res][fs] - d['15'][res][fs]
            elif fs in d['30'][res].keys():
                d_difference[res][fs] = d['30'][res][fs]
            elif fs in d['15'][res].keys():
                d_difference[res][fs] = -d['15'][res][fs]
            else:
                print('problem')
    d_diff_filtered = {}
    for res in residues:
        d_diff_filtered[res] = {}
        for i in d_difference[res]:
            if abs(d_difference[res][i]) > 10000:
                d_diff_filtered[res][i] = d_difference[res][i]
        d_diff_filtered[res] = {k: v for k, v in sorted(d_diff_filtered[res].items(), key=lambda item: item[1])}


    # idea: make a hypergraph for each condition
        # hyperedge denotes frozensets, with weight equal to fraction of time they exist
        # then 
    
    # idea: hypergraph that combines conditions
        # hyperedge is frozensets, with weight equal to difference in time in high vs low cholesterol conditions

    fig, axs = plt.subplots(2,3,figsize=(16,9), gridspec_kw={'height_ratios': [1, 0.05]})

    # plt1 = fig.add_subplot(1,3,1)
    # plt2 = fig.add_subplot(1,3,2)
    # plt3 = fig.add_subplot(1,3,3)
    # cbar = fig.add_subplot(2,3,4,figsize=(6.4, 0.35))
    x = list(range(len(d_difference[res])))
    y = sorted(d_difference[res].values(),reverse=True)
    plt.plot(x,y)

    for p in [1,2,3]:
        res = residues[p-1]
        ax = axs[0,p-1]
        plt.sca(ax)
        hedges = {}
        hp = {}
        key = 0
        for k,v in sorted(d_diff_filtered[res].items(), key=lambda item: item[1]):
            key += 1
            hedges[key] = frozenset(map(lambda x: MD_to_resid(x,u),k))
            hp[key] = {'weight':v/100000}
        H = hnx.Hypergraph(hedges,edge_properties = hp)
        color = [H.edge_properties[e]['weight'] for e in H.edges()]
        norm = plt.Normalize(-1, 1)
        color = cm.bwr(norm(color))
        alpha = 0.9
        # hnx.draw(H)
        hnx.drawing.draw(H,
            edges_kwargs={
                'facecolors': color*(1, 1, 1, alpha),
                'edgecolors': 'black',
                'linewidths': 2
            },
            with_edge_labels=False
        )
        plt.title(f'residue {MD_to_resid(res,u)}')
    ax = axs[1,0]
    plt.sca(ax)
    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack((gradient, gradient))
    plt.imshow(gradient, aspect='auto', cmap=cm.bwr)
    ax.text(-0.01, 0.5, '15',ha='right', va='center',fontsize=10,transform=ax.transAxes)
    ax.text(1.01, 0.5, '30', ha='left', va='center',fontsize=10,transform=ax.transAxes)

    for ax in axs[1]:
        plt.sca(ax)
        plt.xticks([])
        plt.yticks([])
    fig.delaxes(axs[1,1])
    fig.delaxes(axs[1,2])

    plt.show()
    
def wholenetwork_stable_PCA():
    gdds = []
    for d in ldirs:
        for f in os.listdir(d):
            if f == '.DS_Store': continue
            gdd = graphlet_degree_distribution(os.path.join(d, f))

            if '30' in f: chol = 30
            else: chol = 15
            if 'closed' in f: state = 'closed'
            else: state = 'open'
            gdds.append([f,chol,state] + list(gdd))
    
    df = pd.DataFrame(gdds,columns=["name", "chol", "state"] + list(range(73)))
    # https://builtin.com/machine-learning/pca-in-python
    features = list(range(73))
    x = df.loc[:, features].values
    y = df.loc[:,['chol']].values
    x = StandardScaler().fit_transform(x)
    pca = PCA()
    principalComponents = pca.fit_transform(x)
    evr = pca.explained_variance_ratio_.cumsum()
    for i,j in enumerate(evr):
        if j > 0.99:
            nPCs = i + 1
    pca = PCA(n_components=nPCs)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data = principalComponents
             , columns = [f'PC{x}' for x in range(1,nPCs+1)])
    finalDf = pd.concat([principalDf, df[['chol','name']]], axis = 1)
    finalDf['start'] = finalDf['name'].str.split('-').str[3]
    finalDf['start'] = finalDf['start'].apply(int)

    return finalDf,pca


def node_PCA_windowed(r,ldirs,output=False):
    '''for each node do logistic regression to fit node distance ~ cholesterol. output p values.'''

    if r == 'all':
        ldirs = [f'../orca/output/{r}-{c}-closed-contact' for r in ['R1','R2','R3'] for c in ['15,30']]
    else:
        ldirs = [f'../orca/output/{r}-15-closed-contact',f'../orca/output/{r}-30-closed-contact']

    rows = []
    for d in ldirs:
        for fn in os.listdir(d):
            if fn == '.DS_Store': continue
            chol = f.split('-')[4][1:]
            if 'closed' in f: state = 'closed'
            else: state = 'open'
            if '-15-' in f:
                chol_condition = 15
            else:
                chol_condition = 30

            with open(os.path.join(d, fn)) as f:
                for i,line in enumerate(f): 
                    l = line.split(' ')
                    rows.append([fn,i+1,chol,state] + l)
    df = pd.DataFrame(rows,columns=["name", 'node',"chol",'chol_condition', "state"] + list(range(73)))


    # load mda universe. just for getting residue names
    #u = mda.Universe(f'{r}-30-closed/{r}-0-start-membrane-3JYC.pdb',f'{r}-30-closed/{r}-0-1000-3JYC.xtc')
    u = mda.Universe(f'R1-30-closed/R1-0-start-membrane-3JYC.pdb',f'R1-30-closed/R1-0-1000-3JYC.xtc')

    for n in range(1,1349):
            
        #PCA
        features = list(range(73))
        x = df.loc[:, features].values
        y = df.loc[:,['chol']].values
        x = StandardScaler().fit_transform(x)
        pca = PCA()
        principalComponents = pca.fit_transform(x)
        evr = pca.explained_variance_ratio_.cumsum()
        for i,j in enumerate(evr):
            if j > 0.99:
                nPCs = i + 1
                break
        pca = PCA(n_components=nPCs)
        print(f'using {nPCs} components')
        principalComponents = pca.fit_transform(x)
        principalDf = pd.DataFrame(data = principalComponents,
                                    columns = [f'PC{x}' for x in range(1,nPCs+1)])

        finalDf = pd.concat([principalDf, df[['chol','name','node']]], axis = 1)
        PCs = [f'PC{x}' for x in range(1,nPCs + 1)]

        #now fit logit model
    
    # PCA see which nodes move the most
    df[features] = df[features].astype(int)
    result_group_node= finalDf.groupby(['chol','node'])
    distance_by_node = finalDf.groupby('node').apply(lambda x: spatial.distance.pdist(np.array(x[PCs])).mean())
    if output:
        ax = distance_by_node.plot(kind='bar' ,y='distance_by_node',rot=0,ec='blue')
        thresh = 25
        for p in ax.patches:
            resid = int(p.get_x() + 1.25)
            res = int(MD_to_chimerax(resid)[5:])
            resname = one_letter[u.residues[resid-1].resname]
            label = resname + str(res)
            if p.get_height() > thresh:
                ax.annotate(
                    label, xy=(p.get_x(), p.get_height() + 0.1), fontsize=4
                )
        plt.tight_layout()
        plt.savefig(f'{r}-nodemovement-windowed.png')


    d_s = distance_by_node.sort_values()

    dd=d_s.to_dict()
    if output: dd_to_csv(dd,f'{r}-nodepca',u)


    return dd