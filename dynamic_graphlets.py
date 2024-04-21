from network_statistics import *
from network_generation import *
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy import spatial
from chimera_script import *
from matplotlib.ticker import AutoLocator
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression



one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
'GLY':'G', 'PRO':'P', 'CYS':'C'}



def output_temporal_graph(u,filename,s=0,d=None):
    res = u.select_atoms('not resname CHOL and not resname POPC')

    with open(filename, 'w') as f:
        for ts in u.trajectory[s:d]:
            frame = u.trajectory.frame
            r = res.atoms.center_of_mass(compound='residues')
            mat = distances.contact_matrix(r, cutoff=6)
            np.fill_diagonal(mat, 0)

            G = nx.from_numpy_array(mat)
            for line in nx.generate_edgelist(G,delimiter = ' ',data=False):
                linesp = line.split()
                f.write(f'{linesp[0]}\t{linesp[1]}\t{frame}\n')
    return


def output_uniform_graphs(u,numgraphs,basename,lengraph=1000):
    starts = np.linspace(0,len(u.trajectory)-lengraph,numgraphs)
    starts = np.rint(starts).astype(np.int32)
    for i in tqdm.tqdm(starts):
        fn = f'{basename}-{i}-{i+lengraph}.in'
        output_consensus_graph(u,fn,s=i,d=i+lengraph)
    return

def dyn_graphlet_degree_distribution(filename,get_graphlet_names = False):
    # get graphlet degree distribution of dcounts output
    if get_graphlet_names: x = 0
    else: x = 1
    gdd = []
    with open(filename) as f:
        for line in f:
            l = line.split()
            gdd.append(l[x])
    return gdd


def dyngraphlets_cholesterol(u):
    winlen = 50
    starts = random.sample(range(1000,len(u.trajectory)-winlen), 1000)
    cholesterol = {}
    sites = binding_sites('closed')

    ### CHOLESTEROL
    # # find how many cholesterols are binding in each sampled window
    # for t in tqdm.tqdm(starts):
    #     bound = sterol_occupancy_at_t(u,t,winlen,0.40)
    #     if bound not in cholesterol:
    #         cholesterol[bound] = [t]
    #     else:
    #         cholesterol[bound].append(t)
    # # compute graphlets for a sample of each level of cholesterol binding
    # for n in cholesterol:
    #     for f in cholesterol[n]:
    #         output_temporal_graph(u,f'{basename}-contact-c{n}-f{f}.in',s=f,d=f+winlen)
    ####
    # Just comparing 15 and 30
    
    for i in ['30','15']:
        xtcs = []
        for file in os.listdir(f'R1-{i}-closed'):
            if file.endswith('.xtc'):
                xtcs.append(f'R1-{i}-closed/'+file)
        xtcs.sort(key=lambda x: int(x.split('-')[1]))
        u = mda.Universe(f'R1-{i}-closed/R1-0-start-membrane-3JYC.pdb',*xtcs,continuous=True)
        starts = np.linspace(0,len(u.trajectory)-winlen,50).astype(np.int32)
        basename = f'R1-{i}-closed'
        for j in tqdm.tqdm(starts):
            fn = f'{basename}-{j}-{j+winlen}.in'
            output_temporal_graph(u,fn,s=j,d=j+winlen)


    # run count_tmp_graphlets
    input()
    gotnames = False
    

    gdds = []
    ldirs = ['/Users/Tommy/Desktop/thesis/dynamic_graphlets/output/R1-30-closed-len50/dcounts','/Users/Tommy/Desktop/thesis/dynamic_graphlets/output/R1-15-closed-len50/dcounts']
    for d in ldirs:
        for f in os.listdir(d):
            if f == '.DS_Store': continue
            if not gotnames: 
                graphlet_names = dyn_graphlet_degree_distribution(os.path.join(d, f),get_graphlet_names = True)
                gotnames = True
            gdd = dyn_graphlet_degree_distribution(os.path.join(d, f))

            if '-30-' in f: chol = 30
            else: chol = 15
            if 'closed' in f: state = 'closed'
            else: state = 'open'
            if 'R1' in f: r = 1
            elif 'R2' in f: r = 2
            elif 'R3' in f: r = 3
            gdds.append([f,chol,state,r] + list(gdd))
    
    df = pd.DataFrame(gdds,columns=["name", "chol", "state","replicate"] + graphlet_names)
    # https://builtin.com/machine-learning/pca-in-python
    features = graphlet_names
    x = df.loc[:, features].values
    y = df.loc[:,['chol']].values
    x = StandardScaler().fit_transform(x)
    pca = PCA(n_components = 4)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data = principalComponents
             , columns = [f'PC{x}' for x in range(1,5)])
    finalDf = pd.concat([principalDf, df[['chol','name']]], axis = 1)
    finalDf['start'] = finalDf['name'].str.split('-').str[3]
    finalDf['start'] = finalDf['start'].apply(int)

    #finalDf['distance'] = finalDf[]



def plot_PCA_dyn_gdd(finalDf,pca,remote=False,fn=None,PCs=('PC1','PC2'),colorby='chol'):
    # PLOTTING PCA. only for 15 and 30
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1) 
    # ax.set_xlabel('Principal Component 1', fontsize = 15)

    
    evr = pca.explained_variance_ratio_
    xExp = round(evr[int(PCs[0][2:])-1]*100,1)
    yExp = round(evr[int(PCs[1][2:])-1]*100,1)
    ax.set_xlabel(f'{PCs[0]} ({xExp}% Variance)', fontsize = 15)
    ax.set_ylabel(f'{PCs[1]} ({yExp}% Variance)', fontsize = 15)
    ax.set_title('PCA dynamic graphlet degree distribution', fontsize = 20)

    targets = sorted(finalDf[colorby].unique())
    if colorby == 'chol':
        color = iter(['#18a5ff','#d41159'])
    else:
        color = iter(cm.rainbow(np.linspace(0, 1, len(targets))))


    # ax.scatter(finalDf['start']
    #             , finalDf['PC2']
    #             , c = 'r')
    x,y = PCs
    for target in targets:
        indicesToKeep = finalDf[colorby] == target
        if len(indicesToKeep) > 10:
            c = next(color)
            
            
            ax.scatter(finalDf.loc[indicesToKeep, x]
                    , finalDf.loc[indicesToKeep, y],
                    c=c,s=8)
            
    if colorby == 'chol':
        ax.legend(['15 mol%','30 mol%'])
    else:
        ax.legend(targets)

    ax.grid(False)
    ax.yaxis.set_major_locator(AutoLocator())
    ax.xaxis.set_major_locator(AutoLocator())
    
    if remote:
        plt.savefig(f'{fn}.png')
    else:
        plt.show()
    return

def PCA_logistic_selection(finalDf,pca,nPCs,output_fig=False):
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
    if output_fig:
        feature_importance.plot(x='Feature', y='Importance', kind='barh', figsize=(10, 6))
        plt.savefig(f'dyngraphlet-varimportance.png')
        plt.clf()
    return (feature_importance['Feature'].iloc[-1],feature_importance['Feature'].iloc[-2])
    
def find_pca_loadings(pca,component):
    c = pca.components_[component-1]
    indices = np.argsort(c)
    loadings = [c[i] for i in indices]
    return


def dynamic_PCA_nodes(r):
    outdirs = [f'/Users/Tommy/Desktop/thesis/dynamic_graphlets/output/{r}-15-closed-len50/dgdv',f'/Users/Tommy/Desktop/thesis/dynamic_graphlets/output/{r}-30-closed-len50/dgdv']
    #outdirs = ['/Users/Tommy/Desktop/thesis/orca/output/R1-15-closed-uniform','/Users/Tommy/Desktop/thesis/orca/output/R1-30-closed-uniform']
    nfeatures = 3728 # length of vector
    plt.clf()
    # load mda universe. just for getting residue names and n residues
    if r != 'all':
        u = mda.Universe(f'{r}-30-closed/{r}-0-start-membrane-3JYC.pdb',f'{r}-30-closed/{r}-0-1000-3JYC.xtc')
    else:
        u = mda.Universe('R1-30-closed/R1-0-start-membrane-3JYC.pdb','R1-30-closed/R1-0-1000-3JYC.xtc')
    nres = len(u.select_atoms('not resname CHOL and not resname POPC').residues)

    rows = np.zeros((nres*2,nfeatures+2))
    for i in range(len(rows)):
        rows[i][0] = i % nres + 1 # set res id in [0]
        # set chol state in [1]
        if i //nres == 0:
            rows[i][1] = 15
        else: rows[i][1] = 30
        
    d = '../dynamic_graphlets/output'
    for fn in tqdm.tqdm(os.listdir(d)):
        if fn == '.DS_Store': continue
        if 'dgdv' not in fn: continue
        if r != 'all':
            if r not in fn: continue
        if '-15-' in fn: 
            start = 0
            end = nres
        else: 
            start = nres
            end = None
        try:
            rows[start:end, 2:] += np.loadtxt(os.path.join(d, fn),dtype=int)
        
        # sometimes a node has no edges through the length of the window.
        # so have to insert a row with zero.
        except ValueError:
            notseen = lookformissing(fn)
            x = np.loadtxt(os.path.join(d, fn),dtype=int)
            for i in notseen:
                x = np.insert(x,int(i),np.zeros(nfeatures),0)
            print(notseen)
            rows[start:end, 2:] += x

    

    df = pd.DataFrame(rows,columns=['node',"chol"] + list(range(nfeatures)))

    


    #PCA
    features = list(range(nfeatures))
    x = df.loc[:, features].values
    y = df.loc[:,['chol']].values
    x = StandardScaler().fit_transform(x)

    # pick # components
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
    

    finalDf = pd.concat([principalDf, df[['node','chol']]], axis = 1)
    PCs = [f'PC{x}' for x in range(1,nPCs+1)]

    # PCA see which nodes move the most
    df[features] = df[features].astype(int)
    result_group_node= finalDf.groupby('node')
    distance_by_node = result_group_node[PCs].apply(lambda x: spatial.distance.pdist(x.to_numpy(),metric='euclidean').mean())
    distance_by_node.plot(kind='bar' ,y='distance_by_node',rot=0,ec='blue')
    thresh = 25
    # for p in ax.patches:
    #     resid = int(p.get_x() + 1.25)
    #     res = int(MD_to_chimerax(resid)[5:])
    #     resname = one_letter[u.residues[resid-1].resname]
    #     label = resname + str(res)
    #     if p.get_height() > thresh:
    #         ax.annotate(label, xy=(p.get_x(), p.get_height() + 0.1), fontsize=4)
    plt.tight_layout()
    plt.savefig(f'{r}-dyn-nodemovement.png')
    plt.clf()



    d_s = distance_by_node.sort_values()

    dd=d_s.to_dict()
    dd_to_csv(dd,f'{r}-dyn-nodepca',u)
    color_by_centrality(dd,f'{r}-dyn-nocepca-color')
    



def lookformissing(fn):
    f = fn.split('.')[0]
    fn = f'../dynamic_graphlets/input/{f}.in'
    seen = []
    with open(fn) as f:
        for line in f:
            l = line.split()
            if l[0] not in seen: seen.append(l[0])
            if l[1] not in seen: seen.append(l[1])
    notseen = []
    for i in range(1348):
        if str(i) not in seen:
            notseen.append(str(i))
    return notseen