from network_statistics import *
import pandas as pd


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
        u = mda.Universe(f'R1-{i}-closed/R1-0-start-membrane-3JYC.pdb',*xtcs)
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
            gdds.append([f,chol,state] + list(gdd))
    
    df = pd.DataFrame(gdds,columns=["name", "chol", "state"] + graphlet_names)
    # https://builtin.com/machine-learning/pca-in-python
    features = graphlet_names
    x = df.loc[:, features].values
    y = df.loc[:,['chol']].values
    x = StandardScaler().fit_transform(x)
    pca = PCA(n_components = 3)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data = principalComponents
             , columns = [f'PC{x}' for x in range(1,4)])
    finalDf = pd.concat([principalDf, df[['chol','name']]], axis = 1)
    finalDf['start'] = finalDf['name'].str.split('-').str[3]
    finalDf['start'] = finalDf['start'].apply(int)


def plot_PCA_dyn_gdd(finalDf,pca):
    # PLOTTING PCA. only for 15 and 30
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1) 
    # ax.set_xlabel('Principal Component 1', fontsize = 15)
    ax.set_xlabel('PC 1', fontsize = 15)

    ax.set_ylabel('PC 2', fontsize = 15)
    ax.set_title('PCA dynamic graphlet degree distribution', fontsize = 20)

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
                    c=c)
            
    ax.legend(targets)
    ax.grid()
    ax.yaxis.set_major_locator(AutoLocator())
    ax.xaxis.set_major_locator(AutoLocator())

    plt.show()
