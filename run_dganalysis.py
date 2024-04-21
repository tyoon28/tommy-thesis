from dynamic_graphlets import *
from hotelling.stats import hotelling_t2

def main():
    gotnames = False
    gdds = []
    ldirs = ['../dynamic_graphlets/output/']
    # ldirs = ['../dynamic_graphlets/output/R1-15-closed-len50/dcounts','../dynamic_graphlets/output/R1-30-closed-len50/dcounts']
    print('starting analysis of all')
    for d in ldirs:
        for f in os.listdir(d):
            if f == '.DS_Store': continue
            if 'dcounts' not in f: continue
            if not gotnames: 
                graphlet_names = dyn_graphlet_degree_distribution(os.path.join(d, f),get_graphlet_names = True)
                gotnames = True
            gdd = dyn_graphlet_degree_distribution(os.path.join(d, f))

            if '-30-' in f: chol = 30
            else: chol = 15
            if 'closed' in f: state = 'closed'
            else: state = 'open'
            if 'R1' in f: r = 'R1'
            elif 'R2' in f: r = 'R2'
            elif 'R3' in f: r = 'R3'
            gdds.append([f,chol,state,r] + list(gdd))
    print('making df')
    df = pd.DataFrame(gdds,columns=["name", "chol", "state","replicate"] + graphlet_names)
    # https://builtin.com/machine-learning/pca-in-python
    features = graphlet_names
    x = df.loc[:, features].values
    y = df.loc[:,['chol']].values
    x = StandardScaler().fit_transform(x)
    pca = PCA()
    principalComponents = pca.fit_transform(x)
    evr = pca.explained_variance_ratio_.cumsum()
    for i,j in enumerate(evr):
        if j > 0.99:
            nPCs = i+1
            break
    df.to_csv('dyngraphlets_forR.csv')
    print(f'using {nPCs} components')
    pca = PCA(n_components=nPCs)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data = principalComponents
             , columns = [f'PC{x}' for x in range(1,nPCs+1)])
    finalDf = pd.concat([principalDf, df[['chol','name','replicate']]], axis = 1)
    finalDf['start'] = finalDf['name'].str.split('-').str[3]
    finalDf['start'] = finalDf['start'].apply(int)
    print('done w pca')


    pcpair = PCA_logistic_selection(finalDf,pca,nPCs)

    plot_PCA_dyn_gdd(finalDf,pca,remote=True,fn='dyngraphlets-all',PCs=pcpair)
    plot_PCA_dyn_gdd(finalDf,pca,remote=True,fn='dyngraphlets-1-2-all')
    # plot_PCA_dyn_gdd(finalDf,pca,remote=True,fn='dyngraphlets-colorbyrep',PCs=pcpair,colorby='replicate')


    x = finalDf[finalDf['chol'] == 15].drop(columns = ['name','chol','start','replicate']).to_numpy()
    y = finalDf[finalDf['chol'] == 30].drop(columns = ['name','chol','start','replicate']).to_numpy()

    print(f'hotelling_p for all = {hotelling_t2(x,y)[2]}')



    print('doing nodes')
    nfeatures = 3728 # length of vector
    plt.clf()
    # load mda universe. just for getting residue names and n residues
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
    big= []
    replicates = []
    chol = []
    for fn in tqdm.tqdm(os.listdir(d)):
        if fn == '.DS_Store': continue
        if 'dgdv' not in fn: continue
        if r != 'all':
            if r not in fn: continue
        if '-15-' in fn: 
            c = 15
        else: 
            c = 30
        if 'R3' in fn:
            rep = 'R3'
        else:
            rep = 'R1/R2'
        
        x = np.loadtxt(os.path.join(d, fn),dtype=int).tolist()
        if len(x) == nres:
            b = np.append(x,np.full([len(x),1],c),1)
            big.extend(b)
        else:
            # sometimes a node has no edges through the length of the window.
            # so have to insert a row with zero.
            notseen = lookformissing(fn)
            x = np.loadtxt(os.path.join(d, fn),dtype=int)
            for i in notseen:
                x = np.insert(x,int(i),np.zeros(nfeatures),0)
            print(notseen)
            b = np.append(x,np.full([len(x),1],c),1)
            big.extend(b.tolist())
        
    
    print(f'making df with {len(big)}')
    df = pd.DataFrame(big,columns=list(range(nfeatures)) + ['chol'])
    df['node'] = df.index%nres + 1
    df['chol'] = chol
    df['replicate'] = replicates


    #PCA
    print('doing pca')
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

    # pvals = {}
    # for n in range(1,nres+1):
    #     x = finalDf[(finalDf['chol'] == 15) & (finalDf['node'] == n)].drop(columns = ['node','chol']).to_numpy()
    #     y = finalDf[(finalDf['chol'] == 30) & (finalDf['node'] == n)].drop(columns = ['node','chol']).to_numpy()
    #     try:
    #         pvals[n] = hotelling_t2(x,y)[2]
    #     except np.linalg.LinAlgError:
    #         print(':(',n)
    #         pvals[n] = 1


    # PCA see which nodes move the most
    print('calculating distance')
    df[features] = df[features].astype(int)
    result_group_node= finalDf.groupby(['chol','node'])
    distance_by_node = finalDf.groupby('node').apply(lambda x: spatial.distance.pdist(np.array(x[PCs])).mean())
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
    plt.savefig(f'{r}-dyn-nodemovement-repdiff.png')

    distance_by_node.index += 1
    d_s = distance_by_node.sort_values()
    dd=d_s.to_dict()
    dd_to_csv(dd,f'all-dynnodepca',u)
    for g in ['PC3','PC4']:

        print(f'looking for {g}')
        dic = {}
        for n in range(1,1349):
            d = finalDf.loc[df['node'] == n]
            dic[n] = np.mean(d[g].loc[df['chol'] == 15]) - np.mean(d[g].loc[df['chol'] == 30])
            # dic[n] = np.mean(d[o])
        dd_to_csv(dic,f'{r}-node{g}',u)
        color_by_centrality(dic,f'{r}-node{g}')


    print('done w all')

    for r in ['R1','R2','R3']:
        break
        print(f'running PCA whole for {r}')
        d = df.loc[df['replicate'] == r]
        d = d.reset_index(drop=True)
        
        features = graphlet_names
        x = d.loc[:, features].values
        y = d.loc[:,['chol']].values
        x = StandardScaler().fit_transform(x)
        # pca = PCA()
        # principalComponents = pca.fit_transform(x)
        # evr = pca.explained_variance_ratio_.cumsum()
        # for i,j in enumerate(evr):
        #     if j > 0.99:
        #         nPCs = i + 1
        nPCs = 2
        print(f'using {nPCs} components')
        pca = PCA(n_components=nPCs)
        principalComponents = pca.fit_transform(x)
        print(pca.explained_variance_ratio_.cumsum())
        principalDf = pd.DataFrame(data = principalComponents
                , columns = [f'PC{x}' for x in range(1,nPCs+1)])
        finalDf = pd.concat([principalDf, d[['chol','name']]], axis = 1)
        # print(finalDf[finalDf.isnull().any(axis=1)])
        pcpair = PCA_logistic_selection(finalDf,pca,nPCs,output_fig=False)
        plot_PCA_dyn_gdd(finalDf,pca,remote=True,fn=f'dyngraphlets-{r}',PCs=pcpair)
        print(f'done with PCA whole for {r}')
    
    # for r in ['R1','R2','R3','all']:
    #     print(f'running PCA nodes for {r}')
    #     dynamic_PCA_nodes(r)
    #     print(f'done with PCA nodes for {r}')
        

    
    return

def R12vs3():
    #compare node orbits in R1 and R2 vs R3
    print('comparing node orbits in R1 and R2 vs R3')
    nfeatures = 3728 # length of vector
    plt.clf()
    # load mda universe. just for getting residue names and n residues
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
    big= []
    replicates = []
    chol = []
    for fn in tqdm.tqdm(os.listdir(d)):
        if fn == '.DS_Store': continue
        if 'dgdv' not in fn: continue
        if r != 'all':
            if r not in fn: continue
        if '-15-' in fn: 
            c = 15
            start = 0
            end = nres
        else: 
            c = 30
            start = nres
            end = None
        if 'R3' in fn:
            rep = 'R3'
        else:
            rep = 'R1/R2'
        
        x = np.loadtxt(os.path.join(d, fn),dtype=int).tolist()
        if len(x) == nres:
            big.extend(x)
        else:
            # sometimes a node has no edges through the length of the window.
            # so have to insert a row with zero.
            notseen = lookformissing(fn)
            x = np.loadtxt(os.path.join(d, fn),dtype=int)
            for i in notseen:
                x = np.insert(x,int(i),np.zeros(nfeatures),0)
            print(notseen)
            big.extend(x.tolist())
        
        chol.extend([c for i in range(len(x))])
        replicates.extend([rep for i in range(len(x))])

    
    print(f'making df with {len(big)}')
    df = pd.DataFrame(big,columns=list(range(nfeatures)))
    df['node'] = df.index%nres + 1
    df['chol'] = chol
    df['replicate'] = replicates


    #PCA
    print('doing pca')
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
    

    finalDf = pd.concat([principalDf, df[['node','chol','replicate']]], axis = 1)
    PCs = [f'PC{x}' for x in range(1,nPCs+1)]


    # PCA see which nodes move the most
    print('calculating distance')
    df[features] = df[features].astype(int)
    result_group_node= finalDf.groupby(['replicate','node'])
    distance_by_node = finalDf.groupby('node').apply(lambda x: spatial.distance.pdist(np.array(x[PCs])).mean())
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
    plt.savefig(f'{r}-dyn-nodemovement-repdiff.png')
    plt.clf()


    d_s = distance_by_node.sort_values()

    dd=d_s.to_dict()
    dd_to_csv(dd,f'{r}-dyn-nodepca-repdiff',u)
    color_by_centrality(dd,f'{r}-dyn-nocepca-repdiff-color')
    
if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    main()