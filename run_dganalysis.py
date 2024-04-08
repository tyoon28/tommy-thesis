from dynamic_graphlets import *

def main():
    gotnames = False
    gdds = []
    ldirs = ['../dynamic_graphlets/output/']
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
            nPCs = i + 1
    print(f'using {nPCs} components')
    pca = PCA(n_components=nPCs)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data = principalComponents
             , columns = [f'PC{x}' for x in range(1,nPCs+1)])
    finalDf = pd.concat([principalDf, df[['chol','name']]], axis = 1)
    finalDf['start'] = finalDf['name'].str.split('-').str[3]
    finalDf['start'] = finalDf['start'].apply(int)
    print('done w pca')
    plot_PCA_dyn_gdd(finalDf,pca,remote=True,fn='dyngraphlets-all')
    print('done w all')

    for r in ['R1','R2','R3']:
        print(f'running PCA whole for {r}')
        d = df.loc[df['replicate'] == r]
        print(len(d))
        features = graphlet_names
        x = d.loc[:, features].values
        y = d.loc[:,['chol']].values
        x = StandardScaler().fit_transform(x)
        pca = PCA()
        principalComponents = pca.fit_transform(x)
        evr = pca.explained_variance_ratio_.cumsum()
        for i,j in enumerate(evr):
            if j > 0.99:
                nPCs = i + 1
        print(evr)
        print(f'using {nPCs} components')
        pca = PCA(n_components=nPCs)
        principalComponents = pca.fit_transform(x)
        principalDf = pd.DataFrame(data = principalComponents
                , columns = [f'PC{x}' for x in range(1,nPCs+1)])
        finalDf = pd.concat([principalDf, d[['chol','name']]], axis = 1)
        plot_PCA_dyn_gdd(finalDf,pca,remote=True,fn=f'dyngraphlets-{r}')
        print(f'done with PCA whole for {r}')
    
    for r in ['R1','R2','R3','all']:
        print(f'running PCA nodes for {r}')
        dynamic_PCA_nodes(r)
        print(f'done with PCA nodes for {r}')



if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    main()