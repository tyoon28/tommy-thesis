from dynamic_graphlets import *

def main():
    gotnames = False
    gdds = []
    ldirs = ['../dynamic_graphlets/output/']
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
    plot_PCA_dyn_gdd(finalDf,pca,remote=True,fn='dyngraphlets-all')

    for r in ['R1','R2','R3']:
        print(f'running PCA whole for {r}')
        d = df.loc[df['replicate'] == r]
        features = graphlet_names
        x = d.loc[:, features].values
        y = d.loc[:,['chol']].values
        x = StandardScaler().fit_transform(x)
        pca = PCA(n_components = 4)
        principalComponents = pca.fit_transform(x)
        principalDf = pd.DataFrame(data = principalComponents
                , columns = [f'PC{x}' for x in range(1,5)])
        finalDf = pd.concat([principalDf, d[['chol','name']]], axis = 1)
        finalDf['start'] = finalDf['name'].str.split('-').str[3]
        finalDf['start'] = finalDf['start'].apply(int)
        plot_PCA_dyn_gdd(finalDf,pca,remote=True,fn=f'dyngraphlets-{r}')
        print(f'done with PCA whole for {r}')
    
    for r in ['R1','R2','R3','all']:
        print(f'running PCA nodes for {r}')
        dynamic_PCA_nodes(r)
        print(f'done with PCA nodes for {r}')



if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    main()