# not an ANOVA anymore really.
# fit a logistic model to each node and print out p values
from graphlets import *
import warnings
from scipy.stats import f_oneway
import statsmodels.api as sm
from statsmodels.formula.api import ols





def main():
    ldirs = [f'../orca/output/{r}-{c}-closed' for r in ['R1','R2','R3'] for c in ['15','30']]
    # load data from uniform slices into df from multiple replicates
    # gets graphlet orbits for each node. slices of 250
    rows = []
    for d in tqdm.tqdm(ldirs):
        for fn in os.listdir(d):
            if fn == '.DS_Store': continue
            if '-0-' in fn: continue # skip the one with the burn in time
            r = fn.split('-')[0]
            if '30' in fn: chol = 30
            else: chol = 15
            if 'closed' in fn: state = 'closed'
            else: state = 'open'
            with open(os.path.join(d, fn)) as f:
                for i,line in enumerate(f): # count resid. does it start at 0...?
                    l = line.split(' ')
                    rows.append([fn,i,chol,state,r] + l)

    gcols = ['g' + str(graphlet) for graphlet in range(73)]
    df = pd.DataFrame(rows,columns=["name", "node","chol", "state",'replicate'] +gcols)
    
    goodnodes = []
    novar = []
    err = []
    for node in df["node"].unique():
        print(node)
        node_df = df.loc[df['node'] == node]
        node_df = node_df.reset_index()

        #doing PCA to avoid correlated variables
        x = node_df.loc[:, gcols].values
        y = node_df.loc[:,['chol']].values
        x = StandardScaler().fit_transform(x)
        pca = PCA()
        principalComponents = pca.fit_transform(x)
        evr = pca.explained_variance_ratio_.cumsum()
        nPCs = 0
        for i,j in enumerate(evr):
            if j > 0.99:
                nPCs = i + 1
                break
        if nPCs == 0:
            print(f'no variance. skipping node {node}')
            novar.append(str(node))
            continue
        pca = PCA(n_components=nPCs)
        principalComponents = pca.fit_transform(x)
        principalDf = pd.DataFrame(data = principalComponents
                , columns = [f'PC{x}' for x in range(1,nPCs+1)])
        finalDf = pd.concat([principalDf, node_df[['chol','name','node']]], axis = 1)

        modelP = build_model(finalDf,node,nPCs)

        if modelP == 0: err.append(str(node))
        elif modelP < 0.05: goodnodes.append(str(node))
        
        print(f'node {node}: P = {modelP}, nPCs = {nPCs}')
    print(f'good nodes: {", ".join(goodnodes)}')
    print(f'no variance nodes: {", ".join(novar)}')
    print(f'err nodes: {", ".join(err)}')

    print(f'none: {1348-len(goodnodes)+len(novar)+len(err)}')

    return



# for each node: is this graphlet vector different from that one?
# build a model for each node. get p value of the model.
    
def build_model(finalDf,node,nPCs):
    columns = [f'PC{x}' for x in range(1,nPCs+1)]
    X = finalDf[columns]
    y = pd.get_dummies(finalDf['chol'])[15] # True means 15 mol%, false means 30.
    
    X_train, X_test, y_train, y_test = train_test_split(X,y , 
                                    random_state=104,  
                                    train_size=0.8,  
                                    shuffle=True) 
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    
    # building the model and fitting the data 
    try:
        log_reg = sm.Logit(y_train, X_train).fit(maxiter=100,) 
    except np.linalg.LinAlgError: # I think this happens when there is perfect correlation or things are constant.
        print(f'node {node}: singular matrix error')
        return 0

    pvalue = log_reg.llr_pvalue # p value for the whole model

    return pvalue





if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    main()