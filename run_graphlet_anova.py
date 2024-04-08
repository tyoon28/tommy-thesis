# not an ANOVA anymore really.
# fit a logistic model to each node and print out p values
from graphlets import *
import warnings
from scipy.stats import f_oneway
import statsmodels.api as sm
from statsmodels.formula.api import ols
from sklearn.metrics import (confusion_matrix,accuracy_score) 
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import roc_curve, roc_auc_score 
from sklearn.model_selection import KFold, cross_val_score





def main(remote=False):
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
                    l = list(map(int,line.split(' ')))
                    rows.append([fn,i,chol,state,r] + l)

    gcols = ['g' + str(graphlet) for graphlet in range(73)]
    df = pd.DataFrame(rows,columns=["name", "node","chol", "state",'replicate'] +gcols)
    
    goodnodes = []
    novar = []
    err = []
    accuracy = {}
    randaccuracy = []
    rocs = {}
    cvs = {}
    for node in tqdm.tqdm(df["node"].unique()):
        node_df = df.loc[df['node'] == node]
        node_df = node_df.reset_index()


        #doing PCA to avoid correlated variables
        x = node_df.loc[:, gcols].values
        y = node_df.loc[:,['chol']].values
        x = StandardScaler().fit_transform(x)
        pca = PCA()
        principalComponents = pca.fit_transform(x)
        evr = pca.explained_variance_ratio_.cumsum()
        print(evr)
        nPCs = 0
        for i,j in enumerate(evr):
            if j > 0.90:
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

        # generate random data
        randDf = pd.DataFrame([np.random.normal(size=nPCs) for x in range(len(node_df))], 
                              columns = [f'PC{x}' for x in range(1,nPCs+1)])
        randDf['chol'] = finalDf['chol']
        # randDf = df.merge(pd.DataFrame(index=np.random.randint(2690608, size=len(finalDf))), left_index=True, right_index=True, how='right')

        acc,roc,cv = build_model(finalDf,node,nPCs)
        rocs[node] = roc
        cvs[node] = cv

        acc_rand,roc, cv = build_model(randDf,node,nPCs)
        accuracy[node] = acc
        randaccuracy.append(acc_rand)
    if not remote:
        print(sorted(accuracy)[-10:])
        print(sorted(rocs)[-10:])
        plt.scatter(accuracy,rocs)
        plt.xlabel('accuracy')
        plt.ylabel('AUC')
        plt.show()
        x = list(range(len(accuracy)))

        plt.plot(x,sorted(accuracy))
        plt.plot(x,sorted(randaccuracy))
        plt.legend()
    # find out whoch nodes have highest accuracy
        
    # output cxc file
    u = mda.Universe('R1-30-closed/R1-0-start-membrane-3JYC.pdb','R1-30-closed/R1-0-1000-3JYC.xtc')
    dd_to_csv(accuracy,'rf_accuracy_allR',u)
    return



# for each node: is this graphlet vector different from that one?
# build a model for each node. get p value of the model.
# look at example node that we know changes a lot
# remove collinear variables
# simulate data as sanity check
# check overfitting
# look at first 10
# check how many PCs
# eigenvalue plot to find ideal number of PCs
# nonlinear predictive models

# dff = pd.DataFrame([[x,2*x,x+15,x-3,np.random.normal()] for x in np.random.normal(size=1996)])
# X = np.concatenate(np.random.normal(size = 1000),np.random.normal(loc=10,size=1000))



def build_model(finalDf,node,nPCs):
    columns = [f'PC{x}' for x in range(1,nPCs+1)]

    X = finalDf[columns]
    y = pd.get_dummies(finalDf['chol'])[15].astype(int) # True means 15 mol%, false means 30.
    
    X_train, X_test, y_train, y_test = train_test_split(X,y , 
                                    random_state=104,  
                                    train_size=0.8,
                                    shuffle=True) 

    scaler = StandardScaler()
    # X_train = scaler.fit_transform(X_train)
    # X_test = scaler.transform(X_test)
    # X_train = sm.add_constant(X_train)
    # X_test = sm.add_constant(X_test)
    X_train = X_train.reset_index(drop=True)
    y_train = y_train.reset_index(drop=True)    


    # building the model and fitting the data 
    try:
        #log_reg = sm.Logit(y_train, X_train).fit_regularized() # nan p values mean they are nothing.

        rf = RandomForestRegressor(n_estimators = 500)
        # Train the model on training data
        rf.fit(X_train, y_train)
        probs = rf.predict(X_test)
        preds = np.round(probs)
        # Calculate the absolute errors
        acc = (y_test == preds).sum() / len(preds)
        TP = ((preds == 1) & (y_test == 1)).sum()
        FP = ((preds == 1) & (y_test == 0)).sum()
        precision = TP / (TP+FP)

        fpr, tpr, thresholds = roc_curve(y_test, probs, pos_label=1)
        roc_auc = roc_auc_score(y_test, probs) 
        roc_auc

        # CV:
        k_folds = KFold(n_splits = 5)

        scores = cross_val_score(rf, X, y, cv = k_folds)


        # plt.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % roc_auc) 
        # # roc curve for tpr = fpr  
        # plt.plot([0, 1], [0, 1], 'k--', label='Random classifier') 
        # plt.xlabel('False Positive Rate') 
        # plt.ylabel('True Positive Rate') 
        # plt.title('ROC Curve') 
        # plt.legend(loc="lower right") 
        # plt.show()



        return acc, roc_auc,scores.mean()

    except np.linalg.LinAlgError: # I think this happens when there is perfect correlation or things are constant.
        print(f'node {node}: singular matrix error')
        return 0
    
    
    
    y_pred = log_reg.predict(X_test)
    preds = list(map(round, y_pred)) 
    return accuracy_score(y_test, preds)



    for i in log_reg.pvalues.dropna()[1:]:
        if i < 0.05:
            return i
    
    pvalue = log_reg.llr_pvalue # p value for the whole model

    return 1
    return pvalue





if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    main()