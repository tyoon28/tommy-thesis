from hypergraphs import *
from network_statistics import *
import pandas as pd
from sklearn.tree import DecisionTreeClassifier # Import Decision Tree Classifier
from sklearn.model_selection import train_test_split # Import train_test_split function
from sklearn import metrics #Import scikit-learn metrics module for accuracy calculation
from sklearn import tree



class uhhgraphs(AnalysisBase):
    '''
    atomgroup: selection of cholesterol(s) and protein or all atoms
    results[id] is a dictionary containing:
            contacts - list of frames that contact the channel
            contact_time - total contact time
            longest_contact - duration (frames) of longest contact with channel
            binding_events - ranges of continuous interaction
            mostcontacts - most residues that were touched in one interaction
            residuecontacts - # residue contacts in each interaction
            need to add start of longest contact.
            
    '''

    def __init__(self, atomgroup, **kwargs):
        super(uhhgraphs, self).__init__(atomgroup.universe.trajectory,
                                          **kwargs)
        self.ag = atomgroup

    def _prepare(self):
        # OPTIONAL
        # Called before iteration on the trajectory has begun.
        # Data structures can be set up at this time
        # get all cholesterol ids
        self.sterols = self.ag.select_atoms(f'(resname CHOL and not (name RC1 or name RC2))').residues
        self.protein = self.ag.select_atoms('not resname CHOL and not resname POPC').residues
        self.lenprotein = len(self.protein)
        self.lensterols = len(self.sterols)


        self.sterol_residuecontacts = {i:0 for i in self.sterols.resids}
        self.currentinteractions = {i:[] for i in self.sterols.resids}


        self.protein_sterols =  self.protein + self.sterols
        self.contact_threshold = 6

        self.pairwise = pd.DataFrame()
        self.hypergraph = pd.DataFrame()
        self.hypergraph_dict = {}
        self.pairwise_dict = {}
        
        

    def _single_frame(self):
        # Called after the trajectory is moved onto each new frame.

        res = self.protein.atoms.center_of_mass(compound='residues')
        contact_mat = distances.contact_matrix(res, cutoff=self.contact_threshold)
        HG = self._hypergraph_1_frame(contact_mat)
        G = np.transpose(np.triu(contact_mat).nonzero())
        f = self.ag.universe.trajectory.frame

        self.pairwise_dict[f] = [set(i) for i in G]

        self.hypergraph_dict[f] = []
        for k in HG:
            for edge in HG[k]:
                self.hypergraph_dict[f].append(set(edge))

        


    def _conclude(self):
        #remove self loops
        pdp = pd.DataFrame.from_dict(self.pairwise_dict, orient='index')
        pdh = pd.DataFrame.from_dict(self.hypergraph_dict, orient='index')

        p = pd.get_dummies(pdp.map(repr), prefix='', prefix_sep='')
        h = pd.get_dummies(pdh.map(repr), prefix='', prefix_sep='')

        self.pairwise = p.T.groupby(p.columns).max().T
        self.hypergraph = h.T.groupby(h.columns).max().T

        self.pairwise = self.pairwise.filter(regex=',',axis=1)
        self.hypergraph = self.hypergraph.filter(regex=',',axis=1)

        self.results['hypergraph'] = self.hypergraph
        self.results['pairwise'] = self.pairwise
        pass

    def _hypergraph_1_frame(self,mat):
    #infer hypergraph from 1 frame. return list of hyperedges

        NITERS = 1000

        G = gt.Graph(directed=False)
        G.add_edge_list(np.transpose(mat.nonzero()))

        ntimes = 1

        # from https://graph-tool.skewed.de/static/doc/autosummary/graph_tool.inference.CliqueState.html
        ## WHAT IS THIS
        state = gt.inference.CliqueState(G)
        output = state.mcmc_sweep(niter=NITERS)

        cliques = {} # indexed by k
        # iterate through factor graph
        for v in state.f.vertices():
            if state.is_fac[v]:
                continue # skip over factors

            if state.x[v] > 0: # activated hyperedge
                k = len(state.c[v])
                if k not in cliques:
                    cliques[k] = []
                cliques[k].append(state.c[v])
        
        return cliques





def main():
    molpcts = ['15','30']

    universes = [uone,utwo]
    print(uone)
    print(utwo)

    pairwise = pd.DataFrame()
    hypergraphs = pd.DataFrame()
    ind = 0
    for u in universes:
        graphs = uhhgraphs(u)
        graphs.run(start=0,stop=3,verbose=True)
        graphs.results['pairwise']['molpct'] = molpcts[ind]
        graphs.results['hypergraph']['molpct'] = molpcts[ind]

        pairwise = graphs.results['pairwise']
        hypergraphs = graphs.results['hypergraph']

        ind += 1

    
    # making decision tree for pairwise
    X = pairwise.drop('molpct',axis=1) # Features
    y = pairwise.molpct # Target variable

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=1) # 70% training and 30% test

    # Create Decision Tree classifer object
    # TODO: only pre-pruned for now
    # TODO: for training, sample randomly from simulation
    clf = DecisionTreeClassifier(criterion="entropy", max_depth=10)
    # Train Decision Tree Classifer
    clf = clf.fit(X_train,y_train)
    #Predict the response for test dataset
    y_pred = clf.predict(X_test)

    with open('out.txt','w') as f:
        print("Accuracy, pairwise:",metrics.accuracy_score(y_test, y_pred),file=f)

    # again for hypergraph

     # making decision tree for pairwise
    X2 = hypergraphs.drop('molpct',axis=1) # Features
    y2 = hypergraphs.molpct # Target variable
    X_train2, X_test2, y_train2, y_test2 = train_test_split(X2, y2, test_size=0.3, random_state=1) # 70% training and 30% test
    clf2 = DecisionTreeClassifier(criterion="entropy", max_depth=10)
    # Train Decision Tree Classifer
    clf2 = clf2.fit(X_train2,y_train2)
    #Predict the response for test dataset
    y_pred2 = clf2.predict(X_test2)


    with open('out.txt','w') as f:
        print("Accuracy, hypergraph:",metrics.accuracy_score(y_test2, y_pred2))


    tree.plot_tree(clf2,feature_names=list(X2.columns),filled=True,class_names=['15','30'])
    plt.savefig('hypergraph_tree.png')

    tree.plot_tree(clf,feature_names=list(X.columns),filled=True,class_names=['15','30'])
    plt.savefig('pairwise_tree.png')











        
        
if __name__ == "__main__":
    warnings.filterwarnings('ignore')

    xtcs = []
    for file in os.listdir('R1-15-closed'):
        if file.endswith('.xtc'):
            xtcs.append('R1-15-closed/'+file)
    xtcs.sort(key=lambda x: int(x.split('-')[1]))
    uone = mda.Universe('R1-15-closed/R1-0-start-membrane-3JYC.pdb',*xtcs)

    xtcs = []
    for file in os.listdir('R1-30-closed'):
        if file.endswith('.xtc'):
            xtcs.append('R1-30-closed/'+file)
    xtcs.sort(key=lambda x: int(x.split('-')[1]))
    utwo = mda.Universe('R1-30-closed/R1-0-start-membrane-3JYC.pdb',*xtcs)

    main()