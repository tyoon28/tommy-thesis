from graphlets import *
from itertools import combinations
import seaborn as sns
import time


def main():
    # inefficient to do all of these in series but who cares.
    # print('calculating persistence in RIN')
    # mid_all = persistance()
    # np.save('mid_all.npy', mid_all)    # .npy extension is added if not given

    mid_all = np.load('mid_all.npy')
    print('calculating correlations between contacts')
    correlation(mid_all,read=True)

    # print('calculating cholesterol contact landscape')
    # chol_contact()

    # print('calculating cholesterol contact landscape - binding sites')
    # chol_interactionlength()

    # print('calculating centralities')
    # centralities()

    # print('cholesterol threshold graph')
    # chol_thresh()


def correlation(mid_all,read=False):
    # perhaps change this to only middle stable contacts.
    # i'm only doing this for 30 mol%.
    
    c = '30'
    xtcs = []
    for file in os.listdir(f'R1-30-closed'):
        if file.endswith('.xtc'):
            xtcs.append(f'R1-30-closed/'+file)
    xtcs.sort(key=lambda x: int(x.split('-')[1]))
    u = mda.Universe('R1-30-closed/R1-0-start-membrane-3JYC.pdb',*xtcs,continuous=True)

    lenwin = len(u.trajectory)

    # store all contacts in mm
    nedges = len(mid_all[0])
    comp = np.zeros((len(u.trajectory[:lenwin])*3,nedges),dtype=bool) # comparison matrix for all replicates

    if not read:
        i = 0
        for r in ['R1','R2','R3']:
            xtcs = []
            for file in os.listdir(f'{r}-{c}-closed'):
                if file.endswith('.xtc'):
                    xtcs.append(f'{r}-{c}-closed/'+file)
            xtcs.sort(key=lambda x: int(x.split('-')[1]))
            u = mda.Universe('R1-30-closed/R1-0-start-membrane-3JYC.pdb',*xtcs,continuous=True)

            contact_threshold = 6
            sterols = u.select_atoms(f'(resname CHOL and not (name RC1 or name RC2))').residues
            protein = u.select_atoms('not resname CHOL and not resname POPC').residues

            for ts in tqdm.tqdm(u.trajectory[:lenwin]):
                frame = u.trajectory.frame
                r = protein.atoms.center_of_mass(compound='residues')
                mat = distances.contact_matrix(r, cutoff=6)
                np.fill_diagonal(mat, 0)
                # only look at middling edges.
                comp[i] = mat[mid_all[0],mid_all[1]]
                i += 1
        
        df = pd.DataFrame(comp)
        df.to_csv('biggo.csv')
    else:   
        t = time.time()
        corrcoef_xy = np.load('cormat-30.npy')

    print('startgraph')
    G = nx.from_numpy_array(corrcoef_xy)
    selected_nodes = [n for n,v in G.nodes(data=True) if v['weight'] > 0.25]
    H = G.subgraph(selected_nodes) 
    nx.write_gexf(H,'correlationgraph.gexf')
    print('donegraph')

    print('starclus')
    p = (corrcoef_xy > 0.25)
    pp = np.where(p.any(axis=1))[0]
    plot = sns.clustermap(corrcoef_xy[pp,:][:,pp])
    fig = plot.get_figure()
    fig.savefig("clustermap-30.png")

    print('startplot')
    plot = sns.heatmap(corrcoef_xy, annot=True)
    fig = plot.get_figure()
    fig.savefig("cormat-30.png")
    print('doneplot')

    plt.clf()


    # print('done making df. calculating cormat')
    # a =df.sample(100).corr(method='pearson')
    # print(f'100 done in {time.time()-t}')
    # a =df.sample(1000).corr(method='pearson')
    # print(f'1000 done in {time.time()-t}')
    # a =df.sample(10000).corr(method='pearson')
    # print(f'10000 done in {time.time()-t}')
    # a =df.sample(100000).corr(method='pearson')
    # print(f'100000 done in {time.time()-t}')

    return
    

    cor = df.corr(method='pearson')\
    
    cor.to_csv('cormat-30.csv',index=False)
    plot = sns.heatmap(cor, annot=True)
    fig = plot.get_figure()
    fig.savefig("cormat-30.png")

    plt.clf()

    G = nx.from_pandas_adjacency(cor)
    selected_nodes = [n for n,v in G.nodes(data=True) if v['weight'] > 0.5]
    H = G.subgraph(selected_nodes) 
    nx.write_gexf(H,'correlationgraph.gexf')


    

    return


def correlation2(mid_all):
    # perhaps change this to only middle stable contacts.
    # i'm only doing this for 30 mol%.
    c = '30'
    

    xtcs = []
    for file in os.listdir(f'R1-30-closed'):
        if file.endswith('.xtc'):
            xtcs.append(f'R1-30-closed/'+file)
    xtcs.sort(key=lambda x: int(x.split('-')[1]))
    u = mda.Universe('R1-30-closed/R1-0-start-membrane-3JYC.pdb',*xtcs,continuous=True)

    lenwin = len(u.trajectory)

    # store all contacts in mm
    nedges = len(mid_all[0])
    comp = np.zeros((len(u.trajectory[:lenwin]),nedges),dtype=bool) # comparison matrix for 1 replicate
    store = np.zeros((len(u.trajectory[:lenwin]),nedges),dtype=bool) # average of all replicates
    for r in ['R1','R2','R3']:
        xtcs = []
        for file in os.listdir(f'{r}-{c}-closed'):
            if file.endswith('.xtc'):
                xtcs.append(f'{r}-{c}-closed/'+file)
        xtcs.sort(key=lambda x: int(x.split('-')[1]))
        u = mda.Universe('R1-30-closed/R1-0-start-membrane-3JYC.pdb',*xtcs,continuous=True)

        contact_threshold = 6
        sterols = u.select_atoms(f'(resname CHOL and not (name RC1 or name RC2))').residues
        protein = u.select_atoms('not resname CHOL and not resname POPC').residues

        i = 0
        for ts in tqdm.tqdm(u.trajectory[:lenwin]):
            frame = u.trajectory.frame
            r = protein.atoms.center_of_mass(compound='residues')
            mat = distances.contact_matrix(r, cutoff=6)
            np.fill_diagonal(mat, 0)
            # only look at middling edges.
            comp[i] = mat[mid_all]
            i += 1
        store += comp
    
    df = pd.DataFrame(store/3)
    cor = df.corr(method='pearson')
    cor.to_csv('cormat-30.csv',index=False)
    plot = sns.heatmap(cor, annot=True)
    fig = plot.get_figure()
    fig.savefig("cormat-30.png")

    plt.clf()

    G = nx.from_pandas_adjacency(cor)
    selected_nodes = [n for n,v in G.nodes(data=True) if v['weight'] > 0.5]
    H = G.subgraph(selected_nodes) 
    nx.write_gexf(H,'correlationgraph.gexf')


    

    return


def persistance():
    # did i spell it right?
    
    for c in ['15','30']:
        if c == '15': continue # decided this fig is only for 30%
        store = np.zeros((1348,1348))
        for r in ['R1','R2','R3']:
            contact_threshold = 6
            xtcs = []
            for file in os.listdir(f'{r}-{c}-closed'):
                if file.endswith('.xtc'):
                    xtcs.append(f'{r}-{c}-closed/'+file)
            xtcs.sort(key=lambda x: int(x.split('-')[1]))
            u = mda.Universe(f'{r}-{c}-closed/{r}-0-start-membrane-3JYC.pdb',*xtcs,continuous=True)
            sterols = u.select_atoms(f'(resname CHOL and not (name RC1 or name RC2))').residues
            protein = u.select_atoms('not resname CHOL and not resname POPC').residues
            lenp = len(protein)
            for ts in tqdm.tqdm(u.trajectory):
                r = protein.atoms.center_of_mass(compound='residues')
                mat = distances.contact_matrix(r, cutoff=6)
                np.fill_diagonal(mat, 0)
                store += mat
        
        store = store / (len(u.trajectory) * 3)
        # get all edges in a 1D np array
        no0 = store[np.where(store!=0)]
        if c == '30':
            # mid_all stores the edges that appear between 10 and 90 percent of the time.
            mid_all = (np.where(abs(abs(store)-0.5) < 0.4))
        mid = no0[np.where(no0 > 0.1)]
        mid = mid[np.where(mid < 0.9)]

        plt.hist(no0)
        print(f'persistence {c} mol%: {len(mid)}')
        print(f'total coltacts: {len(no0)}')
        print(f'middle percent total contact time in {c} mol%: {np.sum(mid)/np.sum(no0)}')
        

    
    plt.xlabel('average % persistance')
    plt.ylabel('count')
    # plt.legend(['15 mol%','30 mol%'])

    plt.savefig('figure_edgepersistence.png')
    
    return mid_all



def chol_contact():
    contact_threshold = 6
    u = mda.Universe('R1-30-closed/R1-0-start-membrane-3JYC.pdb','R1-30-closed/R1-0-1000-3JYC.xtc')
    for c in ['15','30']:
        if c == '15': 
            color = 'cyan'
            darkcolor = 'blue'
        else: 
            color = 'yellow'
            darkcolor = 'red'
        avg = np.zeros(len(u.trajectory))
        for r in ['R1','R2','R3']:
            xtcs = []
            for file in os.listdir(f'{r}-{c}-closed'):
                if file.endswith('.xtc'):
                    xtcs.append(f'{r}-{c}-closed/'+file)
            xtcs.sort(key=lambda x: int(x.split('-')[1]))
            u = mda.Universe(f'{r}-{c}-closed/{r}-0-start-membrane-3JYC.pdb',*xtcs,continuous= True)

            protein = u.select_atoms('not resname CHOL and not resname POPC').residues
            chol = u.select_atoms(f'(resname CHOL and not (name RC1 or name RC2))').residues

            protein_sterols = protein + chol

            t = np.zeros(len(u.trajectory))
            y = np.zeros(len(u.trajectory))
            y2=np.zeros(len(u.trajectory))
            # count cholesterols occupying
            # have lines in background, thick line for average of the replicates, and both conditions in same plot.
            for ts in tqdm.tqdm(u.trajectory[:]):
                tim = u.trajectory.time
                t[u.trajectory.frame] = tim


                r = protein_sterols.atoms.center_of_mass(compound='residues')
                contact_mat = distances.contact_matrix(r, cutoff=contact_threshold)

                # take only items representing chol-protein contacts
                z = contact_mat[:len(protein),len(protein):]
                occupancy = np.sum(np.any(z,0))
                occupancy2 = np.sum(np.any(z,1))
                y[u.trajectory.frame] = occupancy
                y2[u.trajectory.frame] = occupancy2

            y2_subsample = y2[::100]
            t_subsample = t[::100]
            
            plt.plot(t,y,color)
            avg += y
        avg = avg / 3
        plt.plot(t,avg,darkcolor)
        print(f'{c} mol%: average {np.mean(avg)}')
    plt.savefig('figure_cholcontact.png', dpi=300, bbox_inches='tight')


def chol_interactionlength():
    # cholesterol contacts tend to be short, with most lasitng X ns. longest contact is X ns.
    # lengths of interactions of binding sites vs with rest of protein
    data = []
    for c in ['15','30']:
        for r in ['R1','R2','R3']:
            xtcs = []
            for file in os.listdir(f'{r}-{c}-closed'):
                if file.endswith('.xtc'):
                    xtcs.append(f'{r}-{c}-closed/'+file)
            xtcs.sort(key=lambda x: int(x.split('-')[1]))
            u = mda.Universe('R1-30-closed/R1-0-start-membrane-3JYC.pdb',*xtcs,continuous=True)


            rog = Cholesterol_contact(u)
            rog.run(start=0,verbose=True)

            all_interactions = [] # list of durations of all interactions
            bs_vs_rest = [[],[]] # 0: durations of interactions that touch a binding site; 1: others
            for ch in rog.results:
                for i in range(len(rog.results[ch]['binding_events'])):
                    start,end = rog.results[ch]['binding_events'][i]
                    dur = end-start

                    bound = rog.results[ch]['binding_events_actual'][i]
                    if bound: ind = 0
                    else: ind = 1
                    bs_vs_rest[ind].append(dur)
                    all_interactions.append(dur)
            
        avg_bs = np.mean(bs_vs_rest[0])
        avg_rest = np.mean(bs_vs_rest[1])

        data.append([avg_bs,avg_rest])


    # grouped bar chart with binding site/non bs and chol concentration.
    index = ['15 mol%', '30 mol%']
    df = pd.DataFrame({'Binding site': data[0], 'Non-binding site': data[1]}, index=index)
    ax = df.plot.bar(rot=0)
    plt.xlabel('Cholesterol concentration')
    plt.ylabel('Average Interaction length (frames (change this to ns))')
    plt.savefig('figure_cholesterol_interactionlength', dpi=300, bbox_inches='tight')
    plt.clear()


    plt.hist([bs_vs_rest[0],bs_vs_rest[1]], stacked=True, density=True)
    plt.legend(['Binding sites','All other residues'])
    plt.savefig('figure_cholesterol_interactionlength2', dpi=300, bbox_inches='tight')
    plt.clear()

    return



        # plot with color corresponding to replicate/cholesterol condition

    # A) a histogram?
    # B) a violinish plot with 2 x values - binding site and not        

def graphs():

    xtcs = []
    for file in os.listdir('R1-30-closed'):
        if file.endswith('.xtc'):
            xtcs.append('R1-30-closed/'+file)
    xtcs.sort(key=lambda x: int(x.split('-')[1]))
    u = mda.Universe('R1-30-closed/R1-0-start-membrane-3JYC.pdb',*xtcs,continuous= True)

    viz_consensus_graph(u,start=0,end=None,threshold = 0.9,outname='',**kwargs)

    # generate contact graphs with increasingly smaller resolutions, also generate colors for protein
    pass

def triangle():
    pass

def centralities():
    # take average centrality.
    betweenesses = {}
    closenesses = {}
    eigenvectors = {}
    i = 0
    for r in ['R1','R2','R3']:
        for c in ['15','30']:
            xtcs = []
            for file in os.listdir(f'{r}-{c}-closed'):
                if file.endswith('.xtc'):
                    xtcs.append(f'{r}-{c}-closed/'+file)
            xtcs.sort(key=lambda x: int(x.split('-')[1]))
            u = mda.Universe(f'{r}-{c}-closed/{r}-0-start-membrane-3JYC.pdb',*xtcs,continuous=True)
            G = get_consensus_graph()
            betweenness = nx.betweenness_centrality(G)
            closeness = nx.closeness_centrality(G)
            eigenvector = nx.eigenvector_centrality_numpy(G)

            if i == 0:
                betweenesses = betweenness
                closenesses = closeness
                eigenvectors = eigenvector
            else:
                for d,s in [(betweenesses,betweenness),(closenesses,closenesses),(eigenvectors,eigenvector)]:
                    for item in d:
                        d[item] += s[item]
            i += 1
    
    for d in [betweenesses,closenesses,eigenvectors]:
        for item in d:
            d[item] = d[item] / i    
        

    color_by_centrality(betweenness,'betweenness-all')
    color_by_centrality(closeness,'closeness-all')
    color_by_centrality(eigenvector,'eigenvector-all')

    dd_to_csv(betweenness,'betweenness',u)
    dd_to_csv(closeness,'closeness',u)
    dd_to_csv(eigenvector,'eigenvector',u)
    return

def chol_thresh():
    data = []
    rog = Cholesterol_contact(u)
    rog.run(start=0,verbose=True)
    for t in np.linspace(0,1,30):
        pass



if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    main()