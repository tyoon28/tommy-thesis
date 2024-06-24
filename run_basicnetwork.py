'''
helper functions for calculating centralities
'''
from graphlets import *

def main():
    centralities()
    
    return


def centralities():
    # get average centralities
    
    firstrun = True
    for r in ['R1','R2','R3']:
        for c in ['15','30']:
            xtcs = []
            for file in os.listdir(f'{r}-{c}-closed'):
                if file.endswith('.xtc'):
                    xtcs.append(f'{r}-{c}-closed/'+file)
            xtcs.sort(key=lambda x: int(x.split('-')[1]))
            u = mda.Universe(f'{r}-{c}-closed/{r}-0-start-membrane-3JYC.pdb',*xtcs,continuous=True)
            G = get_consensus_graph(u)
            bt = nx.betweenness_centrality(G)
            cl = nx.closeness_centrality(G)
            ei = nx.eigenvector_centrality_numpy(G)
            if firstrun:
                betweenness = {k:[bt[k]] for k in bt}
                closeness = {k:[cl[k]] for k in cl}
                eigenvector = {k:[ei[k]] for k in ei}
                firstrun = False
            else:
                for i in bt: 
                    betweenness[i].append(bt[i])
                    closeness[i].append(cl[i])
                    eigenvector[i].append(ei[i])
            
    betweenness = {k:np.mean(betweenness[k]) for k in betweenness}
    closeness = {k:np.mean(closeness[k]) for k in closeness}
    eigenvector = {k:np.mean(eigenvector[k]) for k in eigenvector}

    color_by_centrality(betweenness,'betweenness-avg')
    color_by_centrality(closeness,'closeness-avg')
    color_by_centrality(eigenvector,'eigenvector-avg')

    # write file
    
    with open(f'centralities.csv', 'w', newline='') as csvfile:
        
        wr = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        sh = wr.writerow(['subunit','residue ID','MD simulation ID','residue name','node PCA distance','betweenness','closeness','eigenvector'])
        for i in betweenness:
            resname = u.atoms.select_atoms(f'resid {int(i)}').residues.resnames[0]
            s = MD_to_chimerax(i)
            subunit = s[3]
            residue = s[5:]
            btw = betweenness[i]
            clo = closeness[i]
            eig = eigenvector[i]
            
            mdsimid = i
            ind+=1
    
            sh = wr.writerow([subunit,residue,mdsimid,resname, btw,clo,eig])
    return

def gephi_networks():
    r = 'R1'
    c = '15'
    xtcs = []
    for file in os.listdir(f'{r}-{c}-closed'):
        if file.endswith('.xtc'):
            xtcs.append(f'{r}-{c}-closed/'+file)
    xtcs.sort(key=lambda x: int(x.split('-')[1]))
    u = mda.Universe(f'{r}-{c}-closed/{r}-0-start-membrane-3JYC.pdb',*xtcs,continuous=True)
    for t in [0.8,0.7]:
        na = f'closed_15_{t*100}'
        viz_consensus_graph(u,threshold=t,outname = na)


def graphlet_windowgraph_figure():
    threshold = 0.8
    res = u.select_atoms('not resname CHOL and not resname POPC')
    lenr = len(res.residues)
    edgesmat = np.zeros(shape=(lenr,lenr))

    for ts in u.trajectory[50000:50000+500]:
        frame = u.trajectory.frame
        r = res.atoms.center_of_mass(compound='residues')
        mat = distances.contact_matrix(r, cutoff=6)
        np.fill_diagonal(mat, 0)
        edgesmat += mat

    
    t = len(u.trajectory[50000:50000+500]) * threshold
    G = nx.from_numpy_array((edgesmat >= t))
    nx.write_gexf(G,f'windowgraph.gexf')


if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    main()