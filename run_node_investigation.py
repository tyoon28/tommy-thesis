'''
This analysis was a little bit too contrived to be useful, not included in thesis
'''
from graphlets import *
import multiprocessing as mp
from itertools import repeat
from functools import partial




def func(d,cond):
    residues = [265,264,178,179,182] # DON"T RUN WITHOUT SETTING RESIDEUSE
    resnames = ['V','L','G','A','A']

    i,r = cond

    print(f'doing {r} {i}')
    xtcs = []
    for file in os.listdir(f'{r}-{i}-closed'):
        if file.endswith('.xtc'):
            xtcs.append(f'{r}-{i}-closed/'+file)
    xtcs.sort(key=lambda x: int(x.split('-')[1]))
    u = mda.Universe(f'{r}-{i}-closed/{r}-0-start-membrane-3JYC.pdb',*xtcs,continuous=True)

    protein = u.select_atoms('not resname CHOL and not resname POPC')
    for res in (residues):
        resi = resid_to_md_subunits(res)
        print(f'{r} {i} chol',u.select_atoms(f'resid {resi[0]}').residues[0],'...')
        # record frequency of each unique set of contacts

        for ts in tqdm.tqdm(u.trajectory):
            frame = u.trajectory.frame
            r_compound = protein.atoms.center_of_mass(compound='residues')
            mat = distances.contact_matrix(r_compound, cutoff=6)
            np.fill_diagonal(mat, 0) 
            # +/- 1 from resid for proper indexing
            for inc in range(4):
                fs = frozenset(np.nonzero(mat[resi-1][inc])[0])
                if fs in d[i][res].keys():
                    d[i][res][fs] += 1 
                else: d[i][res][fs] = 1 

    return d


def resid_to_md_subunits(r):
    return np.array([(i * 337) + (r - 35) for i in range(0,4)])


def main():
    runun = False
    d = {} # key: residue; item: dict of concentrations with dict of contacts and frequencies
     
    residues = [265,264,178,179,182] # DON"T RUN WITHOUT SETTING RESIDEUSE
    resnames = ['V','L','G','A','A']



    print('starting')
    if not runun:
        print('making d')
        conditions = [(i,r) for i in ['30','15'] for r in ['R1','R2','R3']]
        d = {}
        d['15'] = {res:{} for res in residues}
        d['30'] = {res:{} for res in residues}

        for c in conditions:
            d = func(d,c)

        with open('node_invest.pickle', 'wb') as f:     
            pickle.dump(d, f)
    else:
        with open('/Users/Tommy/Desktop/thesis/tommy-thesis/node_invest.pickle', 'rb') as f:
            d = pickle.load(f)


    dd = convert_d(d)
    d_difference = {}
    for res in residues:
        d_difference[res] = {}

        s = set(dd['15'][res].keys()).union(set(dd['30'][res].keys()))

        for fs in s:
            if fs in dd['30'][res].keys() and fs in dd['15'][res].keys():
                d_difference[res][fs] = dd['30'][res][fs] - dd['15'][res][fs]
            elif fs in dd['30'][res].keys():
                d_difference[res][fs] = dd['30'][res][fs]
            elif fs in dd['15'][res].keys():
                d_difference[res][fs] = -dd['15'][res][fs]
            else:
                print('problem')
    for k in d_diff_filtered:
        print(k,min(d_diff_filtered[k].values()))

    d_diff_filtered = {}
    for res in residues:
        d_diff_filtered[res] = {}
        for i in d_difference[res]:
            if abs(d_difference[res][i]) > 40000:
                d_diff_filtered[res][i] = d_difference[res][i]
        d_diff_filtered[res] = {k: v for k, v in sorted(d_diff_filtered[res].items(), key=lambda item: item[1])}





    print('making fig')
    fig, axs = plt.subplots(2,5,figsize=(16,9), gridspec_kw={'height_ratios': [1, 1]})
    for p in [1,2,3,4,5]:
        res = residues[p-1]
        for s in [0,1]:
            if s == 0:
                plt.title(f'residue {resnames[p-1]}{res}')
                if max(d_diff_filtered[res].values()) < 0: continue
            elif min(d_diff_filtered[res].values()) > 0: continue
            ax = axs[s,p-1]
            plt.sca(ax)
            hedges = {}
            hp = {}
            key = 0
            for k,v in sorted(d_diff_filtered[res].items(), key=lambda item: len(item[0])):
                key += 1
                if (s == 0 and v > 0) or (s == 1 and v < 0):
                    hedges[key] = frozenset(map(lambda x: MD_to_resid(int(x),u),k))
                    hp[key] = {'weight':v/300000}

            
            H = hnx.Hypergraph(hedges,edge_properties = hp)
            color = [H.edge_properties[e]['weight'] for e in H.edges()]
            norm = plt.Normalize(-1, 1)
            color = cm.bwr(norm(color))
            alpha = 0.9
            # hnx.draw(H)
            hnx.drawing.draw(H, edges_kwargs={'facecolors': color*(1, 1, 1, alpha),'edgecolors': 'black','linewidths': 2},with_edge_labels=False)


    # ONE FROM EACH
    d_diff_filtered = {}
    for res in residues:
        d_diff_filtered[res] = {}
        ma=0
        mi = 1000000
        for i in d_difference[res]:
            ma = max(d_difference[res], key=d_difference[res].get)
            mi = min(d_difference[res], key=d_difference[res].get)
            d_diff_filtered[res][ma] = d_difference[res][ma]
            d_diff_filtered[res][mi] = d_difference[res][mi]

    fig, axs = plt.subplots(1,5,figsize=(16,9))
    for p in [1,2,3,4,5]:
        res = residues[p-1]
        s=0
        ax = axs[p-1]
        plt.sca(ax)
        plt.title(f'residue {resnames[p-1]}{res}')

        hedges = {}
        hp = {}
        key = 0
        for k,v in sorted(d_diff_filtered[res].items(), key=lambda item: len(item[0]),reverse=True):
            key += 1
            hedges[key] = frozenset(map(lambda x: int(x)+1,k))
            hp[key] = {'weight':v/300000}
        
        H = hnx.Hypergraph(hedges,edge_properties = hp)
        color = [H.edge_properties[e]['weight'] for e in H.edges()]
        norm = plt.Normalize(-1, 1)
        color = cm.bwr(norm(color))
        alpha = 0.9
        # hnx.draw(H)
        hnx.drawing.draw(H, edges_kwargs={'facecolors': color*(1, 1, 1, alpha),'edgecolors': 'black','linewidths': 2},with_edge_labels=False)

        
    # ax = axs[1,0]
    # plt.sca(ax)
    # gradient = np.linspace(0, 1, 256)
    # gradient = np.vstack((gradient, gradient))
    # plt.imshow(gradient, aspect='auto', cmap=cm.bwr)
    # ax.text(-0.01, 0.5, '15',ha='right', va='center',fontsize=10,transform=ax.transAxes)
    # ax.text(1.01, 0.5, '30', ha='left', va='center',fontsize=10,transform=ax.transAxes)

    # for ax in axs[1]:
    #     plt.sca(ax)
    #     plt.xticks([])
    #     plt.yticks([])
    # fig.delaxes(axs[1,1])
    # fig.delaxes(axs[1,2])

    plt.savefig('hgraphs.png')
    print('done')

def convert_d(d):
    dd={}
    dd['15'] = {res:{} for res in residues}
    dd['30'] = {res:{} for res in residues}
    for c in d:
        for res in d[c]:
            for k in d[c][res]:
                a = [MD_to_resid(p) for p in k]
                fs = frozenset(a)
                if '88' in a: 
                    print(k)
                    break
                if fs in dd[c][res]:
                    dd[c][res][fs] += d[c][res][k]
                else:
                    dd[c][res][fs] = d[c][res][k]
    return dd

def uh():
    r = 'R1'
    i ='15'
    xtcs = []
    for file in os.listdir(f'{r}-{i}-closed'):
        if file.endswith('.xtc'):
            xtcs.append(f'{r}-{i}-closed/'+file)
    xtcs.sort(key=lambda x: int(x.split('-')[1]))
    u = mda.Universe(f'{r}-{i}-closed/{r}-0-start-membrane-3JYC.pdb',*xtcs,continuous=True)

    protein = u.select_atoms('not resname CHOL and not resname POPC')
    res = 182
    resi = resid_to_md_subunits(res)

    targ = resid_to_md_subunits(185)
    print(f'{r} {i} chol',u.select_atoms(f'resid {resi[0]}').residues[0],'...')
    # record frequency of each unique set of contacts
    k=0
    for ts in tqdm.tqdm(u.trajectory):
        frame = u.trajectory.frame
        r_compound = protein.atoms.center_of_mass(compound='residues')
        mat = distances.contact_matrix(r_compound, cutoff=6)
        np.fill_diagonal(mat, 0) 
        if np.any(mat[resi-1,:][:,targ-1][0]):
            print(frame)
            break
        

        # +/- 1 from resid for proper indexing
        for inc in range(4):
            fs = frozenset(np.nonzero(mat[resi-1][inc])[0])



# just calculate pairwise changes in contact probability...
if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    main()