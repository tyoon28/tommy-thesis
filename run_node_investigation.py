from graphlets import *
import multiprocessing as mp
from itertools import repeat


def func(package):
    d = package[0]
    residues = [265,264,178,179,182] # DON"T RUN WITHOUT SETTING RESIDEUSE
    resnames = ['V','L','G','A','A']

    i,r = package[1]

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

        for ts in u.trajectory[:20]:
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
    return


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
        with mp.Manager() as manager:
            d = manager.dict()
            d['15'] = {res:{} for res in residues}
            d['30'] = {res:{} for res in residues}
        
            with manager.Pool(4) as pool:
                s = pool.map(func, zip(repeat(d,len(conditions)),conditions))
            d = dict(d)

        with open('node_invest.pickle', 'ab') as f:     
            pickle.dump(d, f)
    else:
        with open('node_invest.pickle', 'rb') as f:
            d = pickle.load(f)

    print(d)
    return
    d_difference = {}
    for res in residues:
        d_difference[res] = {}

        s = set(d['15'][res].keys()).union(set(d['30'][res].keys()))

        for fs in s:
            if fs in d['30'][res].keys() and fs in d['15'][res].keys():
                d_difference[res][fs] = d['30'][res][fs] - d['15'][res][fs]
            elif fs in d['30'][res].keys():
                d_difference[res][fs] = d['30'][res][fs]
            elif fs in d['15'][res].keys():
                d_difference[res][fs] = -d['15'][res][fs]
            else:
                print('problem')
    d_diff_filtered = {}
    for res in residues:
        d_diff_filtered[res] = {}
        for i in d_difference[res]:
            if abs(d_difference[res][i]) > 60000:
                d_diff_filtered[res][i] = d_difference[res][i]
        d_diff_filtered[res] = {k: v for k, v in sorted(d_diff_filtered[res].items(), key=lambda item: item[1])}



    print('making fig')
    fig, axs = plt.subplots(2,3,figsize=(16,9), gridspec_kw={'height_ratios': [1, 0.05]})

    # plt1 = fig.add_subplot(1,3,1)
    # plt2 = fig.add_subplot(1,3,2)
    # plt3 = fig.add_subplot(1,3,3)
    # cbar = fig.add_subplot(2,3,4,figsize=(6.4, 0.35))
    # x = list(range(len(d_difference[res])))
    # y = sorted(d_difference[res].values(),reverse=True)
    # plt.plot(x,y)

    for p in [1,2,3,4,5,6]:
        res = residues[p-1]
        ax = axs[p//3,(p-1)%3]
        plt.sca(ax)
        hedges = {}
        hp = {}
        key = 0
        for k,v in sorted(d_diff_filtered[res].items(), key=lambda item: item[1]):
            key += 1
            hedges[key] = frozenset(map(lambda x: MD_to_resid(x,u),k))
            hp[key] = {'weight':v/600000}
        H = hnx.Hypergraph(hedges,edge_properties = hp)
        color = [H.edge_properties[e]['weight'] for e in H.edges()]
        norm = plt.Normalize(-1, 1)
        color = cm.bwr(norm(color))
        alpha = 0.9
        # hnx.draw(H)
        hnx.drawing.draw(H,
            edges_kwargs={
                'facecolors': color*(1, 1, 1, alpha),
                'edgecolors': 'black',
                'linewidths': 2
            },
            with_edge_labels=False
        )
        plt.title(f'residue {resnames[p-1]}{res}')
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

if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    main()