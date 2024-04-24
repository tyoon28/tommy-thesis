from graphlets import *
import seaborn as sns
from scipy.spatial.distance import cdist


# calculate average pairwise absolute cosine distance for all 8 I177 and M181 in protein.

def resid_to_md_subunits(r):
    return np.array([(i * 337) + (r - 35) for i in range(0,4)])

def main():
    for i in ['15','30']:
        residues = [177,181]    
        selstring = 'resid '
        for res in residues:
            a = resid_to_md_subunits(res)
            selstring += ' or resid '.join(list(map(str,a)))
            selstring += ' or resid '
        selstring = selstring[:-10]

        angles = []
        iters = 0
        for r in ['R1','R2','R3']:
            xtcs = []
            for file in os.listdir(f'{r}-{i}-closed'):
                if file.endswith('.xtc'):
                    xtcs.append(f'{r}-{i}-closed/'+file)
            xtcs.sort(key=lambda x: int(x.split('-')[1]))
            u = mda.Universe(f'{r}-{i}-closed/{r}-0-start-membrane-3JYC.pdb',*xtcs,continuous=True)
            print(i,r,len(u.atoms))
    #         protein = u.select_atoms(selstring)
    #         for ts in tqdm.tqdm(u.trajectory[500:1000]):
    #             frame = u.trajectory.frame
    #             poss = protein.positions
    #             vectors = np.zeros((int(len(poss)/2),2))
                
    #             for j,p in enumerate(poss.reshape((int(len(poss)/2), 2,3))):
    #                 vectors[j] = (p[0] - p[1])[:2]
    #             distmat = 1 - cdist(vectors, vectors, metric='cosine')
    #             flattened = distmat[np.triu_indices_from(distmat, k=1)]
    #             angle = np.mean(abs(flattened))
    #             angles.append(angle)
    #             iters += 1

    #     if i == '15':
    #         np.save('storeangles15_flat.npy', np.copy(angles))
    #     else:
    #         np.save('storeangles30_flat.npy', np.copy(angles))

    # ####
    # a15 = np.load('storeangles15_flat.npy')
    # a30 = np.load('storeangles30_flat.npy')
    # np.mean(a15)
    # np.mean(a30)

    # t = range(len(a15[:100001]))
    # plt.plot(t,savgol_filter(a15[:100001], 10000, 2),color='red')
    # plt.plot(t,savgol_filter(a15[100001:200002], 10000, 2),color='red')
    # plt.plot(t,savgol_filter(a15[200002:], 10000, 2),color='red')
    # plt.plot(t,savgol_filter(a30[:100001], 10000, 2),color='blue')
    # plt.plot(t,savgol_filter(a30[100001:200002], 10000, 2),color='blue')
    # plt.plot(t,savgol_filter(a30[200002:], 10000, 2),color='blue')

    # fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    # ax.plot(t, a15,'o')
    # ax.set_rmax(2)
    # ax.set_rticks([0.5, 1, 1.5, 2])  # Less radial ticks
    # ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
    # ax.grid(True)



if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    main()