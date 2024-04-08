from graphlets import *
from itertools import combinations


def main():
    contact_threshold = 6

    for r,c in [(i,j) for i in ['R1','R2','R3'] for j in ['15','30']]:
        xtcs = []
        for file in os.listdir(f'{r}-{c}-closed'):
            if file.endswith('.xtc'):
                xtcs.append(f'{r}-{c}-closed/'+file)
        xtcs.sort(key=lambda x: int(x.split('-')[1]))
        u = mda.Universe('R1-30-closed/R1-0-start-membrane-3JYC.pdb',*xtcs)

        protein = u.select_atoms('not resname CHOL and not resname POPC').residues
        chol = u.select_atoms(f'(resname CHOL and not (name RC1 or name RC2))').residues

        protein_sterols = protein + chol

        t = np.zeros(len(u.trajectory)-20)
        y = np.zeros(len(u.trajectory)-20)
        y2=np.zeros(len(u.trajectory)-20)
        # count cholesterols occupying
        # have lines in background, thick line for average of the replicates, and both conditions in same plot.
        for ts in tqdm.tqdm(u.trajectory[20:]):
            tim = u.trajectory.time
            t[u.trajectory.frame-20] = tim


            r = protein_sterols.atoms.center_of_mass(compound='residues')
            contact_mat = distances.contact_matrix(r, cutoff=contact_threshold)

            # take only items representing chol-protein contacts
            z = contact_mat[:len(protein),len(protein):]
            occupancy = np.sum(np.any(z,0))
            occupancy2 = np.sum(np.any(z,1))
            y[u.trajectory.frame-20] = occupancy
            y2[u.trajectory.frame-20] = occupancy2

        y2_subsample = y2[::100]
        t_subsample = t[::100]
        plt.plot(t,y)
        plt.plot(t,y2)
        plt.plot(t_subsample,y2_subsample)
        plt.show()

def persistence():
    # how long do edges last?
    for r,c in [(i,j) for i in ['R1','R2','R3'] for j in ['15','30']]:
        xtcs = []
        for file in os.listdir(f'{r}-{c}-closed'):
            if file.endswith('.xtc'):
                xtcs.append(f'{r}-{c}-closed/'+file)
        xtcs.sort(key=lambda x: int(x.split('-')[1]))
        u = mda.Universe('R1-30-closed/R1-0-start-membrane-3JYC.pdb',*xtcs)

        contact_threshold = 6
        sterols = u.select_atoms(f'(resname CHOL and not (name RC1 or name RC2))').residues
        protein = u.select_atoms('not resname CHOL and not resname POPC').residues
        lenp = len(protein)

        # for correlation values, try going through time and recording how many frames edges have in common.
        numedges = int(((lenp**2)-lenp)/2)

        # matrix storing how many times edge i and edge j are there together or absent together
        r = protein.atoms.center_of_mass(compound='residues')
        mat = distances.contact_matrix(r, cutoff=6)
        mm = mat[np.triu_indices(len(mat), k = 1)]
        comp = np.zeros((len(mm),len(mm)))

        results = np.zeros((numedges,numedges))
        for ts in tqdm.tqdm(u.trajectory[20:]):
            frame = u.trajectory.frame
            r = protein.atoms.center_of_mass(compound='residues')
            mat = distances.contact_matrix(r, cutoff=6)
            np.fill_diagonal(mat, 0)

            mm = mat[np.triu_indices(len(mat), k = 1)] # take upper triangle of contact array

            comp += (mm == mm.reshape(len(mm),1))

        # comp contains how many times each pair of edges is co-occurent
        
        
        print('setup done')
        edges = []

        # how many edges exist for less than 90%?
        for ts in tqdm.tqdm(u.trajectory[20:]):
            frame = u.trajectory.frame
            r = protein.atoms.center_of_mass(compound='residues')
            mat = distances.contact_matrix(r, cutoff=6)
            np.fill_diagonal(mat, 0)
            edgesmat += mat
        t = len(u.trajectory[20:]) * 0.9
        mask = (edgesmat < t)
        num_non_persistent = mask.sum()/2

        # histogram of edges versus how long they last




def chol_contact():
    #cholesterol contacts tend to be short, with most lasitng X ns. longest contact is X ns.
    # lengths of interactions of binding sites vs with rest of protein

    for r,c in [(i,j) for i in ['R1','R2','R3'] for j in ['15','30']]:
        xtcs = []
        for file in os.listdir(f'{r}-{c}-closed'):
            if file.endswith('.xtc'):
                xtcs.append(f'{r}-{c}-closed/'+file)
        xtcs.sort(key=lambda x: int(x.split('-')[1]))
        u = mda.Universe('R1-30-closed/R1-0-start-membrane-3JYC.pdb',*xtcs)


        rog = Cholesterol_contact(u)
        rog.run(start=0,verbose=True)

        all_interactions = [] # list of durations of all interactions
        bs_vs_rest = [[],[]] # 0: durations of interactions that touch a binding site; 1: others
        for c in rog.results:
            for i in range(len(rog.results[c]['binding_events'])):
                start,end = rog.results[c]['binding_events'][i]
                dur = end-start

                bound = rog.results[c]['binding_events_actual'][i]
                if bound: ind = 0
                else: ind = 1
                bs_vs_rest[ind].append(dur)
                all_interactions.append(dur)

        # plot with color corresponding to replicate/cholesterol condition

    # A) a histogram?
    # B) a violinish plot with 2 x values - binding site and not
    




        



    pass
        



    

if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    main()