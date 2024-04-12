from graphlets import *
from itertools import combinations


def main():
    persistance()
    chol_contact()
    chol_interactionlength()

def correlation():
    # how long do edges last?
    for r,c in [(i,j) for i in ['R1','R2','R3'] for j in ['15','30']]:
        xtcs = []
        for file in os.listdir(f'{r}-{c}-closed'):
            if file.endswith('.xtc'):
                xtcs.append(f'{r}-{c}-closed/'+file)
        xtcs.sort(key=lambda x: int(x.split('-')[1]))
        u = mda.Universe('R1-30-closed/R1-0-start-membrane-3JYC.pdb',*xtcs,continuous=True)

        contact_threshold = 6
        sterols = u.select_atoms(f'(resname CHOL and not (name RC1 or name RC2))').residues
        protein = u.select_atoms('not resname CHOL and not resname POPC').residues
        lenp = len(protein)

        # for correlation values, try going through time and recording how many frames edges have in common.
        numedges = int(((lenp**2)-lenp)/2)

        # matrix storing how many times edge i and edge j are there together or absent together
        r = protein.atoms.center_of_mass(compound='residues')
        mat = distances.contact_matrix(r, cutoff=6)
        mm = np.zeros((lenp,lenp))

        results = np.zeros((numedges,numedges))
        for ts in tqdm.tqdm(u.trajectory):
            frame = u.trajectory.frame
            r = protein.atoms.center_of_mass(compound='residues')
            mat = distances.contact_matrix(r, cutoff=6)
            np.fill_diagonal(mat, 0)
            mm += (mat != 0)
        
        numedges = np.count_nonzero(mm)
        mask = (mm>0)
        mask.nonzero()

        # comp contains how many times each pair of edges is co-occurent.
        comp = np.zeros((len(u.trajectory),numedges),dtype=bool)
        i = 0
        for ts in tqdm.tqdm(u.trajectory):
            r = protein.atoms.center_of_mass(compound='residues')
            mat = distances.contact_matrix(r, cutoff=6)
            comp[i] = mat[mask]
            i += 1

        df = pd.DataFrame(comp)
        df.corr(method='pearson')




        


        
        
        
        print('setup done')
        edges = []

        # how many edges exist for less than 90%?
        for ts in tqdm.tqdm(u.trajectory):
            frame = u.trajectory.frame
            r = protein.atoms.center_of_mass(compound='residues')
            mat = distances.contact_matrix(r, cutoff=6)
            np.fill_diagonal(mat, 0)
            edgesmat += mat
        t = len(u.trajectory[20:]) * 0.9
        mask = (edgesmat < t)
        num_non_persistent = mask.sum()/2
    return

def persistance():

    # histogram of edges versus how long they last
    for c in ['15','30']:
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
            store = np.zeros((lenp,lenp))
            for ts in tqdm.tqdm(u.trajectory):
                r = protein.atoms.center_of_mass(compound='residues')
                mat = distances.contact_matrix(r, cutoff=6)
                np.fill_diagonal(mat, 0)
                store += mat
        
        store = store / (len(u.trajectory) * 3)
        # get all edges in a 1D np array
        no0 = store[np.where(store!=0)]
        mid = no0[np.where(no0 > 0.1)]
        mid = mid[np.where(mid > 0.9)]

        plt.hist(no0)
        print(f'persistence {c}: len(mid)')
    plt.xlabel('average % persistance')
    plt.ylabel('count')
    plt.legend(['15 mol%','30 mol%'])

    plt.savefig('figure_edgepersistence.png')



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




if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    main()