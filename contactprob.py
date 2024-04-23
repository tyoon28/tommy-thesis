from graphlets import *
import seaborn as sns


def resid_to_md_subunits(r):
    return np.array([(i * 337) + (r - 35) for i in range(0,4)])


def main(residues):
    # calculate contact probabilities between residues, identify contacts between subunits
    # make sure to record ones in subunit 1.
    residues = [75,175,176,177,178,179,180,181,182,183]
    selstring = 'resid '
    for res in residues:
        a = resid_to_md_subunits(res)
        selstring += ' or resid '.join(list(map(str,a)))
        selstring += ' or resid '
    selstring = selstring[:-10]
    
    storemat15 = np.zeros((len(residues)*4,len(residues)*4))
    dmat = np.zeros((len(residues)*4,len(residues)*4))
    for i in ['15','30']:
        for r in ['R1','R2','R3']:
            xtcs = []
            for file in os.listdir(f'{r}-{i}-closed'):
                if file.endswith('.xtc'):
                    xtcs.append(f'{r}-{i}-closed/'+file)
            xtcs.sort(key=lambda x: int(x.split('-')[1]))
            u = mda.Universe(f'{r}-{i}-closed/{r}-0-start-membrane-3JYC.pdb',*xtcs,continuous=True)
            protein = u.select_atoms(selstring)
            for ts in tqdm.tqdm(u.trajectory):
                frame = u.trajectory.frame
                r_compound = protein.atoms.center_of_mass(compound='residues')
                mat = distances.contact_matrix(r_compound, cutoff=6)
                np.fill_diagonal(mat, 0) 
                storemat += mat
        storemat = storemat / (len(u.trajectory)*3)
        if i == '15':
            storemat15 = np.copy(storemat)
        else:
            dmat = storemat - storemat15

    # dmat stores difference in contact prob between pairs of residues.
    np.save('dmat.npy', dmat)


    ax = sns.heatmap(dmat, linewidth=0.5)
    plt.savefig('dmat.png')
    return

if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    main()