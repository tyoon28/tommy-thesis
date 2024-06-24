'''
for select residues (rnames) calculate how likely they are to be in contact. This is figure 10B.
'''
from graphlets import *
import seaborn as sns


def resid_to_md_subunits(r):
    return np.array([(i * 337) + (r - 35) for i in range(0,4)])


def main():
    # calculate contact probabilities between residues, identify contacts between subunits
    # make sure to record ones in subunit 1.
    rnames = ['V75','F140','G143','I136''F175','M176','I177','G178','A179','I180','M181','A182','K183','M184','A185']

    residues = [75,140,143,136,175,176,177,178,179,180,181,182,183,184,185]
    selstring = 'resid '
    for res in residues:
        a = resid_to_md_subunits(res)
        selstring += ' or resid '.join(list(map(str,a)))
        selstring += ' or resid '
    selstring = selstring[:-10]
    
    dmat = np.zeros((len(residues)*4,len(residues)*4))
    for i in ['15','30']:
        storemat = np.zeros((len(residues)*4,len(residues)*4))

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
    np.save('dmat3.npy', dmat)


    ax = sns.heatmap(dmat, linewidth=0.5,xticklabels=rnames*4, yticklabels=rnames*4)
    plt.savefig('dmat3.png')

    
    shif = len(residues)
    rnames = ['V75','F175','M176','I177','G178','A179','I180','M181','A182','K183','M184','A185']
    dmat2 = np.load('dmat3.npy')
    a = np.zeros((shif,shif))

    for row in range(len(a)):
        ro = np.mean(dmat2[[row,row+shif,row+(2*shif),row+(3*shif)]],axis=0)
        for j in range(shif):
            a[row][j] = np.mean(ro[[j,j+shif,j+(2*shif),j+(3*shif)]])
    ax = sns.heatmap(a, linewidth=0.5,xticklabels=rnames, yticklabels=rnames)
    plt.savefig('a3.png')



    contactpairs = [(181,177),(181,182),(181,181),(181,75),(184,75),(177,177)]
    selstring = 'resid '
    residues = 181,177,182,75,184
    for res in residues:
        a = resid_to_md_subunits(res)
        selstring += ' or resid '.join(list(map(str,a)))
        selstring += ' or resid '
    selstring = selstring[:-10]

    contactpairs = [(2,1),(2,3),(2,2),(2,0),(4,0),(1,1)]

    shif = 5
    for i in ['15','30']:
        cond = np.zeros((len(contactpairs),100001))
        for r in ['R1','R2','R3']:
            xtcs = []
            for file in os.listdir(f'{r}-{i}-closed'):
                if file.endswith('.xtc'):
                    xtcs.append(f'{r}-{i}-closed/'+file)
            xtcs.sort(key=lambda x: int(x.split('-')[1]))
            u = mda.Universe(f'{r}-{i}-closed/{r}-0-start-membrane-3JYC.pdb',*xtcs,continuous=True)
            protein = u.select_atoms(selstring)
            a = np.zeros((len(contactpairs),len(u.trajectory)))
            for ts in tqdm.tqdm(u.trajectory):
                frame = u.trajectory.frame
                r_compound = protein.atoms.center_of_mass(compound='residues')
                mat = distances.contact_matrix(r_compound, cutoff=6)
                np.fill_diagonal(mat, 0) 
                j = 0
                for p1,p2 in contactpairs:
                    resids1 = np.array([p1,p1+5,p1+10,p1+15])
                    resids2 = np.array([p2,p2+5,p2+10,p2+15])
                    k = np.count_nonzero(mat[resids1,:][:,resids2])
                    a[j][frame] = k
                    j+=1
            a = a/2 # since matrix is symmetrical
            cond += a
        
    t = list(range(100001))
    for k in a:
        plt.plot(t,savgol_filter(k, 1000, 2))
    # plt.legend([('M181','I177'),('M181','A182'),('M181','M181'),('M181','V75'),('M184','V75'),('I177','I177')])
    plt.ylabel('# contacts (subunits)')
    plt.xlabel('time')
    plt.ylim(0,3.4)
    plt.show()

    return

if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    main()