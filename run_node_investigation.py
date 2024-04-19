from graphlets import *

def main():
    d = {} # key: residue; item: dict of concentrations with dict of contacts and frequencies
     
    residues = None # DON"T RUN WITHOUT SETTING RESIDEUSE
    for i in ['30','15']:
        for r in ['R1','R2','R3']:
            d[i] = {}
            xtcs = []
            for file in os.listdir(f'{r}-{i}-closed'):
                if file.endswith('.xtc'):
                    xtcs.append(f'{r}-{i}-closed/'+file)
            xtcs.sort(key=lambda x: int(x.split('-')[1]))
            u = mda.Universe(f'{r}-{i}-closed/{r}-0-start-membrane-3JYC.pdb',*xtcs,continuous=True)

            protein = u.select_atoms('not resname CHOL and not resname POPC')
            for res in residues:
                # record frequency of each unique set of contacts
                d[i][res] = {}
                print(res,i)
                for ts in tqdm.tqdm(u.trajectory):
                    frame = u.trajectory.frame
                    r_compound = protein.atoms.center_of_mass(compound='residues')
                    mat = distances.contact_matrix(r_compound, cutoff=6)
                    np.fill_diagonal(mat, 0) 
                    # +/- 1 from resid for proper indexing
                    nz = np.nonzero(mat[res-1])[0] +1
                    fs = frozenset(nz)
                    if fs in d[i][res].keys():
                        d[i][res][fs] += 1 
                    else: d[i][res][fs] = 1 
    
        
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
            if abs(d_difference[res][i]) > 10000:
                d_diff_filtered[res][i] = d_difference[res][i]
        d_diff_filtered[res] = {k: v for k, v in sorted(d_diff_filtered[res].items(), key=lambda item: item[1])}


def resid_to_md_subunits(r):
    r-35