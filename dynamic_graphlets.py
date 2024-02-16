from network_statistics import *


def output_temporal_graph(u,filename,d=None):
    res = u.select_atoms('not resname CHOL and not resname POPC')

    with open(filename, 'w') as f:
        for ts in u.trajectory[:d]:
            frame = u.trajectory.frame
            r = res.atoms.center_of_mass(compound='residues')
            mat = distances.contact_matrix(r, cutoff=6)
            np.fill_diagonal(mat, 0)

            G = nx.from_numpy_array(mat)
            for line in nx.generate_edgelist(G,delimiter = ' ',data=False):
                linesp = line.split()
                f.write(f'{linesp[0]}\t{linesp[1]}\t{frame}\n')



