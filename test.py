from network_generation import *
from network_statistics import *
import os

xtcs = []
for file in os.listdir('R1-30-closed'):
    if file.endswith('.xtc'):
        xtcs.append('R1-30-closed/'+file)
xtcs.sort(key=lambda x: int(x.split('-')[1]))


u = mda.Universe('R1-30-closed/R1-0-start-membrane-3JYC.pdb','R1-30-closed/R1-0-1000-3JYC.xtc')
u = mda.Universe('R1-30-closed/R1-0-start-membrane-3JYC.pdb',*xtcs)

for ts in u.trajectory[:20]:
    time = u.trajectory.time
    rgyr = u.atoms.radius_of_gyration()
    print(f"Frame: {ts.frame:3d}, Time: {time:4.0f} ps, Rgyr: {rgyr:.4f} A")


def testspeed():


    # using contact matrix
    ress = u.select_atoms(f'(not resname CHOL and not resname POPC) or \
                                    (resname CHOL and resid {sterol_id} and not (name RC1 or name RC2))')
    start = time.time()
    for ts in u.trajectory[:10]:
        
        # all residues and sterol
        res = u.select_atoms(f'(not resname CHOL and not resname POPC) or \
                                    (resname CHOL and resid {sterol_id} and not (name RC1 or name RC2))').atoms.center_of_mass(compound='residues')
        contact_mat = distances.contact_matrix(res, cutoff=contact_threshold)


    end = time.time()
    print(end-start)


    # using distance matrix
    start = time.time()
    for ts in u.trajectory[:1000]:
        

        res = (u.select_atoms(f'(not resname CHOL and not resname POPC) or \
                                    (resname CHOL and resid {sterol_id} and not (name RC1 or name RC2))')
                                    .atoms.center_of_mass(compound='residues'))
        # worried about efficiency using distance array but can't find built in contact function for two selections...
        #distance array between sterol and all other residues
        dist_arr = distances.distance_array(res, 
                                    res, 
                                    box=u.dimensions)

    end = time.time()
    print(end-start)


rog = Cholesterol_contact(u)
rog.run(start=0,verbose=True)

rog._prepare()
rog._single_frame()

u.trajectory[1].time

for i in rog.results:
    if rog.results[i]['longest_contact'] > 30:
        print(i,rog.results[i]['longest_contact'])

def anybinding(r):
    n = []
    for i in r:
        n.append(max(r[i]['residuecontacts']))
    print(max(n))


r = u.select_atoms('not resname CHOL and not resname POPC').residues


for i in range(len(r.resnames)-1):
    if r.resnames[i] == 'LEU' and r.resnames[i+4] == 'SER'  and r.resnames[i+84] == 'ILE':
        print('aaa',i)
    

a = np.array([[0,0,0,0,0,0,0,0,0],
              [0,0,0,0,0,0,0,0,0],
              [0,0,2,3,4,5,6,0,0],
              [0,1,4,5,6,7,8,0,0],
              [0,1,4,5,6,7,8,0,0]])




# plot distribution of # frames of contact vs # sterols



degree_sequence = sorted((d for n, d in G.degree()), reverse=True)
dmax = max(degree_sequence)
plt.bar(*np.unique(degree_sequence, return_counts=True))
plt.title("Degree histogram")
plt.xlabel("Degree")
plt.ylabel("# of Nodes")
plt.show()