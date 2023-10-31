from network_generation import *
import os

xtcs = []
for file in os.listdir('R1-15-closed'):
    if file.endswith('.xtc'):
        xtcs.append('R1-15-closed/'+file)
xtcs.sort(key=lambda x: int(x.split('-')[1]))


u = mda.Universe('R1-15-closed/R1-0-start-membrane-3JYC.pdb','R1-15-closed/R1-0-1000-3JYC.xtc')
u = mda.Universe('R1-15-closed/R1-0-start-membrane-3JYC.pdb',*xtcs)




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
rog.run(start=0,stop=20000,verbose=True)


res = (u.select_atoms('resid 2602 or not resname CHOL and not resname POPC').atoms.center_of_mass(compound='residues'))
x = distances.contact_matrix(res, cutoff=contact_threshold)

