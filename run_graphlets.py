from graphlets import *

# this file does nothing. don't use it
def main():
    for r in ['R1','R2,R3']:
        do_30_15(r) # output graphs - these go to r-c-closed/

        # ldirs = [f'../orca/output/{r}-15-closed-uniform',f'../orca/output/{r}-30-closed-uniform']
        # finalDf,pca = PCA_gdd(ldirs)
        # plot_PCA_gdd(finalDf,pca,f'{r}_PCA_gdd')
        # xtcs = []
        # for i in ['30','15']:
        #     for file in os.listdir(f'{r}-{i}-closed'):
        #         if file.endswith('.xtc'):
        #             xtcs.append(f'{r}-{i}-closed/'+file)
        #     xtcs.sort(key=lambda x: int(x.split('-')[1]))
        #     u = mda.Universe(f'{r}-{i}-closed/{r}-0-start-membrane-3JYC.pdb',*xtcs)
        # outputs ../orca/input/{basename}-contact/{basename}-contact-c{n}-f{f}.in
        # where n is n cholesterol
        output_graphs_graphlets_cholesterol(r) # output graphs - these go to R-c-closed-contact/

        
        # graphlets_cholesterol_pca(r)
