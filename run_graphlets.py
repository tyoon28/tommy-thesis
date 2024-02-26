from graphlets import *


def main():
    for r in ['R2,R3']:
        do_30_15(r)

        graphlets_cholesterol(r)
        

        ldirs = [f'../orca/output/{r}-15-closed-uniform',f'../orca/output/{r}-30-closed-uniform']
        finalDf,pca = PCA_gdd(ldirs)
        plot_PCA_gdd(finalDf,pca,f'{r}_PCA_gdd')
        xtcs = []
        for i in ['30','15']:
            for file in os.listdir(f'{r}-{i}-closed'):
                if file.endswith('.xtc'):
                    xtcs.append(f'{r}-{i}-closed/'+file)
            xtcs.sort(key=lambda x: int(x.split('-')[1]))
            u = mda.Universe(f'{r}-{i}-closed/{r}-0-start-membrane-3JYC.pdb',*xtcs)
        output_graphs_graphlets_cholesterol(r)

        
        graphlets_cholesterol_pca(r)