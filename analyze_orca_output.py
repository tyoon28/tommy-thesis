# to be run after make_orca_input and running orca
from graphlets import *
import warnings



def main():
    for r in ['R2,R3']:
        ldirs = [f'../orca/output/{r}-15-closed-uniform',f'../orca/output/{r}-30-closed-uniform']
        finalDf,pca = PCA_gdd(ldirs)
        plot_PCA_gdd(finalDf,pca,f'{r}_PCA_gdd')
        graphlets_cholesterol_pca(r)
        node_pca_analysis(r)





if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    main()