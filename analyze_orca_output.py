# to be run after make_orca_input and running orca
from graphlets import *
import warnings



def main():
    for r in ['R2','R3']:
        print(f'working on {r}')
        ldirs = [f'../orca/output/{r}-15-closed-uniform',f'../orca/output/{r}-30-closed-uniform']

        finalDf,pca = PCA_gdd(ldirs)
        print(f'making PCA plots for {r}')
        plot_PCA_gdd(finalDf,pca,f'{r}_PCA_gdd')

        print(f'doing cholesterol plots for {r}')
        graphlets_cholesterol_pca(r)

        print(f'doing individual nodes for {r}')
        node_pca_analysis(r)

        print(f'done with {r}')




if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    main()