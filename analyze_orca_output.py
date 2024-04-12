# to be run after make_orca_input and running orca
from graphlets import *
import warnings



def main():
    for r in ['R1','R2','R3']:
        print(f'working on {r}')
        ldirs = [f'../orca/output/{r}-15-closed',f'../orca/output/{r}-30-closed']

        finalDf,pca = PCA_gdd(ldirs)
        print(f'making PCA plots for {r}')
        plot_PCA_gdd(finalDf,f'{r}_PCA_gdd')
        # DO VARIABLE IMPORTANCE

        print(f'doing cholesterol plots for {r}')
        graphlets_cholesterol_pca(r)

        print(f'doing individual nodes for {r}')
        node_pca_analysis(r,output=True)

        print(f'doing individual nodes windowed for {r}')
        node_PCA_windowed(r,ldirs,output=False)


        

        # print(f'doing individual nodes with cholesterol for {r}')
        # node_pca_analysis(r)

        print(f'done with {r}')
    print('running combined')
    ldirs = [f'../orca/output/{r}-15-closed',f'../orca/output/{r}-30-closed' for r in ['R1','R2','R3']]

    finalDf,pca = PCA_gdd(ldirs)
    print(f'making PCA plots for all')
    plot_PCA_gdd(finalDf,f'all_PCA_gdd')

    print(f'doing cholesterol plots for all')
    graphlets_cholesterol_pca('all')

    print(f'doing individual nodes for all')
    node_pca_analysis(r,output=True)

    print(f'doing individual nodes windowed for all')
    node_PCA_windowed(r,ldirs,output=False)

    # print(f'doing individual nodes with cholesterol for {r}')
    # node_pca_analysis(r)





if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    main()