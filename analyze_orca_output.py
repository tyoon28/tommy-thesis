'''
This script is run after make_orca_input and running orca on the output.
Generates PCA plots and analyses for reach replicate separately, then all replicates combined
'''

from graphlets import *
import warnings


def main():
    rc('font',**{'family':'serif','serif':['Times']})
    for r in ['R1','R2','R3']:
        print(f'working on {r}')
        ldirs = [f'../orca/output/{r}-15-closed',f'../orca/output/{r}-30-closed']

        df,finalDf,pca = PCA_gdd(ldirs)
    
        print(f'making PCA plots for {r}')
        plot_PCA_gdd(finalDf,f'{r}_PCA_gdd',evr = pca.explained_variance_ratio_[:2])

        print(f'calculating variable importance for {r}')
        logistic_selection(df,r)

        print(f'doing cholesterol plots for {r}')
        graphlets_cholesterol_pca(r)

        print(f'doing individual nodes for {r}')
        node_pca_analysis(r,output=True)

        print(f'doing individual nodes windowed for {r}')
        node_PCA_windowed(r,output=False,to_csv=True)


        

        print(f'doing individual nodes with cholesterol for {r}')
        node_pca_analysis(r)

        print(f'done with {r}')

    print('running combined')
    r = 'all'
    ldirs = [f'../orca/output/{r}-{c}-closed' for r in ['R1','R2','R3'] for c in ['15','30']]
    replicate_investigation(ldirs)

    df, finalDf,pca,nPCs = PCA_gdd(ldirs,to_csv=True)
    pcpair = find_PCpair(finalDf,pca,nPCs)
    print(f'making PCA plots for all')
    plot_PCA_gdd(finalDf,f'all_PCA_gdd',evr = pca.explained_variance_ratio_,pcpair = pcpair)

    print(f'calculating variable importance for all')
    logistic_selection(df,r)

    print(f'doing cholesterol plots for all')
    graphlets_cholesterol_pca('all',to_csv=True)

    print(f'doing individual nodes for all')
    node_pca_analysis(r,output=True)

    print(f'doing individual nodes windowed for all')
    node_PCA_windowed(r,output=True,to_csv=True)

    # print(f'doing individual nodes with cholesterol for {r}')
    # node_chol_analysis(r,node)

    print('done')





if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    main()