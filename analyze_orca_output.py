# to be run after make_orca_input and running orca
from graphlets import *
import warnings


#TODO: after inspection look for G14 in response to choleseterol
def main():
    # for r in ['R1','R2','R3']:
    #     print(f'working on {r}')
    #     ldirs = [f'../orca/output/{r}-15-closed',f'../orca/output/{r}-30-closed']

    #     df,finalDf,pca = PCA_gdd(ldirs)
    #     print(f'making PCA plots for {r}')
    #     plot_PCA_gdd(finalDf,f'{r}_PCA_gdd')

    #     print(f'calculating variable importance for {r}')
    #     logistic_selection(df,r)

    #     print(f'doing cholesterol plots for {r}')
    #     graphlets_cholesterol_pca(r)

    #     print(f'doing individual nodes for {r}')
    #     node_pca_analysis(r,output=True)

    #     print(f'doing individual nodes windowed for {r}')
    #     node_PCA_windowed(r,output=False,to_csv=True)


        

        # print(f'doing individual nodes with cholesterol for {r}')
        # node_pca_analysis(r)

    #     print(f'done with {r}')

    print('running combined')
    r = 'all'
    ldirs = [f'../orca/output/{r}-{c}-closed' for r in ['R1','R2','R3'] for c in ['15','30']]
    replicate_investigation(ldirs)

    df, finalDf,pca = PCA_gdd(ldirs,to_csv=True)
    print(f'making PCA plots for all')
    plot_PCA_gdd(finalDf,f'all_PCA_gdd')

    print(f'calculating variable importance for all')
    logistic_selection(df,r)

    print(f'doing cholesterol plots for all')
    graphlets_cholesterol_pca('all',to_csv=True)

    # print(f'doing individual nodes for all')
    # node_pca_analysis(r,output=True)

    print(f'doing individual nodes windowed for all')
    node_PCA_windowed(r,output=True,to_csv=True)


    # TODO: which nodes are highly predictive of cholesterol condition?

    # print(f'doing individual nodes with cholesterol for {r}')
    # node_chol_analysis(r,node)

    print('done')





if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    main()