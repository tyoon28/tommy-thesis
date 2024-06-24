'''
process MD simulation data into output for orca
'''
from graphlets import *
import warnings





def main():
    for r in ['R1','R2','R3']:
        print(f'working on {r}')
        # make make uniform distributed snapshots for 30 and 15 conditions --> full network analysis
        print(f'making uniform snapshots for {r}')
        do_30_15(r)
        

        # make graphs for different cholesterol conditions --> cholesterol analysis
        print(f'making graphs for cholesterol occupancy for {r} in closed-contact')
        output_graphs_graphlets_cholesterol(r,0.4)

        # make full graphs for 15 and 30 conditions --> node analysis
        print('making full graphs for 15 and 30 conditions')
        output_for_node_pca(r)





if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    main()