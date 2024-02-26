from graphlets import *


def main():
    for r in ['R2,R3']:
        # make make uniform distributed snapshots for 30 and 15 conditions --> full network analysis
        do_30_15(r)

        # make graphs for different cholesterol conditions --> cholesterol analysis
        output_graphs_graphlets_cholesterol(r)

        # make full graphs for 15 and 30 conditions --> node analysis
        output_for_node_pca(r)




if __name__ == "__main__":
    main()