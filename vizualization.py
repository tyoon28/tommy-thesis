'''
vizualizing networks
'''
from pyvis.network import Network
from network_generation import *

def viz_graph(G):
    # output an image file of a graph. use this to make animation of dynamic graph?
    color_map = []
    for node in G.nodes:
        if node in Ia: color_map.append('green')
        elif node > 1348: color_map.append('orange')
        else: color_map.append('yellow')

    G.remove_edges_from(nx.selfloop_edges(G))


    nx.draw(G, node_color=color_map, font_weight='bold',with_labels=True,pos=nx.circular_layout(G))
    pass


def cholesterol_occupancy(chol):
    #chol is cholesterol analysisbase
    pass


def cholesterol_network_series(u,sterol_id,rog):
    '''make a series of networks of given cholesterol and the residues it touches
    '''
    allsites = binding_sites('all')
    if not os.path.isdir(str(sterol_id)):
        os.makedirs(str(sterol_id))

    sterol_id = 2723
    frames = rog.results[sterol_id]['contacts']

    # TODO: set trajectory range
    for ts in u.trajectory[1497:1753]:
        plt.clf()
        G = network_1_cholesterol(u,sterol_id)
        color_map = []
        for node in G.nodes:
            if node in allsites: color_map.append('green')
            elif node > 1348: color_map.append('orange')
            else: color_map.append('yellow')

        G.remove_edges_from(nx.selfloop_edges(G))




        center_node = sterol_id  # Or any other node to be in the center
        edge_nodes = set(G) - {center_node}
        # Ensures the nodes around the circle are evenly distributed
        if len(G.nodes) == 2:
            pos = {}
            pos[center_node] = np.array([0, 0])  
            pos[list(edge_nodes)[0]] = np.array([0, 0.5]) 
        else:
            pos = nx.circular_layout(G.subgraph(edge_nodes))
            pos[center_node] = np.array([0, 0])  # manually specify node position
        plt.title(f'frame: {ts.frame}')

            
        nx.draw(G, pos = pos, node_color=color_map, font_weight='bold',with_labels=True)
        plt.xlim([-2, 2])
        plt.ylim([-2, 2])


        plt.savefig(f"{sterol_id}/{ts.frame}.png", format="PNG")




def longest(a):
    # return start/end of longest contact
    if not a: 
        return 0
    longest = 0
    last = a[0]-1
    s = 0
    start = 0
    best = (0,0)
    for i in a:
        if i-1 == last:
            s += 1
            if s > longest: longest = s
            end = i
        else:
            s = 1
            if end-start > best[1]-best[0]:
                best = (start,end)
            start = i

        last = i
    print(best)