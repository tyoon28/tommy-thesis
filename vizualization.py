'''
vizualizing networks
'''

# random forest to predict open/closed or cholesterol percent based on which edges exist in protein?
# compare with same model but with hyperedges?


from pyvis.network import Network
from network_generation import *

def viz_graph(G):
    color_map = []
    for node in G.nodes:
        if node in binding_sites('all'): color_map.append('green')
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


def graph_cholesterol_occupancy(rog):
    '''
    graph the number of cholesterols in contact with binding sites over time
    '''

    y = np.zeros(100020)
    for ch in rog.results:
        for i in rog.results[ch]['contacts']:
            y[i]+=1
    plt.plot(y,linewidth=0.1)
    plt.show()

    counts, bins = np.histogram(y)
    plt.stairs(counts, bins)


    pass

def graph_cholesterol_contact(rog):
    '''
    graph the number of cholesterols in contact with binding sites over time
    '''

    y = []
    x = []
    for ch in rog.results:
        y.append(len(rog.results[ch]['contacts']))

    for i in y:
        if i  > 10000:
            x.append(i)
    counts, bins = np.histogram(x)
    plt.stairs(counts, bins)
    counts, bins = np.histogram(y)
    plt.stairs(counts, bins)

    plt.xlabel("total # frames in contact with channel")
    plt.ylabel("# cholesterols")

    x = []
    for ch in rog.results:
        x.append(rog.results[ch]['longest_contact'])
    counts, bins = np.histogram(x)
    plt.stairs(counts, bins)
    plt.xlabel("longest contact")
    plt.ylabel("# cholesterols")

    plt.show()


    pass


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

def residue_network_series(u,r_ids,highlight=[]):
    '''
    make a series of contact networks of given residues - networks only include those residues.
    highlight selected residues in highlight (not implemented)
    '''
    allsites = binding_sites('all')
    if not os.path.isdir(str(r_ids[0])):
        os.makedirs(str(r_ids[0]))


    # TODO: set trajectory range
    for ts in u.trajectory:
        plt.clf()
        G = network_selected(u,r_ids)
        color_map = []
        for node in G.nodes:
            if node in allsites: color_map.append('green')
            elif node > 1348: color_map.append('orange')
            else: color_map.append('yellow')

        G.remove_edges_from(nx.selfloop_edges(G))

        if len(G.edges) == 0: continue

        center_node = r_ids[0]  # Or any other node to be in the center
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


        plt.savefig(f"{r_ids[0]}/{ts.frame}.png", format="PNG")

def sterol_contact_hist(rog):
    '''
    y axis: num sterols
    x axis: num contacts
    '''
    y = []
    for ch in rog.results:
        if len(rog.results[ch]['contacts']) > 10000 and len(rog.results[ch]['contacts']) <15000:
            y.append(len(rog.results[ch]['contacts']))

    counts, bins = np.histogram(y)
    plt.stairs(counts, bins)

    plt.ylabel("total # frames in contact with channel")
    plt.xlabel("# cholesterols")
    plt.show()

    return

