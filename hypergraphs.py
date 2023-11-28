import graph_tool as gt
import sys
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PSF, DCD, GRO, XTC
from MDAnalysis.analysis import contacts, distances
from MDAnalysis.analysis.base import AnalysisBase
import numpy as np



def hypergraph_1_frame(u,res):
    '''
    infer hypergraph from 1 frame. return list of hyperedges
    '''
    NITERS = 10000

    G = gt.Graph(directed=False)
    mat = distances.contact_matrix(res, cutoff=6)
    G.add_edge_list(np.transpose(mat.nonzero()))

    ntimes = 1

    # from https://graph-tool.skewed.de/static/doc/autosummary/graph_tool.inference.CliqueState.html
    ## WHAT IS THIS
    state = gt.inference.CliqueState(G)
    output = state.mcmc_sweep(niter=NITERS)

    cliques = {} # indexed by k
    # iterate through factor graph
    for v in state.f.vertices():
        if state.is_fac[v]:
            continue # skip over factors

        if state.x[v] > 0: # activated hyperedge
            k = len(state.c[v])
            if k not in cliques:
                cliques[k] = []
            cliques[k].append(state.c[v])
    
    return cliques

def hypergraph_1_frame_simple(u,res):
    '''
    infer hypergraph from 1 frame. return list of hyperedges
    '''
    NITERS = 10000

    G = gt.Graph(directed=False)
    mat = distances.contact_matrix(res, cutoff=6)
    G.add_edge_list(np.transpose(mat.nonzero()))

    ntimes = 1

    # from https://graph-tool.skewed.de/static/doc/autosummary/graph_tool.inference.CliqueState.html
    ## WHAT IS THIS
    state = gt.inference.CliqueState(G)
    output = state.mcmc_sweep(niter=NITERS)

    cliques = {} # indexed by k
    # iterate through factor graph
    for v in state.f.vertices():
        if state.is_fac[v]:
            continue # skip over factors

        if state.x[v] > 0: # activated hyperedge
            k = len(state.c[v])
            if k not in cliques:
                cliques[k] = []
            cliques[k].append(state.c[v])
    
    return cliques


#look for residues which move as one unit --> hypergraph edge?

