
from rdflib import Graph as RDFGraph
from rdflib.extras.external_graph_libs import rdflib_to_networkx_graph
import networkx as nx
from networkx import Graph as NXGraph
import matplotlib.pyplot as plt
import statistics
import collections
from rdflib.namespace import RDF
import rdflib
import random
import gzip
import stringdb


def load_graph(datapath):
    with gzip.open(datapath,'rb') as inputdata:
        rg = RDFGraph()
        rg.parse(inputdata.read(), format='turtle')
        print("rdflib Graph loaded successfully with {} triples".format(len(rg)))
        return rg
    
#recursive function to get some neighbours
def grab_neighbours(seed,rg, limit = 10, recur = 0 , maxrecur = 3 , sample_run = 10 ):
    nodes = []
    #sample the neighbourhood randomly
    obj_gen =  rg.objects(subject = seed )
    objs = []
    subgraph = None
    for obj in obj_gen:
        triples = [ t for t in rg.triples((obj,None,seed))]
        nodes.append(obj)
        if subgraph is None:
            subgraph = triples
        else:
            subgraph += triples
        objs.append(obj)
        if len( objs) > sample_run-1:
            break
    random.shuffle(objs)
    for i,obj in enumerate(objs):
        if i < limit and recur < maxrecur:
            subnodes, subg = grab_neighbours(obj,rg,limit = limit, recur = recur+1 )    
            nodes += subnodes
            subgraph += subg
        elif i < limit:
            nodes.append( obj )
            subgraph += [ t for t in rg.triples((obj,None,seed))]
        if i > limit:
            break
    return nodes, subgraph

def sample( rg  , seed , sample_run= 10, limit = 5,  layers = 3 , layer_limit = 3 , verbose = True):
    nodes , subgraph = grab_neighbours(seed,rg, limit = layer_limit , recur = 0 , sample_run=sample_run, maxrecur = layers)
    print(len(nodes) , len(subgraph))
    subG = RDFGraph()
    [   subG.add( (s,p,o) ) for s,p,o in subgraph ]
    return subG

def add_xrefs(rg_info, subg):
    print(subg)
    cross_ref = rdflib.term.URIRef('http://purl.org/lscr#xrefUniprot')
    string_ids_dict = {}
    # put in dictionary the IDs of original proteins (for which we got orthologs and paralogs)
    #iterate over all the subjects and objects in subg
    prots = set([ o for o in subg.objects() ] ).union(
        set([ s for s in subg.subjects() ] )
    )
    print(len(prots))
    triples = [(s,p,o)  for p in prots for s,p,o in rg_info.triples((p, cross_ref, None)) ]
    print(len(triples))
    [ subg.add(t) for t in triples ]
    return subg
