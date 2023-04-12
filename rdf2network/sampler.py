
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
def grab_neighbours(seed,rg, limit = 10, recur = 0 , maxrecur = 3 , sample_run = 100 ):
    nodes = []
    #sample the neighbourhood randomly
    objs =  rg.objects(subject = seed ) 
    objs = [  next(objs) for i in range(sample_run) ]
    random.shuffle(objs)
    for i,obj in enumerate(objs):
        if i < limit and recur < maxrecur:
            nodes += grab_neighbours(obj , rg,limit = limit, recur = recur+1 )    
        elif i < limit:
            nodes.append( obj )
        if i > limit:
            break
    return nodes

def sample( rg , subject_iter = None , seed = None, sample_run= 10,  layers = 3 , layer_limit = 3 , verbose = True):
    if not seed:
        seed = next(subject_iter)
    if verbose == True:
        print('seed', seed)
    nodes = grab_neighbours(seed,rg, limit = layer_limit, recur = 0 , sample_run=sample_run, maxrecur = layers)
    subG = RDFGraph()
    triples = [ rg.triples((n1,None,n2)) for n1 in nodes for n2 in nodes if n1 != n2 ]
    for t in triples:
        subG+=t
    return subG

def add_xrefs(rg_info, subg):
    subg += rg_info
    print(subg)
    cross_ref = rdflib.term.URIRef("http://purl.org/lscr#xrefUniprot")
    string_ids_dict = {}
    # put in dictionary the IDs of original proteins (for which we got orthologs and paralogs)
    cross_ref = rdflib.term.URIRef("http://purl.org/lscr#xrefUniprot")
    for (s,p,o) in subg.triples((None, cross_ref, None)):
        # the triple is in the form (STRING_ID, xrefUniprot, UNIPROT_ID)
        uniprot_url = str(o)
        string_url = str(s)
        if(uniprot_url not in string_ids_dict):
            string_ids_dict[uniprot_url] = string_url
    return string_ids_dict , subg
