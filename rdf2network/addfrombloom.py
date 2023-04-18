
import itertools
import rdflib

#load bloom filters for string data
import pickle
import glob

def check_filters(val,filters = filters):
    for f in filters:
        if val in f[0][0]:
            return True
    return False

def load_filters(filters = './filters/bloomfinal_big*.pkl' ):
    filters = glob.glob(filters)
    filters = [ pickle.load(open(f, 'rb')) for f in filters ]
    #check bloom filters for interactors in diff species
    return filters

def check_allvall( objects , predicate = rdflib.term.URIRef('https://string-db.org/rdf/high-confidence-cutoff') ):
    #replace a call to a server with a check on bloom filters
    triples = [ (o1, predicate, o2) for o1,o2 in itertools.combinations(objects,2)
    if check_filters (o1 + '_' + o2  , filters=filters) ]
    return triples


    
    

