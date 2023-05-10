
import itertools
import rdflib
import functools
#load bloom filters for string data
import pickle
import glob
import multiprocessing as mp

import numpy as np

def check_filters(val,filters ):
    for f in filters:

        if val in f[0][0]:
            return True
    return False

def checkone(val,f):
    if val in f[0][0]:
        return True
    return False

def check_all(vals,f):
    #check all filters for each val
    #return true if any are true
    #return false if all are false
    return [ checkone(v,f) for v in vals ]


def check_filters_mp(vals,filters , pool , verbose = False ):
    #use map async to check all filters
    #return true if any are true
    #return false if all are false
    check = functools.partial(check_all,vals)
    if verbose == True:
        print(vals[0:10])
    res = pool.map(check, filters , chunksize = 1)

    # boolean or on list of lists
    res = np.array(res)
    res = np.any(res, axis=0)
    return list(res)

def load_filters(filters = './filters/bloomfinal_big*.pkl' ):
    filters = glob.glob(filters)
    filters = [ pickle.load(open(f, 'rb')) for f in filters ]
    #check bloom filters for interactors in diff species
    return filters

def URIfmt(object , urlstring = 'https://string-db.org/network/' ):
    return rdflib.term.URIRef( urlstring + object)

def check_allvall( objects , predicates = [rdflib.term.URIRef('https://string-db.org/rdf/high-confidence-cutoff') , rdflib.term.URIRef('https://string-db.org/rdf/any-confidence') ]
                   , urlstring= 'https://string-db.org/network/', filters = None ):
    #replace a call to a server with a check on bloom filters
    if filters is None:
        filters = load_filters()
    triples = []
    for predicate in predicates:
        triples += [ (URIfmt(o1 , urlstring) ,predicate , URIfmt(o2, urlstring)) for o1,o2 in itertools.combinations(objects,2)
        if check_filters(o1 + '_'+o2 , filters=filters) ]
    return triples

def check_allvall_mp( objects , predicate = rdflib.term.URIRef('https://string-db.org/rdf/high-confidence-cutoff') , urlstring= 'https://string-db.org/network/', ncpu = 10, filters = None ):
    #replace a call to a server with a check on bloom filters
    if filters is None:
        filters = load_filters()
    combos = [  o1+'_'+o2 for o1,o2 in itertools.combinations(objects,2) ]
    if ncpu is None:
        pool = mp.Pool(mp.cpu_count())
    else:
        pool = mp.Pool(ncpu)
    #configure the check filers function with the filters using partial 
    #check filters for all combos
    results = [ check_filters_mp(c,filters, pool) for c in combos ]
    print(results[0:10])
    triples = [ ( URIfmt(o1 , urlstring) ,predicate , URIfmt(o2, urlstring) ) for o1,o2 in itertools.combinations(objects,2) if results.pop(0) ]
    return triples