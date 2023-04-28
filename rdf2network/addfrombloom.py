
import itertools
import rdflib

#load bloom filters for string data
import pickle
import glob
import multiprocessing as mp


def check_filters(val,filters ):
    for f in filters:
        if val in f[0][0]:
            return True
    return False

def load_filters(filters = './filters/bloomfinal_big*.pkl' ):
    filters = glob.glob(filters)
    filters = [ pickle.load(open(f, 'rb')) for f in filters ]
    #check bloom filters for interactors in diff species
    return filters

def URIfmt(object , urlstring = 'https://string-db.org/network/' ):
    return rdflib.term.URIRef( urlstring + object)

def check_allvall( objects , predicate = rdflib.term.URIRef('https://string-db.org/rdf/high-confidence-cutoff') , urlstring= 'https://string-db.org/network/', filters = None ):
    #replace a call to a server with a check on bloom filters
    if filters is None:
        filters = load_filters()
    triples = [ (URIfmt(o1 , urlstring) ,predicate , URIfmt(o2, urlstring)) for o1,o2 in itertools.combinations(objects,2)
    if check_filters(o1 + '_'+o2 , filters=filters) ]
    return triples



def check_allvall_mp( objects , predicate = rdflib.term.URIRef('https://string-db.org/rdf/high-confidence-cutoff') , urlstring= 'https://string-db.org/network/', ncpu = 10, filters = None ):
    #replace a call to a server with a check on bloom filters
    if filters is None:
        filters = load_filters()
    combos = [  o1+'_'+o2 for o1,o2 in itertools.combinations(objects,2) ]
    if cpu is None:
        pool = mp.Pool(mp.cpu_count())
    else:
        pool = mp.Pool(ncpu)
    #configure the check filers function with the filters
    check_filters_set = lambda x: check_filters(x , filters=filters)
    results = pool.map(check_filters_set, combos)
    triples = [ (URIfmt(o1 , urlstring) ,predicate , URIfmt(o2, urlstring)) for o1,o2 in itertools.combinations(objects,2) if results.pop(0) ]
    return triples


    
    

