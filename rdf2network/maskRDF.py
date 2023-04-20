import rdflib

#some functions to remove particular triples


#select a particular predicate obeject pair and remove all triples with that predicate object pair
def remove_by_predicate_object(rg, predicate, object):
    for s,p,o in rg.triples((None,predicate,object)):
        rg.remove((s,p,o))

#sample all the predicte object pairs for a given subject
def sample_by_subject(rg, subject):
    for s,p,o in rg.triples((subject,None,None)):
        yield p,o

