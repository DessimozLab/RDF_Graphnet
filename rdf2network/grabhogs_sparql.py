import time, traceback
import sys
from SPARQLWrapper import SPARQLWrapper, JSON , JSONLD
from rdflib import Graph, URIRef
import json
from rdflib.plugins.sparql.results.jsonresults import JSONResultParser
import io
import rdflib
import sys, os, time
import pandas as pd
import stringdb

# always display full column results (don't truncate output)
pd.set_option('display.max_colwidth', None)
# just for reference, check how many of the proteins require multiple batch calls
over_limits = 0

# this is the batch size we want for SPARQL queries
LIMIT = 10000


def results_to_pandas(results):

    # how to transform SPARQL results into Pandas dataframes
    if(len(results["results"]["bindings"]) == 0):
        # this might actually correspond to an error since
        # at least the header (row 0) should be returned
        return pd.DataFrame()
    # get header (column names) from results
    header = results["results"]["bindings"][0].keys()
    # display table of results:
    table = []
    # the SPARQL JSON results to the query are available in the "results", "bindings" entry:
    for entry in results["results"]["bindings"]:
        # append entries from the results to a regular Python list of rows, which we can then transform to a Pandas DF
        row = [entry[column]["value"] if entry.get(column, None) != None else None for column in header]
        table.append(row)
        #add the entry to rdflib graph
    df = pd.DataFrame(table, columns=list(header))
    return df 

    
def get_string_id(uniprot_url, taxon_url):
    accession = uniprot_url.rsplit('/', 1)[1]
    taxon_id = taxon_url.rsplit('/', 1)[1]
    if(uniprot_url not in string_ids_dict):
    # NOTE! Here if no species mentioned, defaults to HUMAN
    # also, if protein not available in STRING, it will through a ValueException error
        try:
            string_ids_dict[uniprot_url] = stringdb.get_string_ids([accession], int(taxon_id))["stringId"]
        except:
            pass



# this is the batch size we want for SPARQL queries
LIMIT = 10000
def results_to_pandas(results):
    # how to transform SPARQL results into Pandas dataframes
    if(len(results["results"]["bindings"]) == 0):
        # this might actually correspond to an error since
        # at least the header (row 0) should be returned
        return pd.DataFrame()
    # get header (column names) from results
    header = results["results"]["bindings"][0].keys()
    # display table of results:
    table = []
    # the SPARQL JSON results to the query are available in the "results", "bindings" entry:
    for entry in results["results"]["bindings"]:
        # append entries from the results to a regular Python list of rows, which we can then transform to a Pandas DF
        row = [entry[column]["value"] if entry.get(column, None) != None else None for column in header]
        table.append(row)
        #add the entry to rdflib graph
    df = pd.DataFrame(table, columns=list(header))
    return df 

def query_SPARQL_with_limit_offset(sparql_query, limit, offset, sparql_endpoint , fmt = JSON):
    global over_limits
    query_with_limit_offset = f'{sparql_query} LIMIT {limit} OFFSET {offset}'
    #print(query_with_limit_offset)
    if(sparql_endpoint is not None):
        sparql_endpoint.setQuery(query_with_limit_offset)
        sparql_endpoint.setReturnFormat(fmt)
        results_OMA = sparql_endpoint.query().convert()
    else:
        print("Error: Configure SPARQL Endpoint")
        return None    
    import pdb ; pdb.set_trace()

    if fmt == JSON  and (len(results_OMA["results"]["bindings"]) == LIMIT):
        # Request next set of results until end, append to previous
        res = results_OMA["results"]["bindings"]
        while len(res) == LIMIT and offset < 5*LIMIT:
            #print("Over 10K results for {}".format(results_OMA["results"]["bindings"][1]))    
            additional_bindings =  query_SPARQL_with_limit_offset(sparql_query, offset , offset + LIMIT, sparql_endpoint )
            offset += LIMIT
            results_OMA["results"]["bindings"] += additional_bindings["results"]["bindings"]
            res = additional_bindings["results"]["bindings"]
    return results_OMA

def query_orthologs(uniprot_entry, sparql_endpoint , return_graph = False):

    offset = 0
    
    query_OMA_orthologs = """
    PREFIX obo: <http://purl.obolibrary.org/obo/>
    PREFIX orth: <http://purl.org/net/orth#>
    PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
    PREFIX up: <http://purl.uniprot.org/core/>
    PREFIX lscr: <http://purl.org/lscr#>
    
    SELECT ?protein1 ?protein1_uniprot ?taxon ?orth_protein ?orth_protein_uniprot ?taxon_orth {
        VALUES ?protein1_uniprot { """ + uniprot_entry + "}" + """    
        ?protein1 a orth:Protein.
        ?protein1 lscr:xrefUniprot ?protein1_uniprot.
        ?protein1 orth:organism/obo:RO_0002162 ?taxon.
        
        ?orth_cluster a orth:OrthologsCluster.
        ?orth_cluster orth:hasHomologousMember ?orth_node1.
        ?orth_cluster orth:hasHomologousMember ?orth_node2.
        ?orth_node1 orth:hasHomologousMember* ?protein1.
        ?orth_node2 orth:hasHomologousMember* ?orth_protein.
        ?orth_protein a orth:Protein.
        ?orth_protein orth:organism/obo:RO_0002162 ?taxon_orth.
        ?orth_protein lscr:xrefUniprot ?orth_protein_uniprot.
        FILTER(?orth_node1 != ?orth_node2)
    }  
    """
    
    if return_graph == True:
        fmt = JSONLD
        data = query_SPARQL_with_limit_offset(query_OMA_orthologs, LIMIT, offset, sparql_endpoint , fmt= fmt)
        graph = Graph()
        graph.parse(data=json.dumps(data), format='json-ld')
        return  graph
    else:       
        fmt = JSON
        results_OMA = query_SPARQL_with_limit_offset(query_OMA_orthologs, LIMIT, offset, sparql_endpoint)
        results_OMA = results_to_pandas(results_OMA)
        return results_OMA
        
def query_paralogs(uniprot_entry, sparql_endpoint , return_graph = False):
    
    offset = 0
    
    query_OMA_paralogs = """
    PREFIX obo: <http://purl.obolibrary.org/obo/>
    PREFIX orth: <http://purl.org/net/orth#>
    PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
    PREFIX up: <http://purl.uniprot.org/core/>
    PREFIX lscr: <http://purl.org/lscr#>
    
    SELECT ?protein1 ?protein1_uniprot ?taxon ?para_protein ?para_protein_uniprot ?taxon_para {
        VALUES ?protein1_uniprot { """ + uniprot_entry + "}" + """    
        ?protein1 a orth:Protein.
        ?protein1 lscr:xrefUniprot ?protein1_uniprot.
        ?protein1 orth:organism/obo:RO_0002162 ?taxon.
        
        ?para_cluster a orth:ParalogsCluster.
        ?para_cluster orth:hasHomologousMember ?para_node1.
        ?para_cluster orth:hasHomologousMember ?para_node2.
        ?para_node1 orth:hasHomologousMember* ?protein1.
        ?para_node2 orth:hasHomologousMember* ?para_protein.
        ?para_protein a orth:Protein.
        ?para_protein orth:organism/obo:RO_0002162 ?taxon_para.
        ?para_protein lscr:xrefUniprot ?para_protein_uniprot.
        FILTER(?para_node1 != ?para_node2)
    }
    """

    if return_graph == True:
        fmt = JSONLD
        data = query_SPARQL_with_limit_offset(query_OMA_paralogs, LIMIT, offset, sparql_endpoint , fmt= fmt)
        graph = Graph()
        graph.parse(data=json.dumps(data), format='json-ld')
        return  graph
    else:       
        fmt = JSON
        results_OMA = query_SPARQL_with_limit_offset(query_OMA_paralogs, LIMIT, offset, sparql_endpoint)
        results_OMA = results_to_pandas(results_OMA)
        return results_OMA
"""
def grab_hogs_graph( subg , cross_ref , sparql_endpoint= None , USE_CASE = 1 , verbose = True , cross_ref_prop = rdflib.term.URIRef("http://purl.org/lscr#xrefUniprot")):
    start_time = time.time()
    num_queries = 0
    num_errors = 0
    global over_limits 
    over_limits = 0
    proteins_by_species = {}
    # FIRST, need to get UniProt IDs for all string prots needed
    subjs = [[o for (s,p,o) in subg.triples((None, cross_ref, None))]]
    subjs = [ s[0] for s in subjs if len(s) > 0 ]
    subjs = set(subjs)
    print(len(subjs))
    orthograph = Graph()
    for cross_ref_subj in subjs:
        # get cross reference for interacting protein 1, cross_ref_subj                
        # options: 1. UNIL OMA endpoint; 2. official OMA endpoint;
        if(USE_CASE == 1):
            sparql_OMA = SPARQLWrapper("https://unil.omabrowser.org/sparql/")
        elif(USE_CASE == 2):
            sparql_OMA = SPARQLWrapper("https://sparql.omabrowser.org/sparql/")
        else:
            sparql_OMA = None
        # get HOG of interacting protein 1, cross_ref_subj (orthologous proteins + species)
        cross_ref_subj = "<" + cross_ref_subj + ">"  
        #print("Getting HOGs of " + cross_ref_subj)
        # try querying the remote server but skip over errors if any
        # NOTE: it seems like errors happen when number of results is too big,
        # e.g. <http://purl.uniprot.org/uniprot/B6JXK8> seems to have 55K paralogs!?
        try:
            g1 = query_orthologs(cross_ref_subj, sparql_OMA , return_graph=True )
            orthograph += g1
        except Exception:
            traceback.print_exc()
            num_errors += 1
            pass
        try:
            results_subj_para , g1 = query_paralogs(cross_ref_subj, sparql_OMA , return_graph = True)
            if("para_protein_uniprot" in results_subj_para.columns):
                df1_grouped = results_subj_para[["para_protein_uniprot","taxon_para"]].groupby('taxon_para')
                orthograph += g1
        except Exception:
            #traceback.print_exc()
            num_errors += 1
            pass
        #results_to_pandas(results_obj)
        num_queries += 2
    end_time = time.time()
    if verbose == True:
        print("Total time for {} SPARQL queries: {} seconds (multiple batch calls in: {} cases)".format(num_queries, end_time - start_time, over_limits))
        print("Num errors: {}".format(num_errors)) 
    return orthograph
"""


def grab_hogs_graph( subg , cross_ref , sparql_endpoint= None , USE_CASE = 1 , verbose = True , cross_ref_prop = rdflib.term.URIRef("http://purl.org/lscr#xrefUniprot")):
    start_time = time.time()
    num_queries = 0
    num_errors = 0
    global over_limits 
    over_limits = 0
    proteins_by_species = {}
    # FIRST, need to get UniProt IDs for all string prots needed
    subjs = [[o for (s,p,o) in subg.triples((None, cross_ref, None))]]
    subjs = [ s[0] for s in subjs if len(s) > 0 ]
    subjs = set(subjs)
    print(len(subjs))
    orthograph = Graph()
    for cross_ref_subj in subjs:
        # get cross reference for interacting protein 1, cross_ref_subj                
        # options: 1. UNIL OMA endpoint; 2. official OMA endpoint;
        if(USE_CASE == 1):
            sparql_OMA = SPARQLWrapper("https://unil.omabrowser.org/sparql/")
        elif(USE_CASE == 2):
            sparql_OMA = SPARQLWrapper("https://sparql.omabrowser.org/sparql/")
        else:
            sparql_OMA = None
        # get HOG of interacting protein 1, cross_ref_subj (orthologous proteins + species)
        cross_ref_subj = "<" + cross_ref_subj + ">"  
        #print("Getting HOGs of " + cross_ref_subj)
        # try querying the remote server but skip over errors if any
        # NOTE: it seems like errors happen when number of results is too big,
        # e.g. <http://purl.uniprot.org/uniprot/B6JXK8> seems to have 55K paralogs!?
        g1 = query_orthologs(cross_ref_subj, sparql_OMA , return_graph=True )
        orthograph += g1
        g1 = query_paralogs(cross_ref_subj, sparql_OMA , return_graph = True)
        orthograph += g1
        #results_to_pandas(results_obj)
        num_queries += 2
    end_time = time.time()
    if verbose == True:
        print("Total time for {} SPARQL queries: {} seconds (multiple batch calls in: {} cases)".format(num_queries, end_time - start_time, over_limits))
        print("Num errors: {}".format(num_errors)) 
    return orthograph

#get the OMA orthology info for the string subgraph
def grab_hogs( subg , cross_ref , sparql_endpoint= None , USE_CASE = 1 , verbose = True , cross_ref_prop = rdflib.term.URIRef("http://purl.org/lscr#xrefUniprot")):
    start_time = time.time()
    num_queries = 0
    num_errors = 0
    global over_limits 
    over_limits = 0
    proteins_by_species = {}
    # FIRST, need to get UniProt IDs for all string prots needed
    subjs = [ o for (s,p,o) in subg.triples((None, cross_ref, None))]
    subjs = set(subjs)
    print(len(subjs))
    for cross_ref_subj in subjs:
        # get cross reference for interacting protein 1, cross_ref_subj                
        # options: 1. UNIL OMA endpoint; 2. official OMA endpoint;
        
        if(USE_CASE == 1):
            sparql_OMA = SPARQLWrapper("https://unil.omabrowser.org/sparql/")
        elif(USE_CASE == 2):
            sparql_OMA = SPARQLWrapper("https://sparql.omabrowser.org/sparql/")
        else:
            sparql_OMA = None
        # get HOG of interacting protein 1, cross_ref_subj (orthologous proteins + species)
        cross_ref_subj = "<" + cross_ref_subj + ">"  
        #print("Getting HOGs of " + cross_ref_subj)
        # try querying the remote server but skip over errors if any
        # NOTE: it seems like errors happen when number of results is too big,
        # e.g. <http://purl.uniprot.org/uniprot/B6JXK8> seems to have 55K paralogs!?
        try:
            results_subj_ortho  = query_orthologs(cross_ref_subj, sparql_OMA, return_graph = False  )
            # fill dictionary with accessions to look for by taxon_id
            # this will be used to query STRING API (batched by species)
            if("orth_protein_uniprot" in results_subj_ortho.columns):
                df1_grouped = results_subj_ortho[["orth_protein_uniprot","taxon_orth"]].groupby('taxon_orth')
                # iterate over each group
                for group_name, df_group in df1_grouped:
                    taxon_id = group_name.rsplit('/', 1)[1]
                    if(proteins_by_species.get(taxon_id) == None):
                        proteins_by_species[taxon_id] = set()
                    #print('group {}'.format(taxon_id))
                    for row_index, row in df_group.iterrows():
                        protein_url = row["orth_protein_uniprot"]
                        accession = protein_url.rsplit('/', 1)[1]
                        proteins_by_species[taxon_id].add(accession)
        except Exception:
            traceback.print_exc()
            num_errors += 1
            pass
        try:
            results_subj_para  = query_paralogs(cross_ref_subj, sparql_OMA , return_graph = False)
            if("para_protein_uniprot" in results_subj_para.columns):
                df1_grouped = results_subj_para[["para_protein_uniprot","taxon_para"]].groupby('taxon_para')
                # iterate over each group
                for group_name, df_group in df1_grouped:
                    taxon_id = group_name.rsplit('/', 1)[1]
                    if(proteins_by_species.get(taxon_id) == None):
                        proteins_by_species[taxon_id] = set() 
                    #print('group {}'.format(taxon_id))
            
                    for row_index, row in df_group.iterrows():
                        protein_url = row["para_protein_uniprot"]
                        accession = protein_url.rsplit('/', 1)[1]
                        proteins_by_species[taxon_id].add(accession)
        except Exception:
            #traceback.print_exc()
            num_errors += 1
            pass
        #results_to_pandas(results_obj)
        num_queries += 2
    end_time = time.time()
    if verbose == True:
        print("Total time for {} SPARQL queries: {} seconds (multiple batch calls in: {} cases)".format(num_queries, end_time - start_time, over_limits))
        print("Num errors: {}".format(num_errors)) 
    return proteins_by_species , results_subj_para , results_subj_ortho 
