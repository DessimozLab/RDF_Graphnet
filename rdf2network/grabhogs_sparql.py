import time, traceback
import sys
from SPARQLWrapper import SPARQLWrapper, JSON
import rdflib
import sys, os, time
import pandas as pd
import stringdb

# always display full column results (don't truncate output)
pd.set_option('display.max_colwidth', None)

# just for reference, check how many of the proteins require multiple batch calls
over_limits = 0

# function to print in a table results of a SPARQL query
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
    df = pd.DataFrame(table, columns=list(header))
    return df

# this is the batch size we want for SPARQL queries
LIMIT = 10000

def query_SPARQL_with_limit_offset(sparql_query, limit, offset, sparql_endpoint):
    global over_limits
    
    query_with_limit_offset = f'{sparql_query} LIMIT {limit} OFFSET {offset}'
    
    #print(query_with_limit_offset)
    
    if(sparql_endpoint is not None):
        sparql_endpoint.setQuery(query_with_limit_offset)
        sparql_endpoint.setReturnFormat(JSON)

        results_OMA = sparql_endpoint.query().convert()

    else:
        print("Error: Configure SPARQL Endpoint")
        return None
    
    if(len(results_OMA["results"]["bindings"]) == LIMIT):
        # Request next set of results until end, append to previous
        if(offset == 0):
            over_limits += 1
            #print("Over 10K results for {}".format(results_OMA["results"]["bindings"][1]))
        return pd.concat([results_to_pandas(results_OMA), query_SPARQL_with_limit_offset(sparql_query, LIMIT, offset + LIMIT, sparql_endpoint)])
   
    return results_to_pandas(results_OMA)

def query_orthologs(uniprot_entry, sparql_endpoint):

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
    
    results_OMA = query_SPARQL_with_limit_offset(query_OMA_orthologs, LIMIT, offset, sparql_endpoint)
    
    return results_OMA
        
def query_paralogs(uniprot_entry, sparql_endpoint):
    
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
    
    results_OMA = query_SPARQL_with_limit_offset(query_OMA_paralogs, LIMIT, offset, sparql_endpoint)
    
    return results_OMA
    
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

# HERE, get HOGS for all vertices in G from oma sparql endpoint

def grab_hogs( G , subg , cross_ref , sparql_endpoint= None):
    start_time = time.time()

    num_queries = 0
    num_errors = 0

    global over_limits 
    over_limits = 0

    proteins_by_species = {}
    # FIRST, need to get UniProt IDs for all string prots needed

    for i,e in enumerate(G.edges(data = True , keys=True)):
        uri_subj = e[0]
        uri_obj = e[1]
        
        cross_ref_prop = rdflib.term.URIRef("http://purl.org/lscr#xrefUniprot")
        
        cross_ref_subj = [o for (s,p,o) in subg.triples((uri_subj, cross_ref, None))][0]
        cross_ref_obj = [o for (s,p,o) in subg.triples((uri_obj, cross_ref, None))][0]
        
        # options: 1. UNIL OMA endpoint; 2. official OMA endpoint;
        
        USE_CASE = 1
        
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
            results_subj_ortho = query_orthologs(cross_ref_subj, sparql_OMA)
            
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
            results_subj_para = query_paralogs(cross_ref_subj, sparql_OMA)
            
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
        
        # get HOG of interacting protein 2, cross_ref_obj (orthologous proteins + species )
        
        #print("Getting HOGs of " + cross_ref_obj)
        cross_ref_obj = "<" + cross_ref_obj + ">"    
        try:
            results_obj_ortho = query_orthologs(cross_ref_obj, sparql_OMA)
            if("orth_protein_uniprot" in results_obj_ortho.columns):
                df1_grouped = results_obj_ortho[["orth_protein_uniprot","taxon_orth"]].groupby('taxon_orth')

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
            #traceback.print_exc()
            num_errors += 1
            pass
        
        try:
            results_obj_para = query_paralogs(cross_ref_obj, sparql_OMA)

            if("para_protein_uniprot" in results_obj_para.columns):
                df1_grouped = results_obj_para[["para_protein_uniprot","taxon_para"]].groupby('taxon_para')

                # iterate over each group
                for group_name, df_group in df1_grouped:
                    taxon_id = group_name.rsplit('/', 1)[1]
                    if(proteins_by_species.get(taxon_id) == None):
                        proteins_by_species[taxon_id] = set()
            
                    for row_index, row in df_group.iterrows():
                        protein_url = row["para_protein_uniprot"]
                        accession = protein_url.rsplit('/', 1)[1]
                        proteins_by_species[taxon_id].add(accession)
                
        except Exception:
            #traceback.print_exc()
            num_errors += 1
            pass
        #results_to_pandas(results_obj)
        num_queries += 4
    end_time = time.time()
    print("Total time for {} SPARQL queries: {} seconds (multiple batch calls in: {} cases)".format(num_queries, end_time - start_time, over_limits))
    print("Num errors: {}".format(num_errors)) 
    return proteins_by_species , results_subj_para , results_subj_ortho , results_obj_para , results_obj_ortho
