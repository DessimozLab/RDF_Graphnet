# first we need to define the SPARQL endpoints of each source, to use later in the protocols
import sys
from SPARQLWrapper import SPARQLWrapper, JSON
import sys, os, time
import pandas as pd
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
    