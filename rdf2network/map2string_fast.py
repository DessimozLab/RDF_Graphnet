# SPEED UP for cross-refs Uniprot - STRING:
# Using Fuseki server with local TDB2 database for cross-refs only
import time
from SPARQLWrapper import SPARQLWrapper, JSON , JSONLD

error_taxa = set()
error_proteins = set()
def map2string_SPARQL(uniprot_id , serverurl , fmt = JSON):
    sparql_STRING = SPARQLWrapper(serverurl)
    get_string_id_for_gene = """
    PREFIX lscr: <http://purl.org/lscr#>
    SELECT DISTINCT ?string_id where {
        ?string_id lscr:xrefUniprot \"""" + uniprot_id  + """\"
    } LIMIT 1
    
    """
    
    sparql_STRING.setQuery(get_string_id_for_gene)
    sparql_STRING.setReturnFormat(JSON)
    json_results = sparql_STRING.query().convert()
    
    if(len(json_results["results"]["bindings"]) > 0):
        return json_results["results"]["bindings"][0]["string_id"]["value"]
    error_proteins.add(uniprot_id)
    return None


def mapall(proteins_by_species, serverurl= "http://dna081:3030/string_fuseki/sparql" , verbose = False):
    results_string_per_protein = {}
    num_api_errors = 0
    
    start_time = time.time()
    found_in_api = set()
    for taxon_id in proteins_by_species.keys():
        for protein in proteins_by_species[taxon_id]:
            try:
                string_id = map2string_SPARQL("http://purl.uniprot.org/uniprot/" + protein , serverurl= serverurl)
                if(string_id is not None):
                    results_string_per_protein[protein] = string_id
                    continue
                # invoking the STRING API is costly, only do this for validation:
                #try:
                #    string_id = stringdb.get_string_ids([protein], int(taxon_id))
                #    if not string_id.empty:
                #        results_string_per_protein[protein] = string_id
                #        found_in_api.add(protein)
                #except Exception:
                #    pass
            except Exception:
                #print(proteins_by_species[taxon_id])
                #print("Error computing IDs for taxon {}".format(taxon_id))
                traceback.print_exc()
                error_taxa.add(taxon_id)
                error_proteins.add(protein)
                num_api_errors += 1
                pass
    end_time = time.time()
    if verbose == True:
        print("Total time in STRING SPARQL calls: {} seconds. Number of results {} with {} errors".format(end_time - start_time,len(results_string_per_protein), num_api_errors))
    return results_string_per_protein

def mapall_graph(proteins_by_species, serverurl= "http://dna081:3030/string_fuseki/sparql" , verbose = False):
    results_string_per_protein = {}
    num_api_errors = 0
    
    start_time = time.time()
    found_in_api = set()
    for taxon_id in proteins_by_species.keys():
        for protein in proteins_by_species[taxon_id]:
            try:
                string_id = map2string_SPARQL("http://purl.uniprot.org/uniprot/" + protein , serverurl= serverurl)
                if(string_id is not None):
                    results_string_per_protein[protein] = string_id
                    continue
                # invoking the STRING API is costly, only do this for validation:
                #try:
                #    string_id = stringdb.get_string_ids([protein], int(taxon_id))
                #    if not string_id.empty:
                #        results_string_per_protein[protein] = string_id
                #        found_in_api.add(protein)
                #except Exception:
                #    pass
            except Exception:
                #print(proteins_by_species[taxon_id])
                #print("Error computing IDs for taxon {}".format(taxon_id))
                traceback.print_exc()
                error_taxa.add(taxon_id)
                error_proteins.add(protein)
                num_api_errors += 1
                pass
    end_time = time.time()
    if verbose == True:
        print("Total time in STRING SPARQL calls: {} seconds. Number of results {} with {} errors".format(end_time - start_time,len(results_string_per_protein), num_api_errors))
    return results_string_per_protein