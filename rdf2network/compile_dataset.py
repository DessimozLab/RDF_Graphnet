import sampler
import grabhogs_sparql
import map2string_fast
import addfrombloom
import rdflib
import networkx as nx
import glob
from matplotlib import pyplot as plt
from rdflib import Graph, URIRef


if __name__ == "__main__":        
    

    #parse the command line arguments
    import argparse
    parser = argparse.ArgumentParser(description='Compile a dataset from STRING RDF')
    parser.add_argument('--taxa', type=str, nargs='+', help='taxa to include in the dataset')
    parser.add_argument('--proteins', type=str, nargs='+', help='proteins to include in the dataset')
    parser.add_argument('--serverurl', type=str, help='URL of the SPARQL endpoint to use for mapping Uniprot IDs to STRING IDs')
    parser.add_argument('--layer_limit', type=int, help='maximum number of layers to include in the dataset')
    parser.add_argument('--sample_run', type=int, help='number of times to sample the network')
    parser.add_argument('--sample_size', type=int, help='number of nodes to sample from the network')
    parser.add_argument('--nseeds', type=str, help='seed node to use for sampling the network')
    parser.add_argument('--output', type=str, help='output dir to write the dataset to')
    parser.add_argument('--bloom', type=str, help='path to bloom filter to use for sampling the network')
    
    args = parser.parse_args()
    
    taxa = args.taxa
    proteins = args.proteins
    serverurl = args.serverurl
    layer_limit = args.layer_limit
    sample_run = args.sample_run
    sample_size = args.sample_size
    nseeds = args.nseeds
    output = args.output
    bloom = args.bloom

    #if no taxa or proteins are specified, use all of them
    if taxa is None:
        links = '/work/FAC/FBM/DBC/cdessim2/default/dmoi/datasets/STRING/rdf/protein.links.rdf.v11.5/*.protein.links.rdf.v11.5.txt.gz'
        linkfiles = glob.glob(links)
        linkfiles = { l:{ 'links':l , 'info':l.replace('protein.links' , 'protein.info' ) } for l in linkfiles}
        print(len(linkfiles ))
    else:
        linkfiles = []
        taxa = taxa.split(',')
        for t in taxa:
            links = '/work/FAC/FBM/DBC/cdessim2/default/dmoi/datasets/STRING/rdf/protein.links.rdf.v11.5/{}.protein.links.rdf.v11.5.txt.gz'.format(t)
            linkfiles += glob.glob(links)
        linkfiles = { l:{ 'links':l , 'info':l.replace('protein.links' , 'protein.info' ) } for l in linkfiles}
        print(len(linkfiles ))
    #load string bloom filters
    if bloom:
        filters = addfrombloom.load_filters(bloom)
    else:
        filters = addfrombloom.load_filters(bloom)

    for file in linkfiles:
        #go through all the link files and sample some networks
        datapath = linkfiles[file]['links']
        tax = datapath.split('/')[-1].split('.')[0]
        datapath2 = linkfiles[file]['info']
        rg = sampler.load_graph(datapath)
        rg_info  = sampler.load_graph(datapath2)
        print( 'loaded string graph with {} triples'.format(len(rg)) )
        print( 'loaded xref graph with {} triples'.format(len(rg_info)) )
        subjs = rg.subjects( unique = True)
        for iseed in range(nseeds):
            seed = next(subjs)
            print('seed:',seed)
            subg = sampler.sample( rg = rg , seed = seed,  layer_limit= 2 , sample_run = 20 )
            print(set([p for p in subg.predicates()]))
            print("rdflib Graph sampled successfully with {} triples".format(len(subg)))
            subg = sampler.add_xrefs( rg_info , subg )
            print("sample Graph annotated successfully with {} triples".format(len(subg)))
            print(set([p for p in subg.predicates()]))
            cross_ref = rdflib.term.URIRef("http://purl.org/lscr#xrefUniprot")
            orthograph =  grabhogs_sparql.grab_hogs_graph( subg , cross_ref , sparql_endpoint= None
                        , USE_CASE = 1 , verbose = True , cross_ref_prop = rdflib.term.URIRef("http://purl.org/lscr#xrefUniprot"))
            print("rdflib Graph annotated successfully with {} triples".format(len(orthograph)))
            #get all species
            taxa = [ 'protein1_uniprot_taxon_orth' , 'protein1_uniprot_taxon_para']
            species = set( [ o for t in taxa  for s,p,o in orthograph.triples((None, URIRef(t), None))  ])
            taxa = [ 'orth_protein_uniprot_taxon_orth' , 'para_protein_uniprot_taxon_para']
            #get all proteins for each species
            prots_by_species = { spec: set([s  for s,p,o in orthograph.triples((None, URIRef(t), spec ))]) for t in taxa for spec in species  }
            prots_by_species = { spec:prots_by_species[spec] for spec in prots_by_species if len(prots_by_species[spec])   }
            prots_by_species = { spec:[ p.replace('https://string-db.org/network/' , '' ) for p in prots_by_species[spec] ] for spec in prots_by_species }
            #"http://"+server+":3030/string_fuseki/sparql"
            ortho_xrefgraph = map2string_fast.mapall(prots_by_species , serverurl=serverurl  , retgraph = True)
            orthograph += ortho_xrefgraph
            subg += orthograph
            #lets add the interactions for all using the bloom filters
            print( len(orthograph ) , set([p for p in orthograph.predicates() ]) )
            #get string ids by species
            pred = rdflib.term.URIRef('http://purl.org/lscr#xrefUniprot')
            interactions = []
            for spec in prots_by_species:
                stringids = [ s for prot in prots_by_species[spec] for s,p,o in orthograph.triples((None, pred , prot )) ]
                #stringids = [ s.replace('https://string-db.org/network/' , '' ) for s in stringids ]
                if len(stringids ) > 2 :
                    interactions += addfrombloom.check_allvall( objects = stringids , urlstring = 'https://string-db.org/network/' , filters = filters )
            [subg.add(t) for t in interactions]
            #halelujah we have a graph with everything in it
            #serialize to turtle format
            v = subg.serialize(format="ttl")
            #check for the output dir
            if not os.path.exists(output):
                os.makedirs(output)
            #use the seed to create the filename
            filename = output + '/tax_{}_seed_{}.ttl'.format(seed.split('/')[-1] , tax )
            with open(filename, 'w') as graphout:
                graphout.write(v)
