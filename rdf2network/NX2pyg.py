import torch
from torch_geometric.data import HeteroData , Data
import networkx as nx
from rdflib.namespace import RDF
import rdflib
# create a new heterodata object
from torch_geometric.data import HeteroData , Data
import torch_geometric.transforms as T
import torch_geometric.utils 
from rdflib import RDF , URIRef
import itertools
import numpy as np
import scipy.sparse
import torch


def sparse2pairs(sparsemat, matrows = None):
    '''
    This functino takes a sparse matrix and returns a list of pairs of the non zero entries
    args:
        sparsemat: a sparse matrix
        matrows: a list of the matrix rows to keep
    Returns:    
        a list of pairs of the non zero entries
    '''
    if matrows :
        sparsemat = sparsemat[matrows,:]
        sparsemat = sparsemat[:,matrows]
    sparsemat = scipy.sparse.find(sparsemat)
    return np.vstack([sparsemat[0],sparsemat[1]])

def sparse2pairs(sparsemat, matrows = None):
    '''
    This functino takes a sparse matrix and returns a list of pairs of the non zero entries
    args:
        sparsemat: a sparse matrix
        matrows: a list of the matrix rows to keep
    Returns:    
        a list of pairs of the non zero entries
    '''
    if matrows :
        sparsemat = sparsemat[matrows,:]
        sparsemat = sparsemat[:,matrows]
    sparsemat = scipy.sparse.find(sparsemat)
    return np.vstack([sparsemat[0],sparsemat[1]])

def rdf2hetero(rdf_graph , verbose = True):
    data = HeteroData()
    # assign edge types from the predicate
    edge_types = set([p for s,p,o in rdf_graph])
    edge_types = {e:i for i,e in enumerate(edge_types)}
    uris = {}
    subtypes = {}
    for sub, pred, obj in rdf_graph:
        for o in (sub,obj):
            if isinstance(o, URIRef):
                uri_type = rdf_graph.qname(o).split(":")[0]  # get the namespace prefix
                if uri_type not in uris:
                    uris[uri_type] = []
                    subtypes[uri_type] = set([])
                uris[uri_type].append(o)
                subtypes[uri_type].add( ''.join(o.n3().split('/')[0:-1] )  )
    
    
    allsubtypes = {}
    inputdims = {}

    for t in subtypes:
        #compile the x data matrix for each subtype
        subtype_dict = { ty:i for i,ty in enumerate(subtypes[t])}
        allsubtypes[t] = subtype_dict
        indices = [ subtype_dict[''.join(o.n3().split('/')[0:-1])] for o in uris[t] ]
        x= np.zeros( (len(uris[t]), len(subtype_dict)))
        inputdims[t] = len(subtype_dict)
        x[:,indices] = 1
        #for now this is 1hot for subtype
        #add some more descriptive data here if you have it
        data[t].x = torch.tensor(x , dtype=torch.float )
    
    node_index_by_type = { uritype : { n:i for i,n in enumerate( set(uris[uritype]) ) } for uritype in uris }
    #todo add the diff types of uris within each namespace as 1hot encoded for x attribute of each node
    
    interactions = {}
    for t1,t2 in itertools.product(node_index_by_type,node_index_by_type):
            rows = node_index_by_type[t1]
            columns = node_index_by_type[t2]
            
            if verbose == True:
                print(t1,len(rows),t2, len(columns)) 
            for edge_type in edge_types:
                
                #add the number of input dimension for the rows here
                #for now this is just the subtype 1 hot. theres prob a more elegant way to pass this info along or get it from the data obj
                
                # create a dictionary of nodes
                triples =  [ (s,p,o) for s,p,o in rdf_graph.triples((None, edge_type, None))  if ( s in rows and o in columns )  ]
                if len(triples)>0:

                    interactions[(t1,edge_type.n3(),t2)] = len(allsubtypes[t1])
                    interactions[(t2, 'rev_'+edge_type.n3(),t1)] = len(allsubtypes[t1])


                    adj = scipy.sparse.lil_matrix((len(rows), len(columns)))
                    for s,p,o in triples:
                        adj[rows[s], columns[o]] = 1
                    if verbose == True:
                        print('nadj:', adj.sum())
                        print( t1,len(node_index_by_type[t1]) ,t2, len(node_index_by_type[t2]))
                        print( edge_type , len(triples) )
                    if adj.sum() >0 :
                        #between subgraphs
                        data[ t1 , edge_type.n3() , t2 ].edge_index = torch.tensor( sparse2pairs(adj) ,  dtype=torch.long )
                        if t1 == t2:
                            #within a subgraph of the same namespace
                            torch_geometric.utils.add_self_loops(data[ t1 , edge_type.n3() , t2 ].edge_index )
    
    data = T.ToUndirected()(data)
    data = T.AddSelfLoops()(data)


    
    return data , inputdims , interactions , node_index_by_type



def nx_multigraph_to_heterodata(nx_graph):
    # Create a new HeteroData object
    data = HeteroData()

    # Add nodes to the HeteroData object for each node type
    for node_type in nx_graph.node:
        node_features = torch.tensor([nx_graph.nodes[node]['feature'] for node in nx_graph.nodes if nx_graph.nodes[node]['type'] == node_type], dtype=torch.float)
        data['node_'+node_type].x = node_features

    # Add edges to the HeteroData object for each edge type
    for edge_type in nx_graph.edges:
        edge_list = [(e[0], e[1], d) for e, d in nx_graph.edges[edge_type].items()]
        src, dst, edge_features = zip(*edge_list)
        src = torch.tensor(src)
        dst = torch.tensor(dst)
        edge_index = torch.stack([src, dst])
        edge_index_dst_first = torch.cat([edge_index[1].unsqueeze(0), edge_index[0].unsqueeze(0)], dim=0)
        edge_index = torch.cat([edge_index, edge_index_dst_first], dim=1)
        edge_features = torch.tensor(edge_features, dtype=torch.float)
        data['edge_'+edge_type].edge_index = edge_index
        data['edge_'+edge_type].edge_attr = edge_features
    return data

def hetero_to_multigraph(hetero_data: HeteroData) -> nx.MultiGraph:
    graph = nx.MultiGraph()

    for etype in hetero_data.edge_types:
        src_key, dst_key = hetero_data.edge_index_dict[etype]
        src_idx, dst_idx = src_key.split('_')[0], dst_key.split('_')[0]
        src_node_type, dst_node_type = hetero_data.node_types[src_idx], hetero_data.node_types[dst_idx]

        for i in range(hetero_data[etype].size(1)):
            src, dst = src_key + '_' + str(hetero_data[etype][0, i]), dst_key + '_' + str(hetero_data[etype][1, i])
            edge_data = dict(hetero_data[etype].items())
            edge_data.pop('edge_index')
            graph.add_edge(src, dst, key=etype, attr_dict=edge_data)

    for ntype in hetero_data.node_types:
        num_nodes = hetero_data[ntype].size(0)
        node_data = dict(hetero_data[ntype].items())
        for i in range(num_nodes):
            node_key = ntype + '_' + str(i)
            node_attr = node_data.copy()
            node_attr.pop('x')
            graph.add_node(node_key, attr_dict=node_attr)
    return graph
