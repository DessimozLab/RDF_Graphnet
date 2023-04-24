import torch
from torch_geometric.data import HeteroData , Data
import networkx as nx
from rdflib.namespace import RDF

import rdflib
import sparse




# transform an rdflib graph to a heterodata object
def rdf_to_hetero(rdf_graph , predlinks = None, prednodes = None):
    # create a new heterodata object

    data = HeteroData()
    # assign edge types from the predicate
    edge_types = set([p for s,p,o in rdf_graph])
    edge_types = {e:i for i,e in enumerate(edge_types)}

    # Define the RDF types we want to find
    subject_type = URIRef(RDF.subject)
    object_type = URIRef(RDF.object)

    # Find all subjects and objects and their types
    subject_types = set()
    object_types = set()
    for subject, predicate, obj in g:
        if predicate == subject_type:
            subject_types.add((subject, obj))
        elif predicate == object_type:
            object_types.add((subject, obj))
    for edge_type in edge_types:
        for stype in subject_types:
            for otype in object_types:
                # create a dictionary of nodes
                rows = { tup[0]:i for i,tup in enumerate(rdf_graph.triples(( None , subject_type, stype) ) ) }
                columns = { tup[2]:i for i,tup in enumerate(rdf_graph.triples(( None , object_type, otype) ) ) }
                adj = scipy.sparse.lil_matrix( shape = (len(rows), len(columns)))
                [ adj[rows[s], columns[o]] = 1 for s,p,o in rdf_graph.triples((None, edge_type, None))]
                data[ stype.n3() , p.n3() , otype.n3() ].edge_index = from_scipy_sparse_matrix(adj)
                
    return data


    # create a dictionary of nodes
    
    nodes = {s: {'type': 'subject', 'feature': s} for s in subjects}
    nodes.update({o: {'type': 'object', 'feature': o} for o in objects})


    
    # add edges to the heterodata object for each edge type
        



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
