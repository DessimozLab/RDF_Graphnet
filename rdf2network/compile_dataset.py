#!/usr/bin/env python
# coding: utf-8

# In[2]:


# Install required packages in the current Jupyter kernel
# Uncomment the following lines if you need to install these libraries
# If you run into permission issues, try with the --user option
#import sys
#!pip install -q rdflib networkx matplotlib
#!{sys.executable} -m pip install rdflib networkx matplotlib pandas stringdb --user


# In[7]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')
import sampler
import grabhogs_sparql
import map2string_fast
import addfrombloom
import rdflib
import SPARQLWrapper
import colour
import itertools
import networkx as nx
import glob
from matplotlib import pyplot as plt


# In[4]:


datapath = '/work/FAC/FBM/DBC/cdessim2/default/dmoi/datasets/STRING/rdf/protein.links.rdf.v11.5/402676.protein.links.rdf.v11.5.txt.gz'
# RDF graph loading
rg = sampler.load_graph(datapath)

datapath2 = '/work/FAC/FBM/DBC/cdessim2/default/dmoi/datasets/STRING/rdf/protein.info.rdf.v11.5/402676.protein.info.rdf.v11.5.txt.gz'
# RDF graph loading
rg_info  = sampler.load_graph(datapath2)


# In[10]:


links = '/work/FAC/FBM/DBC/cdessim2/default/dmoi/datasets/STRING/rdf/protein.links.rdf.v11.5/*.protein.links.rdf.v11.5.txt.gz'
linkfiles = glob.glob(links)
linkfiles = { l:{ 'links':l , 'info':l.replace('protein.links' , 'protein.info' ) } for l in linkfiles}
print(len(linkfiles ))


# In[11]:


subjs = rg.subjects( unique = True)
seed = next(subjs)
print(seed)


# In[12]:


subg = sampler.sample( rg = rg , seed = seed,  layer_limit= 2 , sample_run = 20 )
print(set([p for p in subg.predicates()]))
print("rdflib Graph sampled successfully with {} triples".format(len(subg)))


# In[13]:


subg = sampler.add_xrefs( rg_info , subg )
print("rdflib Graph annotated successfully with {} triples".format(len(subg)))
print(set([p for p in subg.predicates()]))
cross_ref = rdflib.term.URIRef("http://purl.org/lscr#xrefUniprot")
i = 0 
for s,p,o in subg.triples((None, cross_ref, None)):
    print(s,p,o)
    i+= 1
    if i > 10:
        break


# In[14]:


#proteins_by_species , results_subj_para , results_subj_ortho  = grabhogs_sparql.grab_hogs( subg , cross_ref = rdflib.term.URIRef("http://purl.org/lscr#xrefUniprot") )


# In[15]:


#print(results_subj_ortho)


# In[16]:


orthograph =  grabhogs_sparql.grab_hogs_graph( subg , cross_ref , sparql_endpoint= None
            , USE_CASE = 1 , verbose = True , cross_ref_prop = rdflib.term.URIRef("http://purl.org/lscr#xrefUniprot"))


# In[17]:


from rdflib import Graph, URIRef
print("rdflib Graph annotated successfully with {} triples".format(len(orthograph)))
#get all species
taxa = [ 'protein1_uniprot_taxon_orth' , 'protein1_uniprot_taxon_para']
species = set( [ o for t in taxa  for s,p,o in orthograph.triples((None, URIRef(t), None))  ])


# In[18]:


#get all the results for that species
taxa = [ 'orth_protein_uniprot_taxon_orth' , 'para_protein_uniprot_taxon_para']
#get all proteins for each species
prots_by_species = { spec: set([s  for s,p,o in orthograph.triples((None, URIRef(t), spec ))]) for t in taxa for spec in species  }
prots_by_species = { spec:prots_by_species[spec] for spec in prots_by_species if len(prots_by_species[spec])   }
prots_by_species = { spec:[ p.replace('https://string-db.org/network/' , '' ) for p in prots_by_species[spec] ] for spec in prots_by_species }


# In[19]:


server = 'dna065'
ortho_xrefgraph = map2string_fast.mapall(prots_by_species , serverurl= "http://"+server+":3030/string_fuseki/sparql" , retgraph = True)


# In[20]:


orthograph += ortho_xrefgraph
subg += orthograph
#we have interactions for one species and ortho info to all others


# In[21]:


filters = addfrombloom.load_filters()


# In[22]:


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


# In[23]:


#halelujah we have a graph with everything in it
#serialize to turtle format
v = subg.serialize(format="ttl")
with open('testgraph.ttl', 'w') as graphout:
    graphout.write(v)


# In[25]:


print(len(subg))
readg = Graph()
readg.parse('testgraph.ttl')
print(len(readg))
#we can save the subgraphs in rdf...


# In[ ]:


G = rdflib.extras.external_graph_libs.rdflib_to_networkx_multidigraph( subg , lambda s, p, o: {'data':{'key':p  , 'weight':1} } )


# In[ ]:


# example usage
uniprot_id = "http://purl.uniprot.org/uniprot/A0A3B5R6M3" 
serverurl = '
results =  map2string_fast.map2string_SPARQL(uniprot_id , serverurl = )
print(results)


# In[ ]:


preds = [ p for p in rg.predicates( unique = True)]
red = colour.Color('red')
blue = colour.Color('blue')
c = [ c.hex_l for c in  list(red.range_to(blue, len(preds))) ]
colors = { p:c[i] for i,p in enumerate(preds)}
delta = 2/len(preds)
curve = { p:delta*i for i,p in enumerate(preds) }
style = itertools.cycle([ '-', '--' ])
line_style = { p:next(style) for p in preds }


# In[ ]:


pos = nx.circular_layout( G )
f = plt.figure()
f.set_figwidth(15)
f.set_figheight(15)
plt.plot()
#plot the whole mess
ax = plt.gca()
for e in G.edges(data = True):
    ax.annotate("",
                xy=pos[e[0]], xycoords='data',
                xytext=pos[e[1]], textcoords='data',
                arrowprops=dict(arrowstyle="-", color=colors[e[2]['data']['key']],
                                shrinkA=5, shrinkB=5, lw = e[2]['data']['weight'], ls = line_style[e[2]['data']['key']],
                                patchA=None, patchB=None, alpha = e[2]['data']['weight'],
                                connectionstyle="arc3,rad=rrr".replace('rrr',str(0.3*curve[e[2]['data']['key']])
                                ),
                                ),
                )

nx.draw_networkx_nodes(G, pos, node_color = 'b', node_size = 150, alpha = .5)
labels=nx.draw_networkx_labels(G , pos = pos , font_size= 15 , font_color='w')

plt.axis('off')
plt.show()


# In[ ]:





# In[ ]:


#map string neighbours to OMA entries

#jump a few steps in HOGs

#fish for STRING interactions in other species



