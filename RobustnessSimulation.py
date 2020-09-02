#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 09:42:59 2019

@author: henriette
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 17:12:11 2019

@author: henriette
"""


import psycopg2
import networkx as nx
import pandas as pd
import numpy as np
import random 

# load networks 

maxdist = {'linear_02' : 250, 'random_02' : 600, 'clustr_02' : 375, 'linear_01' : 400, 'random_01' : 900, 'clustr_01' : 650}    


def loadGraph(style, square, habitat, maxdist, streamnetwork='025025050_'):
    
    
    conn = psycopg2.connect("host=139.14.20.252 port=5432 dbname=DB_PhD user=... password=...")
    cursor = conn.cursor()
    cursor.execute('SELECT ids, ids_org FROM dis_pts_2500_10x10_' + style + '.pts_habitat_red_' + str(square) + '_start_' + str(habitat))
     
    ids = cursor.fetchall()
    ids = [i[0] for i in ids]

    conn = psycopg2.connect("host=139.14.20.252 port=5432 dbname=DB_PhD user=... password=...")
    cursor = conn.cursor()
    cursor.execute('SELECT ids, ids_org FROM dis_pts_2500_10x10_' + style + '.pts_habitat_red_' + str(square) + '_start_' + str(habitat))
    
    ids_org = cursor.fetchall()
    ids_org = [i[0] for i in ids_org]

    cursor = conn.cursor()
    cursor.execute('SELECT start, aim, costs FROM stream_network_' + streamnetwork + style + '.habitats_shortpath_red_nlmr_testarea_50x50_0_'+str(square)+'_start_'+str(habitat))

    arcs =  cursor.fetchall()
    
    arcs = [list(x) for x in arcs]
            
    G = nx.Graph()
    G.add_nodes_from(ids) 

    # populate all habitats
    for v in G.nodes():
        G.node[v]['pop'] = 1    
    
    for e in arcs:
        if e[2] < maxdist:
            G.add_edge(e[0],e[1], weight = e[2])
    
    return G
    



###
# # simulate robustness
###
def extinction(G, extb):
    Cliquesize = nx.node_clique_number(G) 
    
    for v in G.nodes():    
        cl = Cliquesize[v]
        val = random.uniform(0,1)
        if val < extb**(1-cl):   # 2**(-cl/2)
            G.node[v]['pop'] = 0 
    return G

    
def colonisation(G, sigma):
            
    dvw = dict(nx.shortest_path_length(G, weight = 'weight'))
    disp = {}
    Kernel = {}
    for v in G.nodes():
        disp[v] = {}
        if G.node[v]['pop'] == 1:
            Kernel = {}
            kernels = 0
            for w in G.nodes():
                if G.node[w]['pop'] == 0:
                    try:
                        ker = 1/(np.sqrt(2*np.pi) * sigma) *np.exp(-(dvw[v][w])**2/(2*sigma**2))
#                        if ker > 0.001:
                        Kernel[w] = ker
#                        else:
#                            Kernel[w] = 0
                        kernels += ker
                    except:
                        1 == 1
            if kernels > 0:
                for w in Kernel:
                    disp[v][w] = Kernel[w]/kernels
    for v in disp:
        for w in disp[v]:
            val = random.uniform(0,1)
            if val < disp[v][w]:
                G.node[w]['pop'] = 1

    return G
  
def randomPerturbation(G, phi, extb, sigma):
    torem = []
    for v in G.nodes():
        if random.uniform(0,1) < phi: 
            torem.append(v)
    G.remove_nodes_from(torem)
    popsize = G.number_of_nodes()
    if popsize == 0:
        return(0)

    else:
        oldnumbers = 0
        
        while popsize != oldnumbers: 
            oldnumbers = popsize                            

            G = extinction(G, extb)
            G = colonisation(G, sigma)
        
            popsize = 0
            for v in G.nodes():
                if G.node[v]['pop'] == 1:
                    popsize += 1 
    
        return(popsize/float(G.number_of_nodes()))

        
        
streamnetwork = '025025050_'
Wdh = 10
Phi = np.linspace(0,1, 100) 
Summen = {}

Perturbation = 'random'
addEdges = 'Optimal'

Squares = 26
Habitats = 10
for nEdges in [2,5,10,15,20,25,30]:
    for style in ['clustr_01','linear_01', 'clustr_01']:

        for extb in [2,5,9]:
            for sigma in [2,5,9]:
                print(Perturbation, extb, sigma, style, addEdges, nEdges)
                for square in range(1, Squares):
                    for habitat in range(Habitats):
                        print(square, habitat)
                        Gc = loadGraph(style, square, habitat, maxdist[style], streamnetwork)
                  
                        # add edges                           
                        if addEdges == 'LazyGreedy':
                            Edges = pd.read_csv('/LazyGreedy/' + style + '/' +str(nEdges) + 'edges/'
                                                + style + '_' + str(square) + '_' + str(habitat) + '.csv', delimiter = ';')
                           
                        if addEdges == 'Greedy':
                            Edges = pd.read_csv('/Greedy/' +style + '/' + str(nEdges) 
                                                + 'edges/edges_' + str(square) + '_' + str(habitat) + '.csv', delimiter = ';')
                        if addEdges == 'Random':
                            Edges = pd.read_csv('/RandomEdges/' +style + '/' + str(nEdges) 
                                                + 'edges/edges_' + str(square) + '_' + str(habitat) + '.csv', delimiter = ';')
                       
                        Edges = [[Edges['start'][i],Edges['aim'][i]] for i in range(nEdges)]
                        Gc.add_edges_from(Edges)
                        
                        Summen[(square, habitat)] = []
                        for wdh in range(Wdh):
                            summe = 0
                            Popsize = []
                            for phi in Phi:
                                G = Gc.copy()
                                if Perturbation == 'random':
                                    popsizephi = randomPerturbation(G, phi, extb, sigma)
                               
                                Popsize.append(popsizephi)
                            RI = 1.0/100 * (sum(Popsize) - 0.5 * (Popsize[0] + Popsize[99]))    
                            Summen[(square, habitat)].append(RI)
                
                if addEdges == 'LazyGreedy':
                    file = open('/LazyGreedy/' + Perturbation + '_' + str(nEdges) + 'edges_' + style + '_' + str(extb) + str(sigma) + '.csv', 'w')
                if addEdges == 'Greedy':
                    file = open('/Greedy/' + Perturbation + '_' + str(nEdges) + 'edges_' + style + '_' + str(extb) + str(sigma) + '.csv', 'w')
                if addEdges == 'Random':
                    file = open('/RandomEdges/' + Perturbation + '_' + str(nEdges) + 'edges_' + style + '_' + str(extb) + str(sigma) + '.csv', 'w')
            
                file.write('Square;Habitatselection')
            
                for wdh in range(Wdh):
                    file.write(';' + str(wdh))
                file.write('\n')
                
                for square in range(1,Squares):
                    for habitat in range(Habitats):
                
                        file.write(str(square) + ';' + str(habitat))
                        for wdh in range(Wdh):
                            file.write(';' + str(Summen[(square, habitat)][wdh]))
                        file.write('\n')
                file.close()
        