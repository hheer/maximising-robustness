#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 16:51:14 2019

@author: henriette
"""
import networkx as nx
import pandas as pd
import psycopg2
import random


def maxClustering(G):        
    deg = nx.degree(G)
    cl = {}
    cl0 = nx.average_clustering(G) # clustering of original graph
    cl[0,0] = cl0

    # go through each pair of nodes 
    for u in G.nodes():
        for v in G.nodes():
            if int(u) < int(v) and (u,v) not in G.edges() and (v,u) not in G.edges():
                clneu = 0
                cn = list(nx.common_neighbors(G,u,v))
                k = len(cn)
                if k != 0:
                    # add difference of clustering for neighbours
                    for w in cn:
                        dw = deg[w]
                        clneu += 2/float(dw * (dw-1))
                        
                    # add difference of clustering for u,v
                    if deg[u] > 1:
                        clneu += (2*k*(deg[u]-1) - 4*nx.triangles(G,u))/float(deg[u]*(deg[u]**2 - 1))
                    else: 
                        clneu += 1
                    if deg[v] > 1:
                        clneu += (2*k*(deg[v]-1) - 4*nx.triangles(G,v))/float(deg[v]*(deg[v]**2 - 1))
                    else:
                        clneu += 1
                        
                cl[u,v] = clneu/float(G.number_of_nodes()) + cl0
  
    # find biggest clustering 
    
    m = max(cl, key = cl.get)
    U = m[0]
    V = m[1]

    G.add_edge(U,V)
    clmax = nx.average_clustering(G)
    G.remove_edge(U,V)
    
    return U,V,clmax

    
def Greedy(G):
    u1,v1,clmax = maxClustering(G)
    G.add_edge(u1,v1)
    u2,v2,clmax = maxClustering(G)
    G.remove_edge(u1,v1)
    
    return u1,v1,u2,v2, clmax

    
def LazyGreedy(G):        
    deg = nx.degree(G)
    cl = {}
    cl0 = nx.average_clustering(G) # clustering of original graph
    cl[(0,0)] = cl0
    cl[(-1,-1)] = cl0

    # go through each pair of nodes 
    for u in G.nodes():
        for v in G.nodes():
            if int(u) < int(v) and (u,v) not in G.edges() and (v,u) not in G.edges():
                clneu = 0
                cn = list(nx.common_neighbors(G,u,v))
                k = len(cn)
                if k != 0:
                    # add difference of clustering for neighbours
                    for w in cn:
                        dw = deg[w]
                        clneu += 2/float(dw * (dw-1))
                        
                    # add difference of clustering for u,v
                    if deg[u] > 1:
                        clneu += (2.0*k*(deg[u]-1) - 4.0*nx.triangles(G,u))/float(deg[u]*(deg[u]**2 - 1))
                    else: 
                        clneu += 1
                    if deg[v] > 1:
                        clneu += (2.0*k*(deg[v]-1) - 4.0*nx.triangles(G,v))/float(deg[v]*(deg[v]**2 - 1))
                    else:
                        clneu += 1
                        
                cl[u,v] = clneu/float(G.number_of_nodes()) + cl0

    # find biggest clustering 
    m = max(cl, key = cl.get)
    u1 = m[0]
    v1 = m[1]

    del cl[m]
    m = max(cl, key = cl.get)
    u2 = m[0]
    v2 = m[1]
    G.add_edges_from([(u1,v1), (u2,v2)])
    clmax = nx.average_clustering(G)
    G.remove_edges_from([(u1,v1), (u2,v2)])
    
    return u1, v1, u2, v2, clmax


def optimal(G):
    cl = {}
    cl0 = nx.average_clustering(G)
    cl[(0,0,0,0)] = cl0
    deg = dict(nx.degree(G))

    for u1 in G.nodes():
        for v1 in G.nodes():
        
            if int(u1) < int(v1) and (u1,v1) not in G.edges():
                clneu = 0
                cn = list(nx.common_neighbors(G,u1,v1))
                k = len(cn)
                if k != 0:
                    for w in cn:
                        dw = deg[w]
                        clneu += 2/float(dw * (dw-1))
                        
                    # add difference of clustering for u
                    if deg[u1] > 1:
                        clneu += (2*k*(deg[u1]-1) - 4*nx.triangles(G,u1))/float(deg[u1]*(deg[u1]**2 - 1))
                    else: 
                        clneu += 1
                    if deg[v1] > 1:
                        clneu += (2*k*(deg[v1]-1) - 4*nx.triangles(G,v1))/float(deg[v1]*(deg[v1]**2 - 1))
                    else:
                        clneu += 1
                        
                    G.add_edge(u1,v1)
                    deg[u1] += 1
                    deg[v1] += 1
                       
                    for u2 in G.nodes():
                        for v2 in G.nodes():
                            if int(u2) < int(v2) and (u1,v1) != (u2,v2) and (u1,v1) != (v2,u2):
                                if (u2,v2) not in G.edges():
                                    
                                    clneu2 = 0
                                    cn2 = list(nx.common_neighbors(G,u2,v2))
                                    k2 = len(cn2)
                                    if k2 != 0:
                                        for w in cn2:
                                            dw = deg[w]
                                            clneu2 += 2/float(dw * (dw-1))
                                            
                                        # add difference of clustering for u
                                        if deg[u2] > 1:
                                            clneu2 += (2*k2*(deg[u2]-1) - 4*nx.triangles(G,u2))/float(deg[u2]*(deg[u2]**2 - 1))
                                        else: 
                                            clneu2 += 1
                                        if deg[v2] > 1:
                                            clneu2 += (2*k2*(deg[v2]-1) - 4*nx.triangles(G,v2))/float((deg[v2]*(deg[v2]**2 - 1)))
                                        else:
                                            clneu2 += 1
                                        
#                                    print(u1,v1,u2,v2)
                                    cl[u1,v1,u2,v2] = (clneu + clneu2)/float(G.number_of_nodes()) + cl0
       
                    G.remove_edge(u1,v1)
                    deg[u1] -= 1
                    deg[v1] -= 1

    # find maximum
    m = max(cl, key = cl.get)
    u1 = m[0]
    v1 = m[1]
    u2 = m[2]
    v2 = m[3]

    G.add_edges_from([(u1,v1), (u2,v2)])
    clmax = nx.average_clustering(G)
    G.remove_edges_from([(u1,v1), (u2,v2)])    

    return u1,v1,u2,v2, clmax
   
    
def RandomEdges(G): 
    possEdges = list(nx.non_edges(G))
    Edges = random.sample(possEdges,2)
    
    G.add_edges_from(Edges)
    clmax = nx.average_clustering(G) 
    G.remove_edges_from(Edges)
    U1, V1 = Edges[0]
    U2, V2 = Edges[1]
    return(U1, V1, U2, V2, clmax)

    
# networks

# sparse networks

def loadWeightedNetwork(graphtype, i):
    if graphtype == 'ER':
        G = nx.read_weighted_edgelist('/networks/ER_4prozWK_' + str(i) +'.weighted.edgelist')
    elif graphtype == 'NWS':
        G = nx.read_weighted_edgelist('/networks/NWS_2Nachb_50prozWK' + str(i) +'.weighted.edgelist')
    elif graphtype == 'R':
        G = nx.read_weighted_edgelist('/networks/regular_4Nachb' + str(i) +'.weighted.edgelist') 
        
    for v in G.nodes():
        G.node[v]['pop'] = 1

    G.remove_edges_from(G.selfloop_edges())

    return G
    
    
# dense networks

def loadDenseNetworks(graphtype, i):
    dense = 0.75
    if graphtype == 'ER':
        G = nx.read_edgelist('/Networks/ER_' + str(dense) + '_'  + str(i) +'_2.edgelist')
    elif graphtype == 'NWS':
        G = nx.read_edgelist('/Networks/NWS_' + str(dense) + '_'  + str(i) +'_2.edgelist')
    elif graphtype == 'R':
        G = nx.read_edgelist('/Networks/R_' + str(dense) + '_'  + str(i) +'_2.edgelist') 

    return G


# landscape based networks

def loadGraph(style, square, habitat, maxdist, streamnetwork='025025050_'):
    maxdist = {'linear_01' : 400, 'random_01' : 900, 'clustr_01' : 650}    
        
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
    

# calculate and save results
for method in ['Optimal', 'Greedy', 'LazyGreedy', 'RandomEdges']:
    for graphtype in ['NWS', 'R', 'ER', 'NWS_dense', 'R_dense', 'ER_dense', 'linear_01', 'random_01', 'clustered_01']: 
        results = []
        
        for i in range(250):
            if graphtype in ['NWS', 'R', 'ER']:
                G = loadWeightedNetwork(graphtype, i)
            if 'dense' in graphtype:
                G = loadDenseNetworks(graphtype,i)
            if graphtype in ['linear_01', 'random_01', 'clustr_01']:
                G = loadGraph(graphtype, int(i/10), i%10, maxdist[graphtype])
    
            if method == 'Optimal':
                u1,v1,u2,v2, clmax = optimal(G)
            elif method == 'Greedy':
                u1,v1,u2,v2, clmax = Greedy(G) 
            elif method == 'LazyGreedy':
                u1,v1,u2,v2, clmax = LazyGreedy(G) 
            elif method == 'RandomEdges':
                u1,v1,u2,v2, clmax = RandomEdges(G)
            
            results.append([u1,v1,u2,v2, clmax])
        
        file = open('Greedy_2edges/' + method + '_' + graphtype + '.csv', 'w')
        file.write('id;u1;v1;u2;v2;clustering\n')
        
        for i in range(len(results)):
            file.write(str(i)+';')
            for text in results[i]: 
                file.write(str(text) + ';')
            file.write('\n')
        file.close()
 