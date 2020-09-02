import psycopg2
import networkx as nx
import pandas as pd
import os


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
            
    # potential edges
    
    cursor = conn.cursor()
    cursor.execute('SELECT start, aim, distance FROM dis_pts_2500_10x10_' + style + '.dist_pts_2500_' + str(square) + '_start_' + str(habitat))
    
    potedges = cursor.fetchall()
    potedges = [list(x) for x in potedges]
    
    return G, potedges
    

def maxClustering(G, potedges):        
    deg = nx.degree(G)
    cl = {}
    cl0 = nx.average_clustering(G) # clustering of original graph
    clmax = cl0
    n = float(G.number_of_nodes())
    
    for i in range(len(potedges)):
        u = potedges[i][0]
        v = potedges[i][1]
        if (u,v) not in G.edges() and (v,u) not in G.edges():
   
                clneu = 0
                cn = list(nx.common_neighbors(G,u,v))
                k = len(cn)
                if k != 0:
                    # add difference of clustering for neighbours
                    for w in cn:
                        dw = deg[w]
                        clneu += 2.0/(dw * (dw-1))
                        
                    # add difference of clustering for u,v
                    if deg[u] > 1:
                        clneu += (2.0*k*(deg[u]-1) - 4.0*nx.triangles(G,u))/(deg[u]*(deg[u]**2 - 1))
                    else: 
                        clneu += 1
                    if deg[v] > 1:
                        clneu += (2.0*k*(deg[v]-1) - 4.0*nx.triangles(G,v))/(deg[v]*(deg[v]**2 - 1))
                    else:
                        clneu += 1
                        
                cl[u,v] = clneu/ n
  
    # find biggest clustering 
    CL = pd.DataFrame(pd.Series(cl))
    maxcl = CL[CL[0] == CL.max()[0]]
               
    U = maxcl.index[0][0]
    V = maxcl.index[0][1]
    clmax = maxcl[0][U][V]+cl0
    return U,V,clmax
    

maxdist = {'linear_02' : 250, 'random_02' : 600, 'clustr_02' : 375, 'linear_01' : 400, 'random_01' : 900, 'clustr_01' : 650}    

# Greedy 
for nEdges in [5,10,15,20,25,30]:
    for style in ['random_01', 'linear_01', 'clustr_01']:
        for square in range(1,26):
            for habitat in range(10):                                                            
                G, potedges = loadGraph(style, square, habitat, maxdist[style])
                potedges = [[potedges[i][0], potedges[i][1]] for i in range(len(potedges))]
              
                torem = []
                for [u,v] in potedges: 
                    if u > v: 
                        torem.append([u,v])
                for e in torem: 
                    potedges.remove(e)
                    
                results = []
                for i in range(nEdges): 
                    u,v,clmax = maxClustering(G, potedges)  
                    G.add_edge(u,v) 
                    potedges.remove([u,v])
                    results.append([u,v])     
    
                folder = '/Greedy/' + style + '/' + str(nEdges) + 'edges/'
                file = open(folder + 'edges_' + str(square) + '_' + str(habitat) + '.csv', 'w')
                file.write('start;aim\n')
                for e in results: 
                    e = sorted(e)
                    file.write(str(e[0]) + ';' + str(e[1]) + '\n')
            
                file.close()
    

# Lazy Greedy    
for nEdges in [2,5,10,15,20,25,30]:
    for style in ['random_01', 'linear_01', 'clustr_01']:
        for square in range(1,26):
            print(style, 'square', square, nEdges)
            for habitat in range(10):
                
                G, potedges = loadGraph(style, square, habitat, maxdist[style])
                potedges = [[potedges[i][0], potedges[i][1]] for i in range(len(potedges))]
              
                torem = []
                for [u,v] in potedges: 
                    if u > v: 
                        torem.append([u,v])
                for e in torem: 
                    potedges.remove(e)
        
                
                deg = nx.degree(G)
                cl = {}
                CL = pd.DataFrame(index = ['U', 'V', 'Clustering'])
                cl0 = nx.average_clustering(G)
                n = float(G.number_of_nodes())
                count = 0
                
                for i in range(len(potedges)):
                    u = potedges[i][0]
                    v = potedges[i][1]

                    if (u,v) not in G.edges() and (v,u) not in G.edges():             
                        clneu = 0 
                        cn = list(nx.common_neighbors(G,u,v))
                        k = len(cn) 
                        
                        if k != 0: 
                            # add difference of clustering for neighbours
                            for w in cn:
                                dw = deg[w]
                                clneu += 2.0/(dw * (dw-1))
                                    
                            # add difference of clustering for u,v
                            if deg[u] > 1:
                                clneu += (2.0*k*(deg[u]-1) - 4.0*nx.triangles(G,u))/(deg[u]*(deg[u]**2 - 1))
                            else: 
                                clneu += 1
                            if deg[v] > 1:
                                clneu += (2.0*k*(deg[v]-1) - 4.0*nx.triangles(G,v))/(deg[v]*(deg[v]**2 - 1))
                            else:
                                clneu += 1
                                
                        cl[count] = [u,v, clneu/n + cl0] 
                        count += 1
                                
                CL = pd.DataFrame(cl) 
                CL.index = ['U', 'V', 'Clustering']
                CL = CL.T
                CL.sort_values('Clustering', ascending = False, inplace = True)
                CL.index = range(len(CL))
                
                newEdges = [[CL['U'][i], CL['V'][i]] for i in range(nEdges)]
                

                folder = '/LazyGreedy/' + style + '/' + str(nEdges) + 'edges/'
  
                file = open(folder + style +'_' + str(square) + '_' + str(habitat) + '.csv', 'w')
                
                file.write('start;aim\n')
                
                for i in range(nEdges): 
                    file.write(str(newEdges[i][0]) + ';' + str(newEdges[i][1]) + '\n')
                file.close()
    
# random edges

for nEdges in [5,10,15,20,25,30]:
    for style in ['random_01','linear_01', 'clustr_01']:
        for square in range(1,26):
            for habitat in range(10):
                G, potedges = loadGraph(style, square, habitat, maxdist[style])
                potedges = [[potedges[i][0], potedges[i][1]] for i in range(len(potedges))]
                torem = []
                for e in potedges: 
                    if e in G.edges(): 
                        torem.append(e)
                for e in torem: 
                    potedges.remove(e)
                Edges = random.sample(potedges,nEdges)
                
                
                folder = '/RandomEdges/' + style + '/' + str(nEdges) + 'edges/'
                file = open(folder + 'edges_' + str(square) + '_' + str(habitat) + '.csv', 'w')
                file.write('start;aim\n')
                for e in Edges: 
                    file.write(str(e[0]) + ';' + str(e[1]) + '\n')
                file.close()

    
        