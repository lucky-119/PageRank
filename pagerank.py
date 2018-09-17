import networkx as nx
from igraph import *
import numpy as np
from scipy.sparse import csc_matrix

def pageRank(G, s = .85, maxerr = 0.0001): 
    n = G.shape[0]
    for i in range(n):
        x=0
        for j in range(n):
           if(G[i,j]==1):
               x+=1;
        if(x!=0):    
            for j in range(n):
                G[i,j]/=x
        else:
            for j in range(n):
                G[i,j]=1/n;
    print('\n\nStochastic Matrix: \n',G)
    print('\n')
    M = csc_matrix(G,dtype=np.float)
    rsums = np.array(M.sum(1))[:,0]
    sink = rsums==0
    ro, r = np.zeros(n), np.ones(n)
    k=0
    while np.sum(np.abs(r-ro)) > maxerr:
        ro = r.copy()
        k+=1
        for i in range(0,n):
            Ii = np.array(M[:,i].todense())[:,0]
            Si = sink / float(n)
            Ti = np.ones(n) / float(n)
            r[i] = ro.dot( Ii*s + Si*s + Ti*(1-s) )
        print('Page Rank After Iteration ',k,': ',r/sum(r))
    print('\nTotal Number of Iterations: ',k)
    return r/sum(r)

g=nx.DiGraph()
g.add_nodes_from([1,2,3,4,5,6])
g.add_edges_from([(1,2),(1,3),(3,1),(3,2),(3,5),(4,5),(4,6),(5,4),(5,6),(6,4)])
G=nx.to_numpy_matrix(g)
print('Adjancency Matrix: \n',G)
print ('\nFinal Page Rank: ',pageRank(G,0.9))

