'''
spectral_cluster(graph,atrb):

Clusters nodes using a directed spectral method found in Gleich(2006)
https://www.cs.purdue.edu/homes/dgleich/publications/Gleich%202005%20-%20hierarchical%20directed%20spectral.pdf

Input:
    graph(networkx DiGraph)
    
    atrb(string): Attribute to cluster on
    
    min_size(int): Smallest number of nodes allowed in a cluster
    
Output:
    cluster(list of lists): Nodes of each cluster in lists
    
    graph(networkx DiGraph): graph with a cluster label under 'Cluster' added
    
Example:
    (cluster,graph) = spectral_cluster(graph,'ShipTypePleasureCraft')
'''
import numpy as np
import networkx as nx
import scipy.sparse.linalg.eigen.arpack as arpack

 
def spectral_cluster(graph,atrb,min_size):

    #initialize lists
    cluster = [graph.nodes()] #initial cluster
    clstr_temp = [] #stores new partitions
    prts = 1 #number of partitions
    end = 0 #number of partitions not cut

    for i in range(10): #maximum 10 iterations
        
        for x in cluster:
            
            #create temporary graph for calculating eigenvalues
            temp_graph = graph.subgraph(x)

            #Find laplacian and eigenvalues for new node list
            try:
                lap = nx.directed_laplacian_matrix(temp_graph,weight=atrb,walk_type='pagerank')
                (vals,vecs) = np.linalg.eig(lap)

            #If error (generally occurs with small unconnected graph) then make no cut
            except(arpack.ArpackNoConvergence,ValueError):
                print 'No eigenvalue convergence for cluster of size %d, did not paritition' % len(x)
                end+=1
                clstr_temp.append(temp_graph.nodes())
                
            else:
                #get eigenvector associated with largest non-trivial eigenvalue
                vec = [l[0] for l in np.real(vecs[:, np.argsort(np.real(vals))[1]]).tolist()] 
            
                #sort eigenvector and get sorted indices
                P = np.sort(vec)
                order = np.argsort(vec)

                #find where the sorted eigenvector crosses 0. if the cut will parition less than min_size nodes then make no cut
                cut = [a+1 for a in range(len(P)-1) if P[a]<0 and P[a+1]>0 ]
                if len(cut) > 0:
                    cut = cut[0] #debugging, sometimes returned in list format
            
                #don't make cut if partitioning less than min_size nodes
                if cut < min_size or len(x)-cut < min_size:
                    cut = 0;
                    end+=1
                    
                #order nodelist using sorted eigenvector indices, split at cut value
                nodes = np.array(temp_graph.nodes())[order]
                if cut > 0:
                    clstr_temp.append(list(nodes[0:cut]))
                clstr_temp.append(list(nodes[cut:len(nodes)]))
        
        #add new clusters to the final list, reset temp
        cluster = clstr_temp
        clstr_temp = []
        
        #if no cuts were made in all clusters then break from loop,
        #else update the number of clusters, reset counter for clusters with no cuts made
        if end == prts:
            break
        else:
            end = 0
            prts = len(cluster)

    #assign node attribute to signify cluster
    label = 1
    if len(cluster) <= 2:
        denom = 4
    else:
        denom = (len(cluster)+len(cluster)%2)
    clstr_label = {}
    for x in cluster:
        for y in x:
            clstr_label[y] = label*np.pi/denom #pi is used to separate around a unit circle (draw function uses np.sin(clstr_label))
            
        label +=2
        
    nx.set_node_attributes(graph,'Cluster',clstr_label)
    
    print '%d Clusters found' % len(cluster)

    return (cluster,graph) 