import json
import numpy as np
import subprocess
from scipy.spatial.distance import *
from scipy.cluster.hierarchy import *
from matplotlib import pyplot as plt
from scipy.cluster import hierarchy
from fastcluster import linkage as linkage_fc
# @Author: Soubès Franck


def optimal_scores(Z, rd, dists):
"""
# code initially adapted @markak on GitHub, algorithm is from
# Ziv Bar-Joseph et al., Bioinformatics 2001
"""
    # Z - linkage matrix from scipy.cluster.hierarchy
    # rd - ClusterNode dictionary from to_tree
    # dists - distance matrix

    n_nodes = Z.shape[0] + 1
    M = {}

    # iterating through the linkage matrix guarantees
    # we never see a node before its children
    for i in xrange(Z.shape[0]):
        # linkage matrix starts at first non-leaf node
        v = n_nodes + i
        # the left and right nodes
        j,k = int(Z[i, 0]), int(Z[i, 1])

        if Z[i, 3] == 2:
            # both j and k are leaves, so there is no ordering to be done
            M[v, j, k] = M[v, k, j] = Z[i, 2]
        elif rd[j].is_leaf():
            # if j is a leaf, we calculate the distances to all
            # subtrees of k
            kwns = [kwn for kwn in M if kwn[0] == k]
            for k,w,n in kwns:
                M[v, j, n] = M[v, n, j] = M[k, w, n] + dists[j,w]
                M[v, j, w] = M[v, w, j] = M[k, w, n] + dists[j,n]
        elif rd[k].is_leaf():
            # symmetrically if k is a leaf
            jums = [jum for jum in M if jum[0] == j]
            for j,u,m in jums:
                M[v, m, k] = M[v, k, m] = M[j, u, m] + dists[k,u]
                M[v, u, k] = M[v, k, u] = M[j, u, m] + dists[k,m]
        else:
            # neither j nor k are leaves, so we consider combinations of subtrees
            LL,LR = rd[j].left.pre_order(), rd[j].right.pre_order()
            RL,RR = rd[k].left.pre_order(), rd[k].right.pre_order()

            for (this_L,that_L),(this_R,that_R) in itertools.product(((LL,LR), (LR,LL)),
                                                                     ((RL,RR), (RR,RL))):
                for u,w in itertools.product(this_L, this_R):
                    m_order = sorted(that_L, key=lambda m: M[j, u, m])
                    n_order = sorted(that_R, key=lambda n: M[k, w, n])
                    C = dists[np.ix_(m_order, n_order)].min()
                    Cmin = 1e10
                    for m,n in itertools.product(m_order, n_order):
                        if M[j, u, m] + M[k, w, n] + C >= Cmin:
                            break
                        C = M[j, u, m] + M[k, w, n] + dists[m,n]
                        if C < Cmin:
                            Cmin = C

                    M[v, u, w] = M[v, w, u] = Cmin

    return M



def optimal_ordering(Z, dists):
"""
# code initially adapted @markak on GitHub, algorithm is from
# Ziv Bar-Joseph et al., Bioinformatics 2001
"""
    # Z - linkage matrix
    # dists - the distance matrix

    # get the tree and a list of handles to its leaves
    tree,rd = hierarchy.to_tree(Z, True)

    # Generate scores
    M = optimal_scores(Z, rd, dists)
    # re-order the tree accordingly
    order_tree(Z, rd, M)

    # new leaf ordering
    row_reorder = tree.pre_order()

    return row_reorder



def DIANA(database, param, dendogram ):
  
    vertices = [database.keys()]
    for parameter in param:
        result=[]
        print(parameter)
        for vertex in vertices:
            clustered = []
            uprot_id = []
            for number in database:
                if number in vertex:
                    uprot_id.append(number)
                    clustered.append([database[number][parameter]])
            if clustered != []: # avoir empty clusters 
                try:
                    Z = linkage_fc( pdist(clustered, 'euclidean'), 'ward') #first argument is a distance matrix from the previous command line with the method for forming clusters
                    #optimal_Z = optimal_ordering(Z, matrix_distance)
                    #cutree = cluster.hierarchy.cut_tree(Z, n_clusters=[5, 10]) 
                    vertex1, vertex2 = cut_tree(Z)
                    result.append(get_backid(vertex1,uprot_id))
                    result.append(get_backid(vertex2,uprot_id))
                except:
                    pass
        vertices = result
    plot(Z,dendogram)
    return vertices


def plot(Z,flag):
    if flag:
        plt.figure(figsize=(35, 17))
        plt.title('Hierarchical Clustering Dendrogram: DIANA')
        plt.xlabel('number of clusters by leaf')
        plt.ylabel('euclidian distance')
        dendrogram(
            Z,
            leaf_rotation=90.,
            leaf_font_size=8.,  
        )
        plt.savefig('dendogram.png',bbox_inches='tight')
        plt.show()

        

def cut_tree(tree):
    rootnode, nodelist = hierarchy.to_tree(tree, rd=True)
    left = rootnode.get_left().pre_order(lambda x: x.id)
    right = rootnode.get_right().pre_order(lambda x: x.id)
    return left, right


def lecture(fichier):
    f=open(fichier, "r")
    database=json.load(f)
    return(database)
    f.close()

def get_backid(number,uprot_id):
    number_to_id=[]
    for i in number:
        number_to_id.append(uprot_id[i])
    return number_to_id

def number_cluster(cluster):
    size=1
    for dim in np.shape(cluster): size *= dim
    return size

def save(cluster,size):
    header = ["There is", size, "cluster"]
    cluster[0] = header
    np.savetxt('cluster.txt', cluster,fmt='%5s',delimiter=',')
    

#########################################
##################Main###################
#########################################

#look at the PCA for the parameters orders.
param =  ['longueur_sequence','molecularweight','Sheet','Helix','Turn','phi','hydrophobicity','cystein','aromaticity'] 
precluster = lecture("Human.json")
cluster = DIANA(precluster, param,True)
number_of_cluster = number_cluster(cluster)
save(cluster,number_of_cluster)
subprocess.call("./purge.sh", shell=True) #just need to add the number for each cluster.!


######################################################################################################################"
