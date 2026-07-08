import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.spatial import distance


## nearest neighbors that the algorithm will consider (0 = will consider all genes)
NN = 100

#------------------------------ Functions---------------------------------------#
def closest_node(node, nodes, specialGenesIndex):
    dist = 0
    nodesL = nodes.copy()
    #nodesL = np.delete(nodes, specialGenesIndex, axis=0)
    while (dist == 0):
        mat = distance.cdist([node], nodesL)
        dist = mat.min() 
        #print ("dist min ", dist)
        closest_index = mat.argmin()
        if closest_index in specialGenesIndex:
            dist = 0
            nodesL[closest_index][0] = 1000
            nodesL[closest_index][1] = 1000

    return closest_index, dist

def getClosestPCA(N, specialGenesIndex, X_input):
    # Apply PCA
    pca = PCA(n_components=2)
    pca.fit(X_input)
    X = pca.transform(X_input)

    #print ("N = ", N)

    matRank = []
    indexRemove = specialGenesIndex.copy()
    for j in range(0, N):
        #print ("------------------------")
        for i in specialGenesIndex:
            x_specialGene = (X[i,0], X[i,1]); #print ("sg", i); print (indexRemove)
            closest_index, dist = closest_node(x_specialGene, X, indexRemove)
            matRank.append([closest_index, dist] )
            indexRemove.append(closest_index)
    
    matRankArray = np.array(matRank, dtype='float')
    matRankArraySorted= matRankArray[matRankArray[:, 1].argsort()]
    return matRankArraySorted



def saveTop(sorted_genes, top, geneIndex, output_filename="", sep=","):
    info = "GeneId,GeneName,smallest distance\n"
    for geneInd, rank in enumerate(sorted_genes[:top], start=1):
        info += str(int(rank[0])) + sep + geneIndex[int(rank[0])] + sep +  str(rank[1]) + '\n'
    if output_filename != "":
        with open(output_filename, "w") as f:
            f.write(info)
    else:
        return info






