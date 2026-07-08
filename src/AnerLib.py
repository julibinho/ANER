
import networkx as nx
from networkx.algorithms.approximation import clique
from paretoset import paretoset, paretorank
from scipy.spatial import distance
from scipy.spatial.distance import pdist, squareform
import numpy as np

import time

## minimum number of neighbors for the algorithm to consider that gene as "relevant"
MinNB=1
## nearest neighbors that the algorithm will consider (0 = will consider all genes)
NN = 100
## standard deviation, higher the value more relaxed are the clean graph distance procedure, negative values will restric the search
std_dev = 0.5
## number of genes that will be plotted taking into consideration the results in the table
num_genes_table = 20
## number of genes that will be plotted taking into consideration the pareto frontier
num_genes_pareto = 20

def createGraph(indicesKNN, distanceKNN, specialGenesIndex, geneIndex):
  """
  create a graph network based on KNN distances
  input1 indicesKNN: matrix with indices of k closest neighbors
  input2 distanceKNN: matrix with distances of k closest neighbors
  input3 specialGenesIndex: identifier of special genes (red)
  input4 geneIndex: identifier of all genes (red)
  output1 G: graph networkx
  output2 adj: list of blue nodes
  """
  adj = []
  G = nx.Graph()


  for i in specialGenesIndex:
    G.add_node(str(i), label=str(geneIndex[i]), color="red")
    lKNN = indicesKNN[i, :] #taking all neigbors of special gene i
    dKNN = distanceKNN[i, :] #taking all distances  --> previously taking indicesKNN
    for j in range(1,len(lKNN)): #ignore first position it's itself
      if lKNN[j] not in G:
        if lKNN[j] not in specialGenesIndex: #and (dKNN[j] < 0.1):
          G.add_node(str(lKNN[j]), label=str(geneIndex[lKNN[j]]), color="blue")
          adj.append(lKNN[j])
        else:
          G.add_node(str(lKNN[j]), label=str(geneIndex[lKNN[j]]), color="red")
      if(dKNN[j] != 0):
        G.add_edge(str(i), str(lKNN[j]), weight=float(dKNN[j]))

  return G, adj

def cleanGraph(G, nbNeighbors, specialGenesIndex):
  """
  Remove nodes with less than nbNeighbors
  input1 G: graph networkx
  input2 nbNeighbors: number of minimum neighbors
  input3 specialGenesIndex: identifier of special genes (red)
  output1 G: graph networkx (updated)
  """
  toBeRemoved = []
  for node in G:
    if(int(node) not in specialGenesIndex):
      if len(list(G.neighbors(node))) < nbNeighbors:
        if G.neighbors(node) not in specialGenesIndex:
          toBeRemoved.append(node)
  G.remove_nodes_from(toBeRemoved)
  return G


def cleanGraphDistance(G, distancesKNN, std_dev):
  """
  Remove nodes with distance less than threshold, in this case we are using mean of distances of that node + standard deviation of those distances.
  input1 G: graph networkx
  input2 distancesKNN: matrix with distances of k closest neighbors
  input3 std_dev: standard deviation used for the threshold
  output1 G: graph networkx (updated)
  """
  to_remove = []
  for a,b,attrs in G.edges(data=True):
    threshold = np.mean(distancesKNN[int(a)]) + (std_dev) * np.std(distancesKNN[int(a)])
    #print(threshold, attrs["weight"])
    if attrs["weight"] >= threshold and a != b:
      to_remove.append((a,b))
  G.remove_edges_from(to_remove)
  return G

def getBlueNodes(G, specialGenesIndex):
  """
  get all the blue nodes (not special genes) from the graph, used to build the adjacency matrix
  input1 G: graph networkx
  input2 specialGenesIndex: identifier of special genes (red)
  output: genes that are connected to the special genes
  """
  rnodes = []
  for node in G.nodes:
    if int(node) not in specialGenesIndex:
      rnodes.append(int(node))
  return rnodes

def countNeighbors(G, node, specialNg):
  """
  count Neighbors of a given node
  input1 G: graph networkx
  input2 node: a node of G
  input3 specialNg: set of special nodes
  output1 numNg: number of neighbors under conditions
  output2 totalNg: total number of neighbors (without restrictions)
  """
  neighbors = list(G.neighbors(node))

  count = 0
  for i in neighbors:
    if int(i) in specialNg:
      count = count + 1
  return count, len(neighbors)

def listNeighbors(G, node, specialNg):
  neighbors = list(G.neighbors(node))

  list_special = []
  for i in neighbors:
    if int(i) in specialNg:
      list_special.append(i)
  return list_special

def getEdgesData(edgesData, node, specialNg):
  """
  Compute the mean distance of a node
  input1 edgesData: a list of tuple (node1, node2, {'weight': edge cost})
  input2 node: a node of Graph
  input3 specialNg: set of special nodes
  output1 mean distance: the average distance of node and all its neighbors
  output2 mean distance special nodes: the average distance of a node and its special neighbors (TODO)
  """
  dist = 0; count = 0; distSpecialNg = 0; countSp = 0
  for ed in edgesData:
    n1, n2, dicWeigth = ed
    if node==n1 or node==n2:
      dist = dist + dicWeigth['weight']
      count = count + 1
      if (node==n1 or node==n2) and ((int(n1) in specialNg) or (int(n2) in specialNg)):
        distSpecialNg = distSpecialNg + dicWeigth['weight']
        countSp = countSp + 1

  return (dist/(count)), (distSpecialNg/(countSp))

def getCliquesInfo(cliques, node, specialNg, minSizeClique=0):
  """
  Get clique properties of a given node
  input1 cliques: list of all cliques in a graph
  input2 node: a node of G
  input3 specialNg: set of special nodes
  output1 nbCliques: number of cliques contaning the node
  output2 sizeMaxCliques: size of largest clique contaning the node
  output3 specialCliques: number of special nodes in the clique
  """
  nbCliques = 0; largestClique = 0
  for clique in cliques:
    if len(clique) >= minSizeClique:
      if node in clique:
        nbCliques = nbCliques + 1
        if len(clique) > largestClique:
          largestClique = len(clique)
  return nbCliques, largestClique


def rankingNodes(G, specialGenesIndex):
  """
  rank all nodes in a graph
  input1 G: graph networkx
  input2 specialNodes: set of special nodes
  output1 rankMat: numpy array containind nodes and its properties
  nbSpecialNg --> number of neighbors under conditions
  nbTotalNg --> total number of neighbors (without restrictions)
  nbCliques --> number of cliques contaning the node
  largestClique --> size of largest clique contaning the node
  specialCliques --> number of special nodes in the clique
  meanDistance --> the average distance of node and all its neighbors
  dictionary_index[str(node)]
  """
  cliques = nx.enumerate_all_cliques(G)
  lCliques = list(cliques)

  edgesData = G.edges.data()
  rankMat = []
  rankMat2 = []
  count = 0
  for node in G:
    if int(node) not in specialGenesIndex:
      nbSpecialNg, nbTotalNg = countNeighbors(G, node, specialGenesIndex)
      #print (node, nbSpecialNg, nbTotalNg); sys.exit()
      listSpecialNg = listNeighbors(G, node, specialGenesIndex)
      nbCliques, largestClique =  getCliquesInfo(lCliques, node, specialGenesIndex)
      meanDistance, meanDistanceSp = getEdgesData(edgesData, node, specialGenesIndex)
      mat = []
      mat.append(node); mat.append(nbTotalNg); mat.append(nbCliques); mat.append(largestClique); mat.append(meanDistanceSp); mat.append(listSpecialNg)
      rankMat.append(mat)
      mat2 = []
      mat2.append(nbTotalNg); mat2.append(nbCliques); mat2.append(largestClique); mat2.append(meanDistanceSp)
      rankMat2.append(mat2)
      count +=1

  return rankMat, rankMat2



def argsmallest_n(a, n):
  n = n-1
  ret = np.argpartition(a, n)[:n]
  b = np.take(a, ret)
  return np.take(ret, np.argsort(b))


def getCloseNgs(X_normalized, distanceType='euclidean', log=False):
    start_time = time.time()
    
    nbrsTestIdx = np.zeros((len(X_normalized), NN-1), dtype=int)
    nbrsTestDist = np.zeros((len(X_normalized), NN))

    distances = pdist(X_normalized, metric=distanceType)
    dist_matrix = squareform(distances)
    D = dist_matrix

    for i in range(0, len(X_normalized)):
        #idx = np.argsort(D[i])[:NN]
        idx = argsmallest_n(D[i], NN)

        nbrsTestIdx[i] = idx[:(NN)]
        for j in range(NN-1):
            nbrsTestDist[i][j] = D[i][idx[j]]

    if log and (i+1)%1000 == 0 : 
        print(f"LOG:: Distances sorted for {i+1} genes out of {len(X_normalized)}")
        print(f"LOG:: Distances calculate in %.2f seconds ---" % (time.time() - start_time))

    #   distanceKNN = array to store the distance calculated above
    distancesKNN = nbrsTestDist
    # indicesKNN = array to store the indexes of the closest NN neighbours
    indicesKNN = nbrsTestIdx

    return distancesKNN, indicesKNN


def processGraphAner(arr, specialGenesIndex, geneIndex, distance, log=False):
    return processGraph(arr, specialGenesIndex, geneIndex, distance, log=False)

def processGraph(X_normalized, specialGenesIndex, geneIndex, distance, log=False):
    distancesKNN, indicesKNN = getCloseNgs(X_normalized, distance, log)


    start_time = time.time()
    # CreateGraph = function to create the graph with the distances previously calculated, will take into consideration the special genes (red) and other ones (blue)
    G, adj = createGraph(indicesKNN, distancesKNN, specialGenesIndex, geneIndex)
    if log : print(f"\nLOG::Creating Graph %.2f seconds ---" % (time.time() - start_time))

    if log : print("LOG:: Number of nodes = ", G.size())

    start_time = time.time()
    # cleanGraphDistance = function to clear the graph with the distance threshold and the standard deviation, just change the standard deviation if you want to relax the distance considered
    G = cleanGraphDistance(G, distancesKNN, std_dev)
    if log : 
        print(f"LOG:: Cleaning per distance e %.2f seconds ---" % (time.time() - start_time))
        print("LOG:: Number of nodes = ", G.size())

    start_time = time.time()
    # cleanGraph = function to clear the graph with the minimum neighbors threshold, can increase the parameter NB to restrict the number of miminum neighors to be considered (this will apply for the blue genes)
    G = cleanGraph(G, MinNB, specialGenesIndex)
    if log : 
        print(f"LOG:: Cleaning per neigbords %.2f seconds ---" % (time.time() - start_time))
        print("LOG:: Number of nodes = ", G.size())

    start_time = time.time()
    # getBlueNodes = fucntion to get the blue nodes
    adj = getBlueNodes(G, specialGenesIndex)
    if log : print(f"LOG:: Computing adj matrix %.2f seconds ---" % (time.time() - start_time))

    return G, adj


def printRank(matrank, geneIndex, top, sep="\t"):
    
    info = "GeneId"+ sep+ "GeneName" +sep + "neighbors" + sep + "cliques" + sep + "largest clique" + sep + "mean distance to training set\n"
    print (len(matrank))
    
    for i in range(0, min(top, len(matrank))):
        info += matrank[i][0] + sep + geneIndex[int(matrank[i][0])] + sep + str(matrank[i][1]) + sep +str(matrank[i][2]) + sep + str(matrank[i][3]) + sep + str(matrank[i][4]) + "\n"
    return info


def getParetoSet(rankMat2):
    cost = paretorank(rankMat2, sense=["max", "max", "max", "min"])

    pareto_set = []; pareto_final = []
    for j in range(1, max(cost)):
        pareto_set = []
        for i in range(0, len(cost)):
            if cost[i] == j:
                pareto_set.append(rankMat2[i])
        pareto_final.append(pareto_set)
    return pareto_final


def printFrontPareto(lista_completa, geneIndex, trainIndex, top=12):
    ### final list of the genes separated by frontier
    info = ""
    for i in range(0, min(top, len(lista_completa))):
        info += "candidates of frontier " + str(i) + " ->" 
        for j in lista_completa[i]:
            info += " " + str( geneIndex[int(j)] )
        info += "\n"
    return info

def multiObjRanking(pareto_final, rankMat):
    lista_genes = []; lista_completa = []; lista_completa_genes = []; complete_list_pareto = []
    best_result = [rankMat[i][0] for i in range(len(rankMat))]
    
    for k in range(0, len(pareto_final)):
        lista_genes = []; lista_genes_id = []
        for j in range(0, len(pareto_final[k])):
            for i in range(0, len(rankMat)):
                if pareto_final[k][j][0] in rankMat[i] and pareto_final[k][j][1] in rankMat[i] and pareto_final[k][j][2] in rankMat[i] and pareto_final[k][j][3] in rankMat[i]:
                    lista_genes.append(rankMat[i][0])
                    lista_genes_id.append(best_result[i])
                    complete_list_pareto.append(best_result[i])
        lista_completa.append(lista_genes)
        lista_completa_genes.append(lista_genes_id)
    return lista_completa_genes, lista_completa

def saveOutput(topranked, outputFile):
    with open(outputFile, "w") as f:
        f.write(topranked)