import AnerLib
import parameters
import networkx as nx
from networkx.algorithms.approximation import clique
from scipy.spatial import distance
from scipy.spatial.distance import pdist, squareform
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import sys
from tabulate import tabulate
from optparse import OptionParser
import time



## minimum number of neighbors for the algorithm to consider that gene as "relevant"
MinNB=1
## nearest neighbors that the algorithm will consider (0 = will consider all genes)
NN = 100
## standard deviation, higher the value more relaxed are the clean graph distance procedure, negative values will restric the search
std_dev = 0.5
## number of genes that will be plotted taking into consideration the results in the table
num_genes_table = 30
## number of genes that will be plotted taking into consideration the pareto frontier
num_genes_pareto = 20



usage = "python run_Aner.py -d <DATASET> -o <OUTPUTFILE> -t <TOP> - s <TRAININGSET> \n"
parser = OptionParser(usage)
parser.add_option("-d", "--dataset", dest="dataset",  help="dataset in CSV format")
parser.add_option("-o", "--outputFile", dest="outputFile",  help="output file name")
parser.add_option("-t", "--top", dest="top",  help="number of top predictions (default t=20)")
parser.add_option("-s", "--trainingSet", dest="trainingSet",  help="set 1 to use ccnb1,tpx2,aurka,cdc20,ccna2; set 2 to use mos,cdc6,slbp,prc1,btg4,cnot7,cnot8 or type your list of genes separated by commas and without spaces")
parser.add_option("-l", "--log", dest="log",  help="set to 1 to print log messages; set to 0 otherwise")

if len(sys.argv) < 2:
  parser.print_help()
  sys.exit(1)
	
(options, args) = parser.parse_args()

dataset = options.dataset
outputFile = options.outputFile
top = options.top
log = options.log



data = parameters.checkDataset(dataset)
geneIndex = data.index.to_list()

trainingSet = parameters.getTrainingSetFromArgs(sys.argv); print (trainingSet)
trainingSet = parameters.getTrainingSet(trainingSet, geneIndex)

print ("LOG:: Using training set: ", trainingSet) 
    
#Checking parameter format
top = parameters.checkTop(top, num_genes_table)
outputFile = parameters.checkOutputFile(outputFile)
log = parameters.checkLog(log)



if log: parameters.printParameters(dataset, outputFile, top, log, trainingSet)


# Loading data
if log: print ("\nLOG:: Data loaded")
data_aux = pd.read_csv(dataset)
arr = data 


#saving complete data
df = pd.DataFrame(data=arr, index=data.index, columns=data.columns)
df_aux = pd.DataFrame(data_aux)

# convert gene index to list
geneIndex = data.index.to_list()

# Special genes that are relevant
specialGenes = trainingSet

specialGenesIndex = [geneIndex.index(gene) for gene in specialGenes]
specialGenesValue = [arr.iloc[idx] for idx in specialGenesIndex]

# Generate geneIndexValue directly (0, 1, 2, ..., len(geneIndex)-1)
geneIndexValue = list(range(len(geneIndex)))

#geneNames = [ df_aux.loc[int(geneid), 'geneid'] for geneid in geneIndex]

#print (geneNames)

#sys.exit()
# Normalize the data (if not already normalized)
scaler = StandardScaler()
X_normalized = scaler.fit_transform(arr)

X_normalized = arr

if log: print ("\nLOG:: Data normalized")

if NN == 0:
  NN = len(X_normalized)
  print(NN)

nbrsTestIdx = np.zeros((len(X_normalized), NN-1), dtype=int)
nbrsTestDist = np.zeros((len(X_normalized), NN))

# if you want to test many distances and check the distance matrix, just use a list of distances, but the final result will only show the last distance used in the list of distances!!! Work in progress to store all possible distances results!!!
#distances = ['correlation', 'canberra', 'braycurtis', 'chebyshev', 'euclidean']
distance = 'euclidean'

G, adj = AnerLib.processGraph(X_normalized, specialGenesIndex, geneIndex, distance, log)



start_time = time.time()
rankMat, rankMat2 = AnerLib.rankingNodes(G, specialGenesIndex) # The idea is to rank genes and give the top genes
if log: print(f"\nLOG:: Ranking %.2f seconds ---" % (time.time() - start_time))


rankMat.sort(key=lambda k: (k[1], -k[4], k[5], k[2], k[3]), reverse=True) ### k[1] - k[5] descending order, k[6] & k[7] ascending order


infoRank = AnerLib.printRank(rankMat, geneIndex,  top, "\t")
AnerLib.saveOutput(infoRank, outputFile + ".topRanked.txt")


pareto_final = AnerLib.getParetoSet(rankMat2)
lista_completa_genes, lista_completa = AnerLib.multiObjRanking(pareto_final, rankMat)

topRankedPareto = AnerLib.printFrontPareto(lista_completa, geneIndex, specialGenesIndex, top)

AnerLib.saveOutput(topRankedPareto, outputFile + ".pareto.txt")

print ("\nOutputs saved:"); 
print ("Ranking = ", outputFile + ".topRanked.txt"); 
print ("Pareto fontiers =",  outputFile + ".pareto.txt" )

#print (topRankedPareto)
#print(tabulate(rankMat[0:20], headers=['geneid', '#TotalNeighbors', '#Cliques', '#MaximalClique', 'AvgDistance', 'Specific Genes']))


