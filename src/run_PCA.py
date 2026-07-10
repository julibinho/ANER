from collections import defaultdict
from optparse import OptionParser
from sklearn.preprocessing import StandardScaler
import random
import sys
import pandas as pd
import numpy as np
import aner.pca_ranking as pca_ranking
import aner.parameters as parameters

#----------------------------- Args and constants ----------------------
# Define the number of permutations
## nearest neighbors that the algorithm will consider (0 = will consider all genes)
topDefault = 20
logDefault = 0

# Initialize a dictionary to store gene ranks across permutations
gene_ranks = defaultdict(list)

def ts_callback(option, opt, value, parser):
    setattr(parser.values, option.dest, value.split(','))


#===================================================================================
#						Main
#===================================================================================
def main():
    usage = "python run_PCA.py -d <DATASET> -o <OUTPUTFILE> -t <TOP> - s <TRAININGSET> \n"
    parser = OptionParser(usage)
    parser.add_option("-d", "--dataset", dest="dataset",  help="dataset in CSV format")
    parser.add_option("-o", "--outputFile", dest="outputFile",  help="output file name")
    parser.add_option("-t", "--top", dest="top",  help="number of top predictions (default t=20)")
    parser.add_option("-s", "--trainingSet", dest="trainingSet",  help="1 to use ccnb1,tpx2,aurka,cdc20,ccna2; 2 to use mos,cdc6,slbp,prc1,btg4,cnot7,cnot8 or type your list of genes separated by commas and without spaces")
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
    # convert gene index to list
    geneIndex = data.index.to_list()

    trainingSet = parameters.getTrainingSetFromArgs(sys.argv); print (trainingSet)
    trainingSet = parameters.getTrainingSet(trainingSet, geneIndex)
    log = parameters.checkLog(log)


	#Checking parameter format
    top = parameters.checkTop(top, topDefault)
    outputFile = parameters.checkOutputFile(outputFile)

    if log == 1 : parameters.printParameters(dataset, outputFile, top, log, trainingSet)
    
    
    if log == 1 : print ("\nLOG:: Data loaded")
     
    data_aux = pd.read_csv(dataset)
    arr = data

    #specialGenes = specialGenes3

    # Getting the index and value of the special genes
    specialGenesIndex = [geneIndex.index(gene) for gene in trainingSet]
    specialGenesValue = [arr.iloc[idx] for idx in specialGenesIndex]

    #PCA 
    #Data scalling
    if log == 1 :print ("\nLOG:: Data scalling")

    scaler = StandardScaler()
    X_normalized = scaler.fit_transform(arr)
    

    if log == "1":  print ("LOG:: Getting closest neighbor of trainingSet using  PCA coordinates")
    neighborPCA = pca_ranking.getClosestPCA(top, specialGenesIndex, X_normalized)
    pca_ranking.saveTop(neighborPCA, top, geneIndex, outputFile)

    if log == 1 : print ("LOG:: done")
    print ("***********************")
#===================================================================================
if __name__ == "__main__":
	main()
