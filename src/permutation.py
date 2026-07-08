from collections import defaultdict
from sklearn.preprocessing import StandardScaler
import random
import sys
import pandas as pd
import numpy as np
import pca_ranking
import AnerLib
import parameters
from optparse import OptionParser


#----------------------------- Args and constants ----------------------
# Define the number of permutations
num_permutations = 100  # Replace with the desired number of permutations
## nearest neighbors that the algorithm will consider (0 = will consider all genes)
NN = 200

# Initialize a dictionary to store gene ranks across permutations
gene_ranks = defaultdict(list)


#----------------------------- Functions -------------------------------

def select_cats(df1):
    # Initialize a list to store all selected CATS values across subcategories
    all_selected_cats = []
    
    # Iterate over each subcategory and select 50% of its CATS values
    for subcat, cat_list in zip(df1["Subcat"], df1["CATS"]):
        # Calculate 50% of the CATS values for this subcategory
        num_cats_to_select = max(1, len(cat_list) // 2)  # Ensure at least 1 category is selected
        # Randomly select 50% of the CATS values for this subcategory
        selected_cats = random.sample(cat_list, num_cats_to_select)
        # Add the selected CATS values to the global list
        all_selected_cats.extend(selected_cats)
        print(f"  Selected 50% of CATS for {subcat}: {selected_cats}")
    return all_selected_cats


def gene_rank_tool(gene_ranks_local, geneIndex, tool):
    
    for rank, geneInfo in enumerate(tool[:200], start=1):
        #print (int(rank), geneInfo[0], geneInfo[1])
        geneNum = geneInfo[0]#; print ("Gene Num ", geneNum)
        geneName = geneIndex[int(geneNum)]
        gene_ranks_local[geneName].append((perm + 1, rank))  # Store permutation number and rank

    return gene_ranks_local

def permutateRank(df1, arr, specialGenesIndex, geneIndex):

    # Initialize a dictionary to store gene ranks across permutations
    gene_ranks = defaultdict(list)
    for perm in range(num_permutations):
        print(f"Processing Permutation {perm + 1}")
    
        all_selected_cats = select_cats(df1)
        # Convert selected CATS numbers to column names
        cat_cols_to_drop = [f"CAT{cat}" for cat in all_selected_cats if f"CAT{cat}" in arr.columns]
    
        # Create new data frame with exclusions
        df_excluded = arr.drop(columns=cat_cols_to_drop, errors="ignore")

        print(f"  Excluded Columns: {cat_cols_to_drop}")
        X = df_excluded
        # Normalize the data (if not already normalized)
        scaler = StandardScaler()
        X_normalized = scaler.fit_transform(X)
    
        # Getting closest neighbor of specialGenes using  PCA coordinates
        neighborPCA = pca_ranking.getClosestPCA(20, specialGenesIndex, X_normalized)

        gene_ranks = gene_rank_tool(gene_ranks, geneIndex,  neighborPCA)
        
            
    return gene_ranks
    

def sortGenesPermutation(gene_ranks):
    
    # Create a DataFrame to store ranks for each permutation
    rank_data = []

    # Iterate through each gene and its ranks
    for gene, ranks in gene_ranks.items():
        rank_dict = {"Gene": gene}
        for perm, rank in ranks:
            rank_dict[f"Permutation {perm}"] = rank
        rank_data.append(rank_dict)

    # Convert the list of dictionaries to a DataFrame
    gene_ranks_df = pd.DataFrame(rank_data)

    # Fill missing values with a high rank (e.g., 1000) and convert to integers
    gene_ranks_df = gene_ranks_df.fillna(1000)

    # Calculate the average rank for each gene
    gene_ranks_df['AverageRank'] = gene_ranks_df.iloc[:, 1:].mean(axis=1)

    # Sort genes by their average rank
    sorted_genes = gene_ranks_df.sort_values(by='AverageRank')

    return sorted_genes

def saveTop(sorted_genes, top, output_filename):
    top_genes = sorted_genes.head(top)

    # Cap ranks at 20 for rank columns only
    rank_columns = [col for col in top_genes.columns if col.startswith('Permutation') or col == 'AverageRank']
    top_genes[rank_columns] = top_genes[rank_columns].map(lambda x: min(x, 2022))
    top_genes.to_csv(output_filename, index=False)
    
#----------------------------- Main -------------------------------
usage = "python permutation.py -d <DATASET> -o <OUTPUTFILE> -t <TOP> - s <TRAININGSET> \n"
parser = OptionParser(usage)
parser.add_option("-d", "--dataset", dest="dataset",  help="dataset in CSV format")
parser.add_option("-o", "--outputFile", dest="outputFile",  help="output file name")
parser.add_option("-p", "--permutations", dest="p",  help="number of permutations (default p=100)")
parser.add_option("-s", "--trainingSet", dest="trainingSet",  help="set 1 to use ccnb1,tpx2,aurka,cdc20,ccna2; set 2 to use mos,cdc6,slbp,prc1,btg4,cnot7,cnot8 or type your list of genes separated by commas and without spaces")
parser.add_option("-l", "--log", dest="log",  help="set to 1 to print log messages; set to 0 otherwise")

if len(sys.argv) < 2:
  parser.print_help()
  sys.exit(1)

(options, args) = parser.parse_args()

dataset = options.dataset
outputFile = options.outputFile
p = options.p
log = options.log

data = parameters.checkDataset(dataset)
geneIndex = data.index.to_list()

trainingSet = parameters.getTrainingSetFromArgs(sys.argv); print (trainingSet)
trainingSet = parameters.getTrainingSet(trainingSet, geneIndex)

print ("LOG:: Using training set: ", trainingSet) 
    
#Checking parameter format
p = parameters.checkNbPermutation(p, num_permutations)
outputFile = parameters.checkOutputFile(outputFile)
log = parameters.checkLog(log)




if log: parameters.printParameters(dataset, outputFile, p, log, trainingSet)




dataCat = {
    "Subcat": ["Xenopus", "Mouse", "Human", "Cells"],
    "CATS": [
        [7, 10, 12, 13, 25, 34, 38, 39],
        [11, 14, 17, 19, 21, 22, 26, 27, 28, 32, 33, 35],
        [20, 23, 24, 30, 31],
        [1, 2, 3, 4, 5, 6, 8, 9, 16, 36, 37],
    ],
}


#cols = list(pd.read_csv(fileName, nrows=1)) #print(cols)
# Loading data
data = pd.read_csv(dataset, index_col=0, header=0)
# convert gene index to list
geneIndex = data.index.to_list()
data_aux = pd.read_csv(dataset)
arr = data


if log: print ("Data loaded")


# Special genes that are relevant
specialGenes = trainingSet


# Getting the index and value of the special genes
specialGenesIndex = [geneIndex.index(gene) for gene in specialGenes]
specialGenesValue = [arr.iloc[idx] for idx in specialGenesIndex]



#Creating a data frame for categories
df1 = pd.DataFrame(dataCat)
subcat_averages = {}

#gene_ranks = permutateRank(df1, arr, specialGenesIndex, geneIndex)


gene_ranks_pca = defaultdict(list); gene_ranks_aner =  defaultdict(list);

for perm in range(p):
    print(f"Processing Permutation {perm + 1}")
    
    all_selected_cats = select_cats(df1)
    # Convert selected CATS numbers to column names
    cat_cols_to_drop = [f"CAT{cat}" for cat in all_selected_cats if f"CAT{cat}" in arr.columns]
    
    # Create new data frame with exclusions
    df_excluded = arr.drop(columns=cat_cols_to_drop, errors="ignore")

    print(f"  Excluded Columns: {cat_cols_to_drop}")
    X = df_excluded
    # Normalize the data (if not already normalized)
    scaler = StandardScaler()
    X_normalized = scaler.fit_transform(X)

    #PCA 
    neighborPCA = pca_ranking.getClosestPCA(20, specialGenesIndex, X_normalized) # Getting closest neighbor of specialGenes using  PCA coordinates
    gene_ranks_pca = gene_rank_tool(gene_ranks_pca, geneIndex,  neighborPCA)
    
    # Aner
    G, adj = AnerLib.processGraph(X, specialGenesIndex, geneIndex, 'euclidean', False)
    rankMat, rankMat2 = AnerLib.rankingNodes(G, specialGenesIndex) # The idea is to rank genes and give the top genes
    rankMat.sort(key=lambda k: (k[1], -k[4], k[5], k[2], k[3]), reverse=True) ### k[1] - k[5] descending order, k[6] & k[7] ascending order


    #print (rankMat)
    gene_ranks_aner = gene_rank_tool(gene_ranks_aner, geneIndex,  rankMat)
    

#print (gene_ranks_pca)
print ("***********************")
#print (gene_ranks_aner)



sorted_genes_pca = sortGenesPermutation(gene_ranks_pca)
saveTop(sorted_genes_pca, 20, f"{outputFile}_pca_{p}_permutations.csv")

sorted_genes_aner = sortGenesPermutation(gene_ranks_aner)
saveTop(sorted_genes_aner, 20, f"{outputFile}_aner_{p}_permutations.csv")


print ("done")
