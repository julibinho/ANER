import sys
import os
import pandas as pd
import numpy as np


# Special genes or traning set
specialGenes1 = ['ccnb1', 'tpx2', 'aurka', 'cdc20', 'ccna2']
specialGenes2 = ['ccnb1', 'ccna2', 'tpx2', 'cdc20', 'aurka', 'bub1b', 'bora', 'ckap2', 'cep55', 'blm']
specialGenes3 = ['bub1b', 'bora', 'ckap2', 'cep55', 'blm']

specialGenes4 = ["mos", "cdc6", "slbp", "prc1", "btg4", "cnot7", "cnot8"]


def printParameters(dataset, outputFile, top, log, trainingSet):
    print ("LOG:: Parameters checked")
    print ("LOG:: dataset = ", dataset)
    print ("LOG:: outputFile =  ", outputFile)
    print ("LOG:: top =  ", top)
    print ("LOG:: trainingSet =  ", trainingSet)


def checkNbPermutation(permutations, defaultValue=100):
    try:
        if not permutations:
            return defaultValue
        topInt = int(permutations)

    except ValueError:
        print("Error: --permutation parameter must be an integer :", permutations)  
        sys.exit(1)
    else:
        return topInt

def checkLog(log):
    try:
        if not log:
            return False
        topInt = int(log)
        if topInt == 0:
            return False
        return True

    except ValueError:
        print("Error: log parameter must be an integer :", log)  
        sys.exit(1)
    

def checkTop(top, topDefault: int):
    try:
        if not top:
            return topDefault
        topInt = int(top)

    except ValueError:
        print("Error: top parameter must be an integer :", top)  
        sys.exit(1)
    else:
        return topInt

def checkDataset(dataset):
    try:
        if not dataset:
            raise Exception("Error:  parameter -d <dataset> containing CSV file is required")
        data = pd.read_csv(dataset, index_col=0, header=0)
    except UnicodeDecodeError:
        print("Error: provide a csv File for parameter -d")  
        sys.exit(1)
    except IOError:
        print("Error: Couldn’t open file: ", dataset)
        sys.exit(1)
    except FileNotFoundError:
        print("Error: File NOT found, provide a CVS file")
        sys.exit(1)
    except Exception as e:
        print (e)
        sys.exit(1)
    else:
        return data

def checkOutputFile(output_filename):
    try:
        if not output_filename:
            raise Exception("Error:  parameter -o <outputFile> containing output file name is required")
        
        directory = os.path.dirname(output_filename)
        print (directory)
        if directory and not os.path.isdir(directory):
            raise Exception("Error:  The outputFile directory cannot be created " +  output_filename)

    except Exception as e:
        print (e)
        sys.exit(1)
    else:
        return output_filename
    
def getTrainingSet(trainingSet, GeneIds):
    try:
        if not trainingSet or trainingSet == "":
            raise ValueError()
        if len(trainingSet) == 1:
            ts = int(trainingSet)
            if ts not in [1, 2]:
                raise Exception("Error: when parameter -s or --trainingSet is an integer, it should be 1 or 2")
            if ts == 1:
                #print ("LOG:: Using training set 1")
                return specialGenes1
            else:
                #print ("LOG:: Using training set 2")
                return specialGenes4
        else:
            ts = trainingSet.split(',')
            for i in ts:
                if i not in GeneIds:
                    raise Exception("Error: the gene name ", i, "is not valid")
    except ValueError:
        print("Error: parameter -s or --trainingSet must be an integer (1 or 2) or a list of gene names separated by a comma :", trainingSet)  
        sys.exit(1)
    except Exception as e:
        print (e)
        sys.exit(1)

    return ts


def getTrainingSetFromArgs(args):
    for i in range(len(args)):
        if args[i]=='-s' or args[i]=='--trainingSet':
            if (i+1) < len(args):
                return args[i+1]
    return ""