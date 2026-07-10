# ANeR

**Attribute Network-based Ranking (ANeR) identifies conserved translational fingerprints of M-phase.**

Attribute Network-based Ranking (ANeR) model that integrates
diverse large-scale datasets to extract shared properties. This allowed identification of new
proteins whose translational regulation mirrors that of a highly curated training set of proteins
with demonstrated translation activation during cell division in both mitotic and meiotic cells,
across multiple organisms, using low-throughput approaches.

**REFERENCE**
Sarah B. Eivers, Axelle Baumgartner, Marcio Dorn, Catherine Jessus, Juliana Silva Bernardes, and Enrico Maria Daldello. Attribute Network-based Ranking (ANeR) identifies conserved translationalfingerprints of M-phase. To be submitted.

**CONTACT**
  E-mail: 
  enrico.daldello@sorbonne-universite.fr,
  juliana.silva_bernardes@sorbonne-universite.fr 


## Requirements 

  * We strongly recommend [anaconda](https://docs.anaconda.com/anaconda/install/) environment.
  * See all requirements at requirements.txt file.

## Using ANeR 
   
#### Step 1: Preparing data
   
  ```
  $ python src/prepareData.py <INPUT_FILE> <OUTPUT_FILE>
  ```
#### Required arguments 
  * INPUT_FILE is a xlsx file containing the datasets
  * OUTPUT_FILE is the path of the ouptut file

  For instance:
  ```
  $python src/prepareData.py examples/Oocyte_datasets.xlsx examples/Oocyte_datasets.csv
  ```

#### Step 2: Running ANeR
   
  ```
  $ python src/run_Aner.py -d <DATASET> -o <OUTPUTFILE> -t <TOP> -s <TRAININGSET>
  ```

#### Arguments

  * -d DATASET the dataset in CSV format (produced in step 1)
  * -o OUTPUTFILE is the path of the ouptut file
  * -t TOP  number of top predictions (default t=20)
  * -s TRAININGSET set 1 to use [ccnb1,tpx2,aurka,cdc20,ccna2]; set 2 to use [mos,cdc6,slbp,prc1,btg4,cnot7,cnot8] or type your list of genes separated by commas and without spaces
  * -l LOG set to 1 to print log messages (default 0)
  * -h print help



  For instance:
  ```
  $ python src/run_Aner.py -d examples/Oocyte_datasets.csv -o outputs/aner_Oocytes -s 2 -l 1 -t 30

  $ python src/run_Aner.py -d examples/Oocyte_datasets.csv -o outputs/aner_Oocytes -s mos,cdc6,slbp,prc1 -l 1 -t 30
  ```

## Using PCA ranking

#### Running PCA ranking
Use the formatted dataset produced during Step 1 of the 'Using ANeR' section.
  
  ```
  $ python src/run_PCA.py -d examples/Oocyte_datasets.csv -o outputs/pca_Oocytes_topRanked.txt -s 2 -l 1 -t 30

  $ python src/run_PCA.py -d examples/Oocyte_datasets.csv -o outputs/pca_Oocytes_topRanked.txt -s mos,cdc6,slbp,prc1 -l 1 -t 30
  ```

#### Arguments

  * -d DATASET the dataset in CSV format (produced in step 1)
  * -o OUTPUTFILE is the path of the ouptut file
  * -t TOP number of top predictions (default t=20)
  * -s TRAININGSET set 1 to use [ccnb1,tpx2,aurka,cdc20,ccna2]; set 2 to use [mos,cdc6,slbp,prc1,btg4,cnot7,cnot8] or type your list of genes separated by commas and without spaces
  * -l LOG set to 1 to print log messages (default 0)
  * -h print help

## Comparing PCA ranking and ANeR ranking
You can comparing both appraches 

 ```
  $ python src/permutation.py -d <DATASET> -o <OUTPUTFILE> -p <PERMUTATIONS> -s <TRAININGSET> 

  ```
#### Arguments

  * -d DATASET the dataset in CSV format (produced in step 1)
  * -o OUTPUTFILE is the path of the ouptut files 
  * -p PERMUTATIONS number of permutations (default p=100)
  * -s TRAININGSET set 1 to use [ccnb1,tpx2,aurka,cdc20,ccna2]; set 2 to use [mos,cdc6,slbp,prc1,btg4,cnot7,cnot8] or type your list of genes separated by commas and without spaces
  * -l LOG set to 1 to print log messages (default 0)
  * -h print help


  For instance:
  ```
  $ python src/permutation.py -d examples/Oocyte_datasets.csv -o outputs/Oocytes -s 2 -l 1 -p 50
  ```

Two files will be produced:
* Oocytes_aner_50_permutations.csv
* Oocytes_pca_50_permutations.csv

## License, Patches, and Ongoing Developements
[Licence]([other_file.md](https://github.com/julibinho/ANER/blob/src/LICENCE)
