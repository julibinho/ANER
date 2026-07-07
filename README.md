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


## Using ANeR 
   
#### Step 1: Preparing data
   
  ```
  $ python src/prepareData.py <INPUT_FILE> <OUTPUT_FILE>
  ```
#### Required arguments 
  * <INPUT_FILE> is a xlsx file containing the datasets
  * <OUTPUT_FILE> is the path of the ouptut file

  For instance:
  ```
  $python src/prepareData.py examples/Oocyte_datasets.xlsx examples/Oocyte_datasets.csv
  ```

#### Step 2: Running ANeR
   
  ```
  $ python src/run_Aner.py -d <DATASET> -o <OUTPUTFILE> -t <TOP> -s <TRAININGSET>
  ```

#### Arguments

  * DATASET the dataset in CSV format (produced in step 1)
  * OUTPUTFILE is the path of the ouptut file
  * TOP  number of top predictions (default t=20)
  * TRAININGSET set 1 to use [ccnb1,tpx2,aurka,cdc20,ccna2]; set 2 to use [mos,cdc6,slbp,prc1,btg4,cnot7,cnot8] or type your list of genes separated by commas and without spaces



  For instance:
  ```
  $python src/run_Aner.py -d examples/Oocyte_datasets.csv -o outputs/aner_Oocytes -s 2 -l 1 -t 30

  $python src/run_Aner.py -d examples/Oocyte_datasets.csv -o outputs/aner_Oocytes -s mos,cdc6,slbp,prc1 -l 1 -t 30
  ```

