#Loading data through pandas

import pandas as pd
import numpy as np
import sys

import warnings
warnings.filterwarnings("ignore")
#----------------------------- Args and constants ----------------------


# File name of the CSV that you want to calculate the distance/objectives
outputFile   =   sys.argv[2]  #fileName = '3Cluster_Clean10Mar26.xlsx'
catFile    =   sys.argv[1]


#----------------------------- Functions -------------------------------
def process_excel_to_csv(inputFile, outputFile):
    fileName = pd.read_excel(inputFile)
    #print(fileName)

    #def process_excel_to_csv(input_file, fileName):
    # Load the Excel file into a DataFrame
    df = fileName

    # Separate the 'gene' column from the rest of the data
    gene_column = df["Gene"]
    data_columns = df.drop(columns=["Gene"])

    # Calculate the threshold for 60% missing values (excluding the 'gene' column)
    threshold = 0.6 * data_columns.shape[1]

    # Remove rows with 60% or more missing values in the data columns
    valid_rows = data_columns.dropna(thresh=threshold)
    valid_genes = gene_column[valid_rows.index]

    # Combine the valid 'gene' column with the valid data columns
    df_cleaned = pd.concat([valid_genes, valid_rows], axis=1)

    # Iterate over each row to fill missing values with the row's average (excluding the 'gene' column)
    for index, row in df_cleaned.iterrows():
        # Calculate the number of missing values in the row (excluding the 'gene' column)
        missing_count = row[1:].isnull().sum()

        # If more than 60% of values are present, fill missing values with the row's average
        if missing_count > 0 and missing_count / data_columns.shape[1] < 0.6:
            row_mean = row[1:].mean()  # Calculate the mean of the data columns only
            # Use .loc with explicit column labels to update the row
            #df_cleaned.loc[index, df_cleaned.columns[1:]] = row[1:].fillna(row_mean)
            df_cleaned.loc[index, df_cleaned.columns[1:]] = row[1:].fillna(0)
    # Save the processed DataFrame to a CSV file
    df_cleaned.to_csv(outputFile, index=False)

    print(f"Processed CSV saved to {outputFile}")
#----------------------------- Main -------------------------------
# Example usage
#TODO : check if file exists

print("It can take some time :-)")

# Query to search for the file in the specified folder
#query = f"'{folder_id}' in parents and title = '{file_name}' and trashed = false"

process_excel_to_csv(catFile, outputFile)

print ("done")

