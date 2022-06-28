'''    
The intended purpose is to match occupancies of hydrogen bonds accross systems and across replicates

This command is indented to be issued on a directory with hbond-avg.dat files from cpptraj

Objective of this software
Import the Acceptor, Donor, and Frac columns from hbond avg files
Align and combine based on the hydrogen bond identity
Clean and rename the columns to refer to the source files
finally to output a .csv file for further analysis in excel

Author: Makay Murray 2022
Lab: Wetmore @ uLeth
'''

import pandas as pd
import os
import glob

df_final = None # initilize an empty df, this will be filled once in the next loop

#Get all the DAT files in a folder
path = os.getcwd()
dat_files = glob.glob(os.path.join(path, "*.dat"))

# loop over the list of dat files
for f in dat_files: 
   
    # Fill the data frame with the first dat file and adjust columns
    if df_final is None:   # fill the final variable with the first .dat file
        df_final = pd.read_csv(f, delim_whitespace=True, usecols = ['#Acceptor', 'Donor', 'Frac'])
        df_final['hbond'] =  df_final["#Acceptor"].astype(str) + str("---") + df_final["Donor"].astype(str)
        df_final = df_final.sort_values(by=['hbond'])
        df_final.pop('#Acceptor')
        df_final.pop('Donor')
        df_final = df_final.reindex(columns=['hbond','Frac'])
    #append the remaining hbond files in the dir into the df_final data frame. 
    df = pd.read_csv(f, delim_whitespace=True, usecols = ['#Acceptor', 'Donor', 'Frac'])
    df['hbond'] =  df["#Acceptor"].astype(str) + str("---") + df["Donor"].astype(str)
    df = df.sort_values(by=['hbond'])
    df.pop('#Acceptor') #delete the #Acceptor column
    df.pop('Donor') #Delete the Donor Column
    df = df.reindex(columns=['hbond','Frac'])
    df = df.rename(columns = {"Frac":f}) #convert the Frac column header into the file name the data is from. 
    df_final = df_final.merge(df, how='left', on='hbond',) #final merge of the loop to append to the final DF

df_final.drop(df_final.columns[[1]], axis = 1, inplace = True)
df_final.to_csv('final_hbond.csv', header=True)

print('The Hbond Data has been combined! Please check final_hbond.csv for your data')
