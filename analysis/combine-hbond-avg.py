'''    
COMBINE HYDROGEN BOND AVERAGES SCRIPT
"Align occupancies of hydrogen bonds to their name accross similar systems or across replicates"

How to use:
    1. fill a dir with hbond-avg.dat files from cpptraj (for mutated or replicate systems where hbonds will align)
    2. Ensure the names of the hbond files are meaningful (system-1.dat, system-2.dat etc.)
    3. Run this script with "combine-hbond-avg.py" (AFTER installing in your ~/bin) 
    4. take the output .csv file for futher inspection elsewhere
    
Workflow of the script: 
    Import the Acceptor, Donor, and Frac columns from hbond-avg.dat files (from cpptraj)
    Merge Acceptor and Donor columns into a single field "hbond"
    Sort hbond by string
    Clean, reorder and rename the Frac column to ref the source file(s)
    Combine n number of additional .dat files and merge until complete
    Output a .csv file of the 'final' data frame for further analysis in excel

Installation instructions for cedar users:
	1. copy this file "combine-hbond-avg.py" to your ~/bin dir
	2. bash: chmod +x ~/bin/combine-hbond-avg.py
	3. done 
	
Modules for Cedar users:
    module load StdEnv/2020 scipy-stack/2022a

Author: Makay Murray 2022
Lab: Wetmore @ uLeth
'''

import pandas as pd
import os
import glob
#%% 
df_final = None # initilize an empty df, this will be filled in the next loop
#%%
#Get all the .dat file names in current path
path = os.getcwd()
dat_files = glob.glob(os.path.join(path, "*.dat"))
#%%
# loop over the list of dat files
for f in dat_files: 
   
    # Fill the data frame with the first dat file and adjust columns
    if df_final is None:
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
#%%
df_final.drop(df_final.columns[[1]], axis = 1, inplace = True)
df_final.to_csv('final_hbond.csv', header=True)

print('The Hbond Data has been combined! Please check final_hbond.csv for your data')
#%%