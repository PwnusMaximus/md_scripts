'''
This command is indented to be issued on a directory with hbond-avg.dat files from cpptraj hbond analysis. 
AFTER this series of bash commands is issued: 
    for file in *avg.dat; do awk '{ print $1 "---" $3 " " $5}' $file > cut$file; done
    mkdir originals
    find . -type f -not -name '*cut*' -exec mv {} originals/ \;
    
The intended purpose is to match occupancies of hydrogen bonds accross systems and across replicates

The basic workflow is sort the hydrogen bond names
then to append column #2 (Occupancy %) to the right for each additional hbond-avg.dat file 
and finally to output a .csv file for further analysis in excel

Author: Makay Murray 2022
'''

import pandas as pd
import os
import glob

df_final = None # initilize an empty variable, this will be filled once in the next loop

#Get all the DAT files in a folder
path = os.getcwd()
dat_files = glob.glob(os.path.join(path, "*.dat")) #create a list of all .dat files in this dir

# loop over the list of csv files
for f in dat_files: 
   
    if df_final is None:   ## fill the final variable with the first .dat file
        df_final = pd.read_csv(f, delim_whitespace=True, usecols = ['#Acceptor', 'Donor', 'Frac'])
        df_final['hbond'] =  df_final["#Acceptor"].astype(str) + str("---") + df_final["Donor"].astype(str)
        df_final = df_final.sort_values(by=['hbond'])
        df_final.pop('#Acceptor')
        df_final.pop('Donor')
        df_final = df_final.reindex(columns=['hbond','Frac'])
        
    # itterate through each dat file and merge it into an ever-growing "df_final"
    # the first column is the key used for merging 
    # files are sorted by the hbonds column before being merged
    df = pd.read_csv(f, delim_whitespace=True, usecols = ['#Acceptor', 'Donor', 'Frac'])
    df['hbond'] =  df["#Acceptor"].astype(str) + str("---") + df["Donor"].astype(str)
    df = df.sort_values(by=['hbond'])
    df.pop('#Acceptor')
    df.pop('Donor')
    df.sort_values(by=['hbond'])
    df = df.reindex(columns=['hbond','Frac'])
    df_final = df_final.merge(df, how='left', on='hbond',)
        
#df.drop(df.columns[[2]], axis = 1, inplace = True)
df_final.to_csv('final_hbond.csv', header=True)
