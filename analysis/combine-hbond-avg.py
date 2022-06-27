import pandas as pd
import os
import glob

#Get all the DAT files in a folder
path = os.getcwd()
dat_files = glob.glob(os.path.join(path, "*.dat")) #create a list of all .dat files in this dir

df_final = None # initilize an empty variable, this will be filled once in the next loop

# loop over the list of csv files
for f in dat_files: 
    

    ## fill the final varialbe with the first .dat file
    if df_final is None:
        header_list = ["Acceptor---Donor", "First"]
        df_final = pd.read_csv(f, delimiter=' ', names=header_list)
    # itterate through each dat file and merge it into an ever-growing "df_final"
    # the first column is the key used for merging 
    # files are sorted by the hbonds column before being merged
    header_list = ["Acceptor---Donor", f]
    df = pd.read_csv(f, delimiter=' ', names=header_list)
    df.sort_values(by='Acceptor---Donor', ascending=False)
    df_final = df_final.merge(df, on=['Acceptor---Donor'])
        

df_final.to_csv('final_hbond.csv', header=True)
