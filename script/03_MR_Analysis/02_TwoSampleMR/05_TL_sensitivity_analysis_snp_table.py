#############################################
### Import libraries ########################
import pandas as pd
import os

#############################################
### Set working directory ###################
wd_path = "../"
export_path = "../"

#############################################
### List all files in folder ################
ma_file_names = []
for file in os.listdir(wd_path):
    if file.endswith(".ma"):
        ma_file_names.append(os.path.join(wd_path, file))

print("List of files stored")

#############################################
### Read rsdid file #########################
rsid_file = open(export_path + "TL_IVS_345.txt", "r") ### ! Used export path here
file_content = rsid_file.read()
rsids = file_content.split("\n")
rsid_file.close()
rsids = rsids[0:len(rsids)-1]

print("List of rsids stored")
#############################################
### Generate pvalue file ####################

### Make empty pandas dataframe
df = pd.DataFrame(rsids, columns=['SNP'])

### Loop over all files and extract SNP pvalues
for file in ma_file_names:
    ### Create temporary dataframe
    df_temp = pd.DataFrame(columns=['SNP','p'])
    ### Read by chuncksizes
    chunksize = 10 ** 6
    reader = pd.read_csv(file, chunksize=chunksize, sep = '\t')
    for chunk in reader: ### Merge with SNPs found
        df_temp = pd.concat([df_temp, chunk.loc[chunk['SNP'].isin(df['SNP'])][['SNP','p']]], ignore_index=True)
    ### Merge dataframe
    df = df.merge(df_temp, how='left', on='SNP')
    ### Rename new column
    newcolname = file.split('/')[-1].split('_')[0]
    df.rename(columns={'p':newcolname}, inplace=True)
    print("Added:",newcolname)

#############################################
### Save file ###############################
df.to_csv(export_path + '/TL_IV_snp_pvalues.txt', header=True, index=None, sep=' ', mode='w')
print("Dataframe saved")