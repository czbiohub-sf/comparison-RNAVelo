"""
This script processes and integrates loom files from multiple TDR (fish) samples.
It combines the loom files for a single subset and iteration to create one adata object. 
"""

# Load environment vel38zebra 

import loompy
import os
import scvelo as scv
import pandas as pd

destdir="/velocity_subset_take2_extra/"
savedir="/full_adatas_subset/"
destdir_barcode = "/velocity_subset_take2/"

#only include the barcodes (cells) that have been included in the original full/clean anndata object
whitelist_barcodes=pd.read_csv(destdir_barcode + 'whitelist_barcodes.csv', index_col=0)


fishfolders = []
#find TDR (fish) samples 
for folder in os.listdir(destdir):
    if folder.startswith('TDR'):
        fishfolders.append(folder)



#input the subset of fraction of reads to concat samples together
fraction=0.98

#we ran five iterations for each fraction of reads
for iteration in range(1,6):
    frit=str(fraction)+'-'+str(iteration)+'.loom'
    file1 = [i for i in os.listdir(destdir+fishfolders[0]) if i.split("-",maxsplit=1)[1]==frit][0]
    file2 = [i for i in os.listdir(destdir+fishfolders[1]) if i.split("-",maxsplit=1)[1]==frit][0]
    file3 = [i for i in os.listdir(destdir+fishfolders[2]) if i.split("-",maxsplit=1)[1]==frit][0]
    file4 = [i for i in os.listdir(destdir+fishfolders[3]) if i.split("-",maxsplit=1)[1]==frit][0]
    
    files=[destdir+fishfolders[0]+'/'+file1,destdir+fishfolders[1]+'/'+file2,destdir+fishfolders[2]+'/'+file3,destdir+fishfolders[3]+'/'+file4]
    newfile=destdir+"allfish_"+frit

    #combine loom files for the four fish 
    loompy.combine(files,newfile,key="Accession")
    adata  = scv.read_loom(newfile)

    #include information in adata object about the cell barcode, fish/sample ID
    adata.obs['new_index'] = adata.obs.index.copy()
    adata.obs['fish'] = [i.split("_")[0] for i in adata.obs['new_index']]
    adata.obs['cellbarcode'] = [i.split(":")[1].split("x")[0] for i in adata.obs['new_index']]
    adata.obs['my_ID'] = adata.obs['fish'].astype(str) + "_" + adata.obs['cellbarcode'].astype(str) + "-1"
    adata.obs.set_index('my_ID', inplace=True)

    #subset for barcodes in the full, cleaned adata object
    adata=adata[adata.obs.index.isin(whitelist_barcodes.index.tolist())].copy()
    adata.obs['clusters'] = whitelist_barcodes.clusters.copy()
    adata.write_h5ad(savedir + 'allfish_' + str(fraction) + '-' + str(iteration)+ '.h5ad')