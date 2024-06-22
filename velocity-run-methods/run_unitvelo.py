# conda activate unitvelo2

import scvelo as scv
import unitvelo as utv
import os
import tensorflow as tf
import scanpy as sc
import sys

scv.settings.verbosity = 0

os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = "0"

config = tf.compat.v1.ConfigProto()
sess = tf.compat.v1.Session(config=config)

print ("Num GPUs Available: ", len(tf.config.list_physical_devices('GPU')))

velo = utv.config.Configuration()
velo.R2_ADJUST = True 
velo.IROOT = None
velo.FIT_OPTION = '1'
velo.GPU = 0


adata_path = sys.argv[1]
save_path = sys.argv[2]

#label cluster should be a cluster variable in the adata.obs
label_cluster = sys.argv[3]

adata_fresh = sc.read_h5ad(adata_path)

adata = adata_fresh.copy()

if 'neighbors' in adata.uns:
    del adata.uns['neighbors']
    full_path = '/adata_no_neighbors.h5ad'
    adata.write(full_path)
    dataset = full_path

else:
    dataset = adata_path

label = label_cluster
exp_metrics = {}

velo_config = utv.config.Configuration()
velo_config.R2_ADJUST = True
velo_config.IROOT = None
velo_config.FIT_OPTION = '1'
velo_config.AGENES_R2 = 1


adata = utv.run_model(dataset, label, config_file=velo_config)


adata.write(save_path)