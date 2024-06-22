#load conda environment deepVel6

#prepare counts layer of adata

import scanpy as sc
import torch
import scvelo as scv
import os
import numpy as np
import sys
from deepvelo.utils import velocity, update_dict
from deepvelo.utils.preprocess import autoset_coeff_s
from deepvelo import train, Constants

os.environ['CUDA_VISIBLE_DEVICES'] = '-1'

# fix random seeds for reproducibility
SEED = 123
torch.manual_seed(SEED)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
np.random.seed(SEED)

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.set_figure_params(
    "scvelo", transparent=False
)  # for beautified visualization


adata_path = sys.argv[1]
save_path = sys.argv[2]

adata_fresh = sc.read_h5ad(adata_path)

adata = adata_fresh.copy()


scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_neighbors=30, n_pcs=30)

sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)

# specific configs to overide the default configs
configs = {
    "name": "DeepVelo", # name of the experiment
    "loss": {"args": {"coeff_s": autoset_coeff_s(adata)}},
    "trainer": {"verbosity": 0}, # increase verbosity to show training progress
}
configs = update_dict(Constants.default_configs, configs)


device = torch.device("cuda:0")

# initial velocity
velocity(adata, mask_zero=False)
trainer = train(adata, configs)

scv.tl.velocity_graph(adata)

adata.write(save_path)