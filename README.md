# Challenges and Progress in RNA Velocity: Comparative Analysis Across Multiple Biological Contexts
![Python](https://img.shields.io/badge/python-3.8-blue)

![image](https://github.com/user-attachments/assets/7a087fc7-a4a9-4d86-9967-246b88faf740)

*Figure 1 of the [manuscript](https://www.biorxiv.org/content/10.1101/2024.06.25.600667v1.full), portraying an introduction to RNA velocity, the methods analyzed, and the analyses comparing the methods.*

### Description
This project evaluates five RNA velocity estimation methods across three single-cell RNA sequencing (scRNA-seq) datasets. RNA velocity is a technique that models transcriptional dynamics to predict future cellular states. We analyze five RNA velocity methods: Velocyto, scVelo-stochastic, scVelo-dynamic modes, UniTVelo, and DeepVelo to provide a comprehensive comparison.

Our repository contains a collection of scripts and notebooks to analyze and run the five RNA velocity methods listed above and generate the corresponding figures as published. This repository is part of a [manuscript](https://www.biorxiv.org/content/10.1101/2024.06.25.600667v1.full).

### Organization
The repository is organized as follows, with directories corresponding to the analysis and visualization for each main figure and their corresponding supplementary.
Each folder contains the notebooks and scripts used for the following analyses:

- RNA Velocity methods and UMAP visualization (Figure 1)
- Neighborhood consistency within each method (Figure 2)
- Agreement between methods and with the median vector, as calculated from all methods (Figure 3)
- Overlap of driver genes as identified by CellRank (Figure 4)
- Robustness of each method to simulated sequencing depth, as calculated by subsetting the reads (Figure 5)
- Requirements to run each method and execute the analyses
  
The structure of this repo is illustrated below.

```plaintext
├── README.md
├── driver-genes-overlap #figure_4
│   └── CellRank_driver_genes_identification.ipynb
│   └── overlap_geneplot_drivers_UpSet.ipynb
├── method-agreement #figure_3
│   └── compute_median_vector_cosinesim_df.py
│   └── method_comparison_df_code.ipynb
│   └── visualization_median_vec_cosine_sim_comparison.ipynb
│   └── visualizations_method_comparison_direct_pairs.ipynb
├── neighborhood-consistency #figure_2
│   └── neighborhood_consistency_df.ipynb
│   └── visualization_neighborhood_consistency.ipynb
├── subset-reads-robustness #figure_5
│   └── combine_samples_adata.py
│   └── evaluate_velocity_df.py
│   └── run_allvelocity_methods.sh
│   └── subsample_bam_files_velocyto.sh
│   └── visualization_subsample-reads.ipynb
├── velocity-run-methods #figure_1
│   └── run_deepvelo.py
│   └── run_scvelo_dynamic.py
│   └── run_scvelo_stochastic.py
│   └── run_unitvelo.py
│   └── run_velocyto-scvelo_deterministic.py
│   └── visualization_multiple_umaps.ipynb
└── requirements
```

### Citation
If used please cite: 
```plaintext
@article{ancheta2024challenges,
  title={Challenges and Progress in RNA Velocity: Comparative Analysis Across Multiple Biological Contexts},
  author={Ancheta, Sarah and Dorman, Leah and Le Treut, Guillaume and Gurung, Abel and Royer, Lo{\"i}c A. and Granados, Alejandro and Lange, Merlin},
  journal={bioRxiv},
  pages={2024--06},
  year={2024},
  publisher={Cold Spring Harbor Laboratory}
}
```
