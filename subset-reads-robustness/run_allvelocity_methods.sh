""" 
This script creates scripts to run all the velocity methods for each of the adata objects. 
It will activate each method's conda environment and then runs the method in a screen session.
"""
#written to run on katherine-johnson

SAMPLE_PATH="/subset_adatas/full_adatas_subset/"
DESTDIR="/velocity_subsets_zf24hpf"
WDIR="/scratch/fish24hrs_velocity"
cd $WDIR

CONDA_PATH="/miniconda3/bin/conda"

# Loop through the directory to find samples. Assume the directory contains multiple folders, each of which is one sample. 
for FULLSAMPLE in $SAMPLE_PATH/allfish_*.h5ad; do
    # Check if the file exists to prevent errors when no matching files are found
    [ -e "$FULLSAMPLE" ] || continue

    # Pull out the sample name from the full path
    SAMPLE=${FULLSAMPLE##*/}
    SAMPLENAME=${SAMPLE%".h5ad"}

    # Now you can process your sample with the variable $SAMPLENAME
    echo "Processing $SAMPLENAME"


#write a bash file for each sample that runs all 5 packages sequentially
    cat > $SAMPLENAME.sh <<EOF

#!/bin/bash

# RUN SCRIPT
#paths correspond to running on katherine-johnson

source /miniconda3/etc/profile.d/conda.sh

#scVelo- Stochastic, Dynamical, Deterministic (Velocyto)
conda activate /miniconda3/envs/vel38zebra

#if the directory does not exist (-d), make it recursively (-p)
[ -d "$DESTDIR/scvelo-stochastic" ] || mkdir -p "$DESTDIR/scvelo-stochastic/"
python run_scvelo_stochastic.py $FULLSAMPLE $DESTDIR/scvelo-stochastic/$SAMPLENAME-scvelo-stochastic.h5ad


[ -d "$DESTDIR/scvelo-dynamical" ] || mkdir -p "$DESTDIR/scvelo-dynamical/"
python run_scvelo_dynamic.py $FULLSAMPLE $DESTDIR/scvelo-dynamical/$SAMPLENAME-scvelo-dynamical.h5ad

[ -d "$DESTDIR/scvelo-deterministic" ] || mkdir -p "$DESTDIR/scvelo-deterministic/"
python run_velocyto-scvelo_deterministic.py $FULLSAMPLE $DESTDIR/scvelo-deterministic/$SAMPLENAME-scvelo-deterministic.h5ad

conda deactivate

#DeepVelo
conda activate /miniconda3/envs/deepVel6
#if the directory does not exist (-d), make it recursively (-p)
[ -d "$DESTDIR/deepvelo" ] || mkdir -p "$DESTDIR/deepvelo"
python run_deepvelo.py $FULLSAMPLE $DESTDIR/deepvelo/$SAMPLENAME-deepvelo.h5ad

conda deactivate


#UniTVelo
conda activate /miniconda3/envs/unitvelo2
#if the directory does not exist (-d), make it recursively (-p)
[ -d "$DESTDIR/uniTvelo" ] || mkdir -p "$DESTDIR/uniTvelo"
python run_unitvelo.py $FULLSAMPLE $DESTDIR/uniTvelo/$SAMPLENAME-uniTvelo.h5ad 'clusters'

conda deactivate



EOF

#-dm will starts a new detached session 
#-S means the name of the session
#-c means read some file
screen -dm -S $SAMPLENAME bash -c "bash $SAMPLENAME.sh"



done


