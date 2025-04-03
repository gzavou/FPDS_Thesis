#!/usr/bin/env python3

#%%## score_variants output inspection ###
# run with enrollhd conda environment

import h5py
import pandas as pd

output_path = "/pool01/code/src/genomic_llms/borzoi/tutorials/latest/score_variants/snp_sed/f0c0/sed.h5"

#%%## Open the .h5 file ###
f = h5py.File(output_path,'r')
# Inspect the contents of the .h5 file
print("Keys in the HDF5 file:", list(f.keys()))

# Optionally, inspect the structure of each key
for key in f.keys():
    print(f"Contents of key '{key}':")
    print(f[key])

#%% Extract contribution (from borzoi run_variant_scripts.ipynb)
# Print an example variant effect prediction for a SNP-gene pair (gene-specific expression)

sed_h5 = h5py.File(output_path, 'r')

row_ix = 23
target_ix = 0

print("score: 'logSED', snp: '" + str(sed_h5['snp'][sed_h5['si'][row_ix]].decode()) + "', gene: '" + str(sed_h5['gene'][sed_h5['si'][row_ix]].decode()) + "', track: '" + str(sed_h5['target_labels'][target_ix].decode()) + "' => " + str(round(sed_h5['logSED'][row_ix, target_ix], 4)))
# %%
