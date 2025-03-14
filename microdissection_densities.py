#!/usr/local/bin/python

import numpy as np
import pandas as pd
import os
import PIL
import gc

from skimage import io

import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme(style="darkgrid")

# Increase the limit of allowed images size
PIL.Image.MAX_IMAGE_PIXELS = 10e10

# Define the path to the data
path = 'data_external'

# List the lames
lames = sorted(os.listdir(f'{path}'))

# Create an empty list to store the data
data = []

# Loop over the lames
for lame in lames:
    print(lame)

    # Read the density masks
    densities = sorted([mask for mask in os.listdir(f'{path}/{lame}/results/masks')
                        if 'microdissection' not in mask
                            and 'compressed' not in mask])

    # Read the microdissection masks
    microdissections = sorted([mask for mask in os.listdir(f'{path}/{lame}/results/masks')
                               if 'microdissection' in mask
                                   and 'compressed' not in mask])

    # Loop over the microdissections
    for microdissection in microdissections:
        print(microdissection)
        
        # Load the microdissection mask image
        microdissection_mask = io.imread(f'{path}/{lame}/results/masks/{microdissection}')
        
        # Loop over the densities
        for density in densities:
            print(density)
            
            # Load the density mask image
            density_mask = io.imread(f'{path}/{lame}/results/masks/{density}')
            
            # Compute the density under the microdissection
            density_under_microdissection = np.sum(microdissection_mask * density_mask) / np.sum(microdissection_mask)

            # Append the data
            data.append({'microdissection': '_'.join(microdissection.split('_')[2:-1]),
                         'density': '_'.join(density.split('_')[1:-1]),
                         'density_under_microdissection': density_under_microdissection})
            
            # Clear the memory
            del density_mask
            gc.collect()
        
        # Clear the memory
        del microdissection_mask
        gc.collect()

# Create a dataframe from the data
df = pd.DataFrame(data)

# Pivot the dataframe to have the densities as columns
df_pivot = df.pivot(index='microdissection', columns='density', values='density_under_microdissection').reset_index()

# Remove the 'density' prefix from the column names
df_pivot.columns = [col.replace('density', '') for col in df_pivot.columns]

# Save the dataframe
df_pivot.to_csv('data/MALDI_IHC/microdissection_densities.csv', index=False)