import numpy as np
import pandas as pd
import os

path = 'data_hdd/'
path_results = 'data/MALDI_IHC/correlations/'

# Extract the lames
lames = sorted(os.listdir(path))

# Load and concatenate the datasets
peaks = pd.concat([pd.read_feather(f"{path}{lame}/results/mse_peaks_ref_corr.feather").astype('float32') for lame in lames])
pixels = pd.concat([pd.read_feather(f"{path}{lame}/results/mse_pixels_corr.feather") for lame in lames])

# Reset the index
peaks.reset_index(drop=True, inplace=True)
pixels.reset_index(drop=True, inplace=True)

# Replace the NaN values with 0
pixels = pixels.fillna(0)

# Separate the microdissection columns
microdissections = pixels.filter(regex='Density_microdissection_')
pixels = pixels.drop(columns=microdissections.columns)

# Rename the columns to remove the _mse suffix
microdissections.columns = microdissections.columns.str.replace('Density_microdissection_', '')

# For consistency, change the datatypes of the pixels and microdissection columns to float32
pixels = pixels.astype({col: 'float32' for col in pixels.columns if pixels[col].dtype == 'float64'})
microdissections = microdissections.astype('float32')

# Save the dataframes as pickle files
peaks.to_pickle(f'{path_results}peaks.pkl')
pixels.to_pickle(f'{path_results}pixels.pkl')
microdissections.to_pickle(f'{path_results}microdissections.pkl')