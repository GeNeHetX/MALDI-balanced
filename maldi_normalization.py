import pandas as pd
from utils import maldi_normalize

# Define the path and normalization method
path = 'data/MALDI_IHC/correlations/'
normalization_method = 'TIC'

# Load the MSI data
peaks = pd.read_pickle(f"{path}peaks.pkl")
pixels = pd.read_pickle(f"{path}pixels.pkl")

# Get the number of lames and densities
lames = pixels['run'].unique()

# Noramlize the peaks for each lame using the specified method
peaks_norm = peaks.copy()
for lame in lames:
    print(f"Normalizing {lame}")
    peaks_norm[pixels['run'] == lame] = maldi_normalize(peaks=peaks[pixels['run'] == lame],
                                                        method=normalization_method)

# Save to a pickle file
peaks_norm.to_pickle(f"{path}peaks_normalized_{normalization_method}.pkl")