import numpy as np
import pandas as pd

from utils import maldi_normalize

# Load the MSI data
peaks = pd.read_pickle("data/MALDI_IHC/peaks.pkl")
pixels = pd.read_pickle("data/MALDI_IHC/pixels.pkl")

# Get the number of lames and densities
lames = pixels['run'].unique()

# Noramlize the peaks using PQN
peaks_norm = peaks.copy()
for lame in lames:
    peaks_norm[pixels['run'] == lame] = maldi_normalize(peaks[pixels['run'] == lame])

# Save to a pickle file
peaks_norm.to_pickle("data/MALDI_IHC/peaks_normalized.pkl")