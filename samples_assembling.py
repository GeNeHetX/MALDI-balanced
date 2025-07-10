import pandas as pd
import os
import gc
from tqdm import tqdm

path = 'data_external/new_lames'
path_results = 'data/MALDI_IHC/aligned_peaks'

# Extract the slides
slides = sorted(os.listdir(path))

# Initialize empty dataframes to store the results
peaks, pixels = pd.DataFrame(), pd.DataFrame()

# Iterate over the slides to load and concatenate the datasets
pbar = tqdm(slides, desc="Processing slides")
for slide in pbar:
    pbar.set_description(f"Processing slide {slide}")

    # Load and concatenate the datasets
    peaks_slide = pd.read_pickle(f"{path}/{slide}/results/peaks_aligned.pkl")
    pixels_slide = pd.read_feather(f"{path}/{slide}/results/mse_pixels.feather")

    # Exclude the pixels out of lesion
    peaks_slide = peaks_slide[pixels_slide['Density_Lesion'] > 0.5]
    pixels_slide = pixels_slide[pixels_slide['Density_Lesion'] > 0.5]

    # Remove the defects
    peaks_slide = peaks_slide[pixels_slide['Density_Defects'] < 0.1]
    pixels_slide = pixels_slide[pixels_slide['Density_Defects'] < 0.1]

    # Concatenate the dataframes
    peaks = pd.concat([peaks, peaks_slide], axis=0)
    pixels = pd.concat([pixels, pixels_slide], axis=0)

    # Clear memory
    del peaks_slide, pixels_slide
    gc.collect()

# Reset the index of the final dataframes
peaks.reset_index(drop=True, inplace=True)
pixels.reset_index(drop=True, inplace=True)

# Replace the NaN values with 0
pixels = pixels.fillna(0)

# Separate the microdissection columns
microdissections = pixels.filter(regex='Density_microdissection_')
pixels = pixels.drop(columns=microdissections.columns)

# Rename the columns to remove the _mse suffix
microdissections.columns = microdissections.columns.astype(str).str.replace('Density_microdissection_', '')

# For consistency, change the datatypes of the pixels and microdissection columns to float32
pixels = pixels.astype({col: 'float32' for col in pixels.columns if pixels[col].dtype == 'float64'})
microdissections = microdissections.astype('float32')

# Save the dataframes as pickle files
peaks.to_pickle(f'{path_results}/peaks.pkl')
pixels.to_pickle(f'{path_results}/pixels.pkl')
microdissections.to_pickle(f'{path_results}/microdissections.pkl')