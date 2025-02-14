import pandas as pd
import os

path = 'data_hdd/'

# Extract the lames
lames = sorted(os.listdir(path))

# Extract the lames
lames = sorted(os.listdir(path))

# Load and concatenate the datasets
hes_features = pd.concat([pd.read_pickle(f"{path}{lame}/results/hes_features.pkl") for lame in lames])

# Reset the index
hes_features.reset_index(drop=True, inplace=True)

# Save the dataframes as pickle files
hes_features.to_pickle('data/MALDI_IHC/hes_features.pkl')