import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score, mean_squared_error, mean_poisson_deviance
import os

def maldi_normalize(peaks: pd.DataFrame, method: str='PQN', quantile: float=0.99) -> pd.DataFrame:
    """
    Performs Normalization on a MALDI-MSI dataframe where pixels are rows and M/Z are the columns.

    Args:
    peaks (pd.DataFrame): A pandas DataFrame with m/z values as columns and pixels as rows.
    method (str): A string that specifies the normalization method to use. Can be 'PQN', 'TIC', 'Mean', 'Median', 'RMS', 'MAX' or 'Quantile'. (Default is 'PQN')
    quantile (float): A float that specifies the quantile to use when method is 'Quantile'. (Default is 0.99)

    Returns:
    A normalized DataFrame.
    """
    ## Define the normalization methods

    # Probabilistic Quotient Normalization
    def pqn_normalization(peaks):
        median_spectrum = peaks.median(axis=0)  # Calculate the median spectrum and the median of each pixel
        quotient_ratios = peaks.div(median_spectrum, axis=1)  # Calculate the quotient ratios
        median_quotients = quotient_ratios.median(axis=1)   # Calculate the median of the quotient ratios
        return peaks.div(median_quotients, axis=0)

    # Total Ion Current Normalization
    def tic_normalization(peaks):
        sum_spectrum = peaks.sum(axis=1)  # Calculate the sum spectrum
        return peaks.div(sum_spectrum, axis=0)

    # Mean Normalization
    def mean_normalization(peaks):
        mean_spectrum = peaks.mean(axis=1)  # Calculate the mean spectrum
        return peaks.div(mean_spectrum, axis=0)

    # Median Normalization
    def median_normalization(peaks):
        median_spectrum = peaks.median(axis=1)  # Calculate the median spectrum
        return peaks.div(median_spectrum, axis=0)

    # Root Mean Square Normalization
    def rms_normalization(peaks):
        rms_spectrum = np.sqrt(np.mean(peaks**2, axis=1))  # Calculate the RMS spectrum
        return peaks.div(rms_spectrum, axis=0)

    # Max Normalization
    def max_normalization(peaks):
        max_spectrum = peaks.max(axis=1)  # Calculate the max spectrum
        return peaks.div(max_spectrum, axis=0)

    # Quantile Normalization
    def quantile_normalization(peaks, quantile):
        quantile_spectrum = peaks.quantile(q=quantile, axis=1)  # Calculate the quantile spectrum
        return peaks.div(quantile_spectrum, axis=0)

    normalization_methods = {
        "PQN": pqn_normalization,
        "TIC": tic_normalization,
        "Mean": mean_normalization,
        "Median": median_normalization,
        "RMS": rms_normalization,
        "MAX": max_normalization,
        "Quantile": lambda peaks: quantile_normalization(peaks, quantile)
    }

    if method not in normalization_methods:
        raise ValueError(f"Normalization method '{method}' is not supported.")

    return normalization_methods[method](peaks)


def model_metrics_rms(model, X_train: np.ndarray, y_train: np.ndarray, X_test: np.ndarray, y_test: np.ndarray, y_train_rms: np.ndarray, y_test_rms: np.ndarray, show_plot: bool=True, print_metrics: bool=True) -> tuple:
    """Computes the metrics of a model trained on the RMS values of the data.

    Args:
        model (sklearn model): sklearn model
        X_train (np.ndarray): training data
        y_train (np.ndarray): training labels
        X_test (np.ndarray): testing data
        y_test (np.ndarray): testing labels
        y_train_rms (np.ndarray): training labels RMS
        y_test_rms (np.ndarray): testing labels RMS
        show_plot (bool, optional): plot a histogram of the label distribution. Defaults to True.
        print_metrics (bool, optional): print all the metrics. Defaults to True.

    Returns:
        tuple: r2_train, r2_test, mse_train, mse_test, deviance_train, deviance_test
    """
    # Make predictions
    y_pred_test = model.predict(X_test)
    y_pred_train = model.predict(X_train)

    # Compute R-squared
    r2_train = r2_score(y_train_rms, y_pred_train)
    r2_test = r2_score(y_test_rms, y_pred_test)

    if print_metrics:
        print(f"R-squared (Train): {r2_train:.2f}")
        print(f"R-squared (Test): {r2_test:.2f}")

    # Square the predictions
    y_pred_test = y_pred_test ** 2
    y_pred_train = y_pred_train ** 2

    # Clip the negative values
    y_pred_test = np.clip(a=y_pred_test, a_min=0, a_max=1)
    y_pred_train = np.clip(a=y_pred_train, a_min=0, a_max=1)

    # Compute the mean squared error
    mse_train = mean_squared_error(y_train, y_pred_train)
    mse_test = mean_squared_error(y_test, y_pred_test)

    if print_metrics:
        print(f"Mean Squared Error (Train): {mse_train:.2e}")
        print(f"Mean Squared Error (Test): {mse_test:.2e}")

    # Compute the deviance
    deviance_train = mean_poisson_deviance(y_train, y_pred_train)
    deviance_test = mean_poisson_deviance(y_test, y_pred_test)

    if print_metrics:
        print(f"Deviance (Train): {deviance_train:.2e}")
        print(f"Deviance (Test): {deviance_test:.2e}")

    # Plot a histogram of the predictions
    if show_plot:
        fig, ax = plt.subplots(1, 1, figsize=(5, 3), tight_layout=True)
        bins = np.linspace(min(y_test.min(), y_pred_test.min()), max(y_test.max(), y_pred_test.max()), 500)
        ax.hist(y_test, bins=bins, alpha=0.5, label='Data')
        ax.hist(y_pred_test, bins=bins, alpha=0.5, label='Predictions')
        ax.set_title('Predictions')
        # ax.set_yscale('log')
        ax.set_xlim(0, 0.05)
        ax.legend(fontsize=9)
        plt.show()

    return r2_train, r2_test, mse_train, mse_test, deviance_train, deviance_test


def proteins_share(path_peptides: str, mass_arr: np.ndarray, tolerance: float=0.2) -> dict:
    """Extracts the proteins shared masses with the mass array.

    Args:
        path_peptides (str): path to the proteins peptides masses csv files
        mass_arr (np.ndarray): an array with the masses
        tolerance (float): The tolerance to consider a mass shared between the proteins. (Default is 0.2)

    Returns:
        dict: A dictionary with the proteins names as keys and the shared masses as values.
    """
    # Get the names of the proteins
    protein_names = [name.split('_')[0] for name in os.listdir(path_peptides)]

    # Initialize the dictionary 
    proteins_share = {name: [] for name in protein_names}

    # Iterate over the proteins
    for protein_name in protein_names:

        # Load the masses of the protein
        protein_mass = pd.read_csv(f"{path_peptides}/{protein_name}_HUMAN_peptide_mass.csv")['mass'].values

        # Iterate over the masses of the protein
        for mass in protein_mass:

            # Check if the mass is in the mass array with a tolerance
            if np.min(np.abs(mass - mass_arr)) < tolerance:

                # Append the mass to the list of masses
                proteins_share[protein_name].append(str(mass))

    # Order the proteins by the number of masses
    proteins_share = {k: v for k, v in sorted(proteins_share.items(), key=lambda item: len(item[1]), reverse=True)}
    return proteins_share


def create_protein_peptide_map(path_peptides: str, type=float) -> dict:
    """
    Creates a mapping of protein names to their corresponding peptide masses.

    Args:
        path_peptides (str): Path to the directory containing peptide mass CSV files.

    Returns:
        dict: A dictionary where keys are protein names and values are lists of peptide masses.
    """
    # Extract protein names from the filenames in the specified directory
    protein_names = [name.split('_')[0] for name in os.listdir(path_peptides)]
    
    # Create a dictionary mapping each protein name to its list of peptide masses
    protein_peptide_map = {
        protein_name: pd.read_csv(f"{path_peptides}/{protein_name}_HUMAN_peptide_mass.csv")['mass'].astype(type).to_list()
        for protein_name in protein_names
    }
    
    return protein_peptide_map


def match_protein_peptide_map(protein_peptide_map: dict, match_list: np.ndarray, tolerance: float=0.1) -> dict:
    """
    Updates the values in protein_peptide_map to match the mz values in match_list within a given tolerance.

    Args:
        protein_peptide_map (dict): A dictionary where keys are protein names and values are lists of peptide masses.
        match_list (np.ndarray): An array containing mz values to match.
        tolerance (float): The tolerance within which mz values are considered a match.

    Returns:
        dict: Updated protein_peptide_map with values replaced by matching mz values from match_list.
    """
    updated_map = {}

    mz_values = match_list.astype(float)

    for protein in protein_peptide_map.keys():
        updated_map[protein] = protein_peptide_map[protein][:]
        for i, mz in enumerate(updated_map[protein]):
            if mz in mz_values:
                continue
            elif any(abs(float(mz) - float(mz_value)) <= tolerance for mz_value in mz_values):
                updated_map[protein][i] = str(mz_values.values[np.argmin(np.abs(mz_values - float(mz)))])

    return updated_map