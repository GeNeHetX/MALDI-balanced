import numpy as np
import pandas as pd

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