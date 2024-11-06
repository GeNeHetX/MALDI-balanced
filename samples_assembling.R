library(Cardinal)
library(feather)

# Hyperparameters
lames <- list.files("data_hdd")

# Load the first reference mse
mse_conc <- readMSIData(sprintf("data_hdd/%s/results/mse_ref.imzML", lames[1]))

# Change the centroided attribute to TRUE
centroided(mse_conc) <- TRUE

# Clean the pixels
mse_conc <- mse_conc |>
  subsetPixels(Density_Defects < 0.1)  # Filter the pixels with defects
  # subsetPixels(Density_Lesion > 0.5,  # Filter the pixels out of the lesion
  #              Density_Defects < 0.1)  # Filter the pixels with defects

for (lame in lames[2:length(lames)]) {
  # Load the reference mse
  mse_ref <- readMSIData(sprintf("data_hdd/%s/results/mse_ref.imzML", lame))

  # Change the centroided attribute to TRUE
  centroided(mse_ref) <- TRUE

  # Clean the pixels
  mse_ref <- mse_ref |>
    subsetPixels(Density_Defects < 0.1)  # Filter the pixels with defects
    # subsetPixels(Density_Lesion > 0.5,  # Filter the pixels out of the lesion
    #              Density_Defects < 0.1)  # Filter the pixels with defects

  # Concatenate the reference mse
  mse_conc <- cbind(mse_conc, mse_ref)
}

# Save the concatenated reference mse
writeMSIData(mse_conc, "data/mse_conc_complete.imzML")

# Load the pixels into memory
pixels <- pData(mse_conc) |>
  as.data.frame()

# Save the pixels as a feather file
pData(mse_conc) |>
  as.data.frame() |>
  write_feather("data/mse_conc_pixels_complete.feather")

# Load the peaks with the intensity values to 3 decimal into memeory
peaks <- spectra(mse_conc) |>
  as.matrix() |>
  round(3)

# Change the row names to the m/z values
rownames(peaks) <- mz(mse_conc) |>
  round(4)

# Save the peaks as a feather file
peaks |>
  t() |> # Transpose the matrix
  as.data.frame() |> # Convert to a data frame
  write_feather("data/mse_conc_peaks_complete.feather")