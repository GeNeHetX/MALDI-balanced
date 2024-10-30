library(Cardinal)
library(yaml)

# Read the hyperparameters from the config file
config <- read_yaml("config.yml")

# Hyperparameters
lames <- list.files(config$path_to_data)

# Initialize the peaks
peaks <- c()
for (lame in lames) {
  # Load the detected peaks imzml file
  mse_peaks <- readMSIData(sprintf("%s/%s/results/mse_densities.imzML",
                                   config$path_to_data, lame))

  # Change the data centroided attribute to TRUE
  centroided(mse_peaks) <- TRUE

  # Filter the peaks with low frequency
  mse_peaks <- subsetFeatures(mse_peaks, freq > 0.01)

  # Print the number of peaks
  print(sprintf("%s: %d peaks", lame, length(mz(mse_peaks))))

  # Concatenate the peaks
  peaks <- c(peaks, mz(mse_peaks))
}

# Sort the peaks
peaks <- sort(unique(peaks))

# Print the number of peaks
print(sprintf("Number of peaks: %d", length(peaks)))

# Perform the hierarchical clustering
clusters <- peaks |>
  dist() |>  # Compute the distance matrix
  hclust(method = "complete") |>  # Perform the hierarchical clustering
  cutree(h = config$clustering_tolerance)  # Tree cut with height = tolerance

# Get the mean of the clusters
reference <- tapply(X = peaks, INDEX = clusters, FUN = mean) |>
  as.numeric()

# Print the length of the reference
print(sprintf("Number of the reference peaks: %d", length(reference)))

# Save the reference peaks to a file
write.csv(reference, file = "reference_peaks.csv", row.names = FALSE)

# load the reference peaks
# reference <- read.csv("reference_peaks.csv")$x

for (lame in lames) {
  # print the lame
  print(lame)

  # Load the processed peaks imzml file
  mse_processed <- readMSIData(sprintf("%s/%s/results/mse_processed.imzML",
                                       config$path_to_data, lame))

  # Detect the peaks in the processed data with the reference
  mse_ref <- mse_processed |>
    peakProcess(ref = reference,
                method = config$peak_pick_method,
                SNR = config$signal_to_noise,
                tolerance = config$peak_pick_tolerance,
                units = config$units,
                filterFreq = config$filtered_frequency,
                BPPARAM = MulticoreParam())

  # Change the pixel date of the reference peaks to the densities peaks
  if (file.exists(sprintf("%s/%s/results/mse_ref.imzML",
                          config$path_to_data, lame))) {
    pData(mse_ref) <- readMSIData(sprintf("%s/%s/results/mse_densities.imzML",
                                          config$path_to_data, lame)) |>
      pData()
  }

  # Save the detected peaks with the reference
  writeMSIData(mse_ref, sprintf("%s/%s/results/mse_ref.imzML",
                                config$path_to_data, lame))
}

# Load the the first reference mse
mse_conc <- readMSIData(sprintf("%s/%s/results/mse_ref.imzML",
                                config$path_to_data, lames[1]))

# Change the centroided attribute to TRUE
centroided(mse_conc) <- TRUE

for (lame in lames[2:length(lames)]) {
  # Load the reference mse
  mse_ref <- readMSIData(sprintf("%s/%s/results/mse_ref.imzML",
                                 config$path_to_data, lame))

  # Change the centroided attribute to TRUE
  centroided(mse_ref) <- TRUE

  # Concatenate the reference mse
  mse_conc <- cbind(mse_conc, mse_ref)
}

# Save the concatenated reference mse
writeMSIData(mse_conc, "data/mse_concatenated.imzML")