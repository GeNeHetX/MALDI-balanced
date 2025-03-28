library(Cardinal)
library(feather)
library(yaml)

# Read the hyperparameters from the config file
config <- read_yaml("config.yml")

# Hyperparameters
lames <- list.files(config$path_to_data)

for (lame in lames) {
  print(lame)

  # Load the densities peaks imzml file
  mse_densities <- readMSIData(sprintf("%s/%s/results/mse_densities.imzML",
                               config$path_to_data, lame))

  # Save the pixels as a feather file
  print("Saving the pixels as a feather file")
  mse_densities |>
    pData() |>
    as.data.frame() |>
    write_feather(sprintf("%s/%s/results/mse_pixels.feather",
                          config$path_to_data, lame))

  # Remove the data from the memory
  print("Removing the data from the memory")
  rm(mse_densities)
  gc()
}