##################################
# Export NDVIs from NASA webpage
# ChatGPT, Jelen & Pablo
##################################


rm(list=ls())

library(httr)
library(rvest)
library(xml2)

# Base URL for NDVI GeoTIFFs
base_url <- "https://neo.gsfc.nasa.gov/archive/geotiff/MOD_NDVI_M/"

# Define the years and months we want
years <- 2001:2023
months <- sprintf("%02d", 1:12)  # Create month values with leading zeros (01, 02, ..., 12)

# Create a directory to store downloaded files
dir.create("NDVI_GeoTIFFs", showWarnings = FALSE)

# Function to download a file
download_ndvi_geotiff <- function(year, month) {
  # Create the file URL
  file_name <- paste0("MOD_NDVI_M_", year, "-", month, ".TIFF")
  file_url <- paste0(base_url, file_name)
  
  # Path to save the file locally
  dest_file <- file.path("NDVI_GeoTIFFs", file_name)
  
  # Check if file already exists
  if (!file.exists(dest_file)) {
    cat("Downloading:", file_url, "\n")
    tryCatch({
      # Download the file
      GET(file_url, write_disk(dest_file, overwrite = TRUE))
    }, error = function(e) {
      message(paste("Failed to download:", file_name, " - ", e$message))
    })
  } else {
    cat("File already exists:", dest_file, "\n")
  }
}


# Loop over each year and month to download the GeoTIFFs
for (year in years) {
  for (month in months) {
    download_ndvi_geotiff(year, month)
  }
}

cat("Download process completed!\n")
