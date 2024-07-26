

path_to_ml_labels_file <- 'C:/Users/egrout/Dropbox/calls/Galaxy_labels/machine_labels_all/'
output_folder <- 'C:/Users/egrout/Dropbox/calls/Galaxy_labels/machine_labels_all_cleaned/'

# List all CSV files in the input folder
csv_files <- list.files(path_to_ml_labels_file, pattern = "\\.csv$", full.names = TRUE)

# Loop through each CSV file
for (file in csv_files) {
  # Read the CSV file
  data <- read.csv(file, stringsAsFactors = FALSE, sep = "\t")
  
  # Replace all instances of "cue" with "Cue"
  data <- as.data.frame(lapply(data, function(x) {
    if (is.character(x)) {
      str_replace_all(x, "\\bcue\\b", "Cue")
    } else {
      x
    }
  }))
  
  # Get the base name of the file
  file_name <- basename(file)
  
  # Define the path for the output file
  output_file <- file.path(output_folder, file_name)
  
  # Save the updated data frame to the output folder
  write.csv(data, output_file, row.names = FALSE, quote = FALSE, sep = " ")
}























