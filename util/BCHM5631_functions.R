
#' # setting up the function: all what we just wrote above:
# Get's consensus peak file path, select .broadPeak and then get dbp name from file name
import_peaks <- function(consensus_file_path = broadpeakfilepath) {
  peak_files <- list.files(consensus_file_path, full.names = T, pattern = ".broadPeak")
  dbp_name <- sapply(peak_files, function(x){
    y <-  str_extract(x, "([^\\/]+$)")
    # NOTE ALTERNATIVE Gsub
    # gsub("_peaks.broadPeak", "", y)
    paste(unlist(strsplit(y, "_"))[c(1,2)], collapse = "_") 
    
  })
  
  # NOW FOR LOOP TO IMPORT FILES
  # setting an empty list "peak_list" to be filled by for loop
  peak_list <- c()
  
  # the for loop !
  for(i in 1:length(peak_files)) {
    # Import peaks
    peaks <- rtracklayer::import(peak_files[i])
    # Append this GRanges object to the of the list.
    peak_list <- c(peak_list, peaks)
    # Name the list elements by their TF name.
    names(peak_list)[length(peak_list)] <- dbp_name[i]
  }
  return(peak_list)
}