
#' fun(x * y)
#' A function to multiply two numbers
#'
#' @description 
#' This function will multiply the input values of X and Y
#' 
#' @param x one number you'd like to multiply
#' y the other number you'd like to multiply
fun <- function(x, y) {
  ans <- x * y
  return(ans)
}


#' IMPORT_PEAKS Function
#' set's file path to peak files and extracts DBP name
#' @description 
#' this is the first step in loading in Chip_peak files
#' @param 
#' a file path is needed to the peak files "peak_path"

import_peaks <- function(consensus_file_path = broadpeakfilepath) {
  peak_files <- list.files(consensus_file_path, full.names = T, pattern = ".broadPeak")
  
  dbp_name <- sapply(peak_files, function(x){
    y <-  str_extract(x, "([^\\/]+$)")
    unlist(strsplit(y, "_"))[[1]] })
  
  peak_list <- c()
  for(i in 1:length(peak_files)) {
    peaks <- rtracklayer::import(peak_files[i])
    peak_list <- c(peak_list, peaks)
    names(peak_list)[length(peak_list)] <- dbp_name[i]
  }
  return(peak_list)
}





