
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
    paste(unlist(strsplit(y, "_"))[c(1,2)], collapse = "_") })
  
  peak_list <- c()
  for(i in 1:length(peak_files)) {
    peaks <- rtracklayer::import(peak_files[i])
    peak_list <- c(peak_list, peaks)
    names(peak_list)[length(peak_list)] <- dbp_name[i]
  }
  return(peak_list)
}


#' CREATE CONSENSUS PEAKS
#' this function will take multiple replicate .broadPeak files (also narrow)
#' find peaks that overlap in all the replicates. 
#' @description 
#' input set of chipseq replicate peak files
#' this function then creates one merged file peaks in all samples
#' @param dbp
#' This will be extracted with names(GR_list) in the lapply at end of fun
#' You will need a "dbps" or some object for the lapply that has the 
#' name of each dbp in the named GRanges list
#' 
#' @param peak_list
#' Named list of GRanges for each chipseq replicate
#' peak_list can be generated using import_peaks function above

consensus_from_reduced <- function(dbp, peak_list) {
dbp_peaks <- peak_list[grepl(as.character(dbp), names(peak_list))]
suppressWarnings(all_peaks <- GenomicRanges::reduce(unlist(as(dbp_peaks, "GRangesList"))))
peak_exists <- matrix(NA, nrow = length(all_peaks), ncol = length(dbp_peaks))
for(i in 1:length(dbp_peaks)) {
  suppressWarnings(peak_exists[,i] <- as.numeric(countOverlaps(all_peaks, dbp_peaks[[i]]) > 0))
}

# filter to consensus requiring peaks to be in all replicates
dbp_consensus <- all_peaks[rowSums(peak_exists) == ncol(peak_exists)]
# Required only two replicates == dbp_consensus <- all_peaks[rowSums(peak_exists) > 1]
return(dbp_consensus)
}




