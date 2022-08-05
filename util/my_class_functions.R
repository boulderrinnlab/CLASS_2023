
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



#' PROFILE_TSS 
#' This function takes a promoter window of your choosing
#' Then it will take a peak file and find the overlaps of 
#' all peaks in that window. It will provide a matrix of 
#' 1 (overlapped window) or (0 for no overlap)
#' each peak creates a single vector of the window length
#' these are concatenated, summed and density provided at each pos.
#' This function also takes into consideration of promoter orientation.
#' 
#' @param consensus_peaks can really be any peak file with start-stop - needs to be GRanges
#' @param lncrna_mrna_promoters any promoter GRanges object
#' @param upstream the length upstream of the TSS
#' @param downstream the length downstream of TSS

profile_tss <- function(consensus_peaks, 
                        lncrna_mrna_promoters,
                        upstream = 1e3,
                        downstream = 1e3) {
  
   peak_coverage <- coverage(consensus_peaks)
  coverage_length <- elementNROWS(peak_coverage)
  coverage_gr <- GRanges(seqnames = names(coverage_length),
                         IRanges(start = rep(1, length(coverage_length)), 
                                 end = coverage_length))
  
  promoters_gr <- subsetByOverlaps(lncrna_mrna_promoters, 
                                   coverage_gr, 
                                   type="within", 
                                   ignore.strand=TRUE)
  
  chromosomes <- intersect(names(peak_coverage), 
                           unique(as.character(seqnames(promoters_gr))))
  
  peak_coverage <- peak_coverage[chromosomes]
  promoters_ir <- as(promoters_gr, "IntegerRangesList")[chromosomes]
  promoter_peak_view <- Views(peak_coverage, promoters_ir)
  promoter_peak_view <- lapply(promoter_peak_view, function(x) t(viewApply(x, as.vector)))
  promoter_peak_matrix <- do.call("rbind", promoter_peak_view)
  minus_idx <- which(as.character(strand(promoters_gr)) == "-")
  

  promoter_peak_matrix[minus_idx,] <- promoter_peak_matrix[minus_idx, ncol(promoter_peak_matrix):1]
  promoter_peak_matrix <- promoter_peak_matrix[rowSums(promoter_peak_matrix) > 1,]
  peak_sums <- colSums(promoter_peak_matrix)
  
  peak_dens <- peak_sums/sum(peak_sums)
  metaplot_df <- data.frame(x = -upstream:(downstream-1),
                            dens = peak_dens)
  
  return(metaplot_df)
}



