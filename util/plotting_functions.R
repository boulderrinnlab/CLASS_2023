#' Plot a UCSC style browser track with bigWigs and MACS peak calls
#' 
#' @description 
#' This function will plot browser tracks of the specified DNA binding protein's
#' ChIP signal along with the MACs peak calls and a gene region track from the
#' Gencode annotation. 
#'
#' @param dbps A character vector of the DNA binding proteins to plot
#' @param region A GRanges object specifying the genomic region to plot 
#' @param bw_axis_limits A list axis limits for the bigWig tracks which must
#' be the same length as the dbps vector. This is so that the limits can be set for
#' each DBP separately. ex. list(c(0,0.6), c(0,1.2))
#' @param gencode_gm_path A path to the Gencode gene model csv. This dataframe 
#' contains only the exons from the Gencode annotation. See Gviz vignette for
#' documentation. 
#' @param plot_title Main plot title
plot_replicate_peaks <- function(dbps, region, bw_axis_limits, plot_title,
                                 save_pdf = T, pdf_path, png_path,
                                 height_per_dbp = 5,
                                 fixed_height = 3,
                                 width = 12,
                                 display_plot = F,
                                 gencode_gm_path = "../util/gencode_genemodel_for_Gviz.csv") {
  
  if(length(dbps) != length(bw_axis_limits)) stop("bw_axis_limits must be a list of
                                                  axis limits the same length as
                                                  input dbps.")
  
  if(length(region) > 1) stop("Please plot only one region at a time")
  
  #### Gene region track
  
  # Retrieve gencode gene model
  gencode_gm <- read.csv(gencode_gm_path)
  # Subset to the region of interest
  # TODO: consider what to do if a gene is cut off on one end..
  reg_gm <- gencode_gm %>% filter(chromosome == as.character(seqnames(region)),
                                  start >= start(region),
                                  end <= end(region))
  
  grtrack <- GeneRegionTrack(reg_gm, genome = "hg38",
                             chromosome = as.character(seqnames(region)),
                             name = "Gencode",
                             showId = TRUE,
                             fill = "gray70")
  
  names(dbps) <- dbps
  
  # Peak tracks
  ptracks <- lapply(dbps, create_peak_tracks, region = region)
  
  # Bigwig tracks
  bwtracks <- list()
  for(i in 1:length(dbps)) {
    bwtracks <- c(bwtracks, list(create_bigwig_tracks(dbps[[i]], region = region,
                                                      axis_limits = bw_axis_limits[[i]])))
  }
  names(bwtracks) <- names(dbps)
  
  # Consensus peak tracks
  cptracks <- lapply(dbps, create_consensus_track, region = region)
  
  # Info tracks
  itrack <- IdeogramTrack(genome = "hg38", chromosome = as.character(seqnames(region)))
  gtrack <- GenomeAxisTrack()
  
  # Construct plotting call
  pcall <- "region_plot <- plotTracks(list(itrack, gtrack,"
  
  for(dbp in dbps) {
    for(replicate in 1:length(ptracks[[dbp]])) {
      pcall <- c(pcall, paste0("ptracks[['",dbp,"']][[",replicate,"]],",
                               "bwtracks[['",dbp,"']][[",replicate,"]],"))
    }
    pcall <- c(pcall, paste0("cptracks[['",dbp,"']],"))
  }
  
  pcall <- c(pcall, paste0("grtrack),
                           background.title = 'white',
                           fontcolor = 'black',
                           rotation.title = 0,
                           cex.title = 0.7,
                           col.title = 'black',
                           col.axis = 'black', 
                           main = '", plot_title, "',",
                           " cex.main = 0.8)"))
  
  full_plot_call <- paste(pcall, collapse = " ")
  
  if(save_pdf == F) {
    # This can then be plotted with
    # eval(parse(text = full_plot_call))
    return(full_plot_call)
  } else {
    if(display_plot == T) {
      eval(parse(text = full_plot_call))
    }
    # The pdf height should be about +5 in per dbp
    total_height <- fixed_height + (height_per_dbp * length(dbps))
    # PDF
    pdf(pdf_path, height = total_height, width = width)
    eval(parse(text = full_plot_call))
    dev.off()
    # PNG
    png(png_path, height = total_height, width = width, units = "in",
        res = 300)
    eval(parse(text = full_plot_call))
    dev.off()
  }
}

create_peak_tracks <- function(dbp, region,
                               peak_files_path = "/Shares/rinn_class/data/k562_chip/results/bwa/mergedLibrary/macs/broadPeak/") {
  # Read in peaks 
  peak_files <- list.files(peak_files_path, full.names = T)
  peak_files <- peak_files[grep(dbp, peak_files)]
  peak_files <- peak_files[grep("peaks.broadPeak", peak_files)]
  names(peak_files) <- sapply(peak_files, function(x) {
    unlist(strsplit(unlist(strsplit(x, "//"))[[2]], "_peaks"))[[1]]
  })
  
  peaks <- lapply(peak_files, read_peaks, filter_to_canonical_chr = T)
  peaks <- lapply(peaks, subsetByOverlaps, ranges = region)
  
  peak_tracks <- list()
  for(i in 1:length(peaks)) {
    peak_tracks <- c(peak_tracks, AnnotationTrack(peaks[[i]], 
                                                  col = "#424242", 
                                                  fill = "#424242", 
                                                  name = names(peaks)[[i]]))
    names(peak_tracks)[length(peak_tracks)] <- names(peaks)[[i]]
  }
  return(peak_tracks)
}


create_consensus_track <- function(dbp, region,
                                   consensus_peak_files_path = "../01_consensus_peaks/results/consensus_peaks/filtered_by_peaks/") {
  # Read in peaks 
  peak_files <- list.files(consensus_peak_files_path, full.names = T)
  peak_files <- peak_files[grep(dbp, peak_files)]
  names(peak_files) <- sapply(peak_files, function(x) {
    unlist(strsplit(unlist(strsplit(x, "//"))[[2]], "_consensus"))[[1]]
  })
  
  peaks <- rtracklayer::import(peak_files)
  peaks <- subsetByOverlaps(peaks, region)
  consensus_track <- AnnotationTrack(peaks, 
                                     col = "#424242", 
                                     fill = "#424242", 
                                     name = names(peak_files))
  return(consensus_track)
}


create_bigwig_tracks <- function(dbp, region, axis_limits, bigwig_files_path = "/Shares/rinn_class/data/k562_chip/results/bwa/mergedLibrary/bigwig/") {
  
  # Read in bigwigs 
  bws <- list.files(bigwig_files_path, full.names = T)
  bws <- bws[grep(dbp, bws)]
  names(bws) <- sapply(bws, function(x) {
    unlist(strsplit(unlist(strsplit(x, "//"))[[2]], ".mLb"))[[1]]
  })
  
  bw_tracks <- list()
  for(i in 1:length(bws)) {
    bw_tracks <- c(bw_tracks, DataTrack(range = bws[[i]], 
                                        genome = "hg38", 
                                        type = "h", 
                                        chromosome = as.character(seqnames(region)), 
                                        ylim = axis_limits,
                                        size = 1, 
                                        col = "#212121", 
                                        name = " "))
    names(bw_tracks)[length(bw_tracks)] <- names(bws)[[i]]
  }
  return(bw_tracks)
}






#' Plot genomic ranges from the IRanges object or GRanges object.
#'
#' @param x a IRanges object.
#' @param xlim like in plot.
#' @param main A title of the plot.
#' @param col color. Default "black". e.g. "blue".
#' @param sep Default 0.5.
#' @param ... Other params you can pass to plot function.
#' @return A plot.
#' @export
#' @examples
#' library(IRanges)
#' plotRanges(IRanges(1:10, width=10:1, names=letters[1:10]))

plotRanges <- function(x, xlim = x, main = deparse(substitute(x)),  col = "black", sep = 0.5, ...){
  if (!requireNamespace("IRanges", quietly = TRUE)) {
    stop("IRanges needed for this function to work. Please install it.",
         call. = FALSE)
  }
  height <- 1
  if (is(xlim, "Ranges"))
    xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim, c(0, max(bins) * (height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x) - 0.5, ybottom, end(x) + 0.5, ybottom + height, col = col, ...)
  title(main)
  axis(1)
}



plot_promoter_peak_matrix <- function(promoter_peak_matrix, gene_name, save_pdf = FALSE, show = TRUE,
                                      save_dir = "figures/") {
  # Set the colors for the binary heatmap
  colors <- structure(c("#ffffff","#a8404c"), names = c("0", "1"))
  split <- data.frame(split = c(rep("-3kb",3000), rep("+3kb", 3000)))
  
  ht <- Heatmap(promoter_peak_matrix, cluster_columns = FALSE, col = colors,
                border = "black", show_heatmap_legend = FALSE,
                use_raster = TRUE,
                column_split = split,
                column_gap = unit(0, "mm"),
                row_names_gp = gpar(fontsize = 7))
  
  if(show == TRUE) {
    draw(ht, column_title = paste0(gene_name, " promoter peaks"))
  }
  
  
  if(save_pdf == TRUE) {
    pdf(file.path(save_dir, paste0(gene_name, "_promoter_peaks.pdf")), 
        height = ceiling(0.091*nrow(promoter_peak_matrix)*2)/2, width = 7)
    draw(ht, column_title = paste0(gene_name, " promoter peaks"))
    dev.off()
  }
}

plot_tss_profile <- function(tss_profile_matrix, dbp, save_pdf = FALSE,
                             save_dir = "figures/") {
  show(plot(tss_profile_matrix[dbp,], main = paste0(dbp, " promoter profile"),
            xlab = "", ylab = "", type = "l", lwd = 2))
  
  if(save_pdf == TRUE) {
    pdf(paste0(save_dir, dbp, "_tss_profile.pdf"))
    plot(tss_profile_matrix[dbp,], main = paste0(dbp, " promoter profile"),
         xlab = "", ylab = "", type = "l", lwd = 2)
    dev.off()
  }
  return()
}