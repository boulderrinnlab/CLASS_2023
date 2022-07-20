

get_dbp_name <- function(x) {
  file_name <- str_extract(x, "[\\w-]+\\.broadPeak")
  return(str_extract(file_name, "^[^_]+(?=_)"))
}

#Functions we need

# X features_te_total
# Xfeatures overlapping promoters
# >> seperate R file compare gene expression value
#  X features_promoters
# X peak_feature_matrix


#' import peak .bed files as a list
#' 
#' @description 
#' this function will take consensus peak files and name them by the DNA binding protein and return a list
#' 
#' @param consensus_file_path the path to consensus peak files


import_peaks <- function(consensus_file_path = "scratch/Shares/rinnclass/data/peaks") {
  peak_files <- list.files(consensus_file_path, full.names = T)
  file_names <- str_extract(peak_files, "[\\w-]+\\.bed")
  tf_name <- str_extract(file_names, "^[^_]+(?=_)")
  
  
  peak_list <- c()
  for(i in 1:length(peak_files)) {
    # Import peaks
    peaks <- rtracklayer::import(peak_files[i])
    # Append this GRanges object to the of the list.
    peak_list <- c(peak_list, peaks)
    # Name the list elements by their TF name.
    names(peak_list)[length(peak_list)] <- tf_name[i]
  }
  return(peak_list)
}



#' intersect replicates into a "consensus peak list" 
#' 
#' @description 
#' this function will take the intersect and union of peak widths across replicates for a given DNA binding protein. the function that will take a list of granges objects and return 
# one granges object with merged peaks that are in all replicates
#' 
#' @param 
#'  the path to consensus peak files
#' # We're going to iterate over all the files to make it work. 

create_consensus_peaks <- function(broadpeakfilepath = "data/test_work/all_peak_files") {
  
  
  fl <- list.files(broadpeakfilepath, 
                   full.names=TRUE)
  fl <- fl[grep("peaks.broadPeak", fl)]
  
  tf_name <- sapply(fl, function(x){
    y <-  str_extract(x, "([^\\/]+$)")
    unlist(strsplit(y, "_"))[[1]]
  })
  
  # We don't want to run consensus peak creation
  # on those files that don't have replicates
  tf_df <- data.frame(table(tf_name)) %>%
    # filter those with no replicates
    filter(Freq > 1)
  unique_tf <- as.character(tf_df$tf_name)
  
  consensus_peaks <- list()
  # This for loop will iterate over all dna binding proteins.
  for(i in 1:length(unique_tf)) {
    
    # load all the peak files corresponding to this dna binding proteins.
    tf <- unique_tf[i]
    print(tf)
    tf_index <- grep(tf, tf_name)
    tf_files <- fl[tf_index]
    
    peak_list <- c()
    for(j in 1:length(tf_files)) {
      # See the read peaks function to know what subfunctions are called.
      peak_list <- c(peak_list, read_peaks(tf_files[j]))
    }
    
    
    final_peakset <- intersect_peaks(peak_list = peak_list)
    if(length(final_peakset) > 0) {
      final_peakset$name <- paste0(tf, "_", 1:length(final_peakset))
    }
    
    consensus_peaks <- c(consensus_peaks, list(final_peakset))
    names(consensus_peaks)[length(consensus_peaks)] <- tf
  }
  return(consensus_peaks)
}

# TODO: refactor
read_peaks <- function(broad_peak_file, filter_to_canonical_chr = TRUE) {
  # A broad peak file is just a tab separated file 
  dat <- read.table(broad_peak_file, sep = "\t")
  if(filter_to_canonical_chr == TRUE) {
    dat <- dat[dat$V1 %in% c(paste0("chr", 1:22), "chrM", "chrX", "chrY"),]
  }
  gr <- GRanges(seqnames = dat$V1,
                ranges = IRanges(start=dat$V2,end=dat$V3))
  return(gr)
}
# This is the function that will be doing the core of the work here. 
# When two peaks intercept, we will take their outer boundaries to be the new
# peak -- using the reduce function.
intersect_peaks <- function(peak_list) {
  combined_peaks <- peak_list[[1]]
  for(i in 2:length(peak_list)) {
    suppressWarnings(pl_ov <- findOverlaps(combined_peaks, peak_list[[i]]))
    pl1 <- combined_peaks[unique(pl_ov@from)]
    pl2 <- peak_list[[i]][unique(pl_ov@to)]
    suppressWarnings(combined_peaks <- GenomicRanges::reduce(union(pl1, pl2)))
  }
  return(combined_peaks)
}



#' function finds overlaps between consensus_peaks and genomic features 
#' 
#' @description get_overlapping_peaks
#' this function will intersect the consesus_peaks with gene_body, promoters, mRNA_promoters, lncRNA_promoters, te_family
#' 
#' @param features
#'  set of genomic features as a GRanges object
#'  
#' @param peak_list
#' #list of peaks of dna binding proteins that will be intersected

get_overlapping_peaks <- function(features, peak_list){
  
  overlaps_list <- c()
  for(j in 1:length(peak_list)) {
    suppressWarnings(ov <- findOverlaps(peak_list[[j]], features))
    overlapping_peaks <- peak_list[[j]][unique(ov@from)]
    overlaps_list <- c(overlaps_list, overlapping_peaks)
    names(overlaps_list)[length(overlaps_list)] <- names(peak_list)[j]
  }
  return(overlaps_list)
} 


# subset rmsk

import_repeatmasker <- function(rmsk_file = "/Shares/rinn_class/data/k562_chip/util/rmsk.txt") {
  
  rmsk <- read.table(file = rmsk_file)
  colnames(rmsk) <- c("bin", "swScore", "milliDiv", "milliDel", "milliIns",
                      "genoName", "genoStart", "genoEnd", "genoLeft", "strand",
                      "repName", "repClass", "repFamily", "repStart",	"repEnd",
                      "repLeft",	"id")
  
  rmsk_gr <- GRanges(seqnames = rmsk$genoName,
                     IRanges(rmsk$genoStart, rmsk$genoEnd),
                     strand = rmsk$strand)
  
  # Add metadata colums
  rmsk_gr$rep_class <- rmsk$repClass
  rmsk_gr$rep_family <- rmsk$repFamily
  rmsk_gr$rep_name <- rmsk$repName
  
  return(rmsk_gr)
}

subset_rmsk <- function(rmsk_gr, rep_level = "family") {
  # rep_level needs to be one of "class", "family", or "name"
  if(!(rep_level %in% c("class", "family", "name"))) {
    stop("Repeat level needs to be either: class, family, or name.")
  } 
  
  level_column <- grep(rep_level, names(rmsk_gr@elementMetadata))
  rep_levels <- unique(rmsk_gr@elementMetadata[,level_column])
  
  rmsk_list <- c()
  for(i in 1:length(rep_levels)) {
    rmsk_list <- c(rmsk_list, list(rmsk_gr[rmsk_gr@elementMetadata[,level_column] == rep_levels[i]]))
    names(rmsk_list)[length(rmsk_list)] <- rep_levels[[i]]
  }
  return(rmsk_list)
}


#' function to subset features for promomters. 
#' 
#' @description feature_subset
#' Take a gencode gtf to subset the biotype of promoters we want as a set of GRanges
#' 
#' @param gencode_gr
#'  set of genomic features as a GRanges object
#'  
#' @param biotype
#' this takes "lncRNA" or "protein-coding" as input for promoter type
#'
#' @param upstream
#'To add upstream sequence to feature subset
#'
#' @param downstream
#'To add downstream sequence to feature subset

get_promoter_regions <- function(gencode_gr, biotype, upstream = 3e3, downstream = 3e3) {
  
  genes <- gencode_gr[gencode_gr$type == "gene"]
  genes <- genes[genes$gene_type %in% biotype]
  
  proms <- GenomicRanges::promoters(genes, upstream = upstream, downstream = downstream)
  
  return(proms)
  
}


#' function to subset features for overlapping promomters. 
#' 
#' @description 
#' Take a gencode gtf to subset the biotype of overlapping promoters we want as a set of GRanges
#' 
#' @param gencode_gr
#'  set of genomic features as a GRanges object
overlapping_promoters <- function(gencode_gr, upstream = 1000, downstream = 0) {
  
  genes <- gencode_gr[gencode_gr$type == "gene"]
  proms <- GenomicRanges::promoters(genes, upstream = upstream, downstream = downstream)
  
  reduced_proms <- GenomicRanges::reduce(proms, ignore.strand = T)
  
  reduced_widths <- data.frame(width = width(reduced_proms),
                               index = 1:length(reduced_proms))
  
  ov_proms <- GenomicRanges::findOverlaps(proms, reduced_proms, ignore.strand=TRUE)
  ov_proms <- data.frame("gene_promoters_from" = ov_proms@from,
                         "reduced_promoters_to" = ov_proms@to)
  
  ov_proms_summary <- ov_proms %>% 
    group_by(reduced_promoters_to) %>%
    summarize(count = n())
  ov_proms <- merge(ov_proms, ov_proms_summary)
  overlapped <- ov_proms %>% filter(count > 1)
  
  reduced_promoters_overlapping <- reduced_proms[unique(overlapped$reduced_promoters_to)]
  
  overlapped$gene_id <- proms$gene_id[overlapped$gene_promoters_from]
  overlapped$gene_name <- proms$gene_name[overlapped$gene_promoters_from]
  overlapped$gene_type <- proms$gene_type[overlapped$gene_promoters_from]
  overlapped$strand <- as.character(strand(proms[overlapped$gene_promoters_from]))
  
  overlapped_metadata <- overlapped %>% 
    arrange(strand) %>%
    group_by(reduced_promoters_to, count) %>%
    summarize(nearby_promoter_gene_ids = paste(gene_id, collapse = ";"),
              nearby_promoter_gene_names = paste(gene_name, collapse = ";"),
              nearby_promoter_gene_type = paste(gene_type, collapse = ";"),
              nearby_promoter_gene_strands = paste(strand, collapse = ";"))
  
  # Okay, I think I want to merge this all back into overlapped. 
  overlapped_merge <- merge(overlapped, overlapped_metadata)
  names(overlapped_merge)[2] <- "num_nearby_promoters"
  
  # Let's now filter to just the lncRNA and mRNA genes since that's
  # the subject of our analysis.
  overlapped_merge <- overlapped_merge %>% filter(gene_type %in% c("lncRNA", "protein_coding"))
  
  # Clean up to just the columns we need to add to peak_occurrence_df
  overlapped_merge <- overlapped_merge %>% dplyr::select("gene_id", "strand",
                                                         "num_nearby_promoters",
                                                         "nearby_promoter_gene_ids",
                                                         "nearby_promoter_gene_names",
                                                         "nearby_promoter_gene_type",
                                                         "nearby_promoter_gene_strands")
  
  # Finally let's annotate the bidirectionals
  # So these will be genes that have only two two nearby promoters with
  # opposite direction of transcription.
  # Since we ordered by strand, we only need to account for "-;+" 
  # (and not "+;-")
  overlapped_merge[overlapped_merge$nearby_promoter_gene_strands == "-;+", 
                   "shared_promoter_type"] <- "bidirectional"
  overlapped_merge[overlapped_merge$nearby_promoter_gene_strands == "-;-" |
                     overlapped_merge$nearby_promoter_gene_strands == "+;+"  , 
                   "shared_promoter_type"] <- "nearby_same_strand"
  overlapped_merge[overlapped_merge$num_nearby_promoters > 2, 
                   "shared_promoter_type"] <- "multiple_nearby_promoters"
  
  return(overlapped_merge) 
}


#' function to summarize the number of events in features on each individual promoter. 
#' 
#' @description 
#' Take a gencode gtf to subset the biotype of promoters we want as a set of GRanges
#' 
#' @param features
#' set of genomic features as a GRanges object
#'  
#' @param peak_list
#' #list of peaks of dna binding proteins that will be intersected
#' 
#' @param type
#' Return either a matrix of counts over features or a binary occurrence matrix

count_peaks_per_feature <- function(features, peak_list, type = "counts") {
  
  if(!(type %in% c("counts", "occurrence"))) {
    stop("Type must be either occurrence or counts.")
  }
  
  peak_count <- matrix(numeric(), ncol = length(features), nrow = 0)
  
  for(j in 1:length(peak_list)) {
    suppressWarnings(ov <- countOverlaps(features, peak_list[[j]]))
    peak_count <- rbind(peak_count, ov)
    rownames(peak_count)[nrow(peak_count)] <- names(peak_list)[j]
    colnames(peak_count) <- features$gene_id
  }
  
  peak_matrix <- peak_count
  
  if(type == "occurrence") {
    peak_occurrence <- matrix(as.numeric(peak_count > 0), 
                              nrow = dim(peak_count)[1],
                              ncol = dim(peak_count)[2])
    rownames(peak_occurrence) <- rownames(peak_count)
    colnames(peak_occurrence) <- colnames(peak_count)
    peak_matrix <- peak_occurrence
  }
  
  return(peak_matrix)
  
}

##### These functions are adapted from the ChIPSeeker package.
##### We could not use the ChIPSeeker packages directly because
##### They were hardcoded to use all the available cores -- which
##### is not feasible in a shared compute environment.

get_tag_count <- function(tag_matrix, xlim, conf, ncpus = 6) {
  ss <- colSums(tag_matrix)
  ss <- ss/sum(ss)
  pos <- value <- NULL
  dd <- data.frame(pos=c(xlim[1]:(xlim[2]-1)), value=ss)
  if (!(missingArg(conf) || is.na(conf))){
    tagCiMx <- get_tag_ci_matrix(tag_matrix, conf = conf, ncpus = ncpus)
    dd$Lower <- tagCiMx["Lower", ]
    dd$Upper <- tagCiMx["Upper", ]
  }
  return(dd)
}

get_tag_ci_matrix <- function(tag_matrix, conf = 0.95, resample=500, ncpus = 6){
  RESAMPLE_TIME <- resample
  trackLen <- ncol(tag_matrix)
  if (Sys.info()[1] == "Windows") {
    tagMxBoot <- boot(data = tag_matrix, statistic = get_sgn, R = RESAMPLE_TIME)
  } else {
    tagMxBoot <- boot(data = tag_matrix, statistic = get_sgn, R = RESAMPLE_TIME,
                      parallel = "multicore", ncpus = ncpus)
  }
  cat(">> Running bootstrapping for tag matrix...\t\t",
      format(Sys.time(), "%Y-%m-%d %X"), "\n")
  tagMxBootCi <- sapply(seq_len(trackLen), function(i) {
    bootCiToken <- boot.ci(tagMxBoot, type = "perc", index = i)
    ## parse boot.ci results
    return(parse_boot_ci_perc(bootCiToken))
  }
  )
  row.names(tagMxBootCi) <- c("Lower", "Upper")
  return(tagMxBootCi)
}

get_sgn <- function(data, idx){
  d <- data[idx, ]
  ss <- colSums(d)
  ss <- ss / sum(ss)
  return(ss)
}

parse_boot_ci_perc <- function(bootCiPerc){
  bootCiPerc <- bootCiPerc$percent
  tmp <- length(bootCiPerc)
  ciLo <- bootCiPerc[tmp - 1]
  ciUp <- bootCiPerc[tmp]
  return(c(ciLo, ciUp))
}

plot_profile <- function(tag_count, tf_name) {
  p <- ggplot(tag_count, aes(pos))
  p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper),
                       linetype = 0, alpha = 0.2)
  p <- p + geom_line(aes(y = value))
  
  xlim <-c(-3000, 3000)
  origin_label = "TSS"
  #### Testing below       
  p <- p + geom_vline(xintercept=0,
                      linetype="longdash")
  p <- p + scale_x_continuous(breaks=c(xlim[1], floor(xlim[1]/2),
                                       0,
                                       floor(xlim[2]/2), xlim[2]),
                              labels=c(xlim[1], floor(xlim[1]/2),
                                       origin_label, 
                                       floor(xlim[2]/2), xlim[2]))
  p <- p+xlab("Genomic Region (5'->3')")+ylab("Peak Count Frequency")
  p <- p + theme_bw() + theme(legend.title=element_blank())
  p <- p + theme(legend.position="none")
  p <- p + ggtitle(tf_name)
  return(p)
}

##' calculate the tag matrix
##'
##'
##' @title getTagMatrix
##' @param peak peak file or GRanges object
##' @param weightCol column name of weight, default is NULL
##' @param windows a collection of region with equal size, eg. promoter region.
##' @param flip_minor_strand whether flip the orientation of minor strand
##' @return tagMatrix
##' @export
##' @import BiocGenerics S4Vectors IRanges GenomeInfoDb GenomicRanges
get_tag_matrix <- function(peak.gr, weightCol=NULL, windows, flip_minor_strand=TRUE) {
  
  if (! is(windows, "GRanges")) {
    stop("windows should be a GRanges object...")
  }
  if (length(unique(width(windows))) != 1) {
    stop("width of windows should be equal...")
  }
  
  if (is.null(weightCol)) {
    peak.cov <- coverage(peak.gr)
  } else {
    weight <- mcols(peak.gr)[[weightCol]]
    peak.cov <- coverage(peak.gr, weight=weight)
  }
  cov.len <- elementNROWS(peak.cov)
  cov.width <- GRanges(seqnames=names(cov.len),
                       IRanges(start=rep(1, length(cov.len)),
                               end=cov.len))
  suppressWarnings(windows <- subsetByOverlaps(windows, cov.width,
                                               type="within", ignore.strand=TRUE))
  
  chr.idx <- intersect(names(peak.cov),
                       unique(as.character(seqnames(windows))))
  
  peakView <- Views(peak.cov[chr.idx], as(windows, "IntegerRangesList")[chr.idx])
  tagMatrixList <- lapply(peakView, function(x) t(viewApply(x, as.vector)))
  tagMatrix <- do.call("rbind", tagMatrixList)
  
  ## get the index of windows, that are reorganized by as(windows, "IntegerRangesList")
  idx.list <- split(1:length(windows),  as.factor(seqnames(windows)))
  idx <- do.call("c", idx.list)
  
  rownames(tagMatrix) <- idx
  tagMatrix <- tagMatrix[order(idx),]
  
  ## minus strand
  if (flip_minor_strand) {
    ## should set to FALSE if upstream is not equal to downstream
    ## can set to TRUE if e.g. 3k-TSS-3k
    ## should set to FALSE if e.g. 3k-TSS-100
    minus.idx <- which(as.character(strand(windows)) == "-")
    tagMatrix[minus.idx,] <- tagMatrix[minus.idx, ncol(tagMatrix):1]
  }
  
  tagMatrix <- tagMatrix[rowSums(tagMatrix)!=0,]
  return(tagMatrix)
}


make_promoter_binding_matrix <- function(peak_list, promoter) {
  
  # Filter the peaks to only those overlapping the promoter
  promoter_peaks <- lapply(peak_list, function(x) subsetByOverlaps(x, promoter))
  
  promoter_peaks <- promoter_peaks[sapply(promoter_peaks, length) > 0]
  
  if(length(promoter_peaks) == 0) {
    return(matrix(0, ncol = 6000,nrow = 1))
  }
  
  # Set the seqlevels to only the chromosome that the promoter is on
  for(i in 1:length(promoter_peaks)) {
    seqlevels(promoter_peaks[[i]]) <- as.character(seqnames(promoter))
  }
  
  # Get the Rle coverage values for the promoter
  promoter_coverage <- lapply(promoter_peaks, coverage)
  
  # Change promoter to a GRangesList
  promoter <- GRangesList(promoter)
  
  promoter_peak_view <- lapply(promoter_coverage, extract_peak_view, promoter)
  
  promoter_peak_matrix <- do.call(rbind, promoter_peak_view)
  
  
  if(as.character(strand(promoter[[1]])) == "-") {
    # Then flip the matrix so that downstream is always to the right
    promoter_peak_matrix <- promoter_peak_matrix[,ncol(promoter_peak_matrix):1]
  }
  return(promoter_peak_matrix) 
}


extract_peak_view <- function(peaks, promoter) {
  peak_view <- Views(peaks, promoter)
  peak_view <- lapply(peak_view, function(x) t(viewApply(x, as.vector)))
  
  peak_vector <- as.vector(peak_view[[as.character(seqnames(promoter))]])
  
  # Extend the vector to the promoter window length
  # padding with zeros
  peak_vector <- c(peak_vector, rep(0,width(promoter[[1]]) - length(peak_vector)))
  return(peak_vector)
}


profile_tss <- function(peaks, 
                        promoters_gr,
                        upstream = 3e3,
                        downstream = 3e3) {
  
  
  peak_coverage <- coverage(peaks)
  
  coverage_length <- elementNROWS(peak_coverage)
  coverage_gr <- GRanges(seqnames = names(coverage_length),
                         IRanges(start = rep(1, length(coverage_length)), 
                                 end = coverage_length))
  
  promoters_gr <- subsetByOverlaps(promoters_gr, 
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
  promoter_peak_matrix[minus_idx,] <- promoter_peak_matrix[minus_idx,
                                                           ncol(promoter_peak_matrix):1]
  
  promoter_peak_matrix <- promoter_peak_matrix[rowSums(promoter_peak_matrix) > 1,]
  
  peak_sums <- colSums(promoter_peak_matrix)
  peak_dens <- peak_sums/sum(peak_sums)
  
  metaplot_df <- data.frame(x = -upstream:(downstream-1),
                            dens = peak_dens)
  
  return(metaplot_df)
}

