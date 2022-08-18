---
  title: "01_BIG_KNIT"
author: "JR"
date: "8/17/2022"
output: html_document
editor_options: 
  chunk_output_type: console
---
  ````{r setup, include=FALSE} 
knitr::opts_chunk$set(warning = FALSE, message = FALSE) 
knitr::opts_chunk$set(echo = TRUE)
library(Gviz)
library(ggpubr)
library(ggplot2)
library(tidyverse)
library(GenomicRanges)
library(ComplexHeatmap)
source("../../../../util/my_class_functions.R")
source("../../../../util/_setup.R")

# filepath to import peaks
basepath <- "/scratch/Shares/rinnclass/CLASS_2023/JR"
peak_path <- "CLASS_2023/CLASSES/05_R_analyses/analysis/00_consensus_peaks"
consensusPeakPath <- file.path(basepath, peak_path)

```

# Goal:

Here we aim to download all available DNA binding protein (DBP) profiles in a single cell state (measured by ChIP-seq). This will allow us to investigate the binding properties of hundreds of DBPs in the same cellular context or background. We aim to address several questions: 
  
  (i) What are the number of peaks and genome coverage for each DBP? 
  
  (ii) What are the binding preferences for promoters, gene-bodies and intergenic genomic regions? 
  
  (iii) What are the similarities and differences across DBPs based on their genome-wide binding profiles genome-wide? 
  
  (iv) What properties or preferences do promoters have for binding events. 

(v) How does binding to a promoter affect the transcriptional output of that promoter?
  
  (vi) Are there reservoir promoters in HepG2 as defined in k562 previously? 
  
  
  
  To address these questions we have curated a set of over 1,000 ChIPseq data sets comprised of 48 DBPs in HEPG2 cells from the ENCODE consortium. We required duplicate ChIP-seq experiments for a given DBP and other criterion that can be found here :
  
  <https://www.encodeproject.org/report/?type=Experiment&status=released&assay_slims=DNA+binding&biosample_ontology.term_name=HepG2&assay_title=TF+ChIP-seq&biosample_ontology.classification=cell+line&files.read_length=100&files.read_length=76&files.read_length=75&files.read_length=36&assay_title=Control+ChIP-seq&assay_title=Histone+ChIP-seq&files.run_type=single-ended>
  
  ## These samples were selected on the following criteria:
  
  1)  "chromatin" interaction data, then DNA binding data, cell line HEPG2, "TF-Chip-seq".
2)  We further selected "TF Chip-seq", "Control chip-seq" and "Histone Chip-seq".
3)  We selected several read lengths to get the most DNA binding proteins (DBPs)
4)  Read lengths: 100, 76, 75, 36
5)  ONLY SINGLE END READS (this eliminates 54 samples)

### Experimental data was downloading by (ENCODE report.tsv):

<https://www.encodeproject.org/report.tsv?type=Experiment&status=released&assay_slims=DNA+binding&biosample_ontology.term_name=HepG2&assay_title=TF+ChIP-seq&biosample_ontology.classification=cell+line&files.read_length=100&files.read_length=76&files.read_length=75&files.read_length=36&assay_title=Control+ChIP-seq&assay_title=Histone+ChIP-seq&files.run_type=single-ended>
  
  ### The FASTQ files were downloaded with:
  
  "<https://www.encodeproject.org/metadata/?status=released&assay_slims=DNA+binding&biosample_ontology.term_name=HepG2&assay_title=TF+ChIP-seq&biosample_ontology.classification=cell+line&files.read_length=100&files.read_length=76&files.read_length=75&files.read_length=36&assay_title=Control+ChIP-seq&assay_title=Histone+ChIP-seq&files.run_type=single-ended&type=Experiment>"

MD5sums were checked with all passing (see encode_file_info function to reterive MD5Sum values that are not available from the encode portal (/util)
                                       
                                       ### Processing data:
                                       
                                       We processed all the read alignments and peak calling using the NF_CORE ChIP-seq pipeline: (nfcore/chipseq v1.2.2)
                                       
                                       ## Next we created consensus peaks that overlap in both replicates
                                       
                                       Our strategy was to take peaks in each replicate and find all overlapping peak windows. We then took the union length of the overlapping range in each peak window. We further required that a DBP has at least 1,000 peaks (the min quartile of peaks observed)
                                       
                                       Result: 
                                         (i) we have 430 DBPs that have at least 1,000 peaks in both replicates
                                       
                                       ## We will start with loading filtered consensus peaks (described above)
                                       ## All subequent analyses are contained within this document.
                                       ## To find results, search for "Result:"
                                       
                                       # Loading in all .broadPeak files for all dbps and their replicates
                                       ```{r load .broadPeak files}
                                       
                                       # using import peaks to import .broadPeak files (~10min)
                                       peak_list <- import_peaks(consensus_file_path = broadpeakfilepath)
                                       
                                       ```
                                       
                                       # creating consensus peaks
                                       Using consensus_from_reduced function that requires peak list
                                       and dbp object. This will be done in 2 steps
                                       (1) created dbp object of unique names for dbp
                                       (2) run consensus_from_reduced
                                       
                                       ```{r create consensus peaks}
                                       
                                       # Creating unique DBP object for create_consensus_peaks_from_reduced
                                       dbp <- unique(sapply(names(peak_list), function(x) {
                                         unlist(strsplit(x, "_"))[1]
                                       }))
                                       
                                       # now run our function consensus_from_reduced
                                       consensus_list <- lapply(dbp, consensus_from_reduced, peak_list)
                                       
                                       # adding names to the GRange list
                                       names(consensus_list) <- dbp
                                       
                                       # exploring the number of peaks in the consensus_list
                                       ```
                                       
                                       ```{r exploring number of peaks in consenus peaks}
                                       # creating list of num_peaks per dbp
                                       num_peaks <- sapply(consensus_list, length)
                                       
                                       # plotting
                                       hist(num_peaks, breaks = 1000)
                                       hist(num_peaks, breaks = 1000, xlim = c(0,3000))
                                       
                                       ```
                                       Result:
                                         (i)seems like 1000 peaks should be the min moving forward.
                                       
                                       # filtering consensus_list to dbps with > 1000 peaks
                                       ```{r filtered_consenus_list}
                                       
                                       # filtering to 1000 peaks
                                       filtered_consensus_list <- consensus_list[sapply(consensus_list, length) > 1000]
                                       
                                       # saving 
                                       save(filtered_consensus_list, file = "results/filtered_consensus_list.RData")
                                       
                                       # keeping track of DBPs lost
                                       lost_dbps <- names(consensus_list[sapply(consensus_list, length) < 1000]) %>% as.data.frame()
                                       
                                       # saving 
                                       write.table(lost_dbps, "results/lost_dbps.csv")
                                       
                                       ```
                                       
                                       # exporting filtered_consensus_peaks
                                       ```{r exporting filtered consensus peaks}
                                       
                                       for(i in 1:length(filtered_consensus_list)) {
                                         rtracklayer::export(filtered_consensus_list[[i]], 
                                                             paste0("results/filtered_consensus_peaks/", 
                                                                    names(filtered_consensus_list)[i], 
                                                                    "_filtered_consensus_peaks.bed"))
                                       }
                                       ```
                                       
                                       # loading in genome features
                                       ```{r creating genome feature objects}
                                       
                                       gencode_gr <- rtracklayer::import("/scratch/Shares/rinnclass/CLASS_2023/data/data/genomes/gencode.v32.annotation.gtf")
                                       
                                       # gencode genes
                                       gencode_genes <- gencode_gr[gencode_gr$type == "gene"] 
                                       
                                       # mrna_genes
                                       mrna_genes <- gencode_genes[gencode_genes$gene_type %in% "protein_coding"]
                                       
                                       # lncrna_genes
                                       lncrna_genes <- gencode_genes[gencode_genes$gene_type %in% "lncRNA"] 
                                       
                                       # mrna_lncrna_genes
                                       mrna_lncrna_genes <- gencode_genes[gencode_genes$gene_type %in% c("protein_coding","lncRNA")]
                                       
                                       # lncrna_mrna_promoters
                                       lncrna_mrna_promoters <- promoters(mrna_lncrna_genes, upstream = 1000, downstream = 1000)
                                       
                                       # lncrna_gene_ids
                                       lncrna_gene_ids <- mrna_lncrna_genes$gene_id[mrna_lncrna_genes$gene_type == "lncRNA"]
                                       
                                       # mrna_gene_ids
                                       mrna_gene_ids <-mrna_lncrna_genes$gene_id[mrna_lncrna_genes$gene_type == "protein_coding"]
                                       
                                       ```
                                       
                                       # making data frame of filtered_consensus_peak info
                                       ```{r creating num_peaks_df to track peak properties}
                                       
                                       num_peaks_df <- data.frame("dbp" = names(filtered_consensus_list),
                                                                  "num_peaks" = sapply(filtered_consensus_list, length))
                                       
                                       # total genome covered by peaks
                                       num_peaks_df$total_peak_length <- sapply(filtered_consensus_list, function(x) sum(width(x)))
                                       
                                       # creating number of promoter overlaps entry
                                       promoter_peak_counts <- count_peaks_per_feature(lncrna_mrna_promoters, filtered_consensus_list, type = "counts")
                                       
                                       # creating promoter peak_occurence for clustering - Metaplots later.
                                       promoter_peak_matrix <- count_peaks_per_feature(lncrna_mrna_promoters, filtered_consensus_list, type = "occurrence")
                                       
                                       # saving
                                       write.table(promoter_peak_matrix, "results/promoter_peak_occurrence_matrix.tsv")
                                       
                                       # summing rows to get total number of promoter overlaps
                                       num_peaks_df$peaks_overlapping_promoters <- rowSums(promoter_peak_counts)
                                       
                                       # lncrna promoter overlaps 
                                       num_peaks_df$peaks_overlapping_lncrna_promoters <- rowSums(promoter_peak_counts[,lncrna_gene_ids])
                                       
                                       # mrna promoter overlaps
                                       num_peaks_df$peaks_overlapping_mrna_promoters <- rowSums(promoter_peak_counts[,mrna_gene_ids])
                                       
                                       # Finding overlaps with gene_bodies (will take a few minutes again)
                                       # Note this takes several minutes
                                       genebody_peak_counts <- count_peaks_per_feature(mrna_lncrna_genes, 
                                                                                       filtered_consensus_list, 
                                                                                       type = "counts")
                                       
                                       # All gene bodies overlaps
                                       num_peaks_df$peaks_overlapping_genebody <- rowSums(genebody_peak_counts)
                                       
                                       # lncRNA gene bodies 
                                       num_peaks_df$peaks_overlapping_lncrna_genebody <- rowSums(genebody_peak_counts[,lncrna_gene_ids])
                                       
                                       # mRNA gene bodies
                                       num_peaks_df$peaks_overlapping_mrna_genebody <- rowSums(genebody_peak_counts[,mrna_gene_ids])
                                       
                                       ```
                                       
                                       
                                       # creating promoter peak occurence matrix
                                       This will make a matrix where promoters are cols (30K)
                                       Each will have 1 if overlapped by a given dbp : 0 if no overlaps
                                       
                                       ```{r promoter peak occurence matrix}
                                       
                                       # running count_peaks_per_feature
                                       promoter_peak_occurence_matrix <- count_peaks_per_feature(lncrna_mrna_promoters, filtered_consensus_list, 
                                                                                                 type = "occurrence")
                                       
                                       # Let's double check that all lncrna & mrna genes are accounted for:
                                       stopifnot(all(colnames(promoter_peak_occurence) == lncrna_mrna_promoters$gene_id))
                                       
                                       # saving
                                       write.table(promoter_peak_occurence, "results/lncrna_mrna_promoter_peak_occurence_matrix.tsv")
                                       
                                       # Now let's use the 'data.frame()' fucntion. Set up a bunch of colnames and populate them.
                                       peak_occurence_df <- data.frame("gene_id" = colnames(promoter_peak_occurence),
                                                                       "gene_name" = lncrna_mrna_promoters$gene_name,
                                                                       "gene_type" = lncrna_mrna_promoters$gene_type,
                                                                       "chr" = lncrna_mrna_promoters@seqnames,   
                                                                       "1kb_up_tss_start" = lncrna_mrna_promoters@ranges@start,
                                                                       "strand" = lncrna_mrna_promoters@strand,
                                                                       "number_of_dbp" = colSums(promoter_peak_occurence))
                                       
                                       # saving
                                       write_csv(peak_occurence_df, "results/peak_occurence_dataframe.csv")
                                       
                                       
                                       ```
                                       
                                       # now make a promoter data_frame that tells which dbps are bound
                                       ```{r dbp centric promoter occurence}
                                       
                                       # dbps on promoters object
                                       DBPs_on_promoter <- lncrna_mrna_promoters %>%
                                         as.data.frame() %>%
                                         dplyr::select(gene_id, gene_name)
                                       
                                       # creating promoter dbps by pivot longer of promoter_peak_occurence_matrix
                                       promoter_dbps <- promoter_peak_occurence_matrix %>%
                                         as.data.frame() %>%
                                         rownames_to_column("dbp") %>%
                                         pivot_longer(2:ncol(.), names_to = "gene_id", values_to = "occurrence") %>%
                                         filter(occurrence == 1) %>%
                                         dplyr::select(-occurrence) %>%
                                         left_join(DBPs_on_promoter)
                                       
                                       # checking Firre promoter
                                       firre_promoter <- promoter_dbps %>%
                                         filter(gene_name == "FIRRE")
                                       
                                       # XIST promoter (should be off since male cells)
                                       XIST_promoter <- promoter_dbps %>%
                                         filter(gene_name == "XIST")
                                       
                                       # GAPDH
                                       GAPDH_promoter <- promoter_dbps %>%
                                         filter(gene_name == "GAPDH")
                                       
                                       # saving
                                       promoter_dbps_df <- promoter_dbps %>% as.data.frame()
                                       write.csv(promoter_dbps, "results/promoter_dbps.csv")
                                       
                                       ```
                                       
                                       # Peaks per dbp
                                       ```{r plotting peak features}
                                       
                                       # First let's look at a histogram of peak#/DBP
                                       ggplot(num_peaks_df, aes(x = num_peaks)) + 
                                         geom_histogram(bins = 70)
                                       
                                       # saving
                                       ggsave("figures/num_peaks_hist.pdf")
                                       
                                       ```
                                       Result:
                                         (i) Most DBPs have less than 50,000 peaks
                                       
                                       
                                       # plotting num_peaks versus genome coverage.
                                       ```{r peaks vs coverage}
                                       
                                       ggplot(num_peaks_df, aes(x = num_peaks, y = total_peak_length)) +
                                         geom_point() + 
                                         geom_smooth(method = "gam", se = TRUE, color = "black", lty = 2)+
                                         ylab("BP covered") +
                                         xlab("Number of peaks") +
                                         ggtitle("Peak count vs. total bases covered")
                                       
                                       # saving 
                                       ggsave("figures/peak_num_vs_coverage.pdf")
                                       
                                       
                                       ```
                                       Result: 
                                         (i) there is a linear relationship between number of peaks and total coverage
                                       
                                       
                                       # plotting num peaks on promoters
                                       ```{r number of DBPS on promoters}
                                       
                                       ggplot(num_peaks_df,
                                              aes(x = num_peaks, y = peaks_overlapping_promoters)) +
                                         xlab("Peaks per DBP") +
                                         ylab("Number of peaks overlapping promoters") +
                                         ggtitle("Relationship Between Number of DBP Peaks and Promoter Overlaps")+
                                         geom_point() +
                                         geom_abline(slope = 1, linetype="dashed") +
                                         geom_smooth(method = "lm", se=FALSE, formula = 'y ~ x',
                                                     color = "#a8404c") +
                                         stat_regline_equation(label.x = 35000, label.y = 18000) +
                                         ylim(0,60100) +
                                         xlim(0,60100)
                                       
                                       # saving
                                       ggsave("figures/peak_num_vs_promoter_coverage.pdf")
                                       
                                       ```
                                       Result: 
                                         (i) saturation of binding events -- as you get more peaks 
                                       you stop increasing binding to promoters -- probably saturated.
                                       (ii) inflection of saturation occurs around 20,000 peaks per DBP
                                       
                                       # peak Coverage on gene bodies
                                       ```{r peak coverage on gene bodies}
                                       
                                       ggplot(num_peaks_df,
                                              aes(x = num_peaks, y = peaks_overlapping_genebody)) +
                                         xlab("Peaks per DBP") +
                                         ylab("Number of peaks overlapping genes") +
                                         ggtitle("Relationship Between Number of DBP Peaks and Gene Body Overlaps")+
                                         geom_point() +
                                         geom_abline(slope = 1, linetype="dashed") +
                                         geom_smooth(method = "lm", se=F, formula = 'y ~ x',
                                                     color = "#a8404c") +
                                         stat_regline_equation(label.x = 35000, label.y = 18000) +
                                         ylim(0,60100) +
                                         xlim(0,60100)
                                       
                                       # saving
                                       ggsave("figures/4_peak_num_vs_gene_body_coverage.pdf")
                                       ```
                                       Result: 
                                         (i) Gene bodies explain almost all the places of binding in the genome (.9x slope)
                                       (ii) thus roughly 90% of binding events occur within or overlap gene bodies
                                       
                                       # Density plot of binding events
                                       Let's make a density plot of num DBPs bound per promoter

```{r density plot of DBP localization events}

ggplot(peak_occurence_df, aes(x = number_of_dbp)) +
geom_density(alpha = 0.2, color = "#424242", fill = "#424242") +
  theme_paperwhite() +
  xlab(expression("Number of DBPs")) +
  ylab(expression("Density")) +
  ggtitle("Promoter binding events",
          subtitle = "mRNA and lncRNA genes") 

# saving
ggsave("figures/num_binding_events_per_promoter.pdf")

```
Result: 
(i) very interesting that the promoter binding is bimodal !
(ii) most promoters have upto 100 dbps then a lag and 2nd dist at ~200dbs
(iii) There are two types of promoter binding - (1) normal (2) super-binders


# promoters with out binding events
Lets find how many promoters don't have any DBPs bound
                                       ```{r prmoters with out binding events}
                                       
                                       unbound_promoters <- peak_occurence_df %>% 
                                         filter(peak_occurence_df$number_of_dbp < 1)
                                       
                                       # how many are there?
                                       nrow(unbound_promoters)
                                       # so there are only a few 6,720 promoters that don't have binding evetns (~10%)
                                       
                                       #  let's put it in a folder called results. We will always use this folder structure
                                       write_csv(unbound_promoters, "results/unbound_promoters.csv")
                                       
                                       ```
                                       
                                       # lncRNA versus mRNA promoter binding
                                       Let's compare the binding patterns of lncRNA vs mRNA promoters.
```{r lncrna vs mrna promoter binding}

num_peaks_dfl <- num_peaks_df %>%
  dplyr::select(-peaks_overlapping_promoters) %>%
  pivot_longer(cols = peaks_overlapping_lncrna_promoters:peaks_overlapping_mrna_promoters,
               names_to = "gene_type",
               values_to = "peaks_overlapping_promoters") %>%
  mutate(gene_type = gsub("peaks_overlapping_", "", gene_type))

# plotting
ggplot(num_peaks_dfl, aes(x = num_peaks, y = peaks_overlapping_promoters, 
                         col = gene_type)) +
         geom_point() +
         geom_abline(slope = 1, linetype="dashed") +
  geom_smooth(method = "lm", se = FALSE, formula = "y ~ x") +
  stat_regline_equation() +
  scale_color_manual(values = c("#a8404c", "#424242"))+
  xlab("Peaks per DBP") +
  ylab("Peaks Overlapping Promoters") +
  ggtitle("Number of DBP Peaks and Promoter Overlaps")

# saving
ggsave("figures/peaks_overlaps_relationship_by_gene_type.pdf", height = 5, width = 8)
```
Result:
(i) in general mRNA promoters have more overlaps than lncRNAs
(ii) both have linear trend with more peaks more overlaps
(iii) mRNA promoters are 4x more likely to overlaps a dbp than lncRNA promoters

# creating distance matrix & dendrogram
```{r distance matrix and dendrogram}

# creating distance matrix
peak_occurence_dist <- dist(promoter_peak_matrix, method = "binary")

# clustering distance values
bin_hier <- hclust(peak_occurence_dist, method = "complete")

# Dendrogram of binding profiles by promoter (not binding profile - below)
ggdendro::ggdendrogram(bin_hier, rotate = FALSE,  size = 3,
                       theme_dendro = TRUE) +
   coord_flip() +
   scale_y_continuous() +
   scale_x_continuous(position = "top") +
   scale_x_continuous(breaks = seq_along(bin_hier$labels[bin_hier$order]),
             labels = bin_hier$labels[bin_hier$order], position = "top",
             expand = c(0,0)) +
   theme(axis.text.x = element_text(angle = 90, hjust  = 1)) +
   theme(axis.text.y = element_text(angle = 0,hjust = 1)) +
   scale_y_reverse(expand = c(0.01, 0)) +
   theme(
     plot.background = element_blank(),
     panel.grid.major = element_blank(),
   panel.grid.minor = element_blank(),
     panel.border = element_blank()
   )


ggsave("figures/promoter_overlap_dendrogram.pdf", height = 49, width = 12)
```
Result:
(i) histone mods cluster together and seperate from others
(ii) Pol2 and S5-Pol2 don't cluster together as expected
                                       
                                       
                                       # Using profile_tss for all 430 DBPs
                                       # ! this takes ~45 min !
                                       ```{r metaplot DF of binding profiles over promoter}
                                       
                                       # establishing DF
                                       metaplot_df <- data.frame(x = integer(), dens = numeric(), dbp = character())
                                       
                                       # for loop to populate DF 
                                       for(i in 1:length(filtered_consensus_list)) {
                                         print(names(filtered_consensus_list)[[i]])
                                         tmp_df <- profile_tss(filtered_consensus_list[[i]], lncrna_mrna_promoters)
                                         tmp_df$dbp <- names(filtered_consensus_list)[[i]]
                                         metaplot_df <- bind_rows(metaplot_df, tmp_df)
                                         
                                       }
                                       
                                       # saving
                                       write_rds(metaplot_df, "results/metaplot_df_final.rds")
                                       
                                       ```
                                       
                                       
                                       # creating distance matrix of binding profile correlations
                                       ```{r scaling and plotting dendrogram of binding similarity by promoter}
                                       
                                       metaplot_filtered_matrix <- metaplot_df %>% 
                                         pivot_wider(names_from = x, values_from = dens) %>%
                                         column_to_rownames("dbp") %>%
                                         as.matrix()
                                       mm_scaled <- t(scale(t(metaplot_filtered_matrix)))
                                       metaplot_hclust <- hclust(dist(mm_scaled), method = "complete")
                                       
                                       # plotting relationship between binding profiles
                                       plot(metaplot_hclust)
                                       pdf("figures/tss_profile_dendrogram.pdf", height = 10, width = 27)
                                       par(cex=0.3)
                                       plot(metaplot_hclust)
                                       dev.off()
                                       ```
                                       
                                       
                                       # heat map of binding profile over +/- 1kb TSS window
                                       ```{r heatmap of profile over TSS (+/- 1Kb)}
                                       
                                       pdf("figures/tss_profile_heatmap.pdf", height = 35, width = 10)
                                       Heatmap(mm_scaled, cluster_columns = FALSE, col = col_fun, border = TRUE, 
                                               show_column_names = FALSE,
                                               use_raster = TRUE,
                                               column_split = split,
                                               column_gap = unit(0, "mm"),row_names_gp = gpar(fontsize = 7))
                                       dev.off()
                                       ```
                                       Result:
                                         (i) most binding profiles are over the TSS
                                       (ii) some binding profiles are more narrow directly over TSS
                                       (iii) histone mods tend to be depleted at TSS (probably NFR)
                                       
                                       # looking more into the patterns of ChIP binding profiles relative to TSS
                                       # breaking this down into several clusters and plotting each
                                       ```{r how many groups are there?}
                                       
                                       # Let's see how many DBPs have different patterns.
                                       # manually checking at 6 we get a new group of 10 - after that not much changes
                                       clusters <- cutree(metaplot_hclust, k=6)
                                       table(clusters)
                                       
                                       # plot cluter 1
                                       dev.new()
                                       pdf("figures/cluster_1_heatmap.pdf", height = 35, width = 10)
                                       
                                       Heatmap(mm_scaled[clusters == 1,], cluster_columns = FALSE, col = col_fun, border = TRUE, 
                                               show_column_names = FALSE,
                                               use_raster = TRUE,
                                               column_split = split,
                                               column_gap = unit(0, "mm"),row_names_gp = gpar(fontsize = 7))
                                       graphics.off()
                                       
                                       # cluster 2
                                       dev.new()
                                       pdf("figures/cluster_2_heatmap.pdf", height = 35, width = 10)
                                       
                                       Heatmap(mm_scaled[clusters == 2,], cluster_columns = FALSE, col = col_fun, border = TRUE, 
                                               show_column_names = FALSE,
                                               use_raster = TRUE,
                                               column_split = split,
                                               column_gap = unit(0, "mm"),row_names_gp = gpar(fontsize = 7))
                                       graphics.off()
                                       
                                       # cluster 3
                                       dev.new()
                                       pdf("figures/cluster_3_heatmap.pdf", height = 35, width = 10)
                                       Heatmap(mm_scaled[clusters == 3,], cluster_columns = FALSE, col = col_fun, border = TRUE, 
                                               show_column_names = FALSE,
                                               use_raster = TRUE,
                                               column_split = split,
                                               column_gap = unit(0, "mm"),row_names_gp = gpar(fontsize = 7))
                                       graphics.off()
                                       # very narrow binding over tss !
                                       
                                       # cluster 4
                                       dev.new()
                                       pdf("figures/cluster_4_heatmap.pdf", height = 35, width = 10)
                                       Heatmap(mm_scaled[clusters == 4,], cluster_columns = FALSE, col = col_fun, border = TRUE, 
                                               show_column_names = FALSE,
                                               use_raster = TRUE,
                                               column_split = split,
                                               column_gap = unit(0, "mm"),row_names_gp = gpar(fontsize = 7))
                                       graphics.off()
                                       
                                       # looks like TSS depletion and histone mods
                                       
                                       # cluster 5 only 1 DBP
                                       names(clusters[5])
                                       
                                       # cluster 6
                                       dev.new()
                                       pdf("figures/cluster_6_heatmap.pdf", height = 35, width = 10)
                                       Heatmap(mm_scaled[clusters == 6,], cluster_columns = FALSE, col = col_fun, border = TRUE, 
                                               show_column_names = FALSE,
                                               use_raster = TRUE,
                                               column_split = split,
                                               column_gap = unit(0, "mm"),row_names_gp = gpar(fontsize = 7))
                                       graphics.off()
                                       
                                       # TSS depletion and histone mods again.
                                       
                                       
                                       ```
                                       Result:
                                         (i) Cluster 1 and 2 have broad TSS binding (cluster 2 is more narrow than cluster 1)
                                       (ii) Cluster 3 has very narrow binding right over tss
                                       (iii) cluster 4 and 6 have depletion over TSS (also mostly histone mods and EZH2)
                                       (iv) cluster 5 was not plotted as only has one DBP.
                                       
                                       
                                       # metaplots for each DBP by lncRNA and mRNA promoters
                                       ```{r}
                                       
                                       #TODO test by loading in env variable
                                       #setting up lncrna DF.
                                       lncrna_metaplot_df <- data.frame(x = integer(), dens = numeric(), dbp = character())
                                       
                                       # for loop to populate DF with overlap density in lncrna promoters
                                       for(i in 1:length(filtered_consensus_list)) {
                                         print(names(filtered_consensus_list)[[i]])
                                         tmp_df <- profile_tss(filtered_consensus_list[[i]], lncrna_mrna_promoters = lncrna_promoters)
                                         tmp_df$dbp <- names(filtered_consensus_list)[[i]]
                                         lncrna_metaplot_df <- bind_rows(lncrna_metaplot_df, tmp_df)
                                         
                                       }
                                       
                                       # saving
                                       write_rds(lncrna_metaplot_df, "results/lncRNA_metaplot_df_final.rds")
                                       
                                       # now for mRNAs 
                                       mrna_metaplot_df <- data.frame(x = integer(), dens = numeric(), dbp = character())
                                       
                                       # for loop to populate mRNA_metaplot
                                       for(i in 1:length(filtered_consensus_list)) {
                                         print(names(filtered_consensus_list)[[i]])
                                         tmp_df <- profile_tss(filtered_consensus_list[[i]], lncrna_mrna_promoters = mrna_promoters)
                                         tmp_df$dbp <- names(filtered_consensus_list)[[i]]
                                         mrna_metaplot_df <- bind_rows(mrna_metaplot_df, tmp_df)
                                         
                                       }
                                       
                                       # saving mRNA metaplots
                                       write_rds(mrna_metaplot_df, "results/mrna_metaplot_df_final.rds")
                                       
                                       
                                       # now adding the information of gene type
                                       mrna_metaplot_df$gene_type <- "mRNA"
                                       lncrna_metaplot_df$gene_type <- "lncRNA"
                                       combined_metaplot_profile <- bind_rows(mrna_metaplot_df, lncrna_metaplot_df)
                                       
                                       # saving
                                       write_rds(mrna_metaplot_df, "results/metaplot_df_final.rds")
                                       ```
                                       
                                       # profile plots for lncRNA and mRNA promoters seperated for each DBP
                                       ```{r lncRNA and mRNA profile plots for each dbp}
                                       
                                       ggplot(combined_metaplot_profile, 
                                              aes(x = x, y = dens, color = gene_type )) +
                                         geom_vline(xintercept = 0, lty = 2) + 
                                         geom_line(size = 1.5) + 
                                         facet_wrap(dbp ~ ., scales = "free_y") +
                                         ggtitle("Promoter Metaplot") + 
                                         scale_x_continuous(breaks = c(-1000, 0, 1000),
                                                            labels = c("-1kb", "TSS", "+1kb"),
                                                            name = "") + 
                                         ylab("Peak frequency") +
                                         scale_color_manual(values = c("#424242","#a8404c"))
                                       
                                       # saving
                                       ggsave("figures/mega_meta_plot_lncRNA-mRNA.pdf", width = 49, height = 12)
                                       
                                       ```
                                       Result: 
                                         