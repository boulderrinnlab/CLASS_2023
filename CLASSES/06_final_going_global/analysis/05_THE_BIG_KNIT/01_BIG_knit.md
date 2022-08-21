01\_BIG\_KNIT
================
JR
8/17/2022

# Goal:

Here we aim to download all available DNA binding protein (DBP) profiles
in a single cell state (measured by ChIP-seq). This will allow us to
investigate the binding properties of hundreds of DBPs in the same
cellular context or background. We aim to address several questions:

1.  What are the number of peaks and genome coverage for each DBP?

2.  What are the binding preferences for promoters, gene-bodies and
    intergenic genomic regions?

3.  What are the similarities and differences across DBPs based on their
    genome-wide binding profiles genome-wide?

4.  What properties or preferences do promoters have for binding events.

5.  How does binding to a promoter affect the transcriptional output of
    that promoter?

6.  Are there reservoir promoters in HepG2 as defined in k562
    previously?

To address these questions we have curated a set of over 1,000 ChIPseq
data sets comprised of 48 DBPs in HEPG2 cells from the ENCODE
consortium. We required duplicate ChIP-seq experiments for a given DBP
and other criterion that can be found here :

<https://www.encodeproject.org/report/?type=Experiment&status=released&assay_slims=DNA+binding&biosample_ontology.term_name=HepG2&assay_title=TF+ChIP-seq&biosample_ontology.classification=cell+line&files.read_length=100&files.read_length=76&files.read_length=75&files.read_length=36&assay_title=Control+ChIP-seq&assay_title=Histone+ChIP-seq&files.run_type=single-ended>

## These samples were selected on the following criteria:

1.  “chromatin” interaction data, then DNA binding data, cell line
    HEPG2, “TF-Chip-seq”.
2.  We further selected “TF Chip-seq”, “Control chip-seq” and “Histone
    Chip-seq”.
3.  We selected several read lengths to get the most DNA binding
    proteins (DBPs)
4.  Read lengths: 100, 76, 75, 36
5.  ONLY SINGLE END READS (this eliminates 54 samples)

### Experimental data was downloading by (ENCODE report.tsv):

<https://www.encodeproject.org/report.tsv?type=Experiment&status=released&assay_slims=DNA+binding&biosample_ontology.term_name=HepG2&assay_title=TF+ChIP-seq&biosample_ontology.classification=cell+line&files.read_length=100&files.read_length=76&files.read_length=75&files.read_length=36&assay_title=Control+ChIP-seq&assay_title=Histone+ChIP-seq&files.run_type=single-ended>

### The FASTQ files were downloaded with:

“<https://www.encodeproject.org/metadata/?status=released&assay_slims=DNA+binding&biosample_ontology.term_name=HepG2&assay_title=TF+ChIP-seq&biosample_ontology.classification=cell+line&files.read_length=100&files.read_length=76&files.read_length=75&files.read_length=36&assay_title=Control+ChIP-seq&assay_title=Histone+ChIP-seq&files.run_type=single-ended&type=Experiment>”

MD5sums were checked with all passing (see encode\_file\_info function
to reterive MD5Sum values that are not available from the encode portal
(/util)

### Processing data:

We processed all the read alignments and peak calling using the NF\_CORE
ChIP-seq pipeline: (nfcore/chipseq v1.2.2)

## Next we created consensus peaks that overlap in both replicates

Our strategy was to take peaks in each replicate and find all
overlapping peak windows. We then took the union length of the
overlapping range in each peak window. We further required that a DBP
has at least 1,000 peaks (the min quartile of peaks observed)

Result: (i) we have 430 DBPs that have at least 1,000 peaks in both
replicates

## We will start with loading filtered consensus peaks (described above)

## All subequent analyses are contained within this document.

## To find results, search for “Result:”

# Loading in all .broadPeak files for all dbps and their replicates

``` r
broadpeakfilepath <- "/scratch/Shares/rinnclass/CLASS_2023/data/data/peaks"

# importing peaks as granges with "import_peaks" function
# can take 10+ min !
peak_list <- import_peaks(consensus_file_path = broadpeakfilepath)
```

# creating consensus peaks

Using consensus\_from\_reduced function that requires peak list and dbp
object. This will be done in 2 steps (1) created dbp object of unique
names for dbp (2) run consensus\_from\_reduced

``` r
  # Creating unique DBP object for create_consensus_peaks_from_reduced
dbp <- unique(sapply(names(peak_list), function(x) {
   unlist(strsplit(x, "_"))[1]
}))

# now run our function consensus_from_reduced
consensus_list <- lapply(dbp, consensus_from_reduced, peak_list)

# adding names to the GRange list
names(consensus_list) <- dbp
```

# exploring number of peaks per DBP for filtering

``` r
# creating list of num_peaks per dbp
num_peaks <- sapply(consensus_list, length)

# plotting
hist(num_peaks, breaks = 1000)
```

![](01_BIG_knit_files/figure-gfm/exploring%20number%20of%20peaks%20in%20consenus%20peaks-1.png)<!-- -->

``` r
hist(num_peaks, breaks = 1000, xlim = c(0,3000))
```

![](01_BIG_knit_files/figure-gfm/exploring%20number%20of%20peaks%20in%20consenus%20peaks-2.png)<!-- -->

``` r
ggsave("figures/hist_num_peaks.pdf")
```

Result: (i)seems like 1000 peaks should be the min moving forward.

# filtering consensus\_list to dbps with &gt; 1000 peaks

``` r
# filtering to 1000 peaks
filtered_consensus_list <- consensus_list[sapply(consensus_list, length) > 1000]

# saving 
save(filtered_consensus_list, file = "results/filtered_consensus_list.RData")

# keeping track of DBPs lost
lost_dbps <- names(consensus_list[sapply(consensus_list, length) < 1000]) %>% as.data.frame()

# saving 
write.table(lost_dbps, "results/lost_dbps.csv")
```

# exporting filtered\_consensus\_peaks

``` r
for(i in 1:length(filtered_consensus_list)) {
  rtracklayer::export(filtered_consensus_list[[i]], 
                      paste0("results/filtered_consensus_peaks/", 
                             names(filtered_consensus_list)[i], 
                             "_filtered_consensus_peaks.bed"))
}
```

# loading in genome features

``` r
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

# making data frame of filtered\_consensus\_peak info

``` r
num_peaks_df <- data.frame("dbp" = names(filtered_consensus_list),
                           "num_peaks" = sapply(filtered_consensus_list, length))

# total genome covered by peaks
num_peaks_df$total_peak_length <- sapply(filtered_consensus_list, function(x) sum(width(x)))

# creating number of promoter overlaps entry
promoter_peak_counts <- count_peaks_per_feature(lncrna_mrna_promoters, filtered_consensus_list, type = "counts")

# creating promoter peak_occurrence for clustering - Metaplots later.
promoter_peak_occurrence_matrix <- count_peaks_per_feature(lncrna_mrna_promoters, filtered_consensus_list, type = "occurrence")

# saving
write.table(promoter_peak_occurrence_matrix, "results/promoter_peak_occurrence_matrix.tsv")

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

# creating promoter peak occurrence matrix &gt; DF

This will make a matrix where promoters are cols (30K) Each will have 1
if overlapped by a given dbp : 0 if no overlaps Converting to dataframe

``` r
# Now let's use the 'data.frame()' fucntion. Set up a bunch of colnames and populate them.
peak_occurrence_df <- data.frame("gene_id" = colnames(promoter_peak_occurrence_matrix),
                                "gene_name" = lncrna_mrna_promoters$gene_name,
                                "gene_type" = lncrna_mrna_promoters$gene_type,
                                "chr" = lncrna_mrna_promoters@seqnames,   
                                "1kb_up_tss_start" = lncrna_mrna_promoters@ranges@start,
                                "strand" = lncrna_mrna_promoters@strand,
                                "number_of_dbp" = colSums(promoter_peak_occurrence_matrix))

# saving
write_csv(peak_occurrence_df, "results/peak_occurrence_dataframe.csv")
```

# now make a promoter data\_frame that tells which dbps are bound

``` r
# dbps on promoters object
DBPs_on_promoter <- lncrna_mrna_promoters %>%
                    as.data.frame() %>%
  dplyr::select(gene_id, gene_name)

# creating promoter dbps by pivot longer of promoter_peak_occurrence_matrix
promoter_dbps <- promoter_peak_occurrence_matrix %>%
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

``` r
# First let's look at a histogram of peak#/DBP
 ggplot(num_peaks_df, aes(x = num_peaks)) + 
  geom_histogram(bins = 70)
```

![](01_BIG_knit_files/figure-gfm/plotting%20peak%20features-1.png)<!-- -->

``` r
# saving
ggsave("figures/num_peaks_hist.pdf")
```

Result: (i) Most DBPs have less than 50,000 peaks

# plotting num\_peaks versus genome coverage.

``` r
ggplot(num_peaks_df, aes(x = num_peaks, y = total_peak_length)) +
  geom_point() + 
  geom_smooth(method = "gam", se = TRUE, color = "black", lty = 2)+
  ylab("BP covered") +
  xlab("Number of peaks") +
  ggtitle("Peak count vs. total bases covered")
```

![](01_BIG_knit_files/figure-gfm/peaks%20vs%20coverage-1.png)<!-- -->

``` r
# saving 
ggsave("figures/peak_num_vs_coverage.pdf")
```

Result: (i) there is a linear relationship between number of peaks and
total coverage

# plotting num peaks on promoters

``` r
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
```

![](01_BIG_knit_files/figure-gfm/number%20of%20DBPS%20on%20promoters-1.png)<!-- -->

``` r
# saving
ggsave("figures/peak_num_vs_promoter_coverage.pdf")
```

Result: (i) saturation of binding events – as you get more peaks you
stop increasing binding to promoters – probably saturated. (ii)
inflection of saturation occurs around 20,000 peaks per DBP

# peak Coverage on gene bodies

``` r
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
```

![](01_BIG_knit_files/figure-gfm/peak%20coverage%20on%20gene%20bodies-1.png)<!-- -->

``` r
# saving
ggsave("figures/4_peak_num_vs_gene_body_coverage.pdf")
```

Result: (i) Gene bodies explain almost all the places of binding in the
genome (.9x slope) (ii) thus roughly 90% of binding events occur within
or overlap gene bodies

# Density plot of binding events

Let’s make a density plot of num DBPs bound per promoter

``` r
ggplot(peak_occurrence_df, aes(x = number_of_dbp)) +
geom_density(alpha = 0.2, color = "#424242", fill = "#424242") +
  theme_paperwhite() +
  xlab(expression("Number of DBPs")) +
  ylab(expression("Density")) +
  ggtitle("Promoter binding events",
          subtitle = "mRNA and lncRNA genes") 
```

![](01_BIG_knit_files/figure-gfm/density%20plot%20of%20DBP%20localization%20events-1.png)<!-- -->

``` r
# saving
ggsave("figures/num_binding_events_per_promoter.pdf")
```

Result: (i) very interesting that the promoter binding is bimodal ! (ii)
most promoters have upto 100 dbps then a lag and 2nd dist at \~200dbs
(iii) There are two types of promoter binding - (1) normal (2)
super-binders

# promoters with out binding events

Lets find how many promoters don’t have any DBPs bound

``` r
unbound_promoters <- peak_occurrence_df %>% 
  filter(peak_occurrence_df$number_of_dbp < 1)

# how many are there?
nrow(unbound_promoters)
```

    ## [1] 9448

``` r
# so there are only a few 6,720 promoters that don't have binding evetns (~10%)

#  let's put it in a folder called results. We will always use this folder structure
write_csv(unbound_promoters, "results/unbound_promoters.csv")
```

# lncRNA versus mRNA promoter binding

Let’s compare the binding patterns of lncRNA vs mRNA promoters.

``` r
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
```

![](01_BIG_knit_files/figure-gfm/lncrna%20vs%20mrna%20promoter%20binding-1.png)<!-- -->

``` r
# saving
ggsave("figures/peaks_overlaps_relationship_by_gene_type.pdf", height = 5, width = 8)
```

Result: (i) in general mRNA promoters have more overlaps than lncRNAs
(ii) both have linear trend with more peaks more overlaps (iii) mRNA
promoters are 4x more likely to overlaps a dbp than lncRNA promoters

# creating distance matrix & dendrogram

``` r
# creating distance matrix
peak_occurrence_dist <- dist(promoter_peak_occurrence_matrix, method = "binary")

# clustering distance values
bin_hier <- hclust(peak_occurrence_dist, method = "complete")

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
```

![](01_BIG_knit_files/figure-gfm/distance%20matrix%20and%20dendrogram-1.png)<!-- -->

``` r
ggsave("figures/promoter_overlap_dendrogram.pdf", height = 49, width = 12)
```

Result: (i) histone mods cluster together and seperate from others (ii)
Pol2 and S5-Pol2 don’t cluster together as expected

# Using profile\_tss for all 430 DBPs

# ! this takes \~45 min !

``` r
# establishing DF
metaplot_df <- data.frame(x = integer(), dens = numeric(), dbp = character())

# for loop to populate DF 
for(i in 1:length(filtered_consensus_list)) {
  # print(names(filtered_consensus_list)[[i]])
  tmp_df <- profile_tss(filtered_consensus_list[[i]], lncrna_mrna_promoters)
  tmp_df$dbp <- names(filtered_consensus_list)[[i]]
  metaplot_df <- bind_rows(metaplot_df, tmp_df)
  
}

# saving
write_rds(metaplot_df, "results/metaplot_df_final.rds")
metaplot_df <- read_rds("results/metaplot_df_final.rds")
```

# creating distance matrix of binding profile correlations

``` r
metaplot_filtered_matrix <- metaplot_df %>% 
  pivot_wider(names_from = x, values_from = dens) %>%
  column_to_rownames("dbp") %>%
  as.matrix()
mm_scaled <- t(scale(t(metaplot_filtered_matrix)))
metaplot_hclust <- hclust(dist(mm_scaled), method = "complete")

# plotting relationship between binding profiles
plot(metaplot_hclust)
```

![](01_BIG_knit_files/figure-gfm/scaling%20and%20plotting%20dendrogram%20of%20binding%20similarity%20by%20promoter-1.png)<!-- -->

``` r
pdf("figures/tss_profile_dendrogram.pdf", height = 10, width = 27)
par(cex=0.3)
plot(metaplot_hclust)
dev.off()
```

    ## png 
    ##   2

# heat map of binding profile over +/- 1kb TSS window

``` r
pdf("figures/tss_profile_heatmap.pdf", height = 35, width = 10)
Heatmap(mm_scaled, cluster_columns = FALSE,  border = TRUE, 
        show_column_names = FALSE,
        use_raster = TRUE,
        column_gap = unit(0, "mm"),row_names_gp = gpar(fontsize = 7))
dev.off()
```

    ## png 
    ##   2

Result: (i) most binding profiles are over the TSS (ii) some binding
profiles are more narrow directly over TSS (iii) histone mods tend to be
depleted at TSS (probably NFR)

# looking more into the patterns of ChIP binding profiles relative to TSS

# breaking this down into several clusters and plotting each

``` r
# Let's see how many DBPs have different patterns.
# manually checking at 6 we get a new group of 10 - after that not much changes
clusters <- cutree(metaplot_hclust, k=6)
table(clusters)
```

    ## clusters
    ##   1   2   3   4   5   6 
    ## 281 130  10   5   1   3

``` r
# plot cluter 1
dev.new()
pdf("figures/cluster_1_heatmap.pdf", height = 35, width = 10)

Heatmap(mm_scaled[clusters == 1,], cluster_columns = FALSE,  border = TRUE, 
        show_column_names = FALSE,
        use_raster = TRUE,
        column_gap = unit(0, "mm"),row_names_gp = gpar(fontsize = 7))
graphics.off()

# cluster 2
dev.new()
pdf("figures/cluster_2_heatmap.pdf", height = 35, width = 10)

Heatmap(mm_scaled[clusters == 2,], cluster_columns = FALSE,  border = TRUE, 
        show_column_names = FALSE,
        use_raster = TRUE,
        column_gap = unit(0, "mm"),row_names_gp = gpar(fontsize = 7))
graphics.off()

# cluster 3
dev.new()
pdf("figures/cluster_3_heatmap.pdf", height = 35, width = 10)
Heatmap(mm_scaled[clusters == 3,], cluster_columns = FALSE,  border = TRUE, 
        show_column_names = FALSE,
        use_raster = TRUE,
        column_gap = unit(0, "mm"),row_names_gp = gpar(fontsize = 7))
graphics.off()
# very narrow binding over tss !

# cluster 4
dev.new()
pdf("figures/cluster_4_heatmap.pdf", height = 35, width = 10)
Heatmap(mm_scaled[clusters == 4,], cluster_columns = FALSE,  border = TRUE, 
        show_column_names = FALSE,
        use_raster = TRUE,
        column_gap = unit(0, "mm"),row_names_gp = gpar(fontsize = 7))
graphics.off()

# looks like TSS depletion and histone mods

# cluster 5 only 1 DBP
names(clusters[5])
```

    ## [1] "AKAP8"

``` r
# cluster 6
dev.new()
pdf("figures/cluster_6_heatmap.pdf", height = 35, width = 10)
Heatmap(mm_scaled[clusters == 6,], cluster_columns = FALSE, border = TRUE, 
        show_column_names = FALSE,
        use_raster = TRUE,
        column_gap = unit(0, "mm"),row_names_gp = gpar(fontsize = 7))
graphics.off()

# TSS depletion and histone mods again.
```

Result: (i) Cluster 1 and 2 have broad TSS binding (cluster 2 is more
narrow than cluster 1) (ii) Cluster 3 has very narrow binding right over
tss (iii) cluster 4 and 6 have depletion over TSS (also mostly histone
mods and EZH2) (iv) cluster 5 was not plotted as only has one DBP.

# establishing lncRNA and mRNA promoters (+/- 1kb)

``` r
# creating promoters just in case:
lncrna_promoters <- lncrna_mrna_promoters[lncrna_mrna_promoters$gene_type == "lncRNA"]


mrna_promoters <- lncrna_mrna_promoters[lncrna_mrna_promoters$gene_type == "protein_coding"]
```

# metaplots for each DBP by lncRNA and mRNA promoters

``` r
#setting up lncrna DF.
lncrna_metaplot_df <- data.frame(x = integer(), dens = numeric(), dbp = character())

# for loop to populate DF with overlap density in lncrna promoters
for(i in 1:length(filtered_consensus_list)) {
  #print(names(filtered_consensus_list)[[i]])
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
 # print(names(filtered_consensus_list)[[i]])
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

``` r
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
```

![](01_BIG_knit_files/figure-gfm/lncRNA%20and%20mRNA%20profile%20plots%20for%20each%20dbp-1.png)<!-- -->

``` r
# saving
ggsave("figures/mega_meta_plot_lncRNA-mRNA.pdf", width = 49, height = 12)
```

Result: (i) Most DBPs show similar binding profiles between lncRNA/mRNA
promoters (ii) a few interesting ones EZH2 is more bound over TSS of
lncRNA promoters (iii) H2AFZ is the most depleted over TSS (iv) H3K4me2
is more enriched over lncRNA TSS (v) H3K4me1 is slightly more over
lncRNA promoters (vi) ZNF460 has different pattern between mRNA and
lncRNA promoters

# DBP Binding versus RNA expression output

Now we are going to explore how binding to a promoter relates to gene
expression output. We will use RNAseq data from ENCODE in HEPG2. The
data is RNA from differnet fractions – it was run throguh NF\_CORE
RNAseq pipeline v1.4.2

The data files are here:

``` bash
# ENCODE report on experimental accessions used (fractionated and total RNA HEPG2)

# wget -O samples.txt "https://www.encodeproject.org/report.tsv?type=Experiment&status=released&assay_slims=Transcription&assay_slims=Transcription&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&biosample_ontology.term_name=HepG2&biosample_ontology.classification=cell+line&assay_title=total+RNA-seq&files.read_length=50&limit=all&advancedQuery=date_released:[2009-01-01%20TO%202021-12-31]"

# This will give us a text file of the file names.  

# wget -O files.txt "https://www.encodeproject.org/batch_download/?type=Experiment&status=released&assay_slims=Transcription&assay_slims=Transcription&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&biosample_ontology.term_name=HepG2&biosample_ontology.classification=cell+line&assay_title=total+RNA-seq&files.read_length=50&limit=all&advancedQuery=date_released:[2009-01-01%20TO%202021-12-31]"
```

First we want to load in our final\_samplesheet from
01\_DESEQ\_counts.Rmd \# Reading in sample sheet

``` r
# First let's read in the sample sheet to know what is what
samplesheet <- read_rds("../../../05_R_analyses/05_RNAseq/01_differential_expression/results/final_samplesheet.rds")
```

# reading in TPM values from Salmon for our analyses

``` r
# reading in salmon tpm
salmon_tpm <- read.csv("../../../05_R_analyses/05_RNAseq/00_RNAseq_download_NF_core_pipeline/00_NF_CORE_RNAseq_Pipeline_run/results/salmon/salmon_merged_gene_tpm.csv")

# TPM table is in same order as samplesheet
tpm <- salmon_tpm %>% 
  pivot_longer(cols = 2:ncol(.), names_to = "sample_id", values_to = "tpm") %>%
  merge(samplesheet) %>%
  group_by(gene_id, condition) %>%
  summarize(tpm = mean(tpm, na.rm = T)) %>%
  pivot_wider(names_from = condition, values_from = tpm, names_prefix = "tpm_")
```

# determining how many DBPs are bound at each promoter

``` r
# peak_occurrence_df is loaded in our environment and contains how many DBPs bind to each promoter (rows)
promoter_features_df <- merge(peak_occurrence_df, tpm)

# DBP density plot
ggplot(promoter_features_df, aes(x = number_of_dbp)) +
  geom_density() 
```

![](01_BIG_knit_files/figure-gfm/loading%20in%20peak%20features%20data%20frame-1.png)<!-- -->

``` r
# saving
  ggsave("figures/DBP_binding_density_plot.pdf")
```

result: There appears to be a bimodal distribution of binding at
promoters. those with less than 100dbps and those with more than 250

# Abundance of genes in each cellular fraction

``` r
# First we need to the tpm DF into a matrix


tpm_matrix <- tpm %>% 
  column_to_rownames("gene_id") %>%
  as.matrix()
tpm_scaled <- t(scale(t(tpm_matrix)))
tpm_scaled <- tpm_scaled[complete.cases(tpm_scaled),]


# plotting
new.env()
```

    ## <environment: 0x26f54a8a8>

``` r
pdf("figures/heatmap_expression.pdf", height =49, width = 12)
pheatmap::pheatmap(tpm_scaled, show_rownames = FALSE)
graphics.off()
```

RESULT: (1) Most RNAs are abundant in the nucleus (2) Some RNAs
expressed in total that are not in other fractions

# Plotting binding versus expression

Now let’s examine how binding effects expression. We have published
previously that the more binding events at a promoter the more abundant
the expression is. Let’s see if this holds for our current subset of
DBPs

``` r
# plotting binding vs total RNA expression
ggplot(promoter_features_df, 
            aes(y = log2(tpm_homo_sapiens_hepg2 + 0.001), x = number_of_dbp, color = gene_type)) + 
geom_point(data = promoter_features_df %>% filter(tpm_homo_sapiens_hepg2 > 0.001),
             shape = 17, alpha = 0.7) +
  geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs")) +
  stat_cor() +
  geom_smooth(method = "lm") +
  scale_x_continuous(expand = c(0,0)) +
  scale_color_manual(values = c("#a8404c", "#424242"), name = "Gene type") + 
  ggtitle("Expression vs. promoter binding events") + 
  xlab(expression('Number of TFs')) +
  ylab(expression(log[2](TPM))) 
```

![](01_BIG_knit_files/figure-gfm/DBP%20promoter%20binding%20versus%20total%20RNA%20expression-1.png)<!-- -->

``` r
ggsave("figures/binding_vs_expression_total_rna.pdf")
```

Result: (1) There is a linear trend with number of DBPS and expression
levels (2) There is a population of genes that have numerous DBPs with
low expression

# Binding versus nuclear expression

Let’s see if the binding versus expression holds in the nuclear fraction

``` r
# Now let's make a similar plot for nuclear RNA abundance versus #DBPs bound to their promoter
ggplot(promoter_features_df, 
            aes(y = log2(tpm_homo_sapiens_nuclear_fraction + 0.001), x = number_of_dbp, color = gene_type)) + 
  geom_point(data = promoter_features_df %>% filter(tpm_homo_sapiens_nuclear_fraction > 0.001),
             shape = 17, alpha = 0.7) +
  geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs")) +
  stat_cor() +
  scale_x_continuous(expand = c(0,0)) +
  scale_color_manual(values = c("#a8404c", "#424242"), name = "Gene type") + 
  ggtitle("Nuclear Expression vs. promoter binding events") + 
  xlab(expression('Number of TFs')) +
  ylab(expression(log[2](TPM))) 
```

![](01_BIG_knit_files/figure-gfm/binding%20versus%20nuclear%20expression-1.png)<!-- -->

``` r
  # saving figure
  ggsave("figures/nuclear_expression-vs-promoter_binding.pdf")
```

RESULT: (i) looks very similar to total RNA binding versus expression

# Binding versus cytoplasmic expression

Next we will determine the DBP binding versus cytoplasmic expression

``` binding
# Same thing just seeing if there is a difference of cyto RNAs versus DBPs on promoter
ggplot(promoter_features_df, 
            aes(y = log2(tpm_homo_sapiens_cytosolic_fraction + 0.001), x = number_of_dbp, color = gene_type)) + 
  geom_point(data = promoter_features_df %>% filter(tpm_homo_sapiens_cytosolic_fraction > 0.001),
             shape = 17, alpha = 0.7) +
  geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs")) +
  stat_cor() +
  scale_x_continuous(expand = c(0,0)) +
  scale_color_manual(values = c("#a8404c", "#424242"), name = "Gene type") + 
  ggtitle("Cytoplasmic Expression vs. promoter binding events") + 
  xlab(expression('Number of TFs')) +
  ylab(expression(log[2](TPM))) 
  # saving figure
  ggsave("figures/cytoplasmic_expression-vs-promoter_binding.pdf")
  
```

Result: (1) everything seems to be low abundance (2) Some mRNAs are
expressed in the nucleus – we could look at this more later. (3) The
same linear trend holds but is driven by mostly low expression events.

# lncRNA versus mRNA expression in total RNA

Next we will directly test the lncRNA vs mRNA expression levels in total
RNA.

``` r
# plotting
ggplot(promoter_features_df, aes(x = log2(tpm_homo_sapiens_hepg2 + 0.01), color = gene_type))+
  geom_density()
```

![](01_BIG_knit_files/figure-gfm/determining%20lncRNA%20and%20mRNA%20expression%20levels%20in%20total%20RNA-1.png)<!-- -->

``` r
# saving figure
ggsave("figures/mrna_lncrna_tpm_total_rna.pdf")

# let's also do the same for nuclear since lncRNAs are typically more nuclear
ggplot(promoter_features_df, aes(x = log2(tpm_homo_sapiens_nuclear_fraction + 0.01), color = gene_type))+
  geom_density()
```

![](01_BIG_knit_files/figure-gfm/determining%20lncRNA%20and%20mRNA%20expression%20levels%20in%20total%20RNA-2.png)<!-- -->

``` r
# saving figure
ggsave("figures/mrna_lncrna_tpm_nuclear.pdf")
```

Result: (i) This yet again confirms lncRNAs have lower expression levels
than mRNAs. (ii) In the nuclear fraction it shift’s to closer. (ii)
lot’s of mRNA with nuclear expression – that seems odd

We have previously observed that k562 cells also exhibit high binding
promoters that are not expressed. We termed them ‘reservoirs’ as they
are a reservoir for many Dna-protein interaction sites. Based on the
plot above we observed that this phenomena also exists in hepG2 cells as
well.

Based on this we next wanted to identify the ‘reservoirs’ in hepG2
cells.

``` r
# first we will use a cutoff of 100 DBPs.
promoter_features_df$hepg2_reservoir <- 
  as.numeric(promoter_features_df$number_of_dbp > 5 & 
               promoter_features_df$tpm_homo_sapiens_hepg2 < 0.001)

# seeing what we got with table
table(promoter_features_df$hepg2_reservoir)
```

    ## 
    ##     0     1 
    ## 34352  2360

Result: (i) There are 2,360 reservoirs in this largert data set in HEPG2

Now that we have defined reservoirs in hepG2 cells, we next want to
determine how many are similar genomic regions in k562 and hepG2.

``` r
# reading in k562 promoter_features_DF
k562_df <- read_csv("/scratch/Shares/rinnclass/CLASS_2023/data/data/2020_k562_promoter_peak_df.csv")

# next we want to merge the k562 adn Hepg2 DFs 
k562_df <- k562_df %>% 
  dplyr::select(gene_id, reservoir, conservative_reservoir, tpm, expression, tf_binding, promoter_mean_tpm, promoter_median_tpm, promoter_max_tpm) %>%
  dplyr::rename(k562_reservoir = reservoir, 
                k562_conservative_reservoir = conservative_reservoir,
                k562_expression = expression,
                k562_tpm = tpm,
                k562_tf_binding = tf_binding,
                k562_promoter_mean_tpm =  promoter_mean_tpm,
                k562_promoter_median_tpm = promoter_median_tpm,
                k562_promoter_median_tpm = promoter_median_tpm,
                k562_promoter_max_tpm = promoter_max_tpm)

# save this file in new format
write_csv(k562_df,"results/k562_df.csv")

# renaming promoter_features_df to hepg2_df
hepg2_df <- promoter_features_df %>%
  dplyr::select(gene_id, gene_name, tpm_homo_sapiens_hepg2, tpm_homo_sapiens_cytosolic_fraction, tpm_homo_sapiens_nuclear_fraction, tpm_homo_sapiens_insoluble_cytoplasmic_fraction, tpm_homo_sapiens_membrane_fraction, number_of_dbp, hepg2_reservoir) %>%
   dplyr::rename( tpm_total = tpm_homo_sapiens_hepg2,
                 tpm_cytosolic_fraction =  tpm_homo_sapiens_cytosolic_fraction,
                 tpm_nuclear_fraction = tpm_homo_sapiens_nuclear_fraction ,
                 tpm_insoluble_cytoplasmic_fraction = tpm_homo_sapiens_insoluble_cytoplasmic_fraction ,
                 tpm_membrane_fraction = tpm_homo_sapiens_membrane_fraction)

# let's save this handy file
write_csv(hepg2_df,"results/hepg2_df.csv")
  
# Let's merge the k562 reservoirs in with HEPG2_df
# Merges on Gene_id
hepg2_k562_promoter_features_df <- merge(hepg2_df, k562_df)

# Now saving
write_csv(hepg2_k562_promoter_features_df, "results/hepg2_k562_promoter_features_df.csv")

# Make a table of reservoir status
res_status <- hepg2_k562_promoter_features_df %>% 
  group_by(hepg2_reservoir, k562_reservoir, k562_conservative_reservoir) %>%
  summarize(count = n())

# saving for future
write_csv2(res_status, "results/reservoir_overlap_stats.csv")
```

Result: (i) There are 345 reservoirs in both K562 and HEPG2 (ii) There
are 80 Hepg2 reservoirs that overlap “conservative” k562 reservoirs
