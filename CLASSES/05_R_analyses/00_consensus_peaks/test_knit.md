01\_test\_knit
================
JR
3/17/2022

Acompanying youtube video: <https://youtu.be/BPHQAMAh5GE>

Goal: to create a consensus\_peak file for each dbp, export and format
for UCSC. The first step in analyses is looking at raw data :) Here we
will take our functions to create a consensus peak file and see how it
lines up with raw data on the UCSC genome browser. We will do this in
these steps:

1.  Import peaks
2.  create consensus peaks
3.  export named GRanges list to .bed files
4.  format the file to be uploaded into UCSC
5.  inspect raw data to see if consensus peaks make sense.

# STEP1: import peaks

``` r
# import peaks
peak_list <- import_peaks(consensus_file_path = broadpeakfilepath)

# let's get a list of how many peaks are in each file before we create consensus peaks.
peak_num <- sapply(peak_list, length) %>% 
  as.data.frame(row.names = T)
```

    ## Warning in as.data.frame.integer(., row.names = T): 'row.names' is not a
    ## character vector of length 27 -- omitting it. Will be an error!

``` r
# label column
names(peak_num) <- c("num_peaks")

# make dbp name a col.

peak_num <- peak_num %>%
  rownames_to_column(var = "dbp") %>%
  
  # This is pretty similar to text to columns in excel
  separate(col = dbp,  into = c('dbp', 'replicate'), sep = "_")
  # peak_num <- separate(peak_num, col = dbp,  into = c('dbp', 'replicate'), sep = "_")
peak_num
```

    ##         dbp replicate num_peaks
    ## 1     EP300        R1     62260
    ## 2     EP300        R2     66416
    ## 3     EP300        R3     14319
    ## 4     EP300        R4     16872
    ## 5     EP300        R5     82743
    ## 6     EP300        R6     79479
    ## 7  H3K27me3        R2     78041
    ## 8  H3K27me3        R4    143510
    ## 9  H3K36me3        R1    173838
    ## 10 H3K36me3        R2     99485
    ## 11 H3K36me3        R3      4803
    ## 12 H3K36me3        R4    147910
    ## 13  H3K4me3        R1     30903
    ## 14  H3K4me3        R2     27563
    ## 15  H3K4me3        R3     34650
    ## 16  H3K4me3        R4     40305
    ## 17  H3K4me3        R5     14863
    ## 18  H3K4me3        R6     32403
    ## 19    HDAC2        R1     44226
    ## 20    HDAC2        R2     96342
    ## 21    HDAC2        R3     52310
    ## 22    HDAC2        R4     97967
    ## 23    RCOR1        R1     11781
    ## 24    RCOR1        R2     35131
    ## 25     REST        R1     39702
    ## 26     REST        R2     38225
    ## 27     REST        R4     23298
