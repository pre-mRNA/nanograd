# aim: to produce plots of intron-exon coverage using pileup-data relevant to reference transcriptome
# this analysis uses biomaRt to fetch data about each transcript and analyse splice junction information accordingly

# define global variables
import_filter_threhold <- 20 # the median number of reads needed on the transcript to use it in analysis
threadcount <- 8 # perhaps control this from the shell launcher?

# parse the arguments
args <- commandArgs()
outdir <- args[7]
infile <- args[8]

args <- commandArgs()
outdir <- args[7] # to confirm
infile <- args[8] # to confirm

print(args)

# libs
library(tidyverse)
library(runner)
library(zoo)
library(multidplyr)
library(data.table)
library(scales) # used to add commas to the output
library(ggpubr) # used for qqplot

### part 1. data preprocessing
{
  # import transcriptome data LOCAL
  # import_pileup_data <- read_tsv("../localGadiData/2021-02-22_desplice_filtLib_trimmed.txt.gz", col_names = F, col_types = "cddd") %>%
  #  rename(transcript = 1, coord = 2, wt = 3, mut = 4) %>%
  #  dplyr::select(transcript, coord, wt, mut)

  # import transcriptome data Gadi
    import_pileup_data <- read_tsv(infile, col_names = F, col_types = "cdcd") %>%
      rename(transcript = 1, coord = 2, base = 3, wt = 4) %>%
      dplyr::select(transcript, coord, base, wt)


  imporx_tx_count <- import_pileup_data %>% group_by(transcript) %>% n_groups()
  print(paste("on import, there are", imporx_tx_count, "transcripts identified", sep = " "))

  # filter for abundant genes
  # filtering needs to be careful as not to break the downstream steps

  # we identify genes that have more than 10 counts for both samples as the median
  filter_transcript_targets <- import_pileup_data %>%
    group_by(transcript) %>%
    summarise(transcript=first(transcript), wt_median = median(wt), mut_median = median(mut)) %>%
    filter(wt_median > import_filter_threhold, mut_median > import_filter_threhold)

  # plot median coverage
  import_pileup_data %>%
    group_by(transcript) %>%
    summarise(transcript=first(transcript), wt_median = median(wt), mut_median = median(mut)) %>%
    pivot_longer(cols = c(wt_median, mut_median), names_to = "lib", values_to = "median_coverage") %>%
    ggplot(aes(x=(median_coverage), colour=lib, fill=lib)) + geom_histogram(alpha=0.2, bins=100) + scale_color_brewer(palette="Dark2")+
    scale_fill_brewer(palette="Dark2")+xlim(0,250)+
    facet_grid(lib ~ .)

  # tell user how many transcripts remain
  filt_tx_count <- filter_transcript_targets %>% group_by(transcript) %>% n_groups()
  print(paste("after filtering for transcripts with median coverage of", import_filter_threhold, "there are", filt_tx_count, "transcripts remaining in the analysis", sep = " "))

  # now select the pleup data that correspond to transcripts which meet our filtering criteria
  filtered_pileup_data <- import_pileup_data %>% filter(transcript %in% filter_transcript_targets)
}

### part 2. generate transcriptome annotations for our target genes using biomArt
{

  # make biomaRt query to get transcript information

  # tell ensembl we want information for ensembl mmus build
  ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")

  # fetch attributes for our transcripts of interest
  transcript_attributes <- getBM(attributes=c('ensembl_transcript_id_version', 'transcript_length', 'ensembl_exon_id', 'external_gene_name', 'strand', 'exon_chrom_start', 'exon_chrom_end', 'rank', '5_utr_start', '5_utr_end', '3_utr_start', '3_utr_end'),
                                 filters = 'ensembl_transcript_id_version',
                                 values = filter_transcript_targets$transcript,
                                 mart = ensembl)

  # from the biomart-mined annotation, calculate exon & UTR lengths and split annotation by strand
  # we split by strand because exon "ranks" (refer to ensembl for full explaination) need to be reversed for exons on the negative strand
  processed_transcript_attributes <- transcript_attributes %>%
    mutate(exon_length = (abs(exon_chrom_start - exon_chrom_end) + 1), five_utr_len = (abs(`5_utr_start` - `5_utr_end`) + 1), three_utr_len = (abs(`3_utr_start` - `3_utr_end`) + 1)) %>%
    group_by(ensembl_transcript_id_version) %>%
    summarise(transcript = first(ensembl_transcript_id_version), transcript_length = mean(transcript_length), num_exons=n(), cum_exon_length = sum(exon_length), strand=first(strand)) %>%
    mutate(calculated_transcript_length = cum_exon_length, calculated_diff = calculated_transcript_length - transcript_length) %>%
    group_by(strand) %>%
    group_split
  # splitting yields a list of tibbles
  # now we strictly define these tibbles in list as individual tibbles
  plus_strand_annotation <- as_tibble(processed_transcript_attributes[[2]])
  minus_strand_annotation <- as_tibble(processed_transcript_attributes[[1]])

  # calculate exon positions for plus strand
  # generates a large list where each tibble contains the exons in an ensembl transcript
  plus_strand_exons <- plus_strand_annotation %>%
    group_by(transcript) %>%
    group_split()

  # calculate exon positions for minus strand
  # like above code
  minus_strand_exons <- minus_strand_annotation %>%
    group_by(transcript) %>%
    group_split()

  # calculate the coordinate of each exon relative to the transcript FASTA, and classify exons by position in gene
  # add information about exon types; 0 = single exon, 1 = first exon, -1 = last exon, >1 = nth exon)
  plus_final_annotation <- NULL

  for (p in plus_strand_exons) {
    transcript_start <- min(p$exon_chrom_start)
    p <- p %>%
      arrange(rank) %>%
      mutate(
        exon_end = cumsum(exon_length),
        exon_start = cumsum(exon_length) - exon_length + 1,
      )
  }
  # add information about exon types; 0 = single exon, 1 = first exon, -1 = last exon, >1 = nth exon)
  if (dim(p)[1] == 1) {
    p <- p %>% mutate(exon_type = 0)
  } else if (dim(p)[1] == 2) {
    p <- p %>% mutate(exon_type = c(1,-1))
  } else {
    p <- p %>% mutate(exon_type = c(1,seq(from = 2, to = ((dim(p)[1]-1))), -1))
  }
  plus_final_annotation <- merge(plus_final_annotation, p, all = TRUE)
}

plus_internal_ss <- plus_final_annotation %>%
  filter(exon_type > 1) %>%
  mutate(ss_5_minus_40 = exon_end - 39,
         ss_5 = exon_end,
         ss_3 = exon_start,
         ss_3_plus_40 = exon_start + 40)

##### select the key regions for internal 5'ss
input_internal_5ss <- plus_internal_ss %>% dplyr::select(transcript, ss_5)

##### select the key regions for internal 3'ss
input_internal_3ss <- plus_internal_ss %>% dplyr::select(transcript, ss_3)

# write a function which returns a vector of diff from ss_5 - 39 to 33_5
plus_3_diff <- NULL
plus_3_vector <- function(xa) {
  tryCatch({
    #print(xa[1])
    #print(xa[2])
    #print(paste("input is", xa, "which separately is", xa[1], xa[2], sep = " "))
    #a <- NULL
    #b <- NULL
    a <- normalised_output %>% mutate(absdiff = abs(diff)) %>% filter(absdiff > 0.000000000001) %>% filter(y>20) %>%
      filter(chr == as.character(xa[1])) %>%
      arrange(coord)
    #print(a)
    b <- a[which(a$coord == as.character(xa[2])) + c(-74:100), ]
    #print(a[which(a$coord == as.character(xa[,2])) + c(-39:0), ])
    rownames(b) <- NULL
    if (dim(b)[1] < 20) {
      #print(paste(xa[1], " is empty", ",or equal to", dim(b)[1], sep = ""))
    } else {
      #print("proceeding")
      ca <<- b %>% dplyr::select(diff)
      plus_3_diff <<- bind_cols(plus_3_diff, ca)
    }
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


# test the function
plus_3_vector(input_internal_3ss[1,])

plus_3_diff <- NULL
for (i in seq(dim(input_internal_3ss)[1])) {
  plus_3_vector(input_internal_3ss[i,])
}

plus_3_diff <- plus_3_diff %>% mutate(ave = average)
plus_3_diff$mean <- rowMeans(plus_3_diff, na.rm=TRUE)
plot(plus_3_diff$mean)

}






# now, split the pileup data by strand
# plus_pileup <- import_data %>% filter(chr %in% plus_strand_annotation$ensembl_transcript_id_version)
# minus_pileup <- import_data %>% filter(chr %in% minus_strand_annotation$ensembl_transcript_id_version)

# verify that all positions are accounted for
#dim(plus_pileup)[1] + dim(minus_pileup)[1] - dim(import_data)[1]

# make a funcion to visualize coverage over specific genes
plotCov <- function(targGene) {
  import_data %>% filter(chr == targGene) %>% pivot_longer(cols=c("wt", "mut"), values_to = "count") %>% ggplot(aes(x=coord, y=count, col=count)) + geom_point(aes(shape = name)) + theme_minimal()
}

### part 4. define the 3' feather; that is the number of bases required before coverage at the 3' of the read reaches a threshold compared to the maximum

# define the coverage threhold where the feather initiates
threshold=0.8

plus_summary <-
  filtered_pileup_data %>%
  mutate(index = row_number()) %>%
  group_by(chr) %>%
  summarise(chr = first(chr),
            strand="plus",
            max_cov_wt = max(wt),
            max_cov_mut = max(mut),
            length = n(),
            wt_last_coord_over_thresh = last(coord[wt > (max_cov_wt * threshold)]),
            mut_last_coord_over_thresh = last(coord[mut > (max_cov_mut * threshold)]),
            wt_feather = length((index[coord > wt_last_coord_over_thresh])),
            mut_feather = length((index[coord > mut_last_coord_over_thresh]))
  )

minus_summary <- minus_pileup %>%
  mutate(index = row_number()) %>%
  group_by(chr) %>%
  summarise(chr = first(chr),
            strand="minus",
            max_cov_wt = max(wt),
            max_cov_mut = max(mut),
            length = n(),
            wt_last_coord_over_thresh = last(coord[wt > (max_cov_wt * threshold)]),
            mut_last_coord_over_thresh = last(coord[mut > (max_cov_mut * threshold)]),
            wt_feather = length((index[coord > wt_last_coord_over_thresh])),
            mut_feather = length((index[coord > mut_last_coord_over_thresh]))
  )

# merge the feathering data
merged_feather_data <- add_row(plus_summary, minus_summary) %>%
  dplyr::select(chr, strand, length, wt_feather, mut_feather) %>%
  pivot_longer(c("wt_feather", "mut_feather"), names_to = "feather")

# plot a histogram of the feathering data
ggplot(merged_feather_data, aes(x=value, color=feather)) +
  geom_histogram(fill="white", alpha=0.5, position="identity", bins=60)

### part 5. calculate the gradients

# split the data by transcript
import_by_transcript <- import_data %>% group_by(chr) %>% group_split()

# make an empty data table to which we will write the results
gradient_data <- data.table(chr=numeric(), coord=numeric(), x=numeric(), y=numeric(), exon=numeric(), xgrad=numeric(), ygrad=numeric())

# define a function which takes the interval gradient of a dataset
rollingSlope.lm.fit <- function(vector) {
  a <- coef(.lm.fit(cbind(1, seq(vector)), vector))[2]
  return(a)
}

# define width, either independantly or within a loop
#for (a in seq(6,30,2)) {
a <- 12
width <- a

gradient_data <- NULL
exon_gradients <- NULL
# loop over each "exon data table" in "add_index" and calculate the gradient for y (modified) and x (wild-type) datasets
# the output is appebded to export_data
for (p in add_index) {
  p <- as.data.table(p) %>%
    rename(y = mut, x = wt)

  exon_gradients <- (head(p, -16))
  exon_gradients[, ':=' (xgrad = rollapply(x, width=a, FUN=rollingSlope.lm.fit, fill=NA))]
  exon_gradients[, ':=' (ygrad = rollapply(y, width=a, FUN=rollingSlope.lm.fit, fill=NA))]
  gradient_data <- merge(gradient_data, exon_gradients, all = TRUE)
}


# normalise the new gradients against the depth at that position (adjust for poor coverage in mutant)
normalised_output <- gradient_data %>%
  mutate(ygrad = ygrad / y, xgrad = xgrad / x) %>%
  mutate(diff = ygrad - xgrad)

# plot a histogram of gradient "differences" befire masking runs
ggplot(normalised_output %>% filter(y>25), aes(x=diff)) +
  geom_histogram(binwidth=0.00001, color="black", fill="white") +
  labs(title="Coverage gradient difference, crosslink vs WT cells",x="coverage difference (normalised Δreads/Δnt, crosslink - wt)", y = "number of filtered genome 12mers")+
  theme_classic() +
  xlim(-0.007, 0.007)

# remove rows where diff = 0
normalised_output %>% mutate(absdiff = abs(diff)) %>% filter(absdiff > 0.000000000001) %>% filter(y>20) %>% ggplot(aes(x=diff)) +
  geom_histogram(binwidth=0.0002, color="black", fill="white") +
  labs(title="Coverage gradient difference, crosslink vs WT cells",x="coverage difference (normalised Δreads/Δnt, crosslink - wt)", y = "number of filtered genome 12mers")+
  theme_classic()+
  xlim(c(-0.02,0.02))

# check for repeats
repeat_vector <- as_tibble(c(rle(normalised_output$diff)[1], rle(normalised_output$diff)[2])) %>% uncount(values)

### part 5. gather splice junction coverage data based on annotations

# filter for genes which have coverage
high_coverage_output <- normalised_output %>% mutate(absdiff = abs(diff)) %>% filter(absdiff > 0.000000000001) %>% filter(y>20)
high_coverage_genes <- high_coverage_output %>% group_by(chr) %>% summarize(chr = first(chr), coord = mean(coord, na.rm=TRUE), diff = mean(diff), coverage=median(y))

plotOut <- function(targGene) {
  normalised_output %>% filter(chr == targGene) %>% pivot_longer(cols=c("x", "y"), values_to = "count") %>% ggplot(aes(x=coord, y=count, col=count)) + geom_point(aes(shape = name)) + theme_minimal()
}

plotOut2 <- function(targGene) {
  normalised_output %>% filter(chr == targGene) %>% pivot_longer(cols=c("xgrad", "ygrad"), values_to = "count") %>% ggplot(aes(x=coord, y=count, col=name)) + geom_point(aes(shape = name)) + theme_minimal()
}


### working
plus_final_annotation <- NULL
for (p in plus_strand_exons) {
  transcript_start <- min(p$exon_chrom_start)
  p <- p %>%
    arrange(rank) %>%
    mutate(
      exon_start = exon_chrom_start - transcript_start + 1,
      exon_end = exon_chrom_end - transcript_start + 1)

  # add information about exon types; 0 = single exon, 1 = first exon, -1 = last exon, >1 = nth exon)
  if (dim(p)[1] == 1) {
    p <- p %>% mutate(exon_type = 0)
  } else if (dim(p)[1] == 2) {
    p <- p %>% mutate(exon_type = c(1,-1))
  } else {
    p <- p %>% mutate(exon_type = c(1,seq(from = 2, to = ((dim(p)[1]-1))), -1))
  }
  plus_final_annotation <- merge(plus_final_annotation, p, all = TRUE)
}


# write a function which returns a vector of diff from ss_5 - 39 to 33_5
plus_5_diff <- NULL
plus_5_vector <- function(xa) {
  tryCatch({
    #print(xa[1])
    #print(xa[2])
    #print(paste("input is", xa, "which separately is", xa[1], xa[2], sep = " "))
    #a <- NULL
    #b <- NULL
    a <- normalised_output %>% mutate(absdiff = abs(diff)) %>% filter(absdiff > 0.000000000001) %>% filter(y>20) %>%
      filter(chr == as.character(xa[1])) %>%
      arrange(coord)
    #print(a)
    b <- a[which(a$coord == as.character(xa[2])) + c(-74:100), ]
    #print(a[which(a$coord == as.character(xa[,2])) + c(-39:0), ])
    rownames(b) <- NULL
    if (dim(b)[1] < 20) {
      #print(paste(xa[1], " is empty", ",or equal to", dim(b)[1], sep = ""))
    } else {
      #print("proceeding")
      ca <<- b %>% dplyr::select(diff)
      plus_5_diff <<- bind_cols(plus_5_diff, ca)
    }
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


# test the function
plus_5_vector(input_internal_5ss[1,])

plus_5_diff <- NULL
for (i in seq(dim(input_internal_5ss)[1])) {
  plus_5_vector(input_internal_5ss[i,])
}

plus_5_diff$mean <- rowMeans(plus_5_diff, na.rm=TRUE)
plot(plus_5_diff$mean)

df_for_line <- plus_5_diff %>% rowid_to_column("position") %>% pivot_longer(cols = starts_with("diff"), values_to = "differ", names_to = "transcript") %>% drop_na() %>% group_by(position) %>% summarise(differ = mean(differ));
plus_5_diff %>%
  dplyr::select( -mean) %>%
  rowid_to_column("position") %>%
  pivot_longer(cols = starts_with("diff"), values_to = "differ", names_to = "transcript") %>%
  ggplot(mapping = aes(x = position, y = differ, col = "red")) +
  geom_point(alpha = 0.06) +
  geom_path(data = df_for_line, aes(x = position, y = differ, group = 1)) +
  scale_y_continuous(limits = c(0.0003,0.0015))+
  ylab("gradient_difference")+
  geom_vline(xintercept=75, color="grey", size=1) +
  geom_hline(yintercept=0.0006, color="black", size=0.25)
#stat_summary(fun.y=mean, geom="line", aes(group=1))  +
#stat_summary(fun.y=mean, geom="point") +


plot(plus_5_diff$mean)



df_for_line_plus_3 <- plus_3_diff %>% rowid_to_column("position") %>% pivot_longer(cols = starts_with("diff"), values_to = "differ", names_to = "transcript") %>% drop_na() %>% group_by(position) %>% summarise(differ = mean(differ));
plus_3_diff %>%
  dplyr::select( -mean) %>%
  rowid_to_column("position") %>%
  pivot_longer(cols = starts_with("diff"), values_to = "differ", names_to = "transcript") %>%
  ggplot(mapping = aes(x = position, y = differ, col = "red")) +
  geom_point(alpha = 0.06) +
  geom_path(data = df_for_line_plus_3, aes(x = position, y = differ, group = 1)) +
  scale_y_continuous(limits = c(0.0003,0.0015))+
  ylab("gradient_difference")+
  geom_vline(xintercept=75, color="grey", size=1) +
  geom_hline(yintercept=0.0006, color="black", size=0.25)
#stat_summary(fun.y=mean, geom="line", aes(group=1))  +
#stat_summary(fun.y=mean, geom="point") +


### testing for minus strand

minus_final_annotation <- NULL
for (p in minus_strand_exons) {
  transcript_start <- min(p$exon_chrom_end)
  p <- p %>%
    arrange(desc(rank)) %>%
    mutate(
      exon_end = cumsum(exon_length),
      exon_start = cumsum(exon_length) - exon_length + 1,
    )

  # add information about exon types; 0 = single exon, 1 = first exon, -1 = last exon, >1 = nth exon)
  if (dim(p)[1] == 1) {
    p <- p %>% mutate(exon_type = 0)
  } else if (dim(p)[1] == 2) {
    p <- p %>% mutate(exon_type = c(1,-1))
  } else {
    p <- p %>% mutate(exon_type = c(1,seq(from = 2, to = ((dim(p)[1]-1))), -1))
  }
  minus_final_annotation <- merge(minus_final_annotation, p, all = TRUE)
}

minus_internal_ss <- minus_final_annotation %>%
  filter(exon_type > 1) %>%
  mutate(ss_5_minus_40 = exon_end - 39,
         ss_5 = exon_end,
         ss_3 = exon_start,
         ss_3_plus_40 = exon_start + 40)

##### select the key regions for internal 5'ss
minus_input_internal_5ss <- minus_internal_ss %>% dplyr::select(transcript, ss_5)

# write a function which returns a vector of diff from ss_5 - 39 to 33_5
minus_5_diff <- NULL
minus_5_vector <- function(xa) {
  tryCatch({
    #print(xa[1])
    #print(xa[2])
    #print(paste("input is", xa, "which separately is", xa[1], xa[2], sep = " "))
    #a <- NULL
    #b <- NULL
    a <- normalised_output %>% mutate(absdiff = abs(diff)) %>% filter(absdiff > 0.000000000001) %>% filter(y>20) %>%
      filter(chr == as.character(xa[1])) %>%
      arrange(coord)
    #print(a)
    b <- a[which(a$coord == as.character(xa[2])) + c(-74:50), ]
    #print(a[which(a$coord == as.character(xa[,2])) + c(-39:0), ])
    rownames(b) <- NULL
    if (dim(b)[1] < 20) {
      #print(paste(xa[1], " is empty", ",or equal to", dim(b)[1], sep = ""))
    } else {
      #print("proceeding")
      ca <<- b %>% dplyr::select(diff)
      minus_5_diff <<- bind_cols(minus_5_diff, ca)
    }
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

minus_5_diff <- NULL
for (i in seq(dim(minus_input_internal_5ss)[1])) {
  minus_5_vector(minus_input_internal_5ss[i,])
}

# to do a quick check, add the mean as a row and plot it
minus_5_diff$mean <- rowMeans(minus_5_diff, na.rm=TRUE)
plot(minus_5_diff$mean)
# remove mean as not required for further analysis
minus_5_diff <- minus_5_diff %>% dplyr::select( -mean)

# write the means to a new df
minus_5_df_for_line <- minus_5_diff %>% rowid_to_column("position") %>% pivot_longer(cols = starts_with("diff"), values_to = "differ", names_to = "transcript") %>% drop_na() %>% group_by(position) %>% summarise(differ = mean(differ));
# make position diff plot with means
minus_5_diff %>%
  rowid_to_column("position") %>%
  pivot_longer(cols = starts_with("diff"), values_to = "differ", names_to = "transcript") %>%
  ggplot(mapping = aes(x = position, y = differ, col = "red")) +
  geom_point(alpha = 0.06) +
  geom_path(data = minus_5_df_for_line, aes(x = position, y = differ, group = 1)) +
  scale_y_continuous(limits = c(0.0003,0.0015))+
  ylab("gradient_difference")+
  geom_vline(xintercept=75, color="grey", size=1) +
  geom_hline(yintercept=0.0006, color="black", size=0.25)

#stat_summary(fun.y=mean, geom="line", aes(group=1))  +
#stat_summary(fun.y=mean, geom="point") +

##### select the key regions for internal 5'ss
minus_input_internal_3ss <- minus_internal_ss %>% dplyr::select(transcript, ss_3)

# write a function which returns a vector of diff from ss_5 - 39 to 33_5

minus_3_vector <- function(xa) {
  tryCatch({
    #print(xa[1])
    #print(xa[2])
    #print(paste("input is", xa, "which separately is", xa[1], xa[2], sep = " "))
    #a <- NULL
    #b <- NULL
    a <- normalised_output %>% mutate(absdiff = abs(diff)) %>% filter(absdiff > 0.000000000001) %>% filter(y>20) %>%
      filter(chr == as.character(xa[1])) %>%
      arrange(coord)
    #print(a)
    b <- a[which(a$coord == as.character(xa[2])) + c(-74:50), ]
    #print(a[which(a$coord == as.character(xa[,2])) + c(-39:0), ])
    rownames(b) <- NULL

    if (dim(b)[1] < 20 || dim(b)[1] == 124) {
      print(paste(xa[1], "is empty", ",o'r equal twith an actual value of", dim(b)[1], sep = " "))
    } else {
      #print("proceeding")
      ca <<- b %>% dplyr::select(diff)
      minus_3_diff <<- bind_cols(minus_3_diff, ca)
      #print(dim(ca))
    }
    #return(dim(ca))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

#minus_3_vector(minus_input_internal_3ss[2,])

minus_3_diff <- NULL
for (i in seq(dim(minus_input_internal_3ss)[1])) {
  minus_3_vector(minus_input_internal_3ss[i,])
}

# to do a quick check, add the mean as a row and plot it
minus_3_diff$mean <- rowMeans(minus_3_diff, na.rm=TRUE)
plot(minus_3_diff$mean)
# remove mean as not required for further analysis
minus_3_diff <- minus_3_diff %>% dplyr::select( -mean)

# write the means to a new df
minus_3_df_for_line <- minus_3_diff %>% rowid_to_column("position") %>% pivot_longer(cols = starts_with("diff"), values_to = "differ", names_to = "transcript") %>% drop_na() %>% group_by(position) %>% summarise(differ = mean(differ));
# make position diff plot with means
minus_3_diff %>%
  rowid_to_column("position") %>%
  pivot_longer(cols = starts_with("diff"), values_to = "differ", names_to = "transcript") %>%
  ggplot(mapping = aes(x = position, y = differ, col = "red")) +
  geom_point(alpha = 0.06) +
  geom_path(data = minus_3_df_for_line, aes(x = position, y = differ, group = 1)) +
  scale_y_continuous(limits = c(0.0003,0.0015))+
  ylab("gradient_difference")+
  geom_vline(xintercept=75, color="grey", size=1) +
  geom_hline(yintercept=0.0006, color="black", size=0.25)


##### select the key regions for internal 3'ss
input_internal_3ss <- plus_internal_ss %>% dplyr::select(transcript, ss_3)

# write a function which returns a vector of diff from ss_5 - 39 to 33_5
plus_3_diff <- NULL
plus_3_vector <- function(xa) {
  tryCatch({
    #print(xa[1])
    #print(xa[2])
    #print(paste("input is", xa, "which separately is", xa[1], xa[2], sep = " "))
    #a <- NULL
    #b <- NULL
    a <- normalised_output %>% mutate(absdiff = abs(diff)) %>% filter(absdiff > 0.000000000001) %>% filter(y>20) %>%
      filter(chr == as.character(xa[1])) %>%
      arrange(coord)
    #print(a)
    b <- a[which(a$coord == as.character(xa[2])) + c(-50:50), ]
    #print(a[which(a$coord == as.character(xa[,2])) + c(-39:0), ])
    rownames(b) <- NULL
    if (dim(b)[1] < 20) {
      #print(paste(xa[1], " is empty", ",or equal to", dim(b)[1], sep = ""))
    } else {
      #print("proceeding")
      ca <<- b %>% dplyr::select(diff)
      plus_3_diff <<- bind_cols(plus_3_diff, ca)
    }
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


# test the function
plus_3_vector(input_internal_3ss[1,])

plus_3_diff <- NULL
for (i in seq(dim(input_internal_3ss)[1])) {
  plus_3_vector(input_internal_3ss[i,])
}

plus_3_diff <- plus_3_diff %>% mutate(ave = average)
plus_3_diff$mean <- rowMeans(plus_3_diff, na.rm=TRUE)
plot(plus_3_diff$mean)



dim(bind_rows(plus_5_diff, plus_3_diff))
plot(a)



data1 <- as.data.frame(input_internal_5ss)
xy.list <- split(data1, seq(nrow(data1)))
data_list[["1"]]



sapply(xy.list, plus_5_vector)


# calculate exon positions for minus strand
minus_strand_exons <- minus_strand_annotation %>%
  group_by(transcript) %>%
  group_split()

minus_final_annotation <- NULL
for (p in plus_strand_exons)


  exon_gradients <- (head(p, -16))
exon_gradients <- (tail(exon_gradients, -16))
exon_gradients[, ':=' (xgrad = rollapply(x, width=a, FUN=rollingSlope.lm.fit, fill=NA))]
exon_gradients[, ':=' (ygrad = rollapply(y, width=a, FUN=rollingSlope.lm.fit, fill=NA))]
gradient_data <- merge(gradient_data, exon_gradients, all = TRUE)
}



group_modify(~ {
  .x %>%
    if_else(strand == "1",
            start="arrange(rank)"

            , false, missing = NULL)

  start <- arrange(rank) %>% select(exon_chrom_start)
  mutate(nms = c("min", "Q1", "median", "Q3", "max"))
})
