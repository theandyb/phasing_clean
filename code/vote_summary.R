library(tidyverse)
library(yaml)

source("code/common_functions.R")

config_obj <- yaml::read_yaml("_config.yaml")

vote_switch_dir <- paste0(config_obj$base_dir, "/output/vote_", "15", "/switches/annotated/")
het_loc_dir <- paste0(config_obj$base_dir,"/output/trio_phase_15/het_loc/")
#het_loc_dir <- paste0(config_obj$base_dir,"/output/switch_errors/het_loc/")


gc_content_1kb <- read_tsv("data/chr15_gc1kb_pilot.bed")
colnames(gc_content_1kb) <- c("CHR", "START", "END", "AT", "GC", "A", "C", "G", "T", "TOTAL", "OTHER", "LENGTH")
gc_content_1kb  <- gc_content_1kb %>%
  mutate(bin_id = (START / 1000) + 1)

# sample_info_df <- read_csv("data/1kgp/subject_info.csv") %>%
#   select(SAMPLE_NAME, POPULATION, SUPER)
#
# pair_info_df <- read_delim("data/sample_pairs.csv", col_names = c("POP", "ID1", "ID2"))
# pair_info_df <- left_join(pair_info_df, sample_info_df, by = c("ID1"="SAMPLE_NAME")) %>%
#   rename(SP = SUPER) %>%
#   select(-POPULATION)

df_subj <- read_csv("data/1kgp/subject_info.csv")
ped_df <- read_table("data/1kgp/1kGP.3202_samples.pedigree_info.txt")
child_ids <- ped_df %>% filter(fatherID != "0" & motherID != "0") %>% pull(sampleID)

df_subj_rel <- df_subj %>%
  filter(SAMPLE_NAME %in% child_ids)

df_subj_rel$id2 <- 1:602

vote_summary <- function(pair_id, vote_dir, gc_content_1kb, het_loc_dir, bin_size = 1000){
  switch_err_vote <- read_csv(paste0(vote_dir, "switch_", pair_id, ".csv"), show_col_types = FALSE) %>%
    mutate(bin_id = ceiling(pos_end / bin_size))

  het_df <- read_tsv(paste0(paste0(het_loc_dir, "pair_", pair_id, "_het_loc.txt")),
                     col_names = c("chr", "pos_start", "gt"), show_col_types = FALSE) %>%
    select(chr, pos_start) %>%
    mutate(pos_end = lead(pos_start)) %>%
    drop_na()

  switch_err_vote <- switch_err_vote %>%
    left_join({gc_content_1kb %>% select(bin_id, GC)}, by = "bin_id")

  # get positions of flips
  flip_pos_vote <- get_flip_pos(switch_err_vote)

  # Assign switches flip status
  switch_err_vote$is_flip <- (switch_err_vote$pos_start %in% flip_pos_vote) |
    (switch_err_vote$pos_end %in% flip_pos_vote)


  # Start of flip
  switch_err_vote$start_flip <- (switch_err_vote$pos_end %in% flip_pos_vote)


  # Append switch/flip locations to het_df
  het_df <- het_df %>%
    mutate(vote_switch = pos_start %in% {
      switch_err_vote %>%
        filter(!is_flip) %>%
        pull(pos_start)},
      vote_flip = pos_start %in% {
        switch_err_vote %>%
          filter(start_flip) %>%
          pull(pos_start)
      })
  het_df$id <- 1:length(het_df$chr)

  # stats we want to pull
  n_switch_vote <- length(switch_err_vote$pos_start)
  n_flip_vote <- sum(switch_err_vote$is_flip) / 2
  n_other_vote <- n_switch_vote - 2*n_flip_vote


  n_switch_cpg_vote <- switch_err_vote %>%
    filter(cpg_end == 1) %>%
    filter(is_flip != TRUE | start_flip == TRUE) %>%
    pull(is_flip) %>%
    length()


  n_flip_cpg_vote <- switch_err_vote %>%
    filter(cpg_end == 1) %>%
    pull(start_flip) %>%
    sum()

  n_other_cpg_vote <- switch_err_vote %>%
    filter(is_flip != TRUE) %>%
    pull(cpg_end) %>%
    sum()


  mean_gc_switch_vote <- mean(switch_err_vote$GC)

  # distance metrics
  median_dist_vote <- switch_err_vote %>%
    filter(is_flip == F | start_flip == T) %>%
    pull(pos_end) %>%
    diff(lag = 1) %>%
    median()

  mean_dist_vote <- switch_err_vote %>%
    filter(is_flip == F | start_flip == T) %>%
    pull(pos_end) %>%
    diff(lag = 1) %>%
    mean()

  # number hets between switches
  med_hets_switch_vote <- het_df %>%
    filter(vote_switch) %>%
    mutate(next_id = lead(id)) %>%
    mutate(hets_btw = next_id - id - 1) %>%
    pull(hets_btw) %>%
    median(na.rm = T)

  # number hets between flips
  med_hets_flip_vote <- het_df %>%
    filter(vote_flip) %>%
    mutate(next_id = lead(id)) %>%
    mutate(hets_btw = next_id - id - 1) %>%
    pull(hets_btw) %>%
    median(na.rm = T)

  # mean number of hets between switches and flips
  mean_hets_switch_vote <- het_df %>%
    filter(vote_switch) %>%
    mutate(next_id = lead(id)) %>%
    mutate(hets_btw = next_id - id - 1) %>%
    pull(hets_btw) %>%
    mean(na.rm = T)

  # number hets between flips
  mean_hets_flip_vote <- het_df %>%
    filter(vote_flip) %>%
    mutate(next_id = lead(id)) %>%
    mutate(hets_btw = next_id - id - 1) %>%
    pull(hets_btw) %>%
    mean(na.rm = T)

  return(data.frame(pair_id = pair_id,
                    n_switch_vote = n_switch_vote,
                    n_flip_vote = n_flip_vote,
                    n_other_vote = n_other_vote,
                    n_switch_cpg_vote = n_switch_cpg_vote,
                    n_flip_cpg_vote = n_flip_cpg_vote,
                    n_other_cpg_vote = n_other_cpg_vote,
                    mean_gc_switch_vote = mean_gc_switch_vote,
                    median_dist_vote = median_dist_vote,
                    mean_dist_vote = mean_dist_vote,
                    med_hets_switch_vote = med_hets_switch_vote,
                    med_hets_flip_vote = med_hets_flip_vote,
                    mean_hets_switch_vote = mean_hets_switch_vote,
                    mean_hets_flip_vote = mean_hets_flip_vote,
                    n_hets = length(het_df$pos_start)
  ))
}

df_vcftools <- lapply(c(1:602),
                      function(x){
                        vote_summary(x, vote_switch_dir, gc_content_1kb, het_loc_dir)
                      }) %>%
  bind_rows()

df_vcftools$id2 <- 1:602
df_vcftools <- left_join(df_vcftools, df_subj_rel, by="id2")

write_csv(df_vcftools, paste0(config_obj$base_dir, "/output/vote_", "15", "/summary.csv"))
# df_vcftools$pop <- pair_info_df$SP
# df_vcftools$subpop <- pair_info_df$POP
#
# write_csv(df_vcftools, paste0(config_obj$base_dir, "/output/vote_x/summary.csv"))
