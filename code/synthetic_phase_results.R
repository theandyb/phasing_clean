library(tidyverse)
library(yaml)

source("code/common_functions.R")

config_obj <- yaml::read_yaml("_config.yaml")

eagle_switch_dir <-   paste0(config_obj$base_dir,"/output/switch_errors/switch_errors/eagle/annotated/")
shapeit_switch_dir <- paste0(config_obj$base_dir,"/output/switch_errors/switch_errors/shapeit/annotated/")
beagle_switch_dir <-  paste0(config_obj$base_dir,"/output/switch_errors/switch_errors/beagle/annotated/")
num_sites_dir <-      paste0(config_obj$base_dir,"/output/switch_errors/vcf_n_sites/")
het_loc_dir <-        paste0(config_obj$base_dir,"/output/switch_errors/het_loc/annotated/")

sample_info_df <- read_csv("data/1kgp/subject_info.csv") %>%
  select(SAMPLE_NAME, POPULATION, SUPER)

pair_info_df <- read_delim("data/sample_pairs.csv", col_names = c("POP", "ID1", "ID2"))
pair_info_df <- left_join(pair_info_df, sample_info_df, by = c("ID1"="SAMPLE_NAME")) %>%
  rename(SP = SUPER) %>%
  select(-POPULATION)

n_het_cpg_df <- read_tsv(paste0(config_obj$base_dir,"/output/switch_errors/het_cpg_count.tsv"),
                         col_names = c("pair_id", "n_het_cpg"))

gc_content_1kb <- read_tsv("data/chrX_gc1kb_pilot.bed")
colnames(gc_content_1kb) <- c("CHR", "START", "END", "AT", "GC", "A", "C", "G", "T", "TOTAL", "OTHER", "LENGTH")
gc_content_1kb  <- gc_content_1kb %>%
  mutate(bin_id = (START / 1000) + 1)

df_vcftools <- lapply(c(1:1000),
                      function(x){
                        switch_summary(x, eagle_switch_dir, beagle_switch_dir, shapeit_switch_dir, gc_content_1kb, het_loc_dir)
                      }) %>%
  bind_rows()

df_vcftools$pop <- pair_info_df$SP
df_vcftools$subpop <- pair_info_df$POP

par_size <- (2781479 - 10001) + (156030895 - 155701383)
pd_size <- read_tsv("/net/snowwhite/home/beckandy/research/phasing_clean/output/switch_errors/het_loc/phased_size.txt",
                    col_names = c("id", "size"), show_col_types = FALSE) %>%
  mutate(size_mb = (size - par_size) / 1e6)

df_vcftools$size_mb <- pd_size$size_mb
df_vcftools <- left_join(df_vcftools, n_het_cpg_df)

write_csv(df_vcftools, "output/switch_errors/switch_errors/summary.csv")

# NEW CODE
# annotated het position files for all subjects

get_flip_pos <- function(df){
  flip_list <- c()
  in_progress <- FALSE
  for(i in 1:(length(df$pos_start)-1)){
    if(df$pos_end[i] == df$pos_start[i+1]){
      if(!in_progress){
        in_progress <- TRUE
        flip_list <- c(flip_list, df$pos_end[i])
      } else{
        in_progress <- FALSE
      }
    } else{
      in_progress <- FALSE
    }
  }
  return(flip_list)
}

# function for reading in switch error file for 1 sd for 1 method
# pos_start pos_end cpg_start cpg_end motif_start end_motif
read_switch_err <- function(switch_dir, id){
  f_name <- paste0(switch_dir, "switch_", id, ".csv")
  return(read_csv(f_name, show_col_types = F))
}

# get list of switch errors and flip errors
switch_err_positions <- function(switch_dir, id){
  switch_df <- read_switch_err(switch_dir, id)

  flip_pos <- get_flip_pos(switch_df)

  switch_df$flip_start <- switch_df$pos_end %in% flip_pos
  switch_df$flip_end <- switch_df$pos_start %in% flip_pos
  switch_df$switch <- !(switch_df$flip_start | switch_df$flip_end)

  flip_pos <- switch_df %>% filter(flip_start) %>% pull(pos_end)
  switch_pos <- switch_df %>% filter(switch) %>% pull(pos_end)
  flip_end_pos <- switch_df %>% filter(flip_end) %>% pull(pos_end)

  return(list(flip_pos = flip_pos, switch_pos = switch_pos, flip_end_pos = flip_end_pos))

}

# function for reading in het position file
# pos gt cpg motif
read_het_df <- function(het_loc_dir, id){
  f_name <- paste0(het_loc_dir, "pair_", id, ".csv")
  return(read_csv(f_name, show_col_types = F))
}

# annotate each het position with error status
annotate_het_df <- function(het_loc_dir, beagle_switch_dir, eagle_switch_dir, shapeit_switch_dir, id){
  het_df <- read_het_df(het_loc_dir, id)

  beagle_obj <- switch_err_positions(beagle_switch_dir, id)
  eagle_obj <- switch_err_positions(eagle_switch_dir, id)
  shapeit_obj <- switch_err_positions(shapeit_switch_dir, id)

  het_df$beagle_switch <- as.numeric(het_df$pos %in% beagle_obj$switch_pos)
  het_df$eagle_switch <- as.numeric(het_df$pos %in% eagle_obj$switch_pos)
  het_df$shapeit_switch <- as.numeric(het_df$pos %in% shapeit_obj$switch_pos)

  het_df$beagle_flip <- as.numeric(het_df$pos %in% beagle_obj$flip_pos)
  het_df$eagle_flip <- as.numeric(het_df$pos %in% eagle_obj$flip_pos)
  het_df$shapeit_flip <- as.numeric(het_df$pos %in% shapeit_obj$flip_pos)

  het_df$beagle_flip_end <- as.numeric(het_df$pos %in% beagle_obj$flip_end_pos)
  het_df$eagle_flip_end <- as.numeric(het_df$pos %in% eagle_obj$flip_end_pos)
  het_df$shapeit_flip_end <- as.numeric(het_df$pos %in% shapeit_obj$flip_end_pos)

  return(het_df)
}

# get summary for 1 sd
get_switch_summary <- function(het_loc_dir, beagle_switch_dir, eagle_switch_dir, shapeit_switch_dir, id){
  het_df <- annotate_het_df(het_loc_dir, beagle_switch_dir, eagle_switch_dir, shapeit_switch_dir, id)

  # number of switches
  n_switch_beagle <- sum(het_df$beagle_switch)
  n_switch_eagle <- sum(het_df$eagle_switch)
  n_switch_shapeit <- sum(het_df$shapeit_switch)

  # number of flips
  n_flip_beagle <- sum(het_df$beagle_flip)
  n_flip_eagle <- sum(het_df$eagle_flip)
  n_flip_shapeit <- sum(het_df$shapeit_flip)

  # number hets
  n_het <- length(het_df$pos)

  # hets and errors at CpG
  n_het_cpg <- sum(het_df$cpg)
  n_switch_cpg_beagle <- het_df %>% filter(cpg == 1) %>% pull(beagle_switch) %>% sum()
  n_switch_cpg_eagle <- het_df %>% filter(cpg == 1) %>% pull(eagle_switch) %>% sum()
  n_switch_cpg_shapeit <- het_df %>% filter(cpg == 1) %>% pull(shapeit_switch) %>% sum()
  n_flip_cpg_beagle <- het_df %>% filter(cpg == 1) %>% pull(beagle_flip) %>% sum()
  n_flip_cpg_eagle <- het_df %>% filter(cpg == 1) %>% pull(eagle_flip) %>% sum()
  n_flip_cpg_shapeit <- het_df %>% filter(cpg == 1) %>% pull(shapeit_flip) %>% sum()

  # span of heterozygous sites
  n_bp <- max(het_df$pos) - min(het_df$pos)

  return(data.frame(id = id,
                    n_het = n_het,
                    n_het_cpg = n_het_cpg,
                    n_bp = n_bp,
                    size_MB = n_bp / 1e6,
                    n_switch_beagle = n_switch_beagle,
                    n_switch_eagle = n_switch_eagle,
                    n_switch_shapeit = n_switch_shapeit,
                    n_switch_cpg_beagle = n_switch_cpg_beagle,
                    n_switch_cpg_eagle = n_switch_cpg_eagle,
                    n_switch_cpg_shapeit = n_switch_cpg_shapeit,
                    n_flip_beagle = n_flip_beagle,
                    n_flip_eagle = n_flip_eagle,
                    n_flip_shapeit = n_flip_shapeit,
                    n_flip_cpg_beagle = n_flip_cpg_beagle,
                    n_flip_cpg_eagle = n_flip_cpg_eagle,
                    n_flip_cpg_shapeit = n_flip_cpg_shapeit
  ))
}

# get summary for range of sd
get_all_summary <- function(het_loc_dir, beagle_switch_dir, eagle_switch_dir, shapeit_switch_dir, ids = 1:600){
  results <- vector(mode='list', length=length(ids))

  for(i in 1:length(ids)){
    results[[i]] <- get_switch_summary(het_loc_dir, beagle_switch_dir, eagle_switch_dir, shapeit_switch_dir, ids[i])
  }
  return(bind_rows(results))
}

for(i in 1:1000){
  print(i)
  het_df <- annotate_het_df(het_loc_dir, beagle_switch_dir, eagle_switch_dir, shapeit_switch_dir, i)
  write_csv(het_df, paste0(het_loc_dir, "errors_", i, ".csv"))
}
