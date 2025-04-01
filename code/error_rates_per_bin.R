# get error rates per bin for each subject and method
library(tidyverse)
library(yaml)

config_obj <- yaml::read_yaml("_config.yaml")

source(paste0(config_obj$base_dir,"/code/common_functions.R"))

chrom <- 1

beagle_switch_dir <-   paste0(config_obj$base_dir,"/output/trio_phase_", chrom ,"/switch_errors/beagle/annotated/")
eagle_switch_dir <-   paste0(config_obj$base_dir,"/output/trio_phase_", chrom ,"/switch_errors/eagle/annotated/")
shapeit_switch_dir <-   paste0(config_obj$base_dir,"/output/trio_phase_", chrom ,"/switch_errors/shapeit/annotated/")
het_loc_dir <- paste0(config_obj$base_dir,"/output/trio_phase_", chrom , "/het_loc/annotated/")
out_dir <- paste0(config_obj$base_dir,"/output/trio_phase_", chrom , "/mb_error_rate/")

# beagle_switch_dir <-   paste0(config_obj$base_dir,"/output/switch_errors/switch_errors/beagle/annotated/")
# eagle_switch_dir <-   paste0(config_obj$base_dir,"/output/switch_errors/switch_errors/eagle/annotated/")
# shapeit_switch_dir <-   paste0(config_obj$base_dir,"/output/switch_errors/switch_errors/shapeit/annotated/")
# het_loc_dir <- paste0(config_obj$base_dir, "/output/switch_errors/het_loc/annotated/")


#out_dir <- paste0(config_obj$base_dir, "/", config_obj$trio_result_dir, "/mb_error_rate/")
#out_dir <- paste0(config_obj$base_dir, "/output/switch_errors/mb_error_rate/")


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

# annotated het position files for all subjects
for(i in 1:602){
  print(i)
  het_df <- annotate_het_df(het_loc_dir, beagle_switch_dir, eagle_switch_dir, shapeit_switch_dir, i)
  write_csv(het_df, paste0(het_loc_dir, "errors_", i, ".csv"))
}


# errors per MB bin
for(i in 1:602){
  print(i)
  het_df <- annotate_het_df(het_loc_dir, beagle_switch_dir, eagle_switch_dir, shapeit_switch_dir, i)
  het_df<- het_df %>%
    mutate(bin = ceiling(pos / 1e6)) %>%
    group_by(bin) %>%
    summarise(n_het = n(),
              n_switch_beagle = sum(beagle_switch),
              n_switch_eagle = sum(eagle_switch),
              n_switch_shapeit = sum(shapeit_switch),
              n_flip_beagle = sum(beagle_flip),
              n_flip_eagle = sum(eagle_flip),
              n_flip_shapeit = sum(shapeit_flip)) %>%
    select(bin, n_het,contains("n_switch"), contains("n_flip"))
  write_csv(het_df, paste0(out_dir, "sample_", i, ".csv"))

}
