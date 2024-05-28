library(tidyverse)
library(yaml)
source("code/common_functions.R")

config_obj <- yaml::read_yaml("_config.yaml")

eagle_switch_dir <-   paste0(config_obj$base_dir,"/output/switch_errors/switch_errors/eagle/annotated/")
shapeit_switch_dir <- paste0(config_obj$base_dir,"/output/switch_errors/switch_errors/shapeit/annotated/")
beagle_switch_dir <-  paste0(config_obj$base_dir,"/output/switch_errors/switch_errors/beagle/annotated/")
het_loc_dir <- paste0(config_obj$base_dir, "/output/switch_errors/het_loc/annotated/")
maf_df <- read_tsv(paste0(config_obj$base_dir,"/data/1kgp/chrX_maf.tsv"),
                   col_names = c("chrom", "pos", "maf"))

# single example
df_het <- read_csv(paste0(het_loc_dir, "pair_", 1, ".csv"), show_col_types = FALSE)
df_beagle <- read_csv(paste0(beagle_switch_dir, "switch_", 1, ".csv"), show_col_types = FALSE)
df_eagle <- read_csv(paste0(eagle_switch_dir, "switch_", 1, ".csv"), show_col_types = FALSE)
df_shapeit <- read_csv(paste0(shapeit_switch_dir, "switch_", 1, ".csv"), show_col_types = FALSE)

# Identify flips and switches
flip_pos_beagle <- get_flip_pos(df_beagle)
flip_pos_eagle <- get_flip_pos(df_eagle)
flip_pos_shapeit <- get_flip_pos(df_shapeit)

df_eagle$is_flip <- df_eagle$pos_end %in% flip_pos_eagle
df_beagle$is_flip <- df_beagle$pos_end %in% flip_pos_beagle
df_shapeit$is_flip <- df_shapeit$pos_end %in% flip_pos_shapeit

df_eagle$part_of_flip <- ifelse(df_eagle$is_flip, TRUE, df_eagle$pos_start %in% flip_pos_eagle)
df_beagle$part_of_flip <- ifelse(df_beagle$is_flip, TRUE, df_beagle$pos_start %in% flip_pos_beagle)
df_shapeit$part_of_flip <- ifelse(df_shapeit$is_flip, TRUE, df_shapeit$pos_start %in% flip_pos_shapeit)

df_beagle$is_switch <- !(df_beagle$is_flip | df_beagle$part_of_flip)
df_eagle$is_switch <- !(df_eagle$is_flip | df_eagle$part_of_flip)
df_shapeit$is_switch <- !(df_shapeit$is_flip | df_shapeit$part_of_flip)

# Label hets as switch, flip start, flip end, or no error
df_het$beagle_status <- "correct"
df_het$eagle_status <- "correct"
df_het$shapeit_status <- "correct"

df_het$beagle_status[which(df_het$pos %in%
                             {df_beagle %>%
                                 filter(is_switch) %>%
                                 pull(pos_start)})] <- "switch"
df_het$beagle_status[which(df_het$pos %in%
                             {df_beagle %>%
                                 filter(part_of_flip) %>%
                                 pull(pos_start)})] <- "flip end"
df_het$beagle_status[which(df_het$pos %in%
                             {df_beagle %>%
                                 filter(is_flip) %>%
                                 pull(pos_start)})] <- "flip"

df_het$eagle_status[which(df_het$pos %in%
                             {df_eagle %>%
                                 filter(is_switch) %>%
                                 pull(pos_start)})] <- "switch"
df_het$eagle_status[which(df_het$pos %in%
                             {df_eagle %>%
                                 filter(part_of_flip) %>%
                                 pull(pos_start)})] <- "flip end"
df_het$eagle_status[which(df_het$pos %in%
                             {df_eagle %>%
                                 filter(is_flip) %>%
                                 pull(pos_start)})] <- "flip"

df_het$shapeit_status[which(df_het$pos %in%
                             {df_shapeit %>%
                                 filter(is_switch) %>%
                                 pull(pos_start)})] <- "switch"
df_het$shapeit_status[which(df_het$pos %in%
                             {df_shapeit %>%
                                 filter(part_of_flip) %>%
                                 pull(pos_start)})] <- "flip end"
df_het$shapeit_status[which(df_het$pos %in%
                             {df_shapeit %>%
                                 filter(is_flip) %>%
                                 pull(pos_start)})] <- "flip"
