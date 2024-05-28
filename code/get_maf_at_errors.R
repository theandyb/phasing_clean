library(tidyverse)
library(yaml)
source("code/common_functions.R")

config_obj <- yaml::read_yaml("_config.yaml")

eagle_switch_dir <-   paste0(config_obj$base_dir,"/output/switch_errors/switch_errors/eagle/annotated/")
shapeit_switch_dir <- paste0(config_obj$base_dir,"/output/switch_errors/switch_errors/shapeit/annotated/")
beagle_switch_dir <-  paste0(config_obj$base_dir,"/output/switch_errors/switch_errors/beagle/annotated/")
num_sites_dir <-      paste0(config_obj$base_dir,"/output/switch_errors/vcf_n_sites/")

maf_df <- read_tsv(paste0(config_obj$base_dir,"/data/1kgp/chrX_maf.tsv"),
                   col_names = c("chrom", "pos", "maf"))
sample_info_df <- read_csv(paste0(config_obj$base_dir,"/data/1kgp/subject_info.csv")) %>%
  select(SAMPLE_NAME, POPULATION, SUPER)

pair_info_df <- read_delim(paste0(config_obj$base_dir,"/data/sample_pairs.csv"),
                           col_names = c("POP", "ID1", "ID2"))
pair_info_df <- left_join(pair_info_df, sample_info_df, by = c("ID1"="SAMPLE_NAME")) %>%
  rename(SP = SUPER) %>%
  select(-POPULATION)

# objects to store results
beagle_switches  <- vector(mode = "list", length = 1000)
beagle_flips  <- vector(mode = "list", length = 1000)
eagle_switches  <- vector(mode = "list", length = 1000)
eagle_flips  <- vector(mode = "list", length = 1000)
shapeit_switches  <- vector(mode = "list", length = 1000)
shapeit_flips  <- vector(mode = "list", length = 1000)

for(i in 1:1000){
  print(i)
  # read switch errors
  df_beagle <- read_csv(paste0(beagle_switch_dir, "switch_", i, ".csv"), show_col_types = FALSE)
  df_eagle <- read_csv(paste0(eagle_switch_dir, "switch_", i, ".csv"), show_col_types = FALSE)
  df_shapeit <- read_csv(paste0(shapeit_switch_dir, "switch_", i, ".csv"), show_col_types = FALSE)

  # assign switch/flip status
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

  # append maf to pos start
  df_beagle <- left_join(df_beagle, maf_df, by=c("pos_start" = "pos"))
  df_eagle <- left_join(df_eagle, maf_df, by=c("pos_start" = "pos"))
  df_shapeit <- left_join(df_shapeit, maf_df, by=c("pos_start" = "pos"))

  # store mafs
  beagle_switches[[i]] <- df_beagle %>% filter(is_switch) %>% select(pos_start, maf)
  beagle_flips[[i]] <- df_beagle %>% filter(is_flip) %>% select(pos_start, maf)

  eagle_switches[[i]] <- df_eagle %>% filter(is_switch) %>% select(pos_start, maf)
  eagle_flips[[i]] <- df_eagle %>% filter(is_flip) %>% select(pos_start, maf)

  shapeit_switches[[i]] <- df_shapeit %>% filter(is_switch) %>% select(pos_start, maf)
  shapeit_flips[[i]] <- df_shapeit %>% filter(is_flip) %>% select(pos_start, maf)
}

# assign population labels to each df
for(i in 1:1000){
  pop <- pair_info_df$SP[i]
  beagle_switches[[i]]$pop <- pop
  beagle_flips[[i]]$pop <- pop
  eagle_switches[[i]]$pop <- pop
  eagle_flips[[i]]$pop <- pop
  shapeit_switches[[i]]$pop <- pop
  shapeit_flips[[i]]$pop <- pop
}

# smash into single data frames, get unique rows

beagle_switches_df <- bind_rows(beagle_switches) %>%
  group_by(pop, pos_start) %>%
  summarize(n = n())

beagle_flips_df <- bind_rows(beagle_flips) %>%
  group_by(pop, pos_start) %>%
  summarize(n = n())

eagle_switches_df <- bind_rows(eagle_switches) %>%
  group_by(pop, pos_start) %>%
  summarize(n = n())

eagle_flips_df <- bind_rows(eagle_flips) %>%
  group_by(pop, pos_start) %>%
  summarize(n = n())

shapeit_switches_df <- bind_rows(shapeit_switches)%>%
  group_by(pop, pos_start) %>%
  summarize(n = n())

shapeit_flips_df <- bind_rows(shapeit_flips) %>%
  group_by(pop, pos_start) %>%
  summarize(n = n())

# save for later analysis
write_csv(beagle_switches_df, paste0(config_obj$base_dir, "/output/switch_errors/beagle_switch_maf.csv"))
write_csv(beagle_flips_df, paste0(config_obj$base_dir, "/output/switch_errors/beagle_flip_maf.csv"))

write_csv(eagle_switches_df, paste0(config_obj$base_dir, "/output/switch_errors/eagle_switch_maf.csv"))
write_csv(eagle_flips_df, paste0(config_obj$base_dir, "/output/switch_errors/eagle_flip_maf.csv"))

write_csv(shapeit_switches_df, paste0(config_obj$base_dir, "/output/switch_errors/shapeit_switch_maf.csv"))
write_csv(shapeit_flips_df, paste0(config_obj$base_dir, "/output/switch_errors/shapeit_flip_maf.csv"))
