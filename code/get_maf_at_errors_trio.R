library(tidyverse)
library(yaml)
source("code/common_functions.R")

config_obj <- yaml::read_yaml("_config.yaml")

eagle_switch_dir <-   paste0(config_obj$base_dir,"/output/trio_phase_15/switch_errors/eagle/annotated/")
shapeit_switch_dir <- paste0(config_obj$base_dir,"/output/trio_phase_15/switch_errors/shapeit/annotated/")
beagle_switch_dir <-  paste0(config_obj$base_dir,"/output/trio_phase_15/switch_errors/beagle/annotated/")

maf_df <- read_tsv(paste0(config_obj$base_dir,"/data/1kgp/chr15/chr15_freq.tsv"),
                   col_names = c("chrom", "pos", "maf"))

df_subj <- read_csv(paste0(config_obj$base_dir,"/data/1kgp/subject_info.csv"))
ped_df <- read_table(paste0(config_obj$base_dir,"/data/1kgp/1kGP.3202_samples.pedigree_info.txt"))
child_ids <- ped_df %>% filter(fatherID != "0" & motherID != "0") %>% pull(sampleID)
df_subj_rel <- df_subj %>%
  filter(SAMPLE_NAME %in% child_ids)


# objects to store results
beagle_switches  <- vector(mode = "list", length = 602)
beagle_flips  <- vector(mode = "list", length = 602)
eagle_switches  <- vector(mode = "list", length = 602)
eagle_flips  <- vector(mode = "list", length = 602)
shapeit_switches  <- vector(mode = "list", length = 602)
shapeit_flips  <- vector(mode = "list", length = 602)

for(i in 1:602){
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
for(i in 1:602){
  pop <- df_subj_rel$SUPER[i]
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

shapeit_switches_df <- bind_rows(shapeit_switches) %>%
  group_by(pop, pos_start) %>%
  summarize(n = n())
shapeit_flips_df <- bind_rows(shapeit_flips) %>%
  group_by(pop, pos_start) %>%
  summarize(n = n())

# save for later analysis
write_csv(beagle_switches_df, paste0(config_obj$base_dir, "/output/trio_phase_15/beagle_switch_maf.csv"))
write_csv(beagle_flips_df, paste0(config_obj$base_dir, "/output/trio_phase_15/beagle_flip_maf.csv"))

write_csv(eagle_switches_df, paste0(config_obj$base_dir, "/output/trio_phase_15/eagle_switch_maf.csv"))
write_csv(eagle_flips_df, paste0(config_obj$base_dir, "/output/trio_phase_15/eagle_flip_maf.csv"))

write_csv(shapeit_switches_df, paste0(config_obj$base_dir, "/output/trio_phase_15/shapeit_switch_maf.csv"))
write_csv(shapeit_flips_df, paste0(config_obj$base_dir, "/output/trio_phase_15/shapeit_flip_maf.csv"))
