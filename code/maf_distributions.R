library(tidyverse)
library(yaml)
source("code/common_functions.R")

config_obj <- yaml::read_yaml("_config.yaml")

eagle_switch_dir <-   paste0(config_obj$base_dir,"/output/switch_errors/switch_errors/eagle/annotated/")
shapeit_switch_dir <- paste0(config_obj$base_dir,"/output/switch_errors/switch_errors/shapeit/annotated/")
beagle_switch_dir <-  paste0(config_obj$base_dir,"/output/switch_errors/switch_errors/beagle/annotated/")
het_dir <- paste0(config_obj$base_dir, "/output/switch_errors/het_loc/annotated/")
maf_df <- read_tsv(paste0(config_obj$base_dir,"/data/1kgp/chrX_freq.tsv"),
                   col_names = c("chrom", "pos", "maf", "AC", "AN"))

maf_df <- maf_df %>%
  mutate(maf_cat = cut(maf, breaks = c(-Inf, 0.01, 0.05, 1), labels = c("rare", "uncommon", "common")),
         maf_quin = ntile(maf, 5),
         maf_dec = ntile(maf,10))

get_df_het <- function(id){
  df_beagle <- read_csv(paste0(beagle_switch_dir, "switch_", id, ".csv"), show_col_types = FALSE)
  df_eagle <- read_csv(paste0(eagle_switch_dir, "switch_", id, ".csv"), show_col_types = FALSE)
  df_shapeit <- read_csv(paste0(shapeit_switch_dir, "switch_", id, ".csv"), show_col_types = FALSE)
  df_het <- read_csv(paste0(het_dir, "pair_", id, ".csv"), show_col_types = FALSE)

  df_het <- left_join(df_het, {maf_df %>% select(-chrom)})

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

  df_het$beagle_status <- "correct"
  df_het$eagle_status <- "correct"
  df_het$shapeit_status <- "correct"

  df_het$beagle_status[which(df_het$pos %in%
                               {df_beagle %>%
                                   filter(is_switch) %>%
                                   pull(pos_end)})] <- "switch"
  df_het$beagle_status[which(df_het$pos %in%
                               {df_beagle %>%
                                   filter(part_of_flip) %>%
                                   pull(pos_end)})] <- "flip end"
  df_het$beagle_status[which(df_het$pos %in%
                               {df_beagle %>%
                                   filter(is_flip) %>%
                                   pull(pos_end)})] <- "flip"

  df_het$eagle_status[which(df_het$pos %in%
                              {df_eagle %>%
                                  filter(is_switch) %>%
                                  pull(pos_end)})] <- "switch"
  df_het$eagle_status[which(df_het$pos %in%
                              {df_eagle %>%
                                  filter(part_of_flip) %>%
                                  pull(pos_end)})] <- "flip end"
  df_het$eagle_status[which(df_het$pos %in%
                              {df_eagle %>%
                                  filter(is_flip) %>%
                                  pull(pos_end)})] <- "flip"

  df_het$shapeit_status[which(df_het$pos %in%
                                {df_shapeit %>%
                                    filter(is_switch) %>%
                                    pull(pos_end)})] <- "switch"
  df_het$shapeit_status[which(df_het$pos %in%
                                {df_shapeit %>%
                                    filter(part_of_flip) %>%
                                    pull(pos_end)})] <- "flip end"
  df_het$shapeit_status[which(df_het$pos %in%
                                {df_shapeit %>%
                                    filter(is_flip) %>%
                                    pull(pos_end)})] <- "flip"
  df_het$id <- id
  return(df_het)
}

df_het_summary <- function(id){
  df_het <- get_df_het(id)
  df_res <- data.frame(category = character(), proportion = numeric())
  # all hets categories
  df_res <- bind_rows(df_res, data.frame(category = "het_rare",
                                         proportion = sum(df_het$maf_cat == "rare") / length(df_het$maf_cat)))
  df_res <- bind_rows(df_res, data.frame(category = "het_uncommon",
                                         proportion = sum(df_het$maf_cat == "uncommon") / length(df_het$maf_cat)))
  df_res <- bind_rows(df_res, data.frame(category = "het_common",
                                         proportion = sum(df_het$maf_cat == "common") / length(df_het$maf_cat)))

  # beagle errors
  for(cat_val in unique(df_het$beagle_status)){
    for(maf_cat_val in unique(df_het$maf_cat)){
      cat_name <- paste0("beagle_", cat_val, "_", maf_cat_val)
      prop_val <- sum(df_het$beagle_status == cat_val & df_het$maf_cat == maf_cat_val) / sum(df_het$beagle_status == cat_val)
      df_res <- bind_rows(df_res, data.frame(category = cat_name, proportion = prop_val))}
  }
  # eagle errors
  for(cat_val in unique(df_het$eagle_status)){
    for(maf_cat_val in unique(df_het$maf_cat)){
      cat_name <- paste0("eagle_", cat_val, "_", maf_cat_val)
      prop_val <- sum(df_het$eagle_status == cat_val & df_het$maf_cat == maf_cat_val) / sum(df_het$eagle_status == cat_val)
      df_res <- bind_rows(df_res, data.frame(category = cat_name, proportion = prop_val))}
  }
  #shapeit errors
  for(cat_val in unique(df_het$shapeit_status)){
    for(maf_cat_val in unique(df_het$maf_cat)){
      cat_name <- paste0("shapeit_", cat_val, "_", maf_cat_val)
      prop_val <- sum(df_het$shapeit_status == cat_val & df_het$maf_cat == maf_cat_val) / sum(df_het$shapeit_status == cat_val)
      df_res <- bind_rows(df_res, data.frame(category = cat_name, proportion = prop_val))}
  }

  df_res$id <- id
  return(df_res)
}

df_prop <- df_het_summary(1)
for(i in 2:1000){
  df_prop <- bind_rows(df_prop, df_het_summary(i))
}
write_csv(df_prop, paste0(config_obj$base_dir,"/output/switch_errors/maf_props.csv"))

### MAF distribution at sites flanking switch errors
df_test <- get_df_het(201)
