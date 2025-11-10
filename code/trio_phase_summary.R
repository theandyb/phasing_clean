library(tidyverse)
df_subj <- read_csv("data/1kgp/subject_info.csv")
ped_df <- read_table("data/1kgp/1kGP.3202_samples.pedigree_info.txt")

child_ids <- ped_df %>% filter(fatherID != "0" & motherID != "0") %>% pull(sampleID)
unrel_ids <- scan("data/1kgp/unrelated_subjects.txt", what = character())

df_subj_rel <- df_subj %>%
  filter(SAMPLE_NAME %in% child_ids)
df_subj_rel$id2 <- 1:602

df_subj_unrel <- df_subj %>%
  filter(SAMPLE_NAME %in% unrel_ids)

df_chrom_size <- read_tsv("data/hg38.chrom.sizes", n_max = 24, col_names = c("chrom", "bases"))

load_het_pos <- function(id, het_loc){
  f_name <- paste0(het_loc, "errors_", id, ".csv")
  df <- read_csv(f_name, show_col_types = F)
  df
}

# get error rates for individual subject
single_summary <- function(id, het_loc){
  df <- load_het_pos(id, het_loc)

  df %>%
    summarize(n_het = n(),
              n_het_cpg = sum(cpg),
              switch_b = sum(beagle_switch),
              switch_e = sum(eagle_switch),
              switch_s = sum(shapeit_switch),
              flip_b = sum(beagle_flip),
              flip_e = sum(eagle_flip),
              flip_s = sum(shapeit_flip),
              cpg_switch_b = sum(beagle_switch * cpg),
              cpg_switch_e = sum(eagle_switch * cpg),
              cpg_switch_s = sum(shapeit_switch * cpg),
              cpg_flip_b = sum(beagle_flip * cpg),
              cpg_flip_e = sum(eagle_flip * cpg),
              cpg_flip_s = sum(shapeit_flip * cpg)) %>%
    mutate(r_switch_b = switch_b / n_het,
           r_switch_e = switch_e / n_het,
           r_switch_s = switch_s / n_het,
           r_flip_b = flip_b / n_het,
           r_flip_e = flip_e / n_het,
           r_flip_s = flip_s / n_het,
           r_cpg_switch_b = cpg_switch_b / n_het_cpg,
           r_cpg_switch_e = cpg_switch_e / n_het_cpg,
           r_cpg_switch_s = cpg_switch_s / n_het_cpg,
           r_cpg_flip_b = cpg_flip_b / n_het_cpg,
           r_cpg_flip_e = cpg_flip_e / n_het_cpg,
           r_cpg_flip_s = cpg_flip_s / n_het_cpg
    )
}

all_summary <- function(ids, het_loc){
  results <- vector(mode = "list", length = length(ids))
  for(i in 1:length(ids)){
    results[[i]] <- single_summary(ids[i], het_loc)
    results[[i]]$id <- ids[i]
  }
  return(bind_rows(results))
}

for(chromosome in 2:22){
  print(paste0("Working on chromosome ", chromosome))
  df <- all_summary(1:602, paste0("output/trio_phase_", chromosome, "/no_th/het_loc/"))
  df <- df %>%
    left_join({df_subj_rel %>% select(SUPER, sex, id2) %>% rename(id=id2, pop = SUPER)})
  out_file <- paste0("output/trio_phase_", chromosome, "/no_th/switch_errors/summary.csv")
  write_csv(df, out_file)
  print("Done!")
}
