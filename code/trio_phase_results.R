library(tidyverse)
library(yaml)

config_obj <- yaml::read_yaml("_config.yaml")

source(paste0(config_obj$base_dir,"/code/common_functions.R"))

df_subj <- read_csv("data/1kgp/subject_info.csv")
ped_df <- read_table("data/1kgp/1kGP.3202_samples.pedigree_info.txt")
child_ids <- ped_df %>% filter(fatherID != "0" & motherID != "0") %>% pull(sampleID)

df_subj_rel <- df_subj %>%
  filter(SAMPLE_NAME %in% child_ids)

df_subj_rel$id2 <- 1:602

for(chrom in 1:22){
  beagle_switch_dir <-   paste0(config_obj$base_dir,"/output/trio_phase_", chrom ,"/no_th/switch_errors/beagle/")
  eagle_switch_dir <-   paste0(config_obj$base_dir,"/output/trio_phase_", chrom ,"/no_th/switch_errors/eagle/")
  shapeit_switch_dir <-   paste0(config_obj$base_dir,"/output/trio_phase_", chrom ,"/no_th/switch_errors/shapeit/")
  het_loc_dir <- paste0(config_obj$base_dir,"/output/trio_phase_", chrom , "/no_th/het_loc/")

  gc_content_1kb <- read_tsv(paste0("data/chr", chrom,"_gc1kb_pilot.bed"))
  colnames(gc_content_1kb) <- c("CHR", "START", "END", "AT", "GC", "A", "C", "G", "T", "TOTAL", "OTHER", "LENGTH")
  gc_content_1kb  <- gc_content_1kb %>%
    mutate(bin_id = (START / 1000) + 1)
  n_het_cpg_df <- read_tsv(paste0(config_obj$base_dir,"/output/trio_phase_", chrom, "/no_th/het_cpg_count.tsv"),
                           col_names = c("pair_id", "n_het_cpg"))

  print(chrom)
  df_trio <- lapply(c(1:602),
                      function(x){
                        switch_summary(x, eagle_switch_dir, beagle_switch_dir, shapeit_switch_dir, gc_content_1kb, het_loc_dir)
                        }) %>%
    bind_rows()

  df_trio$id2 <- 1:602
  df_trio <- left_join(df_trio, df_subj_rel, by="id2")
  df_trio <- left_join(df_trio, n_het_cpg_df)

  write_csv(df_trio, paste0(config_obj$base_dir,"/output/trio_phase_", chrom, "/no_th/switch_errors/summary.csv"))
}
