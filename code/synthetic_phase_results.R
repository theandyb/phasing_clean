library(tidyverse)
library(yaml)

source("code/common_functions.R")

config_obj <- yaml::read_yaml("_config.yaml")

eagle_switch_dir <-   paste0(config_obj$base_dir,"/output/switch_errors/switch_errors/eagle/annotated/")
shapeit_switch_dir <- paste0(config_obj$base_dir,"/output/switch_errors/switch_errors/shapeit/annotated/")
beagle_switch_dir <-  paste0(config_obj$base_dir,"/output/switch_errors/switch_errors/beagle/annotated/")
num_sites_dir <-      paste0(config_obj$base_dir,"/output/switch_errors/vcf_n_sites/")
het_loc_dir <-        paste0(config_obj$base_dir,"/output/switch_errors/het_loc/")

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
