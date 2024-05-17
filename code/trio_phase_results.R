library(tidyverse)
library(yaml)

config_obj <- yaml::read_yaml("_config.yaml")

source(paste0(config_obj$base_dir,"/code/common_functions.R"))

beagle_switch_dir <-   paste0(config_obj$base_dir,"/",config_obj$trio_result_dir,"/switch_errors/beagle/annotated/")
eagle_switch_dir <-   paste0(config_obj$base_dir,"/",config_obj$trio_result_dir,"/switch_errors/eagle/annotated/")
shapeit_switch_dir <-   paste0(config_obj$base_dir,"/",config_obj$trio_result_dir,"/switch_errors/shapeit/annotated/")

het_loc_dir <- paste0(config_obj$base_dir,"/",config_obj$trio_result_dir,"/het_pos/")
df_subj <- read_csv("data/1kgp/subject_info.csv")
ped_df <- read_table("data/1kgp/1kGP.3202_samples.pedigree_info.txt")
child_ids <- ped_df %>% filter(fatherID != "0" & motherID != "0") %>% pull(sampleID)

df_subj_rel <- df_subj %>%
  filter(SAMPLE_NAME %in% child_ids)

df_subj_rel$id2 <- 1:602
gc_content_1kb <- read_tsv("data/chr15_gc1kb_pilot.bed")
colnames(gc_content_1kb) <- c("CHR", "START", "END", "AT", "GC", "A", "C", "G", "T", "TOTAL", "OTHER", "LENGTH")
gc_content_1kb  <- gc_content_1kb %>%
  mutate(bin_id = (START / 1000) + 1)


df_trio <- lapply(c(1:602),
                    function(x){
                      switch_summary_trio(x, eagle_switch_dir, beagle_switch_dir, shapeit_switch_dir, gc_content_1kb, het_loc_dir)
                      }) %>%
  bind_rows()

df_trio$id2 <- 1:602
df_trio <- left_join(df_trio, df_subj_rel, by="id2")
