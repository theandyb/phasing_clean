library(tidyverse)
library(yaml)

config_obj <- yaml::read_yaml("_config.yaml")

out_dir <- paste0(config_obj$base_dir, "/data/1kgp/chr15/triple_het/trio_lists/")

# load the subject information
subj_df <- read_csv(paste0(config_obj$base_dir, "/data/1kgp/subject_info.csv"))
unrel_id <- read_csv(paste0(config_obj$base_dir, "/data/1kgp/unrelated_subjects.txt"), col_names = "id")

df_ped <- read_table(paste0(config_obj$base_dir, "/data/1kgp/1kGP.3202_samples.pedigree_info.txt"))

df_trio_ped <- df_ped %>%
  filter(fatherID != "0", motherID != "0")

subj_df_rel <- subj_df %>%
  filter(SAMPLE_NAME %in% df_trio_ped$sampleID)


# Get list of samples to exclude from each sample's specific reference panel
for(i in 1:602){
  child_id <- subj_df_rel$SAMPLE_NAME[i]
  f_id <- subj_df_rel$fatherID[i]
  m_id <- subj_df_rel$motherID[i]

  write_tsv(data.frame(id = c(child_id, f_id, m_id)), paste0(out_dir, "sample_", i, ".txt") ,col_names = F)
}

