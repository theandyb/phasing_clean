library(tidyverse)
library(yaml)

config_obj <- yaml::read_yaml("_config.yaml")

out_dir <- paste0(config_obj$base_dir, "/data/1kgp/chr8/exclude_samples/")

# load the subject information
subj_df <- read_csv(paste0(config_obj$base_dir, "/data/1kgp/subject_info.csv"))
unrel_id <- read_csv(paste0(config_obj$base_dir, "/data/1kgp/unrelated_subjects.txt"), col_names = "id")

df_ped <- read_table(paste0(config_obj$base_dir, "/data/1kgp/1kGP.3202_samples.pedigree_info.txt"))

df_trio_ped <- df_ped %>%
  filter(fatherID != "0", motherID != "0")

data.table::fwrite(list(subj_df_rel$SAMPLE_NAME), paste0(config_obj$base_dir, "/data/1kgp/child_ids.txt"))

subj_df_unrel <- subj_df %>%
  filter(SAMPLE_NAME %in% unrel_id$id)
subj_df_rel <- subj_df %>%
  filter(SAMPLE_NAME %in% df_trio_ped$sampleID)


# Get list of samples to exclude from each sample's specific reference panel
for(i in 1:602){
  child_id <- subj_df_rel$SAMPLE_NAME[i]
  f_id <- subj_df_rel$fatherID[i]
  m_id <- subj_df_rel$motherID[i]

  remove_id <- c()

  if(child_id %in% unrel_id$id){
    remove_id <- c(remove_id,child_id)
  }
  if(m_id %in% unrel_id$id){
    remove_id <- c(m_id, remove_id)
  }
  if(f_id %in% unrel_id$id){
    remove_id <- c(f_id, remove_id)
  }
  write_tsv(data.frame(id = remove_id), paste0(out_dir, "sample_", i, ".txt") ,col_names = F)
}
i <- 33

