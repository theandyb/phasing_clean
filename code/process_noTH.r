library(tidyverse)
library(Rcpp)

base_dir <- "/net/snowwhite/home/beckandy/research/phasing_clean/output"


### Load haplotypes from VCF files
#### Note: we are only concerned with heterozygous positions for these analyses
load_haplotype <- function(chromosome, id, base_dir, method){
  f_name <- paste0(base_dir, "/trio_phase_", chromosome, "/no_th/vcf/",
                   method, "_", id, ".vcf.gz")
  cmd_txt <- paste0('bcftools query -f "%CHROM\t%POS\t[%GT]\n" ', f_name)

  df <- vroom::vroom(pipe(cmd_txt), col_names = c("chrom", "pos", "gt"), delim = "\t", show_col_types = FALSE) %>%
    separate(gt, into = c("h1", "h2"), sep = "\\|") %>%
    filter(h1 != h2)

  return(df)
}

load_switch <- function(chromosome, id, base_dir, method){
  f_name <- paste0(base_dir, "/trio_phase_", chromosome, "/no_th/switch_errors/", method, "/switch_", id ,".csv")
  df <- read_csv(f_name, show_col_types = F)
  #names(df) <- tolower(names(df))
  return(df)
}

load_het <- function(chromosome, id, base_dir){
  f_name <- paste0(base_dir, "/trio_phase_", chromosome, "/no_th/het_loc/annotated/pair_", id, ".csv")
  df <- read_csv(f_name, show_col_types = F)
  return(df)
}

### Aside: check where switches actually occur on the haplotype
df_switch <- load_switch(22, 1, base_dir, "eagle")
df_eagle <- load_haplotype(22, 1, base_dir, "eagle")
df_truth <- load_haplotype(22, 1, base_dir, "truth")
df_het <- load_het(22, 1, base_dir)

df_switch[1,] # chr22  10705519 10930855 NA10847

df_eagle[which(df_eagle$pos == 10705519), ]
df_truth[which(df_truth$pos == 10705519), ]

df_eagle[which(df_eagle$pos == 10930855), ]
df_truth[which(df_truth$pos == 10930855), ]

# Result: the switch occurs after POS_START, but before or at POS_END

### Identify errors as switches or flips
sourceCpp("code/cpp/error_count.cpp")

test <- switch_flip(df_switch$pos_end, df_switch$pos_start)

## Get counts
### all

n_switch <- length(test$switches)
n_flip <- length(test$flips)

### At cpg
n_switch_cpg <- df_switch %>% filter(pos_end %in% test$switches, cpg_end == 1) %>% pull(pos_end) %>% length()
n_flip_cpg <- df_switch %>% filter(pos_end %in% test$flips, cpg_end == 1) %>% pull(pos_end) %>% length()

### Summary of het positions
n_het <- length(df_het$pos)
n_het_cpg <- sum(df_het$cpg)


# Get results for all samples across all autosomes
results <- data.frame(chrom = numeric(),
                      sample = numeric(),
                      n_het = numeric(),
                      n_het_cpg = numeric(),
                      n_switch_b = numeric(),
                      n_flip_b = numeric(),
                      n_switch_cpg_b = numeric(),
                      n_flip_cpg_b = numeric(),
                      n_switch_e = numeric(),
                      n_flip_e = numeric(),
                      n_switch_cpg_e = numeric(),
                      n_flip_cpg_e = numeric(),
                      n_switch_s = numeric(),
                      n_flip_s = numeric(),
                      n_switch_cpg_s = numeric(),
                      n_flip_cpg_s = numeric())

for(chromosome in 1:22){
  for(id in 1:602){
    print(paste0("Chromosome: ", chromosome, "; Sample: ", id))

    df_eagle <- load_switch(chromosome, id, base_dir, "eagle")
    df_beagle <- load_switch(chromosome, id, base_dir, "beagle")
    df_shapeit <- load_switch(chromosome, id, base_dir, "shapeit")
    df_het <- load_het(chromosome, id, base_dir)

    beagle_list <- switch_flip(df_beagle$pos_end, df_beagle$pos_start)
    eagle_list <- switch_flip(df_eagle$pos_end, df_eagle$pos_start)
    shapeit_list <- switch_flip(df_shapeit$pos_end, df_shapeit$pos_start)

    n_switch_b <- length(beagle_list$switches)
    n_flip_b <- length(beagle_list$flips)
    n_switch_e <- length(eagle_list$switches)
    n_flip_e <- length(eagle_list$flips)
    n_switch_s <- length(shapeit_list$switches)
    n_flip_s <- length(shapeit_list$flips)

    ### At cpg
    n_switch_cpg_b <- df_beagle %>% filter(pos_end %in% beagle_list$switches, cpg_end == 1) %>% pull(pos_end) %>% length()
    n_flip_cpg_b <- df_beagle %>% filter(pos_end %in% beagle_list$flips, cpg_end == 1) %>% pull(pos_end) %>% length()
    n_switch_cpg_e <- df_eagle %>% filter(pos_end %in% eagle_list$switches, cpg_end == 1) %>% pull(pos_end) %>% length()
    n_flip_cpg_e <- df_eagle %>% filter(pos_end %in% eagle_list$flips, cpg_end == 1) %>% pull(pos_end) %>% length()
    n_switch_cpg_s <- df_shapeit %>% filter(pos_end %in% shapeit_list$switches, cpg_end == 1) %>% pull(pos_end) %>% length()
    n_flip_cpg_s <- df_shapeit %>% filter(pos_end %in% shapeit_list$flips, cpg_end == 1) %>% pull(pos_end) %>% length()

    ### Summary of het positions
    n_het <- length(df_het$pos)
    n_het_cpg <- sum(df_het$cpg)

    # add to results
    results <- results %>%
      bind_rows(data.frame(chrom = chromosome,
                           sample =id,
                           n_het =n_het,
                           n_het_cpg =n_het_cpg,
                           n_switch_b = n_switch_b,
                           n_flip_b = n_flip_b,
                           n_switch_cpg_b = n_switch_cpg_b,
                           n_flip_cpg_b = n_flip_cpg_b,
                           n_switch_e = n_switch_e,
                           n_flip_e = n_flip_e,
                           n_switch_cpg_e = n_switch_cpg_e,
                           n_flip_cpg_e = n_flip_cpg_e,
                           n_switch_s = n_switch_s,
                           n_flip_s = n_flip_s,
                           n_switch_cpg_s = n_switch_cpg_s,
                           n_flip_cpg_s = n_flip_cpg_s
                           ))
  }
}
