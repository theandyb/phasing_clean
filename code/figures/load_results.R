library(tidyverse)

# Load list of heterozygous positions for an individual synthetic diploid
load_het_sd <- function(id, base_dir = "/net/snowwhite/home/beckandy/research/phasing_clean/output/switch_errors/het_loc/annotated/"){
    f_name <- paste0(base_dir, "errors_", id, ".csv")
    read_csv(f_name, show_col_types = FALSE)
}

# Load list of heterozygous positions for an individual proband
load_het_proband <- function(id, chrom, base_dir = "/net/snowwhite/home/beckandy/research/phasing_clean/output/"){
    res_dir <- paste0(base_dir, "trio_phase_", chrom, "/no_th/het_loc/")
    f_name <- paste0(res_dir, "errors_", id, ".csv")
    read_csv(f_name, show_col_types = FALSE)
}

# Load sample information df for 1kgp
load_1kgp_info <- function(f_name = "/net/snowwhite/home/beckandy/research/phasing_clean/data/1kgp/subject_info.csv"){
    read_csv(f_name, show_col_types = FALSE) %>%
        select(SAMPLE_NAME, POPULATION, SUPER)
}

# Load 1kgp pedigree information
load_1kgp_ped <- function(f_name = "/net/snowwhite/home/beckandy/research/phasing_clean/data/1kgp/1kGP.3202_samples.pedigree_info.txt"){
    read_table(f_name, show_col_types = FALSE)
}

# Get ids of probands in 1kgp pedigree
get_1kgp_probands <- function(ped_df){
    ped_df %>% 
        filter(fatherID != "0" & motherID != "0") %>% 
        pull(sampleID)
}

# Get ids of unrelated individuals in 1kgp
get_1kgp_unrelated <- function(f_name = "/net/snowwhite/home/beckandy/research/phasing_clean/data/1kgp/unrelated_subjects.txt"){
    scan(f_name, what =  character())
}

# Get subject information for probands only
get_1kgp_proband_info <- function(){
    df_subj <- load_1kgp_info()
    child_ids <- get_1kgp_probands(load_1kgp_ped())
    df_subj <- df_subj %>%
        filter(SAMPLE_NAME %in% child_ids)
    df_subj$id2 <- 1:nrow(df_subj)
    df_subj
}

# Get subject information for unrelated individuals only
get_1kgp_unrelated_info <- function(){
    df_subj <- load_1kgp_info()
    unrelated_ids <- get_1kgp_unrelated()
    df_subj %>%
        filter(SAMPLE_NAME %in% unrelated_ids)
}

# Load synthetic diploid pair information
load_sd_info <- function(f_name = "/net/snowwhite/home/beckandy/research/phasing_clean/data/sample_pairs.csv"){
    pair_info_df <- read_delim(f_name, col_names = c("POP", "ID1", "ID2"))
    sample_info_df <- load_1kgp_info()
    pair_info_df <- left_join(pair_info_df, sample_info_df, by = c("ID1"="SAMPLE_NAME")) %>%
        rename(SP = SUPER) %>%
        select(-POPULATION)
    pair_info_df$id <- 1:1000
    pair_info_df
}

# Load X chromosome GC content in 1kb windows
load_gc_content_X <- function(f_name = "/net/snowwhite/home/beckandy/research/phasing_clean/data/chrX_gc1kb_pilot.bed"){
    gc_content_1kb <- read_tsv(f_name)
    colnames(gc_content_1kb) <- c("CHR", "START", "END", "AT", "GC", "A", "C", "G", "T", "TOTAL", "OTHER", "LENGTH")
    gc_content_1kb  <- gc_content_1kb %>%
        mutate(bin_id = (START / 1000) + 1)
    gc_content_1kb
}

# load maf category proportion df for synthetic diploids
load_sd_maf_prop <- function(f_name = "/net/snowwhite/home/beckandy/research/phasing_clean/output/switch_errors/maf_props.csv"){
    maf_prop_df <- read_csv(f_name, show_col_types = FALSE)
    pair_info_df <- load_sd_info()
    maf_prop_df <- left_join(maf_prop_df, {pair_info_df %>% select(id, SP, POP)})
    maf_prop_df
}

# load switch/flip error counts for synthetic diploids
load_sd_errors <- function(f_name = "/net/snowwhite/home/beckandy/research/phasing_clean/output/switch_errors/switch_errors/summary.csv"){
    error_df <- read_csv(f_name, show_col_types = FALSE)
    error_df$n_total_beagle <- (error_df$n_other_beagle + error_df$n_flip_beagle)
    error_df$n_total_shapeit <- (error_df$n_other_shapeit + error_df$n_flip_shapeit)
    error_df$n_total_eagle <- (error_df$n_other_eagle + error_df$n_flip_eagle)
    error_df
}

# load switch/flip error counts for probands
load_proband_errors <- function(f_name = "/net/snowwhite/home/beckandy/research/phasing_clean/output/all_noTH.csv"){
    read_csv(f_name, show_col_types = FALSE)
}

# load switch/flips in MB bins for synthetic diploids
load_sd_errors_mb <- function(f_name = "/net/snowwhite/home/beckandy/research/phasing_clean/output/switch_errors/errors_mb_bin.csv"){
    read_csv(f_name, show_col_types = FALSE)
}
