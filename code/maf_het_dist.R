# For each (proband, chromsome), compute the proportion of heterozygous sites at rare, uncommon, and common variants

library(tidyverse)

load_het_proband <- function(id, chrom, base_dir = "/net/snowwhite/home/beckandy/research/phasing_clean/output/"){
    res_dir <- paste0(base_dir, "trio_phase_", chrom, "/no_th/het_loc/")
    f_name <- paste0(res_dir, "errors_", id, ".csv")
    read_csv(f_name, show_col_types = FALSE)
}

for(chrom in 1:21){
    print(chrom)
    eagle_switch_dir <-   paste0("output/trio_phase_", chrom ,"/no_th/switch_errors/eagle/")
    shapeit_switch_dir <- paste0("output/trio_phase_", chrom, "/no_th/switch_errors/shapeit/")
    beagle_switch_dir <-  paste0("output/trio_phase_", chrom, "/no_th/switch_errors/beagle/")
    maf_df <- read_tsv(paste0("data/1kgp/chr", chrom, "/chr", chrom, "_freq.tsv"), col_names = c("chrom", "pos", "maf"))

    df_subj <- read_csv("data/1kgp/subject_info.csv")
    ped_df <- read_table("data/1kgp/1kGP.3202_samples.pedigree_info.txt")
    child_ids <- ped_df %>% filter(fatherID != "0" & motherID != "0") %>% pull(sampleID)
    df_subj_rel <- df_subj %>%
    filter(SAMPLE_NAME %in% child_ids)
    df_subj_rel$id <- 1:602

    # Lists to store results
    list_het <- vector(mode = "list", length = 602)
    list_error <- vector(mode = "list", length = 602)

    # work through single example

    for(i in 1:602){
    print(i)
        df_het <- load_het_proband(i, chrom)

        prop_res <- df_het %>%
            select(pos) %>%
            left_join(maf_df, by=c("pos" = "pos")) %>%
            mutate(maf_class = case_when(
                maf < 0.01 ~ "rare",
                maf >= 0.01 & maf < 0.05 ~ "uncommon",
                maf >= 0.05 ~ "common"
            )) %>%
            summarize(prop_rare = mean(maf_class == "rare"),
                    prop_uncommon = mean(maf_class == "uncommon"),
                    prop_common = mean(maf_class == "common")) %>%
            mutate(id = i) %>%
            left_join(df_subj_rel %>% select(id, SUPER), by=c("id" = "id"))

        list_het[[i]] <- prop_res

        # Get method switch locations and match to maf
        df_prop <- df_het %>%
            filter(beagle_switch == 1) %>%
            select(pos) %>%
            left_join(maf_df, by=c("pos" = "pos")) %>%
            mutate(maf_class = case_when(
                maf < 0.01 ~ "rare",
                maf >= 0.01 & maf < 0.05 ~ "uncommon",
                maf >= 0.05 ~ "common"
            )) %>%
            summarize(prop_rare = mean(maf_class == "rare"),
                    prop_uncommon = mean(maf_class == "uncommon"),
                    prop_common = mean(maf_class == "common")) %>%
            mutate(id = i, method = "beagle", error_type = "switch") %>%
            left_join(df_subj_rel %>% select(id, SUPER), by=c("id" = "id"))
        df_prop <- bind_rows(df_prop,
                            df_het %>%
                                filter(eagle_switch == 1) %>%
                                select(pos) %>%
                                left_join(maf_df, by=c("pos" = "pos")) %>%
                                mutate(maf_class = case_when(
                                    maf < 0.01 ~ "rare",
                                    maf >= 0.01 & maf < 0.05 ~ "uncommon",
                                    maf >= 0.05 ~ "common"
                                )) %>%
                                summarize(prop_rare = mean(maf_class == "rare"),
                                        prop_uncommon = mean(maf_class == "uncommon"),
                                        prop_common = mean(maf_class == "common")) %>%
                                mutate(id = i, method = "eagle", error_type = "switch") %>%
                                left_join(df_subj_rel %>% select(id, SUPER), by=c("id" = "id"))
        )
        df_prop <- bind_rows(df_prop,
                            df_het %>%
                                filter(shapeit_switch == 1) %>%
                                select(pos) %>%
                                left_join(maf_df, by=c("pos" = "pos")) %>%
                                mutate(maf_class = case_when(
                                    maf < 0.01 ~ "rare",
                                    maf >= 0.01 & maf < 0.05 ~ "uncommon",
                                    maf >= 0.05 ~ "common"
                                )) %>%
                                summarize(prop_rare = mean(maf_class == "rare"),
                                        prop_uncommon = mean(maf_class == "uncommon"),
                                        prop_common = mean(maf_class == "common")) %>%
                                mutate(id = i, method = "shapeit", error_type = "switch") %>%
                                left_join(df_subj_rel %>% select(id, SUPER), by=c("id" = "id"))
        )
        # flips
        df_prop <- bind_rows(df_prop,
                            df_het %>%
                                filter(beagle_flip == 1) %>%
                                select(pos) %>%
                                left_join(maf_df, by=c("pos" = "pos")) %>%
                                mutate(maf_class = case_when(
                                    maf < 0.01 ~ "rare",
                                    maf >= 0.01 & maf < 0.05 ~ "uncommon",
                                    maf >= 0.05 ~ "common"
                                )) %>%
                                summarize(prop_rare = mean(maf_class == "rare"),
                                        prop_uncommon = mean(maf_class == "uncommon"),
                                        prop_common = mean(maf_class == "common")) %>%
                                mutate(id = i, method = "beagle", error_type = "flip") %>%
                                left_join(df_subj_rel %>% select(id, SUPER), by=c("id" = "id"))
        )
        df_prop <- bind_rows(df_prop,
                            df_het %>%
                                filter(eagle_flip == 1) %>%
                                select(pos) %>%
                                left_join(maf_df, by=c("pos" = "pos")) %>%
                                mutate(maf_class = case_when(
                                    maf < 0.01 ~ "rare",
                                    maf >= 0.01 & maf < 0.05 ~ "uncommon",
                                    maf >= 0.05 ~ "common"
                                )) %>%
                                summarize(prop_rare = mean(maf_class == "rare"),
                                        prop_uncommon = mean(maf_class == "uncommon"),
                                        prop_common = mean(maf_class == "common")) %>%
                                mutate(id = i, method = "eagle", error_type = "flip") %>%
                                left_join(df_subj_rel %>% select(id, SUPER), by=c("id" = "id"))
        )
        df_prop <- bind_rows(df_prop,
                            df_het %>%
                                filter(shapeit_flip == 1) %>%
                                select(pos) %>%
                                left_join(maf_df, by=c("pos" = "pos")) %>%
                                mutate(maf_class = case_when(
                                    maf < 0.01 ~ "rare",
                                    maf >= 0.01 & maf < 0.05 ~ "uncommon",
                                    maf >= 0.05 ~ "common"
                                )) %>%
                                summarize(prop_rare = mean(maf_class == "rare"),
                                        prop_uncommon = mean(maf_class == "uncommon"),
                                        prop_common = mean(maf_class == "common")) %>%
                                mutate(id = i, method = "shapeit", error_type = "flip") %>%
                                left_join(df_subj_rel %>% select(id, SUPER), by=c("id" = "id"))
        )
        list_error[[i]] <- df_prop
    }

    df_het_res <- bind_rows(list_het)
    df_error_res <- bind_rows(list_error)

    write_csv(df_het_res, paste0("output/autosome_summaries/het_maf_props_chr", chrom, ".csv"))
    write_csv(df_error_res, paste0("output/autosome_summaries/error_maf_props_chr", chrom, ".csv"))
}