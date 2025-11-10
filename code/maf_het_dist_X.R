library(tidyverse)

load_het_sd <- function(id, base_dir = "/net/snowwhite/home/beckandy/research/phasing_clean/output/switch_errors/het_loc/annotated/"){
    f_name <- paste0(base_dir, "errors_", id, ".csv")
    read_csv(f_name, show_col_types = FALSE)
}
maf_df <- read_tsv(paste0("data/1kgp/chrX_freq.tsv"), col_names = c("chrom", "pos", "maf"))

df_subj <- data.frame(id = 1:1000,
    pop = rep(c("EUR", "AFR", "AMR", "EAS", "SAS"), each = 200))

# lists for results
list_het <- vector(mode = "list", length = 1000)
list_error <- vector(mode = "list", length = 1000)

for(id in 1:1000){
    print(id)
    df_het <- load_het_sd(id)
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
        mutate(id = id) %>% left_join(df_subj %>% select(id, pop), by=c("id" = "id"))

    list_het[[id]] <- prop_res

    prop_error <- df_het %>%
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
        mutate(id = id) %>%
        mutate(method = "beagle", error_type = "switch")
    prop_error <- bind_rows(prop_error,
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
                                mutate(id = id) %>%
                                mutate(method = "eagle", error_type = "switch"))
    prop_error <- bind_rows(prop_error,
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
                                mutate(id = id) %>%
                                mutate(method = "shapeit", error_type = "switch"))
    # flips
    prop_error <- bind_rows(prop_error,
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
                                mutate(id = id) %>%
                                mutate(method = "beagle", error_type = "flip"))
    prop_error <- bind_rows(prop_error,
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
                                mutate(id = id) %>%
                                mutate(method = "eagle", error_type = "flip"))
    prop_error <- bind_rows(prop_error,
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
                                mutate(id = id) %>%
                                mutate(method = "shapeit", error_type = "flip"))
    prop_error <- left_join(prop_error, df_subj %>% select(id, pop), by=c("id" = "id"))
    list_error[[id]] <- prop_error
}

df_het_res <- bind_rows(list_het)
df_error_res <- bind_rows(list_error)
write_csv(df_het_res, "output/autosome_summaries/het_maf_props_chrX.csv")
write_csv(df_error_res, "output/autosome_summaries/error_maf_props_chrX.csv")
