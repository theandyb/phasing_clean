source("code/figures/load_results.R")
source("code/figures/palettes.R")

library(cowplot)
library(rstatix)
filter <- dplyr::filter

# Load switch and flip error counts for synthetic diploids
df_sd_errors <- load_sd_errors()

# Load switch and flip counts in MB bins
df_sd_errors_mb <- load_sd_errors_mb()

# Load load_sd_maf_prop
df_sd_maf_prop <- load_sd_maf_prop()

# Using df_sd_errors, compute fractions of hets at CpG and fraction of errors at CpG so we can do paired T-tests
df_sd_errors_summary <- df_sd_errors %>%
    mutate(beagle_switch =  n_switch_cpg_beagle / n_switch_beagle,
           eagle_switch = n_switch_cpg_eagle / n_switch_eagle,
           shapeit_switch = n_switch_cpg_shapeit / n_switch_shapeit,
           beagle_flip = n_flip_cpg_beagle / n_flip_beagle,
           eagle_flip = n_flip_cpg_eagle / n_flip_eagle,
           shapeit_flip = n_flip_cpg_shapeit / n_flip_shapeit,
           background = n_het_cpg / n_hets) %>%
    select(pair_id, background, beagle_switch:background)

# compute differences between background and each method
df_sd_errors_summary <- df_sd_errors_summary %>%
    mutate(beagle_switch_diff = beagle_switch - background,
           eagle_switch_diff = eagle_switch - background,
           shapeit_switch_diff = shapeit_switch - background,
           beagle_flip_diff = beagle_flip - background,
           eagle_flip_diff = eagle_flip - background,
           shapeit_flip_diff = shapeit_flip - background)

# Reshape to long format to plot densities of *_diff values
df_sd_errors_summary_long <- df_sd_errors_summary %>%
    select(pair_id, ends_with("_diff")) %>%
    pivot_longer(cols = -pair_id, names_to = "method", values_to = "diff") %>%
    mutate(error_type = ifelse(grepl("switch", method), "Switch Error", "Flip Error"),
           method = case_when(
               grepl("beagle", method) ~ "Beagle",
               grepl("eagle", method) ~ "Eagle",
               grepl("shapeit", method) ~ "SHAPEIT"
            ))

# Plot density of differences, facet by error type
ggplot(df_sd_errors_summary_long, aes(x = diff, color = method, fill = method)) +
    geom_density(alpha = 0.4) +
    facet_wrap(~error_type) +
    scale_color_manual(values = ggpubfigs::friendly_pal("muted_nine")) +
    scale_fill_manual(values = ggpubfigs::friendly_pal("muted_nine")) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(x = "Difference in Fraction of Errors at CpG Sites vs. Background",
         y = "Density",
         title = "Density of Differences in Fraction of Errors at CpG Sites vs. Background") +
    theme_bw()

# Perform paired T-tests to see if differences are significant
df_sd_errors_summary_long %>%
    group_by(error_type, method) %>%
    summarize(t_test = list(t.test(x = diff, mu = 0, alternative = "greater") %>% broom::tidy())) %>%
    unnest_wider(t_test, names_sep = "_") %>%
    select(error_type, method, t_test_estimate, t_test_statistic, t_test_p.value, t_test_conf.low, t_test_conf.high) %>%
    knitr::kable()

t.test(x = df_sd_errors_summary$beagle_switch_diff, mu = 0, alternative = "greater")

# In df_sd_maf_prop, for each id I want to compute the difference of the value in the column "proportion" for rows with
# category = "het_rare" to rows with category = "beagle_switch_rare", "eagle_switch_rare", "shapeit_switch_rare", etc.
df_sd_maf_diff <- df_sd_maf_prop %>%
    pivot_wider(names_from = category, values_from = proportion) %>%
    mutate(beagle_switch_rare_diff = beagle_switch_rare - het_rare,
           eagle_switch_rare_diff = eagle_switch_rare - het_rare,
           shapeit_switch_rare_diff = shapeit_switch_rare - het_rare,
           beagle_flip_rare_diff = beagle_flip_rare - het_rare,
           eagle_flip_rare_diff = eagle_flip_rare - het_rare,
           shapeit_flip_rare_diff = shapeit_flip_rare - het_rare) %>%
    select(id, ends_with("_diff")) 

# For all columns ending in _diff, perform a one-sample t-test to see if the mean is significantly greater than 0
df_sd_maf_diff_long <- df_sd_maf_diff %>%
    pivot_longer(cols = -id, names_to = "method", values_to = "diff") %>%
    mutate(error_type = ifelse(grepl("switch", method), "Switch Error", "Flip Error"),
           method = case_when(
               grepl("beagle", method) ~ "Beagle",
               grepl("eagle", method) ~ "Eagle",
               grepl("shapeit", method) ~ "SHAPEIT"
            ))

df_sd_maf_diff_long %>%
    group_by(error_type, method) %>%
    summarize(t_test = list(t.test(x = diff, mu = 0, alternative = "greater") %>% broom::tidy())) %>%
    unnest_wider(t_test, names_sep = "_") %>%
    select(error_type, method, t_test_estimate, t_test_statistic, t_test_p.value, t_test_conf.low, t_test_conf.high) %>%
    knitr::kable()

##############################################################################################################################
##############################################################################################################################
# Autosome Comparisons
df_proband <- load_proband_errors()
df_proband_info <- get_1kgp_proband_info()

df_proband <- df_proband %>%
    left_join(df_proband_info, by = c("sample" = "id2"))

# translate names in df_proband to be easier to use
# n_switch_b -> n_switch_beagle, etc
df_proband <- df_proband %>%
    rename(n_switch_beagle = n_switch_b,
           n_switch_eagle = n_switch_e,
           n_switch_shapeit = n_switch_s,
           n_flip_beagle = n_flip_b,
           n_flip_eagle = n_flip_e,
           n_flip_shapeit = n_flip_s,
           n_switch_cpg_beagle = n_switch_cpg_b,
           n_switch_cpg_eagle = n_switch_cpg_e,
           n_switch_cpg_shapeit = n_switch_cpg_s,
           n_flip_cpg_beagle = n_flip_cpg_b,
           n_flip_cpg_eagle = n_flip_cpg_e,
           n_flip_cpg_shapeit = n_flip_cpg_s,
           pop = SUPER,
           subpop = POPULATION,
           pair_id = sample,
           n_hets = n_het)

df_sd_errors$chrom <- 23

shared_cols <- intersect(colnames(df_sd_errors), colnames(df_proband))
df_errors <- bind_rows(
    df_sd_errors %>% select(all_of(shared_cols)),
    df_proband %>% select(all_of(shared_cols))
)

# Average enrichment of errors at CpG sites by method and dataset and chromosome
df_errors_summary <- df_errors %>%
    mutate(beagle_switch =  n_switch_cpg_beagle / n_switch_beagle,
           eagle_switch = n_switch_cpg_eagle / n_switch_eagle,
           shapeit_switch = n_switch_cpg_shapeit / n_switch_shapeit,
           beagle_flip = n_flip_cpg_beagle / n_flip_beagle,
           eagle_flip = n_flip_cpg_eagle / n_flip_eagle,
           shapeit_flip = n_flip_cpg_shapeit / n_flip_shapeit,
           background = n_het_cpg / n_hets) %>%
    select(pair_id, chrom, background, beagle_switch:background)

# compute differences between background and each method
df_errors_summary <- df_errors_summary %>%
    mutate(beagle_switch_diff = beagle_switch - background,
           eagle_switch_diff = eagle_switch - background,
           shapeit_switch_diff = shapeit_switch - background,
           beagle_flip_diff = beagle_flip - background,
           eagle_flip_diff = eagle_flip - background,
           shapeit_flip_diff = shapeit_flip - background)

# Perform T-tests to see if differences are significant, grouped by chrom, method, and error type
df_errors_summary_long <- df_errors_summary %>%
    select(pair_id, chrom, ends_with("_diff")) %>%
    pivot_longer(cols = -c(pair_id, chrom), names_to = "method", values_to = "diff") %>%
    mutate(error_type = ifelse(grepl("switch", method), "Switch Error", "Flip Error"),
           method = case_when(
               grepl("beagle", method) ~ "Beagle",
               grepl("eagle", method) ~ "Eagle",
               grepl("shapeit", method) ~ "SHAPEIT"
           ))

df_errors_t_test <- df_errors_summary_long %>%
    group_by(chrom, error_type, method) %>%
    summarize(t_test = list(t.test(x = diff, mu = 0, alternative = "greater") %>% broom::tidy())) %>%
    unnest_wider(t_test, names_sep = "_") %>%
    select(chrom, error_type, method, t_test_estimate, t_test_statistic, t_test_p.value, t_test_conf.low, t_test_conf.high)

# Plot t_test_estimate by chromosome, faceted by error_type, colored by method
df_errors_t_test %>% 
    # Change chrom to factor, and relabel 23 as X
    mutate(chrom = factor(ifelse(chrom == 23, "X", as.character(chrom)))) %>%
    ggplot(aes(x = factor(chrom), y = t_test_estimate, color = method, group = method)) +
    geom_point(position = position_dodge(width = 0.5)) +
    # no error bars for now
    # geom_errorbar(aes(ymin = t_test_conf.low, ymax = t_test_conf.high), width = 0.2, position = position_dodge(width = 0.5)) +
    geom_line(position = position_dodge(width = 0.5)) +
    # facet by error_type into rows
    facet_wrap(~error_type, nrow = 2) +
    scale_color_manual(values = ggpubfigs::friendly_pal("muted_nine")) +
    labs(x = "Chromosome",
         y = "Mean Difference in Fraction of Errors at CpG Sites vs. Background") +
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    # Make fonts larger
    theme(text = element_text(size = 16)) +
    # change title of legend from method to Method
    labs(color = "Method")

# Compute mean of diff, first just by chrom and error_type
df_errors_summary_long %>%
    group_by(chrom, error_type) %>%
    summarize(mean_diff = mean(diff, na.rm = TRUE),
              sd_diff = sd(diff, na.rm = TRUE),
              n = n(),
              se_diff = sd_diff / sqrt(n)) %>%
    arrange(chrom, error_type) %>%
    knitr::kable()

# Exclude chromosome 23 and compute mean of diff by error_type (across all autosomes and methods)
df_errors_summary_long %>%
    filter(chrom != 23) %>%
    group_by(error_type) %>%
    summarize(mean_diff = mean(diff, na.rm = TRUE),
              sd_diff = sd(diff, na.rm = TRUE),
              n = n(),
              se_diff = sd_diff / sqrt(n)) %>%
    arrange(error_type) %>%
    knitr::kable()

# Same, but only for chromosome 23
df_errors_summary_long %>%
    filter(chrom == 23) %>%
    group_by(error_type) %>%
    summarize(mean_diff = mean(diff, na.rm = TRUE),
              sd_diff = sd(diff, na.rm = TRUE),
              n = n(),
              se_diff = sd_diff / sqrt(n)) %>%
    arrange(error_type) %>%
    knitr::kable()
