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
# each sample has three rows, one for each MAF category (common, uncommon, rare)

# Scatter plot of total errors, pairwise comparisons of methods,colored by population
# First: Beagle and SHAPEIT

pal_choice <- "vibrant_seven"

df_summary <- df_sd_errors %>%
    mutate(Beagle = (n_other_beagle + n_flip_beagle) / n_hets,
              Eagle = (n_other_eagle + n_flip_eagle) / n_hets,
              SHAPEIT = (n_other_shapeit + n_flip_shapeit) / n_hets)

scatter_panel <- function(df_summary, pal_choice = "vibrant_seven"){
    p1<- df_summary %>%
        ggplot(aes(x= Beagle, y = SHAPEIT, color = pop)) + 
        geom_point() +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
        theme_classic() +
        xlab("Beagle") +
        ylab("SHAPEIT") +
        scale_color_manual(values = friendly_pal(pal_choice)) +
        theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
        theme(legend.text = element_text(size = 14), legend.title = element_text(size = 16)) + 
        labs(color = "Population") + # change legend to have squares instead of circles
        guides(color = guide_legend(override.aes = list(size = 5, shape = 15)))

    # Extract the legend from the first plot
    legend <- cowplot::get_legend(p1)

    # Now we can remove the legend from the first plot
    p1 <- p1 + theme(legend.position = "none") # No legend for this plot

# Second: Beagle and Eagle, no legend
    p2 <- df_summary %>%
        ggplot(aes(x = Beagle, y = Eagle, color = pop)) +
        geom_point() +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
        theme_classic() +
        xlab("Beagle") +
        ylab("Eagle") +
        scale_color_manual(values = friendly_pal(pal_choice)) +
        theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
        theme(legend.text = element_text(size = 14), legend.title = element_text(size = 16)) + 
        labs(color = "Population") +
        theme(legend.position = "none")

# Third: Eagle and SHAPEIT, no legend
    p3 <- df_summary %>%
        ggplot(aes(x = SHAPEIT, y = Eagle, color = pop)) +
        geom_point() +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
        theme_classic() +
        xlab("Eagle") +
        ylab("SHAPEIT") +
        scale_color_manual(values = friendly_pal(pal_choice)) +
        theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
        theme(legend.text = element_text(size = 14), legend.title = element_text(size = 16)) + 
        labs(color = "Population") +
        theme(legend.position = "none")


    combined_plot <- plot_grid(p1, p2, p3, legend, 
                                ncol = 2, nrow = 2, 
                                rel_widths = c(1, 1), 
                                rel_heights = c(1, 1),
                                labels = c("A", "B", "C", ""),
                                label_size = 16)
    combined_plot
}

scatter_panel(df_summary, pal_choice)

# Now just switch errors
df_summary <- df_sd_errors %>%
    mutate(Beagle = n_other_beagle / n_hets,
              Eagle = n_other_eagle / n_hets,
              SHAPEIT = n_other_shapeit / n_hets)
scatter_panel(df_summary, pal_choice)

# same for flips
df_summary <- df_sd_errors %>%
    mutate(Beagle = n_flip_beagle / n_hets,
              Eagle = n_flip_eagle / n_hets,
              SHAPEIT = n_flip_shapeit / n_hets)
scatter_panel(df_summary, pal_choice)

# Density plot of switch rates by phasing method and population
df_sd_errors %>%
    mutate(Beagle = n_other_beagle / n_hets,
              Eagle = n_other_eagle / n_hets,
              SHAPEIT = n_other_shapeit / n_hets) %>%
    select(pair_id, pop, Beagle:SHAPEIT) %>%
    pivot_longer(cols = Beagle:SHAPEIT, names_to = "Method", values_to = "Switch_Rate") %>%
    ggplot(aes(x = Switch_Rate, fill = Method)) +
    geom_density() +
    facet_grid(rows = vars(Method)) + 
    guides(fill = "none") +
    theme_classic() +
    xlab("Switch Error Rate") +
    ylab("Density") +
    scale_fill_manual(values = friendly_pal(pal_choice)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Ensure x axis text is big enough
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
    theme(strip.text = element_text(size = 14))

# Get switch and flip rates by method and plot densities

df_sd_errors %>%
    mutate(Beagle_Switch = n_other_beagle / n_hets,
              Eagle_Switch = n_other_eagle / n_hets,
              SHAPEIT_Switch = n_other_shapeit / n_hets,
              Beagle_Flip = n_flip_beagle / n_hets,
              Eagle_Flip = n_flip_eagle / n_hets,
              SHAPEIT_Flip = n_flip_shapeit / n_hets) %>%
    select(pair_id, pop, Beagle_Switch:SHAPEIT_Flip) %>% 
    pivot_longer(cols = Beagle_Switch:SHAPEIT_Flip, names_to = "Metric", values_to = "Rate") %>%
    separate(Metric, into = c("Method", "Error_Type"), sep = "_") %>% 
    mutate(Error_Type = factor(Error_Type, levels = c("Switch", "Flip"))) %>%
    ggplot(aes(x = Rate, fill = Method)) + 
    geom_density() +
    facet_grid(cols = vars(Error_Type), rows = vars(Method)) +
    theme_classic() +
    xlab("Error Rate") +
    ylab("Density") +
    scale_fill_manual(values = friendly_pal(pal_choice)) +
    guides(fill = "none") + #rotate x axis text
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Ensure x axis text is big enough
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
    theme(strip.text = element_text(size = 14))


# Get switch and flip rates by method and plot boxplots

df_sd_errors %>%
    mutate(Beagle_Switch = n_other_beagle / n_hets,
              Eagle_Switch = n_other_eagle / n_hets,
              SHAPEIT_Switch = n_other_shapeit / n_hets,
              Beagle_Flip = n_flip_beagle / n_hets,
              Eagle_Flip = n_flip_eagle / n_hets,
              SHAPEIT_Flip = n_flip_shapeit / n_hets) %>%
    select(pair_id, pop, Beagle_Switch:SHAPEIT_Flip) %>% 
    pivot_longer(cols = Beagle_Switch:SHAPEIT_Flip, names_to = "Metric", values_to = "Rate") %>%
    separate(Metric, into = c("Method", "Error_Type"), sep = "_") %>% 
    mutate(Error_Type = factor(Error_Type, levels = c("Switch", "Flip"))) %>% 
    # recode Switch to Switch Error Rate and Flip to Flip Error Rate
    mutate(Error_Type = recode(Error_Type, "Switch" = "Switch Error Rate", "Flip" = "Flip Error Rate")) %>%
    ggplot(aes(y = Rate, x = Method , fill = Method)) + 
    geom_boxplot(outlier.size = 0.5) + 
    facet_grid(cols = vars(Error_Type)) + # hide legend
    guides(fill = "none") +
    theme_classic() +
    ylab("Error Rate") + # ensure x axis text is big enough
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
    theme(strip.text = element_text(size = 14)) +
    scale_fill_manual(values = friendly_pal(pal_choice))

# Violin plot of switch rates by phasing method
df_sd_errors %>%
    mutate(Beagle_Switch = n_other_beagle / n_hets,
              Eagle_Switch = n_other_eagle / n_hets,
              SHAPEIT_Switch = n_other_shapeit / n_hets,
              Beagle_Flip = n_flip_beagle / n_hets,
              Eagle_Flip = n_flip_eagle / n_hets,
              SHAPEIT_Flip = n_flip_shapeit / n_hets) %>%
    select(pair_id, pop, Beagle_Switch:SHAPEIT_Flip) %>% 
    pivot_longer(cols = Beagle_Switch:SHAPEIT_Flip, names_to = "Metric", values_to = "Rate") %>%
    separate(Metric, into = c("Method", "Error_Type"), sep = "_") %>% 
    mutate(Error_Type = factor(Error_Type, levels = c("Switch", "Flip"))) %>% 
    # recode Switch to Switch Error Rate and Flip to Flip Error Rate
    mutate(Error_Type = recode(Error_Type, "Switch" = "Switch Error Rate", "Flip" = "Flip Error Rate")) %>%
    ggplot(aes(y = Rate, x = Method , fill = Method)) +
    geom_violin() +
    #geom_boxplot(width = 0.1, colour = "white", alpha = 0.5) +
    facet_grid(cols = vars(Error_Type)) + # hide legend
    guides(fill = "none") +
    theme_classic() +
    xlab("") +
    ylab("Error Rate") + # ensure x axis text is big enough
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
    theme(strip.text = element_text(size = 14)) + 
    scale_fill_manual(values = friendly_pal(pal_choice))

# Violin plot with boxplot

df_sd_errors %>%
    mutate(Beagle_Switch = n_other_beagle / n_hets,
              Eagle_Switch = n_other_eagle / n_hets,
              SHAPEIT_Switch = n_other_shapeit / n_hets,
              Beagle_Flip = n_flip_beagle / n_hets,
              Eagle_Flip = n_flip_eagle / n_hets,
              SHAPEIT_Flip = n_flip_shapeit / n_hets) %>%
    select(pair_id, pop, Beagle_Switch:SHAPEIT_Flip) %>% 
    pivot_longer(cols = Beagle_Switch:SHAPEIT_Flip, names_to = "Metric", values_to = "Rate") %>%
    separate(Metric, into = c("Method", "Error_Type"), sep = "_") %>% 
    mutate(Error_Type = factor(Error_Type, levels = c("Switch", "Flip"))) %>% 
    # recode Switch to Switch Error Rate and Flip to Flip Error Rate
    mutate(Error_Type = recode(Error_Type, "Switch" = "Switch Error Rate", "Flip" = "Flip Error Rate")) %>%
    ggplot(aes(y = Rate, x = Method , fill = Method)) +
    geom_violin(width = 1) +
    geom_boxplot(width = 0.2, colour = "white", alpha = 0.5) +
    facet_grid(cols = vars(Error_Type)) + # hide legend
    guides(fill = "none") +
    theme_bw() + # remove grid lines
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("") +
    ylab("Error Rate") + # ensure x axis text is big enough
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
    theme(strip.text = element_text(size = 14)) + 
    scale_fill_manual(values = friendly_pal(pal_choice))


# Compute enrichment of switch and flip errors at CpG sites
df_sd_errors %>% 
    mutate(Beagle_Switch = (n_other_cpg_beagle / n_het_cpg)/(n_other_beagle / n_hets),
              Eagle_Switch = (n_other_cpg_eagle / n_het_cpg)/(n_other_eagle / n_hets),
              SHAPEIT_Switch = (n_other_cpg_shapeit / n_het_cpg)/(n_other_shapeit / n_hets),
              Beagle_Flip = (n_flip_cpg_beagle / n_het_cpg)/(n_flip_beagle / n_hets),
              Eagle_Flip = (n_flip_cpg_eagle / n_het_cpg)/(n_flip_eagle / n_hets),
              SHAPEIT_Flip = (n_flip_cpg_shapeit / n_het_cpg)/(n_flip_shapeit / n_hets)) %>%
    select(pair_id, pop, Beagle_Switch:SHAPEIT_Flip) %>%
    pivot_longer(cols = Beagle_Switch:SHAPEIT_Flip, names_to = "Metric", values_to = "Enrichment") %>%
    separate(Metric, into = c("Method", "Error_Type"), sep = "_") %>%
    mutate(Error_Type = factor(Error_Type, levels = c("Switch", "Flip"))) %>%
    # recode Switch to Switch Enrichment and Flip to Flip Enrichment
    mutate(Error_Type = recode(Error_Type, "Switch" = "Switch Enrichment", "Flip" = "Flip Enrichment")) %>%
    ggplot(aes(x = Method, y = Enrichment, fill = Method)) +
    geom_boxplot(outlier.size = 0.5) +
    facet_grid(cols = vars(Error_Type)) +
    guides(fill = "none") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ylab("CpG Enrichment") +
    xlab("") + # ensure x axis text is big enough
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
    theme(strip.text = element_text(size = 14)) +
    scale_fill_manual(values = friendly_pal(pal_choice)) +
    guides(fill = "none") + # ensure legend text is big enough
    theme(legend.text = element_text(size = 14)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") # add dashed horizontal line at 1

# Compute mean switch and flip rates in each bin
mb_rate_plot <- df_sd_errors_mb %>%
    mutate(Beagle_Switch = beagle_switch / n_het,
              Eagle_Switch = eagle_switch / n_het,
              SHAPEIT_Switch = shapeit_switch / n_het,
              Beagle_Flip = beagle_flip / n_het,
              Eagle_Flip = eagle_flip / n_het,
              SHAPEIT_Flip = shapeit_flip / n_het) %>%
    select(bin, Beagle_Switch:SHAPEIT_Flip) %>%
    group_by(bin) %>%
    summarise(across(Beagle_Switch:SHAPEIT_Flip, mean)) %>%
    pivot_longer(cols = Beagle_Switch:SHAPEIT_Flip, names_to = "Metric", values_to = "Rate") %>%
    separate(Metric, into = c("Method", "Error_Type"), sep = "_") %>%
    mutate(Error_Type = factor(Error_Type, levels = c("Switch", "Flip"))) %>%
    # recode Switch to Switch Error Rate and Flip to Flip Error Rate
    mutate(Error_Type = recode(Error_Type, "Switch" = "Switch Error Rate", "Flip" = "Flip Error Rate")) %>%
    ggplot(aes(x = bin, y = Rate, color = Method)) +
    geom_point() +
    geom_line() +
    facet_grid(rows = vars(Error_Type)) + 
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("MB Bin") +
    ylab("Error Rate") +
    scale_color_manual(values = friendly_pal(pal_choice)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
    theme(strip.text = element_text(size = 14)) +
    theme(legend.position = "top", legend.box = "horizontal") + # remove legend title
    theme(legend.text = element_text(size = 14)) +
    theme(legend.title = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 5, 
        shape = 22, 
        linetype = "blank",
        fill = friendly_pal(pal_choice)[1:3]),
        title = "")) # make legend points larger
mb_rate_plot

mb_rate_plot + scale_y_log10()

# For each bin plot the ratio of eagle_switch to beagle_switch
df_sd_errors_mb %>%
    group_by(bin) %>%
    summarise(eagle_switch = sum(eagle_switch), 
              beagle_switch = sum(beagle_switch),
              n_het = sum(n_het)) %>%
    mutate(Eagle_rate = eagle_switch / n_het,
              Beagle_rate = beagle_switch / n_het) %>%
    mutate(Eagle_Beagle = Eagle_rate / Beagle_rate) %>%
    ggplot(aes(x = bin, y = Eagle_Beagle)) +
    geom_point() +
    theme_classic() +
    xlab("MB Bin") +
    ylab("Eagle/Beagle Switch Ratio") + # larger axis text
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black")

# Same for SHAPEIT and Eagle
df_sd_errors_mb %>%
    group_by(bin) %>%
    summarize(SHAPEIT_switch = sum(shapeit_switch), 
              eagle_switch = sum(eagle_switch),
              n_het = sum(n_het)) %>%
    mutate(SHAPEIT_rate = SHAPEIT_switch / n_het,
              Eagle_rate = eagle_switch / n_het) %>%
    mutate(Eagle_SHAPEIT =  Eagle_rate / SHAPEIT_rate) %>%
    ggplot(aes(x = bin, y = Eagle_SHAPEIT)) +
    geom_point() +
    theme_classic() +
    xlab("MB Bin") +
    ylab("Eagle/SHAPEIT Switch Ratio") + # larger axis text
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black")

# Ratio of SHAPEIT flip rates to Beagle flip rates in MB bins
df_sd_errors_mb %>%
    group_by(bin) %>%
    summarize(SHAPEIT_flip = sum(shapeit_flip), 
              Beagle_flip = sum(beagle_flip),
              n_het = sum(n_het)) %>%
    mutate(SHAPEIT_rate = SHAPEIT_flip / n_het,
              Beagle_rate = Beagle_flip / n_het) %>%
    mutate(SHAPEIT_Beagle = SHAPEIT_rate / Beagle_rate) %>%
    ggplot(aes(x = bin, y = SHAPEIT_Beagle)) +
    geom_point() +
    theme_classic() +
    xlab("MB Bin") +
    ylab("SHAPEIT/Beagle Flip Ratio") + # larger axis text
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black")

# SHAPEIT Eagle Flip Ratio
df_sd_errors_mb %>%
    group_by(bin) %>%
    summarize(SHAPEIT_flip = sum(shapeit_flip), 
              Eagle_flip = sum(eagle_flip),
              n_het = sum(n_het)) %>%
    mutate(SHAPEIT_rate = SHAPEIT_flip / n_het,
              Eagle_rate = Eagle_flip / n_het) %>%
    mutate(SHAPEIT_Eagle = SHAPEIT_rate / Eagle_rate) %>%
    ggplot(aes(x = bin, y = SHAPEIT_Eagle)) +
    geom_point() +
    theme_classic() +
    xlab("MB Bin") +
    ylab("SHAPEIT/Eagle Flip Ratio") + # larger axis text
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black")


prop_maf_fig <- function(cat_type, title_text, maf_prop_df){
  df <- maf_prop_df %>%
    filter(str_detect(category, paste0(cat_type, "_[rcu]"))) %>%
    mutate(category = factor(category)) %>%
    mutate(category = fct_recode(category, "Rare" = paste0(cat_type, "_rare"),
                               "Uncommon" = paste0(cat_type, "_uncommon"),
                               "Common" = paste0(cat_type, "_common"))) %>%
    mutate(category = fct_relevel(category, c("Rare", "Uncommon", "Common")))

  p <- df %>%
    ggplot(aes(x = category, y = proportion, color = SP)) +
    geom_boxplot() +
    ggtitle(paste0(title_text, " in MAF Bins")) +
    xlab("") +
    ylab("Proportion of Sites") +
    theme_classic() +
    labs(color = "Population") +
    scale_color_manual(values = cbPalette)

#   tab_prop <- df %>%
#     group_by(SP, category) %>%
#     summarize(mean_prop = mean(proportion)) %>%
#     pivot_wider(id_cols = SP, names_from = category, values_from = mean_prop)
  
  p
}
prop_maf_fig("het", "Heterozygous Sites", df_sd_maf_prop)

# For each id, we want to compute the ratio of method_switch_category to het_category
df_sd_maf_prop %>%
    pivot_wider(id_cols = c(id, SP, POP), names_from = category, values_from = proportion) %>%
    mutate(Beagle_Switch_Common_Ratio = beagle_switch_common / het_common,
           Beagle_Switch_Uncommon_Ratio = beagle_switch_uncommon / het_uncommon,
           Beagle_Switch_Rare_Ratio = beagle_switch_rare / het_rare,
           Eagle_Switch_Common_Ratio = eagle_switch_common / het_common,
           Eagle_Switch_Uncommon_Ratio = eagle_switch_uncommon / het_uncommon,
           Eagle_Switch_Rare_Ratio = eagle_switch_rare / het_rare,
           Shapeit_Switch_Common_Ratio = shapeit_switch_common / het_common,
           Shapeit_Switch_Uncommon_Ratio = shapeit_switch_uncommon / het_uncommon,
           Shapeit_Switch_Rare_Ratio = shapeit_switch_rare / het_rare,
           Beagle_Flip_Common_Ratio = beagle_flip_common / het_common,
           Beagle_Flip_Uncommon_Ratio = beagle_flip_uncommon / het_uncommon,
           Beagle_Flip_Rare_Ratio = beagle_flip_rare / het_rare,
           Eagle_Flip_Common_Ratio = eagle_flip_common / het_common,
           Eagle_Flip_Uncommon_Ratio = eagle_flip_uncommon / het_uncommon,
           Eagle_Flip_Rare_Ratio = eagle_flip_rare / het_rare,
           Shapeit_Flip_Common_Ratio = shapeit_flip_common / het_common,
           Shapeit_Flip_Uncommon_Ratio = shapeit_flip_uncommon / het_uncommon,
           Shapeit_Flip_Rare_Ratio = shapeit_flip_rare / het_rare) %>%
    select(id, SP, POP, Beagle_Switch_Common_Ratio:Shapeit_Flip_Rare_Ratio) %>%
    pivot_longer(cols = Beagle_Switch_Common_Ratio:Shapeit_Flip_Rare_Ratio, 
                 names_to = c("Method", "Category", "MAF", "Blank"), 
                    names_sep = "_") %>%
    select(-Blank) %>%
    mutate(Method = factor(Method, levels = c("Beagle", "Eagle", "Shapeit"))) %>%
    mutate(Category = factor(Category, levels = c("Switch", "Flip"))) %>%
    mutate(MAF = factor(MAF, levels = c("Rare", "Uncommon", "Common"))) %>%
    ggplot(aes(x = MAF, y = value, fill = Method)) +
    geom_boxplot() +
    facet_grid(cols = vars(Category)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("MAF Category") +
    ylab("Ratio of MAF Proportions \n Log10 Scale") +
    scale_fill_manual(values = friendly_pal(pal_choice)[1:10]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
    theme(strip.text = element_text(size = 22)) +
    theme(legend.position = "top",
        legend.box = "horizontal",
        legend.text = element_text(size = 18)) + # add dashed horizontal line at 1
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
    guides(fill = guide_legend(override.aes = list(size = 5, 
        shape = 20, 
        linetype = "blank",
        fill = friendly_pal(pal_choice)[1:3]),
        title = "")) +
    scale_y_log10()

# Compute average enrichment in rare MAF category for switch and flip errors
df_sd_maf_prop %>%
    pivot_wider(id_cols = c(id, SP), names_from = category, values_from = proportion) %>%
    mutate(Beagle_Switch_Rare_Ratio = beagle_switch_rare / het_rare,
           Eagle_Switch_Rare_Ratio = eagle_switch_rare / het_rare,
           Shapeit_Switch_Rare_Ratio = shapeit_switch_rare / het_rare,
           Beagle_Flip_Rare_Ratio = beagle_flip_rare / het_rare,
           Eagle_Flip_Rare_Ratio = eagle_flip_rare / het_rare,
           Shapeit_Flip_Rare_Ratio = shapeit_flip_rare / het_rare) %>%
    select(id, SP, Beagle_Switch_Rare_Ratio:Shapeit_Flip_Rare_Ratio) %>%
    pivot_longer(cols = Beagle_Switch_Rare_Ratio:Shapeit_Flip_Rare_Ratio, 
                 names_to = c("Method", "Category", "Blank", "Ratio"), 
                    names_sep = "_") %>%
    select(-Blank, -Ratio) %>%
    group_by(Method, Category) %>%
    summarise(mean_enrichment = mean(value))

# plot enrichment in rare MAF category for switch and flip errors
df_sd_maf_prop %>%
    pivot_wider(id_cols = c(id, SP), names_from = category, values_from = proportion) %>%
    mutate(Beagle_Switch_Rare_Ratio = beagle_switch_rare / het_rare,
           Eagle_Switch_Rare_Ratio = eagle_switch_rare / het_rare,
           Shapeit_Switch_Rare_Ratio = shapeit_switch_rare / het_rare,
           Beagle_Flip_Rare_Ratio = beagle_flip_rare / het_rare,
           Eagle_Flip_Rare_Ratio = eagle_flip_rare / het_rare,
           Shapeit_Flip_Rare_Ratio = shapeit_flip_rare / het_rare) %>%
    select(id, SP, Beagle_Switch_Rare_Ratio:Shapeit_Flip_Rare_Ratio) %>%
    pivot_longer(cols = Beagle_Switch_Rare_Ratio:Shapeit_Flip_Rare_Ratio, 
                 names_to = c("Method", "Category", "Blank", "Ratio"), 
                    names_sep = "_") %>%
    select(-Blank, -Ratio) %>%
    ggplot(aes(x = Method, y = value, fill = Method)) +
    geom_boxplot() +
    facet_grid(cols = vars(Category)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("Phasing Method") +
    ylab("Enrichment in Rare MAF Category") +
    scale_fill_manual(values = friendly_pal(pal_choice)[1:10]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
    theme(strip.text = element_text(size = 22)) +
    theme(legend.position = "top",
        legend.box = "horizontal",
        legend.text = element_text(size = 18)) + # add dashed horizontal line at 1
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
    guides(fill = guide_legend(override.aes = list(size = 5, 
        shape = 20, 
        linetype = "blank",
        fill = friendly_pal(pal_choice)[1:3]),
        title = ""))

# Find values that are infinite after log transform
df_sd_maf_prop %>%
    pivot_wider(id_cols = c(id, SP, POP), names_from = category, values_from = proportion) %>%
    mutate(Beagle_Switch_Common_Ratio = beagle_switch_common / het_common,
           Beagle_Switch_Uncommon_Ratio = beagle_switch_uncommon / het_uncommon,
           Beagle_Switch_Rare_Ratio = beagle_switch_rare / het_rare,
           Eagle_Switch_Common_Ratio = eagle_switch_common / het_common,
           Eagle_Switch_Uncommon_Ratio = eagle_switch_uncommon / het_uncommon,
           Eagle_Switch_Rare_Ratio = eagle_switch_rare / het_rare,
           Shapeit_Switch_Common_Ratio = shapeit_switch_common / het_common,
           Shapeit_Switch_Uncommon_Ratio = shapeit_switch_uncommon / het_uncommon,
           Shapeit_Switch_Rare_Ratio = shapeit_switch_rare / het_rare,
           Beagle_Flip_Common_Ratio = beagle_flip_common / het_common,
           Beagle_Flip_Uncommon_Ratio = beagle_flip_uncommon / het_uncommon,
           Beagle_Flip_Rare_Ratio = beagle_flip_rare / het_rare,
           Eagle_Flip_Common_Ratio = eagle_flip_common / het_common,
           Eagle_Flip_Uncommon_Ratio = eagle_flip_uncommon / het_uncommon,
           Eagle_Flip_Rare_Ratio = eagle_flip_rare / het_rare,
           Shapeit_Flip_Common_Ratio = shapeit_flip_common / het_common,
           Shapeit_Flip_Uncommon_Ratio = shapeit_flip_uncommon / het_uncommon,
           Shapeit_Flip_Rare_Ratio = shapeit_flip_rare / het_rare) %>%
    select(id, SP, POP, Beagle_Switch_Common_Ratio:Shapeit_Flip_Rare_Ratio) %>%
    pivot_longer(cols = Beagle_Switch_Common_Ratio:Shapeit_Flip_Rare_Ratio, 
                 names_to = c("Method", "Category", "MAF", "Blank"), 
                    names_sep = "_") %>%
    select(-Blank) %>%
    mutate(Method = factor(Method, levels = c("Beagle", "Eagle", "Shapeit"))) %>%
    mutate(Category = factor(Category, levels = c("Switch", "Flip"))) %>%
    mutate(MAF = factor(MAF, levels = c("Rare", "Uncommon", "Common"))) %>%
    mutate(value2 = log10(value)) %>%
    filter(is.infinite(value2)) 

# Calculate proportion of heterozygous sites that are CPG sites
df_sd_errors %>%
    mutate(prop_het_cpg = n_het_cpg / n_hets) %>%
    select(pair_id, pop, prop_het_cpg) %>%
    group_by(pop) %>%
    summarise(mean_prop_het_cpg = mean(prop_het_cpg))

df_sd_errors %>%
    mutate(prop_het_cpg = n_het_cpg / n_hets) %>%
    pull(prop_het_cpg) %>%
    summary()

# Recombination Rate Plots
df_rr_1kgp_x <- read_tsv("data/1kgp/chrX_recomb.bed", col_names = c("chrom", "start", "end", "cM_MB"))
df_rr_1kgp_x_mb <- df_rr_1kgp_x %>%
    mutate(binMB = floor(start / 1e6), width = end - start) %>%
    group_by(binMB) %>%
    summarize(rr_MB_1kgp = sum(cM_MB * width) / sum(width)) %>%
    ungroup()
het_dir <- "output/switch_errors/het_loc/annotated/"

load_het_pos <- function(id, base_dir){
    f_name <- paste0(base_dir, "errors_", id, ".csv")
    read_csv(f_name, show_col_types = FALSE)
}
compute_summary_mb <- function(id, base_dir, df_rr_1kgp_x_mb){
    df_het <- load_het_pos(id, base_dir) %>%
        mutate(binMB = floor(pos / 1e6)) %>%
        left_join(df_rr_1kgp_x_mb, by = "binMB")
    # Compute switch and flip counts per MB bin
    df_het %>%
        group_by(binMB) %>%
        summarize(total_hets = n(),
                  cpg_hets = sum(cpg),
                  beagle_switches = sum(beagle_switch),
                  beagle_flips = sum(beagle_flip),
                  eagle_switches = sum(eagle_switch),
                  eagle_flips = sum(eagle_flip),
                  shapeit_switches = sum(shapeit_switch),
                  shapeit_flips = sum(shapeit_flip),
                  avg_rr_MB_1kgp = mean(rr_MB_1kgp, na.rm= TRUE)) %>%
        mutate(id = id) %>%
        ungroup()
}
pop_summary_MB <- function(ids, base_dir, df_rr_1kgp_x_mb, pop_name){
    res_list <- lapply(ids, compute_summary_mb, base_dir = base_dir, df_rr_1kgp_x_mb = df_rr_1kgp_x_mb)
    df_summary <- bind_rows(res_list)

    df_summary %>%
        group_by(binMB) %>%
        summarize(total_hets = sum(total_hets),
                  cpg_hets = sum(cpg_hets),
                  beagle_switches = sum(beagle_switches),
                  beagle_flips = sum(beagle_flips),
                  eagle_switches = sum(eagle_switches),
                  eagle_flips = sum(eagle_flips),
                  shapeit_switches = sum(shapeit_switches),
                  shapeit_flips = sum(shapeit_flips),
                  avg_rr_MB_1kgp = mean(avg_rr_MB_1kgp, na.rm = TRUE)) %>%
        mutate(beagle_switch_rate = beagle_switches / total_hets,
               beagle_flip_rate = beagle_flips / total_hets,
               eagle_switch_rate = eagle_switches / total_hets,
               eagle_flip_rate = eagle_flips / total_hets,
               shapeit_switch_rate = shapeit_switches / total_hets,
               shapeit_flip_rate = shapeit_flips / total_hets) %>%
        ungroup() %>%
        mutate(population = pop_name)
}
df_summary <- vector("list", 5)
# set names. of list elements to population names
names(df_summary) <- c("EUR", "AFR", "AMR", "EAS", "SAS")
df_summary[["EUR"]] <- pop_summary_MB(1:200, het_dir, df_rr_1kgp_x_mb, "EUR")
df_summary[["AFR"]] <- pop_summary_MB(201:400, het_dir, df_rr_1kgp_x_mb, "AFR")
df_summary[["AMR"]] <- pop_summary_MB(401:600, het_dir, df_rr_1kgp_x_mb, "AMR")
df_summary[["EAS"]] <- pop_summary_MB(601:800, het_dir, df_rr_1kgp_x_mb, "EAS")
df_summary[["SAS"]] <- pop_summary_MB(801:1000, het_dir, df_rr_1kgp_x_mb, "SAS")

# Combine all df_summary dataframes for a final comprehensive analysis
df_summary_combined <- bind_rows(df_summary[["EUR"]] %>% mutate(population = "EUR"),
                                 df_summary[["AFR"]] %>% mutate(population = "AFR"),
                                 df_summary[["AMR"]] %>% mutate(population = "AMR"),
                                 df_summary[["EAS"]] %>% mutate(population = "EAS"),
                                 df_summary[["SAS"]] %>% mutate(population = "SAS"))

df_summary_all <- df_summary_combined %>%
    group_by(binMB) %>%
    summarize(total_hets = sum(total_hets),
              cpg_hets = sum(cpg_hets),
              beagle_switches = sum(beagle_switches),
              beagle_flips = sum(beagle_flips),
              eagle_switches = sum(eagle_switches),
              eagle_flips = sum(eagle_flips),
              shapeit_switches = sum(shapeit_switches),
              shapeit_flips = sum(shapeit_flips),
              avg_rr_MB_1kgp = mean(avg_rr_MB_1kgp, na.rm = TRUE)) %>%
    mutate(beagle_switch_rate = beagle_switches / total_hets,
           beagle_flip_rate = beagle_flips / total_hets,
           eagle_switch_rate = eagle_switches / total_hets,
           eagle_flip_rate = eagle_flips / total_hets,
           shapeit_switch_rate = shapeit_switches / total_hets,
           shapeit_flip_rate = shapeit_flips / total_hets,
           cpg_het_rate = cpg_hets / total_hets) %>%
    ungroup() 
cor_beagle_all <- cor(df_summary_all$avg_rr_MB_1kgp, df_summary_all$beagle_switch_rate, use = "complete.obs")
cor_eagle_all <- cor(df_summary_all$avg_rr_MB_1kgp, df_summary_all$eagle_switch_rate, use = "complete.obs")
cor_shapeit_all <- cor(df_summary_all$avg_rr_MB_1kgp, df_summary_all$shapeit_switch_rate, use = "complete.obs")

color_beagle <- friendly_pal(pal_choice)[1]
color_eagle <- friendly_pal(pal_choice)[2]
color_shapeit <- friendly_pal(pal_choice)[3]

df_summary_all %>%
    select(binMB,avg_rr_MB_1kgp, contains("switch_rate")) %>%
    pivot_longer(cols = contains("switch_rate"), names_to = "method", values_to = "switch_rate") %>%
    mutate(method = case_when(
        grepl("beagle", method) ~ "Beagle",
        grepl("eagle", method) ~ "Eagle",
        grepl("shapeit", method) ~ "Shapeit"
    )) %>%
    ggplot(aes(x = avg_rr_MB_1kgp, y = switch_rate, color = method)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    theme_classic() +
    labs(x = "Average Recombination Rate (cM/Mb)",
         y = "Switch Rate",
         color = "Phasing Method") +
    scale_color_manual(values = friendly_pal(pal_choice)) +
    annotate("text", x = 4, y = 0.007, label = paste0("Beagle cor = ", round(cor_beagle_all, 3)), color = color_beagle) +
    annotate("text", x = 4, y = 0.0065, label = paste0("Eagle cor = ", round(cor_eagle_all, 3)), color = color_eagle) +
    annotate("text", x = 4, y = 0.006, label = paste0("Shapeit cor = ", round(cor_shapeit_all, 3)), color = color_shapeit) + # make text larger
    theme(text = element_text(size = 16)) + # larger axis labels
    theme(legend.text = element_text(size = 14), legend.title = element_text(size = 16)) +
    guides(color = guide_legend(override.aes = list(size = 5, shape = 15))) + # move legend to interior
    theme(legend.position = "inside",
          legend.position.inside = c(0.8, 0.2))

# same for flips
cor_beagle_flip_all <- cor(df_summary_all$avg_rr_MB_1kgp, df_summary_all$beagle_flip_rate, use = "complete.obs")
cor_eagle_flip_all <- cor(df_summary_all$avg_rr_MB_1kgp, df_summary_all$eagle_flip_rate, use = "complete.obs")
cor_shapeit_flip_all <- cor(df_summary_all$avg_rr_MB_1kgp, df_summary_all$shapeit_flip_rate, use = "complete.obs")

df_summary_all %>%
    select(binMB,avg_rr_MB_1kgp, contains("flip_rate")) %>%
    pivot_longer(cols = contains("flip_rate"), names_to = "method", values_to = "flip_rate") %>%
    mutate(method = case_when(
        grepl("beagle", method) ~ "Beagle",
        grepl("eagle", method) ~ "Eagle",
        grepl("shapeit", method) ~ "Shapeit"
    )) %>%
    ggplot(aes(x = avg_rr_MB_1kgp, y = flip_rate, color = method)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    theme_classic() +
    labs(x = "Average Recombination Rate (cM/Mb)",
         y = "Flip Rate",
         color = "Phasing Method") +
    scale_color_manual(values = friendly_pal(pal_choice)) +
    annotate("text", x = 4, y = 0.005, label = paste0("Beagle cor = ", round(cor_beagle_flip_all, 3)), color = color_beagle) +
    annotate("text", x = 4, y = 0.0045, label = paste0("Eagle cor = ", round(cor_eagle_flip_all, 3)), color = color_eagle) +
    annotate("text", x = 4, y = 0.004, label = paste0("Shapeit cor = ", round(cor_shapeit_flip_all, 3)), color = color_shapeit) + # make text larger
    theme(text = element_text(size = 16)) + # larger axis labels
    theme(legend.text = element_text(size = 14), legend.title = element_text(size = 16)) +
    guides(color = guide_legend(override.aes = list(size = 5, shape = 15))) + # move legend to interior
    theme(legend.position = "inside",
          legend.position.inside = c(0.8, 0.9))

# RR: enrichment in decile 10
compute_enrichment <- function(id, base_dir, df_rr_1kgp_x){
    df_het <- load_het_pos(id, base_dir)
    df_het_rr <- df_het %>%
        left_join(df_rr_1kgp_x %>% select(start, cM_MB), by = join_by(closest(pos >= start))) %>%
        rename(rr_1kgp = cM_MB) %>%
        mutate(rr_1kgp_decile = ntile(rr_1kgp, 10))
    df_enrich <- df_het_rr %>%
        group_by(rr_1kgp_decile) %>%
        summarize(total_hets = n(),
                  beagle_switches = sum(beagle_switch),
                  eagle_switches = sum(eagle_switch),
                  shapeit_switches = sum(shapeit_switch),
                  beagle_flips = sum(beagle_flip),
                  eagle_flips = sum(eagle_flip),
                  shapeit_flips = sum(shapeit_flip)) %>%
        mutate(id = id) %>%
        ungroup()

    beagle_switch <- df_enrich$beagle_switches / sum(df_enrich$beagle_switches)
    eagle_switch <- df_enrich$eagle_switches / sum(df_enrich$eagle_switches)
    shapeit_switch <- df_enrich$shapeit_switches / sum(df_enrich$shapeit_switches)

    beagle_flip <- df_enrich$beagle_flips / sum(df_enrich$beagle_flips)
    eagle_flip <- df_enrich$eagle_flips / sum(df_enrich$eagle_flips)
    shapeit_flip <- df_enrich$shapeit_flips / sum(df_enrich$shapeit_flips)

    data.frame(id = rep(id,10),
               rr_decile = 1:10,
               beagle_switch_decile = beagle_switch,
               eagle_switch_decile = eagle_switch,
               shapeit_switch_decile = shapeit_switch,
               beagle_flip_decile = beagle_flip,
               eagle_flip_decile = eagle_flip,
               shapeit_flip_decile = shapeit_flip)
}

pop_enrichment <- function(ids, base_dir, df_rr_1kgp_x, pop_name){
    res_list <- lapply(ids, compute_enrichment, base_dir = base_dir, df_rr_1kgp_x = df_rr_1kgp_x)
    df_enrich <- bind_rows(res_list) %>%
        mutate(population = pop_name)
    df_enrich
}

df_enrich <- vector("list", 5)
names(df_enrich) <- c("EUR", "AFR", "AMR", "EAS", "SAS")

df_enrich[["EUR"]] <- pop_enrichment(1:200, het_dir, df_rr_1kgp_x, "EUR")
df_enrich[["AFR"]] <- pop_enrichment(201:400, het_dir, df_rr_1kgp_x, "AFR")
df_enrich[["AMR"]] <- pop_enrichment(401:600, het_dir, df_rr_1kgp_x, "AMR")
df_enrich[["EAS"]] <- pop_enrichment(601:800, het_dir, df_rr_1kgp_x, "EAS")
df_enrich[["SAS"]] <- pop_enrichment(801:1000, het_dir, df_rr_1kgp_x, "SAS")

df_enrich_combined <- bind_rows(df_enrich[["EUR"]] %>% mutate(population = "EUR"), 
    df_enrich[["AFR"]] %>% mutate(population = "AFR"),
    df_enrich[["AMR"]] %>% mutate(population = "AMR"),
    df_enrich[["EAS"]] %>% mutate(population = "EAS"),
    df_enrich[["SAS"]] %>% mutate(population = "SAS"))

df_enrich_combined %>%
    group_by(rr_decile) %>%
    summarize(beagle_switch = mean(beagle_switch_decile),
              eagle_switch = mean(eagle_switch_decile),
              shapeit_switch = mean(shapeit_switch_decile),
              beagle_flip = mean(beagle_flip_decile),
              eagle_flip = mean(eagle_flip_decile),
              shapeit_flip = mean(shapeit_flip_decile)) %>%
    pivot_longer(cols = -rr_decile, names_to = "metric", values_to = "proportion") %>%
    mutate(type = ifelse(grepl("switch", metric), "Switch", "Flip"),
           method = case_when(
               grepl("beagle", metric) ~ "Beagle",
               grepl("eagle", metric) ~ "Eagle",
               grepl("shapeit", metric) ~ "Shapeit"
           )) %>%
    ggplot(aes(x = rr_decile, y = proportion, color = method)) +
    geom_point() +
    geom_line() +
    facet_grid(rows = vars(type)) +
    theme_classic() +
    labs(x = "Recombination Rate Decile",
         y = "Proportion of Switch/Flip Errors",
         color = "Phasing Method") +
    scale_color_manual(values = friendly_pal(pal_choice)) +
    theme(text = element_text(size = 16)) + # larger axis labels
    theme(legend.text = element_text(size = 14), legend.title = element_text(size = 16)) +
    guides(color = guide_legend(override.aes = list(size = 5, shape = 15))) + # move legend to interior
    theme(legend.position = "inside",
          legend.position.inside = c(0.8, 0.2))

# Separate plots for switches and flips
df_enrich_combined %>%
    group_by(rr_decile) %>%
    summarize(beagle_switch = mean(beagle_switch_decile),
              eagle_switch = mean(eagle_switch_decile),
              shapeit_switch = mean(shapeit_switch_decile),
              beagle_flip = mean(beagle_flip_decile),
              eagle_flip = mean(eagle_flip_decile),
              shapeit_flip = mean(shapeit_flip_decile)) %>%
    pivot_longer(cols = -rr_decile, names_to = "metric", values_to = "proportion") %>%
    mutate(type = ifelse(grepl("switch", metric), "Switch", "Flip"),
           method = case_when(
               grepl("beagle", metric) ~ "Beagle",
               grepl("eagle", metric) ~ "Eagle",
               grepl("shapeit", metric) ~ "Shapeit"
           )) %>%
    filter(type == "Switch") %>%
    ggplot(aes(x = rr_decile, y = proportion, color = method)) +
    geom_point() +
    geom_line() +
    theme_classic() +
    labs(x = "Recombination Rate Decile",
         y = "Proportion of Switch/Flip Errors",
         color = "Phasing Method") +
    scale_color_manual(values = friendly_pal(pal_choice)) +
    theme(text = element_text(size = 16)) + # larger axis labels
    theme(legend.text = element_text(size = 14), legend.title = element_text(size = 16)) +
    guides(color = guide_legend(override.aes = list(size = 5, shape = 15))) + # move legend to interior
    theme(legend.position = "inside",
          legend.position.inside = c(0.8, 0.2))

df_enrich_combined %>%
    group_by(rr_decile) %>%
    summarize(beagle_switch = mean(beagle_switch_decile),
              eagle_switch = mean(eagle_switch_decile),
              shapeit_switch = mean(shapeit_switch_decile),
              beagle_flip = mean(beagle_flip_decile),
              eagle_flip = mean(eagle_flip_decile),
              shapeit_flip = mean(shapeit_flip_decile)) %>%
    pivot_longer(cols = -rr_decile, names_to = "metric", values_to = "proportion") %>%
    mutate(type = ifelse(grepl("switch", metric), "Switch", "Flip"),
           method = case_when(
               grepl("beagle", metric) ~ "Beagle",
               grepl("eagle", metric) ~ "Eagle",
               grepl("shapeit", metric) ~ "Shapeit"
           )) %>%
    filter(type == "Flip") %>%
    ggplot(aes(x = rr_decile, y = proportion, color = method)) +
    geom_point() +
    geom_line() +
    theme_classic() +
    labs(x = "Recombination Rate Decile",
         y = "Proportion of Switch/Flip Errors",
         color = "Phasing Method") +
    scale_color_manual(values = friendly_pal(pal_choice)) +
    theme(text = element_text(size = 16)) + # larger axis labels
    theme(legend.text = element_text(size = 14), legend.title = element_text(size = 16)) +
    guides(color = guide_legend(override.aes = list(size = 5, shape = 15))) + # move legend to interior
    theme(legend.position = "inside",
          legend.position.inside = c(0.8, 0.2))

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
           n_flip_cpg_shapeit = n_flip_cpg_s)

# Get average switch and flip rates per method
df_proband_summary <- df_proband %>%
    mutate(beagle_switch_rate = n_switch_beagle / n_het,
           eagle_switch_rate = n_switch_eagle / n_het,
           shapeit_switch_rate = n_switch_shapeit / n_het,
           beagle_flip_rate = n_flip_beagle / n_het,
           eagle_flip_rate = n_flip_eagle / n_het,
           shapeit_flip_rate = n_flip_shapeit / n_het,
           beagle_error_rate = (n_switch_beagle + n_flip_beagle) / n_het,
           eagle_error_rate = (n_switch_eagle + n_flip_eagle) / n_het,
           shapeit_error_rate = (n_switch_shapeit + n_flip_shapeit) / n_het) %>%
    group_by(chrom) %>%
    summarize(mean_beagle_switch = mean(beagle_switch_rate),
              mean_eagle_switch = mean(eagle_switch_rate),
              mean_shapeit_switch = mean(shapeit_switch_rate),
              mean_beagle_flip = mean(beagle_flip_rate),
              mean_eagle_flip = mean(eagle_flip_rate),
              mean_shapeit_flip = mean(shapeit_flip_rate),
              mean_beagle_error = mean(beagle_error_rate),
              mean_eagle_error = mean(eagle_error_rate),
              mean_shapeit_error = mean(shapeit_error_rate)) 

# Get average switch and flip rates from chromosome X
df_proband_summary <- df_sd_errors %>%
    mutate(beagle_switch_rate = n_other_beagle / n_hets,
           eagle_switch_rate = n_other_eagle / n_hets,
           shapeit_switch_rate = n_other_shapeit / n_hets,
           beagle_flip_rate = n_flip_beagle / n_hets,
           eagle_flip_rate = n_flip_eagle / n_hets,
           shapeit_flip_rate = n_flip_shapeit / n_hets,
           beagle_error_rate = (n_other_beagle + n_flip_beagle) / n_hets,
           eagle_error_rate = (n_other_eagle + n_flip_eagle) / n_hets,
           shapeit_error_rate = (n_other_shapeit + n_flip_shapeit) / n_hets) %>%
    summarize(mean_beagle_switch = mean(beagle_switch_rate),
              mean_eagle_switch = mean(eagle_switch_rate),
              mean_shapeit_switch = mean(shapeit_switch_rate),
              mean_beagle_flip = mean(beagle_flip_rate),
              mean_eagle_flip = mean(eagle_flip_rate),
              mean_shapeit_flip = mean(shapeit_flip_rate),
              mean_beagle_error = mean(beagle_error_rate),
              mean_eagle_error = mean(eagle_error_rate),
              mean_shapeit_error = mean(shapeit_error_rate)) %>%
    mutate(chrom = 23) %>%
    select(chrom, everything()) %>%
    bind_rows(df_proband_summary) %>%
    arrange(chrom)


df_proband_summary_pop <- df_proband %>%
    mutate(beagle_switch_rate = n_switch_beagle / n_het,
           eagle_switch_rate = n_switch_eagle / n_het,
           shapeit_switch_rate = n_switch_shapeit / n_het,
           beagle_flip_rate = n_flip_beagle / n_het,
           eagle_flip_rate = n_flip_eagle / n_het,
           shapeit_flip_rate = n_flip_shapeit / n_het,
           beagle_error_rate = (n_switch_beagle + n_flip_beagle) / n_het,
           eagle_error_rate = (n_switch_eagle + n_flip_eagle) / n_het,
           shapeit_error_rate = (n_switch_shapeit + n_flip_shapeit) / n_het) %>%
    group_by(chrom, SUPER) %>%
    summarize(mean_beagle_switch = mean(beagle_switch_rate),
              mean_eagle_switch = mean(eagle_switch_rate),
              mean_shapeit_switch = mean(shapeit_switch_rate),
              mean_beagle_flip = mean(beagle_flip_rate),
              mean_eagle_flip = mean(eagle_flip_rate),
              mean_shapeit_flip = mean(shapeit_flip_rate),
              mean_beagle_error = mean(beagle_error_rate),
              mean_eagle_error = mean(eagle_error_rate),
              mean_shapeit_error = mean(shapeit_error_rate))

# Get average switch and flip rates from chromosome X
df_proband_summary_pop <- df_sd_errors %>%
    mutate(beagle_switch_rate = n_other_beagle / n_hets,
           eagle_switch_rate = n_other_eagle / n_hets,
           shapeit_switch_rate = n_other_shapeit / n_hets,
           beagle_flip_rate = n_flip_beagle / n_hets,
           eagle_flip_rate = n_flip_eagle / n_hets,
           shapeit_flip_rate = n_flip_shapeit / n_hets,
           beagle_error_rate = (n_other_beagle + n_flip_beagle) / n_hets,
           eagle_error_rate = (n_other_eagle + n_flip_eagle) / n_hets,
           shapeit_error_rate = (n_other_shapeit + n_flip_shapeit) / n_hets) %>%
    group_by(pop) %>%
    summarize(mean_beagle_switch = mean(beagle_switch_rate),
              mean_eagle_switch = mean(eagle_switch_rate),
              mean_shapeit_switch = mean(shapeit_switch_rate),
              mean_beagle_flip = mean(beagle_flip_rate),
              mean_eagle_flip = mean(eagle_flip_rate),
              mean_shapeit_flip = mean(shapeit_flip_rate),
              mean_beagle_error = mean(beagle_error_rate),
              mean_eagle_error = mean(eagle_error_rate),
              mean_shapeit_error = mean(shapeit_error_rate)) %>%
    mutate(chrom = 23) %>%
    select(chrom, everything()) %>%
    rename(SUPER = pop) %>%
    bind_rows(df_proband_summary_pop) %>%
    arrange(chrom, SUPER)

# Total Error Rates by Chromosome
df_proband_summary %>%
    select(chrom, ends_with("_error")) %>%
    pivot_longer(cols = contains("mean"), names_to = "Metric", values_to = "Rate") %>%
    separate(Metric, into = c("mean","Method", "Error_Type"), sep = "_") %>%
    select(-mean) %>%
    mutate(Method = factor(Method, levels = c("beagle", "eagle", "shapeit"))) %>%
    mutate(Method = recode(Method, "beagle" = "Beagle", "eagle" = "Eagle", "shapeit" = "SHAPEIT")) %>%
    mutate(chrom = factor(chrom, levels = c(1:22, 23), labels = c(1:22, "X"))) %>%
    ggplot(aes(x = chrom, y = Rate, color = Method)) +
    geom_point(size = 3) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("Chromosome") +
    ylab("Error Rate") +
    scale_color_manual(values = friendly_pal(pal_choice)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
    theme(strip.text = element_text(size = 14)) +
    theme(legend.position = "top", legend.box = "horizontal") + # remove legend title
    theme(legend.text = element_text(size = 14)) +
    theme(legend.title = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 5, 
        shape = 22, 
        linetype = "blank",
        fill = friendly_pal(pal_choice)[1:3]),
        title = ""))

# Plot switch and flip rates by chromosome
df_proband_summary %>%
    select(chrom, ends_with("_flip"), ends_with("_switch")) %>%
    pivot_longer(cols = contains("mean"), names_to = "Metric", values_to = "Rate") %>%
    separate(Metric, into = c("mean","Method", "Error_Type"), sep = "_") %>%
    select(-mean) %>%
    mutate(Error_Type = factor(Error_Type, levels = c("switch", "flip"))) %>%
    # recode switch to Switch  and flip to Flip 
    mutate(Error_Type = recode(Error_Type, "switch" = "Switch", "flip" = "Flip")) %>%
    mutate(Method = factor(Method, levels = c("beagle", "eagle", "shapeit"))) %>%
    mutate(Method = recode(Method, "beagle" = "Beagle", "eagle" = "Eagle", "shapeit" = "SHAPEIT")) %>%
    mutate(chrom = factor(chrom, levels = c(1:22, 23), labels = c(1:22, "X"))) %>%
    ggplot(aes(x = chrom, y = Rate, color = Method)) +
    geom_point(size = 3) +
    facet_grid(rows = vars(Error_Type)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("Chromosome") +
    ylab("Error Rate") +
    scale_color_manual(values = friendly_pal(pal_choice)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
    theme(strip.text = element_text(size = 14)) +
    theme(legend.position = "top", legend.box = "horizontal") + # remove legend title
    theme(legend.text = element_text(size = 14)) +
    theme(legend.title = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 5, 
        shape = 22, 
        linetype = "blank",
        fill = friendly_pal(pal_choice)[1:3]),
        title = ""))

# switch rates by population
df_proband_summary_pop %>%
    select(chrom, SUPER, ends_with("_flip"), ends_with("_switch")) %>%
    pivot_longer(cols = contains("mean"), names_to = "Metric", values_to = "Rate") %>%
    separate(Metric, into = c("mean","Method", "Error_Type"), sep = "_") %>%
    select(-mean) %>%
    mutate(Error_Type = factor(Error_Type, levels = c("switch", "flip"))) %>%
    # recode switch to Switch  and flip to Flip 
    mutate(Error_Type = recode(Error_Type, "switch" = "Switch", "flip" = "Flip")) %>%
    mutate(Method = factor(Method, levels = c("beagle", "eagle", "shapeit"))) %>%
    mutate(Method = recode(Method, "beagle" = "Beagle", "eagle" = "Eagle", "shapeit" = "SHAPEIT")) %>%
    mutate(chrom = factor(chrom, levels = c(1:22, 23), labels = c(1:22, "X"))) %>%
    filter(Error_Type == "Switch") %>%
    ggplot(aes(x = chrom, y = Rate, color = Method)) +
    geom_point(size = 3) +
    facet_grid(rows = vars(SUPER)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("Chromosome") +
    ylab("Switch Rate") +
    scale_color_manual(values = friendly_pal(pal_choice)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
    theme(strip.text = element_text(size = 14)) +
    theme(legend.position = "top", legend.box = "horizontal") + # remove legend title
    theme(legend.text = element_text(size = 14)) +
    theme(legend.title = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 5, 
        shape = 22, 
        linetype = "blank",
        fill = friendly_pal(pal_choice)[1:3]),
        title = ""))

# Flip rates by population
df_proband_summary_pop %>%
    select(chrom, SUPER, ends_with("_flip"), ends_with("_switch")) %>%
    pivot_longer(cols = contains("mean"), names_to = "Metric", values_to = "Rate") %>%
    separate(Metric, into = c("mean","Method", "Error_Type"), sep = "_") %>%
    select(-mean) %>%
    mutate(Error_Type = factor(Error_Type, levels = c("switch", "flip"))) %>%
    # recode switch to Switch  and flip to Flip 
    mutate(Error_Type = recode(Error_Type, "switch" = "Switch", "flip" = "Flip")) %>%
    mutate(Method = factor(Method, levels = c("beagle", "eagle", "shapeit"))) %>%
    mutate(Method = recode(Method, "beagle" = "Beagle", "eagle" = "Eagle", "shapeit" = "SHAPEIT")) %>%
    mutate(chrom = factor(chrom, levels = c(1:22, 23), labels = c(1:22, "X"))) %>%
    filter(Error_Type == "Flip") %>%
    ggplot(aes(x = chrom, y = Rate, color = Method)) +
    geom_point(size = 3) +
    facet_grid(rows = vars(SUPER)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("Chromosome") +
    ylab("Flip Rate") +
    scale_color_manual(values = friendly_pal(pal_choice)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
    theme(strip.text = element_text(size = 14)) +
    theme(legend.position = "top", legend.box = "horizontal") + # remove legend title
    theme(legend.text = element_text(size = 14)) +
    theme(legend.title = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 5, 
        shape = 22, 
        linetype = "blank",
        fill = friendly_pal(pal_choice)[1:3]),
        title = ""))

# Total Error Rates by Population and Method
df_proband_summary_pop %>%
    select(chrom, SUPER, ends_with("_error")) %>%
    pivot_longer(cols = contains("mean"), names_to = "Metric", values_to = "Rate") %>%
    separate(Metric, into = c("mean","Method", "Error_Type"), sep = "_") %>%
    select(-mean) %>%
    mutate(Method = factor(Method, levels = c("beagle", "eagle", "shapeit"))) %>%
    mutate(Method = recode(Method, "beagle" = "Beagle", "eagle" = "Eagle", "shapeit" = "SHAPEIT")) %>%
    mutate(chrom = factor(chrom, levels = c(1:22, 23), labels = c(1:22, "X"))) %>%
    ggplot(aes(x = chrom, y = Rate, shape = Method, color = SUPER)) +
    geom_point(size = 3) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("Chromosome") +
    ylab("Error Rate") +
    scale_color_manual(values = friendly_pal(pal_choice)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
    theme(strip.text = element_text(size = 14)) +
    theme(legend.position = "top", legend.box = "horizontal") + # remove legend title
    theme(legend.text = element_text(size = 14)) +
    theme(legend.title = element_blank())

df_proband_summary_pop %>%
    select(chrom, SUPER, ends_with("_switch")) %>%
    pivot_longer(cols = contains("mean"), names_to = "Metric", values_to = "Rate") %>%
    separate(Metric, into = c("mean","Method", "Error_Type"), sep = "_") %>%
    select(-mean) %>%
    mutate(Method = factor(Method, levels = c("beagle", "eagle", "shapeit"))) %>%
    mutate(Method = recode(Method, "beagle" = "Beagle", "eagle" = "Eagle", "shapeit" = "SHAPEIT")) %>%
    mutate(chrom = factor(chrom, levels = c(1:22, 23), labels = c(1:22, "X"))) %>%
    ggplot(aes(x = chrom, y = Rate, shape = Method, color = SUPER)) +
    geom_point(size = 3) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("Chromosome") +
    ylab("Switch Rate") +
    scale_color_manual(values = friendly_pal(pal_choice)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
    theme(strip.text = element_text(size = 14)) +
    theme(legend.position = "top", legend.box = "horizontal") + # remove legend title
    theme(legend.text = element_text(size = 14)) +
    theme(legend.title = element_blank())

df_proband_summary_pop %>%
    select(chrom, SUPER, ends_with("_flip")) %>%
    pivot_longer(cols = contains("mean"), names_to = "Metric", values_to = "Rate") %>%
    separate(Metric, into = c("mean","Method", "Error_Type"), sep = "_") %>%
    select(-mean) %>%
    mutate(Method = factor(Method, levels = c("beagle", "eagle", "shapeit"))) %>%
    mutate(Method = recode(Method, "beagle" = "Beagle", "eagle" = "Eagle", "shapeit" = "SHAPEIT")) %>%
    mutate(chrom = factor(chrom, levels = c(1:22, 23), labels = c(1:22, "X"))) %>%
    ggplot(aes(x = chrom, y = Rate, shape = Method, color = SUPER)) +
    geom_point(size = 3) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("Chromosome") +
    ylab("Flip Rate") +
    scale_color_manual(values = friendly_pal(pal_choice)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
    theme(strip.text = element_text(size = 14)) +
    theme(legend.position = "top", legend.box = "horizontal") + # remove legend title
    theme(legend.text = element_text(size = 14)) +
    theme(legend.title = element_blank())

# Get rows from df_proband and df_sd_errors 
df_all <- df_proband %>%
    select(chrom, sample, n_het, n_switch_beagle, n_switch_eagle, n_switch_shapeit, n_flip_beagle, n_flip_eagle, n_flip_shapeit,
           n_het_cpg, n_switch_cpg_beagle, n_switch_cpg_eagle, n_switch_cpg_shapeit, n_flip_cpg_beagle, n_flip_cpg_eagle, n_flip_cpg_shapeit, SUPER) %>%
    rename(pair_id = sample,
           pop = SUPER) %>%
    bind_rows(df_sd_errors %>%
        mutate(chrom = 23) %>%
        select(chrom, pair_id, n_hets, n_other_beagle, n_other_eagle, n_other_shapeit, n_flip_beagle, n_flip_eagle, n_flip_shapeit,
            n_het_cpg, n_other_cpg_beagle, n_other_cpg_eagle, n_other_cpg_shapeit, n_flip_cpg_beagle, n_flip_cpg_eagle, n_flip_cpg_shapeit, pop) %>%
        rename(n_het = n_hets,
               n_switch_beagle = n_other_beagle,
               n_switch_eagle = n_other_eagle,
               n_switch_shapeit = n_other_shapeit,
               n_switch_cpg_beagle = n_other_cpg_beagle,
               n_switch_cpg_eagle = n_other_cpg_eagle,
               n_switch_cpg_shapeit = n_other_cpg_shapeit)) %>%
    arrange(chrom, pair_id)

# Test the null hypothesis that mean switch rates are the same across chromosomes
anova_beagle <- aov(beagle_switch_rate ~ as.factor(chrom), data = {df_all %>% mutate(beagle_switch_rate = n_switch_beagle / n_het)})
summary(anova_beagle)
anova_eagle <- aov(eagle_switch_rate ~ as.factor(chrom), data = {df_all %>% mutate(eagle_switch_rate = n_switch_eagle / n_het)})
summary(anova_eagle)
anova_shapeit <- aov(shapeit_switch_rate ~ as.factor(chrom), data = {df_all %>% mutate(shapeit_switch_rate = n_switch_shapeit / n_het)})
summary(anova_shapeit)

# Test the null hypothesis that mean flip rates are the same across chromosomes
anova_beagle_flip <- aov(beagle_flip_rate ~ as.factor(chrom), data = {df_all %>% mutate(beagle_flip_rate = n_flip_beagle / n_het)})
summary(anova_beagle_flip)
anova_eagle_flip <- aov(eagle_flip_rate ~ as.factor(chrom), data = {df_all %>% mutate(eagle_flip_rate = n_flip_eagle / n_het)})
summary(anova_eagle_flip)
anova_shapeit_flip <- aov(shapeit_flip_rate ~ as.factor(chrom), data = {df_all %>% mutate(shapeit_flip_rate = n_flip_shapeit / n_het)})
summary(anova_shapeit_flip) 

# Test whether X is different from autosomes
t_test_beagle <- t.test(beagle_switch_rate ~ is_X, data = {df_all %>% mutate(beagle_switch_rate = n_switch_beagle / n_het, is_X = ifelse(chrom == 23, "X", "Autosome") )})
t_test_beagle
t_test_eagle <- t.test(eagle_switch_rate ~ is_X, data = {df_all %>% mutate(eagle_switch_rate = n_switch_eagle / n_het, is_X = ifelse(chrom == 23, "X", "Autosome") )})
t_test_eagle
t_test_shapeit <- t.test(shapeit_switch_rate ~ is_X, data = {df_all %>% mutate(shapeit_switch_rate = n_switch_shapeit / n_het, is_X = ifelse(chrom == 23, "X", "Autosome") )})
t_test_shapeit

# Flips
t_test_beagle_flip <- t.test(beagle_flip_rate ~ is_X, data = {df_all %>% mutate(beagle_flip_rate = n_flip_beagle / n_het, is_X = ifelse(chrom == 23, "X", "Autosome") )})
t_test_beagle_flip
t_test_eagle_flip <- t.test(eagle_flip_rate ~ is_X, data = {df_all %>% mutate(eagle_flip_rate = n_flip_eagle / n_het, is_X = ifelse(chrom == 23, "X", "Autosome") )})
t_test_eagle_flip
t_test_shapeit_flip <- t.test(shapeit_flip_rate ~ is_X, data = {df_all %>% mutate(shapeit_flip_rate = n_flip_shapeit / n_het, is_X = ifelse(chrom == 23, "X", "Autosome") )})
t_test_shapeit_flip

# Table of switch rates by chromosome
df_all %>%
    mutate(beagle_switch_rate = n_switch_beagle / n_het,
           eagle_switch_rate = n_switch_eagle / n_het,
           shapeit_switch_rate = n_switch_shapeit / n_het,
           beagle_flip_rate = n_flip_beagle / n_het,
           eagle_flip_rate = n_flip_eagle / n_het,
           shapeit_flip_rate = n_flip_shapeit / n_het) %>%
    group_by(chrom) %>%
    summarize(mean_beagle_switch = mean(beagle_switch_rate),
              mean_eagle_switch = mean(eagle_switch_rate),
              mean_shapeit_switch = mean(shapeit_switch_rate),
              mean_beagle_flip = mean(beagle_flip_rate),
              mean_eagle_flip = mean(eagle_flip_rate),
              mean_shapeit_flip = mean(shapeit_flip_rate)) %>%
    mutate(chrom = ifelse(chrom == 23, "X", as.character(chrom))) %>%
    arrange(as.numeric(replace(chrom, chrom == "X", 23))) %>%
    knitr::kable(digits = 4, col.names = c("Chromosome", "Mean Beagle Switch Rate", "Mean Eagle Switch Rate", "Mean SHAPEIT Switch Rate", "Mean Beagle Flip Rate", "Mean Eagle Flip Rate", "Mean SHAPEIT Flip Rate"))

# What is the range of switch rates we observe in the above table?
df_all %>%
    mutate(beagle_switch_rate = n_switch_beagle / n_het,
           eagle_switch_rate = n_switch_eagle / n_het,
           shapeit_switch_rate = n_switch_shapeit / n_het) %>%
    group_by(chrom) %>%
    summarize(mean_beagle_switch = mean(beagle_switch_rate),
              mean_eagle_switch = mean(eagle_switch_rate),
              mean_shapeit_switch = mean(shapeit_switch_rate)) %>%
    summarize(min_beagle = min(mean_beagle_switch),
              max_beagle = max(mean_beagle_switch),
              min_eagle = min(mean_eagle_switch),
              max_eagle = max(mean_eagle_switch),
              min_shapeit = min(mean_shapeit_switch),
              max_shapeit = max(mean_shapeit_switch))

# Plot total error rates by chromosome and method
df_all %>%
    mutate(Beagle = (n_switch_beagle + n_flip_beagle) / n_het,
           Eagle = (n_switch_eagle + n_flip_eagle) / n_het,
           SHAPEIT = (n_switch_shapeit + n_flip_shapeit) / n_het) %>%
    select(chrom, Beagle, Eagle, SHAPEIT) %>%
    pivot_longer(cols = c("Beagle", "Eagle", "SHAPEIT"), names_to = "Method", values_to = "Error_Rate") %>%
    group_by(chrom, Method) %>%
    summarize(Error_Rate = mean(Error_Rate)) %>%
    mutate(chrom = factor(chrom, levels = c(1:22, 23), labels = c(1:22, "X"))) %>%
    ggplot(aes(x = chrom, y = Error_Rate, color = Method)) +
    geom_point(size = 3) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("Chromosome") +
    ylab("Total Error Rate") +
    scale_color_manual(values = friendly_pal(pal_choice)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
    theme(strip.text = element_text(size = 14)) +
    theme(legend.position = "top", legend.box = "horizontal") + # remove legend title
    theme(legend.text = element_text(size = 14)) +
    theme(legend.title = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 5, 
        shape = 22, 
        linetype = "blank",
        fill = friendly_pal(pal_choice)[1:3]),
        title = ""))

# Get table of average error rates by chromsome and method
df_all %>%
    mutate(Beagle = (n_switch_beagle + n_flip_beagle) / n_het,
           Eagle = (n_switch_eagle + n_flip_eagle) / n_het,
           SHAPEIT = (n_switch_shapeit + n_flip_shapeit) / n_het) %>%
    select(chrom, Beagle, Eagle, SHAPEIT) %>%
    group_by(chrom) %>%
    summarize(Mean_Beagle = mean(Beagle),
              Mean_Eagle = mean(Eagle),
              Mean_SHAPEIT = mean(SHAPEIT)) %>%
    mutate(chrom = ifelse(chrom == 23, "X", as.character(chrom))) %>%
    arrange(as.numeric(replace(chrom, chrom == "X", 23))) %>%
    knitr::kable(digits = 4, col.names = c("Chromosome", "Mean Beagle Error Rate", "Mean Eagle Error Rate", "Mean SHAPEIT Error Rate"))

# Get range of total error rates for each method, excluding X
df_all %>%
    filter(chrom != 23) %>%
    mutate(Beagle = (n_switch_beagle + n_flip_beagle) / n_het,
           Eagle = (n_switch_eagle + n_flip_eagle) / n_het,
           SHAPEIT = (n_switch_shapeit + n_flip_shapeit) / n_het) %>%
    select(chrom, Beagle, Eagle, SHAPEIT) %>%
    group_by(chrom) %>%
    summarize(Mean_Beagle = mean(Beagle),
              Mean_Eagle = mean(Eagle),
              Mean_SHAPEIT = mean(SHAPEIT)) %>%
    summarize(Min_Beagle = min(Mean_Beagle),
              Max_Beagle = max(Mean_Beagle),
              Min_Eagle = min(Mean_Eagle),
              Max_Eagle = max(Mean_Eagle),
              Min_SHAPEIT = min(Mean_SHAPEIT),
              Max_SHAPEIT = max(Mean_SHAPEIT))

# Get range of flip rates for each method, excluding X
df_all %>%
    filter(chrom != 23) %>%
    mutate(beagle_flip_rate = n_flip_beagle / n_het,
           eagle_flip_rate = n_flip_eagle / n_het,
           shapeit_flip_rate = n_flip_shapeit / n_het) %>%
    select(chrom, beagle_flip_rate, eagle_flip_rate, shapeit_flip_rate) %>%
    group_by(chrom) %>%
    summarize(mean_beagle_flip = mean(beagle_flip_rate),
              mean_eagle_flip = mean(eagle_flip_rate),
              mean_shapeit_flip = mean(shapeit_flip_rate)) %>%
    summarize(min_beagle = min(mean_beagle_flip),
              max_beagle = max(mean_beagle_flip),
              min_eagle = min(mean_eagle_flip),
              max_eagle = max(mean_eagle_flip),
              min_shapeit = min(mean_shapeit_flip),
              max_shapeit = max(mean_shapeit_flip)) 

# Get the range of flip rates on X
df_all %>%
    filter(chrom == 23) %>%
    mutate(beagle_flip_rate = n_flip_beagle / n_het,
           eagle_flip_rate = n_flip_eagle / n_het,
           shapeit_flip_rate = n_flip_shapeit / n_het) %>%
    select(chrom, beagle_flip_rate, eagle_flip_rate, shapeit_flip_rate) %>%
    summarize(mean_beagle_flip = mean(beagle_flip_rate),
              mean_eagle_flip = mean(eagle_flip_rate),
              mean_shapeit_flip = mean(shapeit_flip_rate))

# Switch rate ranges excluding X
df_all %>%
    filter(chrom != 23) %>%
    mutate(beagle_switch_rate = n_switch_beagle / n_het,
           eagle_switch_rate = n_switch_eagle / n_het,
           shapeit_switch_rate = n_switch_shapeit / n_het) %>%
    select(chrom, beagle_switch_rate, eagle_switch_rate, shapeit_switch_rate) %>%
    group_by(chrom) %>%
    summarize(mean_beagle_switch = mean(beagle_switch_rate),
              mean_eagle_switch = mean(eagle_switch_rate),
              mean_shapeit_switch = mean(shapeit_switch_rate)) %>%
    summarize(min_beagle = min(mean_beagle_switch),
              max_beagle = max(mean_beagle_switch),
              min_eagle = min(mean_eagle_switch),
              max_eagle = max(mean_eagle_switch),
              min_shapeit = min(mean_shapeit_switch),
              max_shapeit = max(mean_shapeit_switch))

# Switch rates on X
df_all %>%
    filter(chrom == 23) %>%
    mutate(beagle_switch_rate = n_switch_beagle / n_het,
           eagle_switch_rate = n_switch_eagle / n_het,
           shapeit_switch_rate = n_switch_shapeit / n_het) %>%
    select(chrom, beagle_switch_rate, eagle_switch_rate, shapeit_switch_rate) %>%
    summarize(mean_beagle_switch = mean(beagle_switch_rate),
              mean_eagle_switch = mean(eagle_switch_rate),
              mean_shapeit_switch = mean(shapeit_switch_rate))

# Compute enrichment of flip and switch errors at CpG
df_all %>%
    mutate(beagle_switch = (n_switch_cpg_beagle / n_switch_beagle)/(n_het_cpg / n_het),
           eagle_switch = (n_switch_cpg_eagle / n_switch_eagle)/(n_het_cpg / n_het),
           shapeit_switch = (n_switch_cpg_shapeit / n_switch_shapeit)/(n_het_cpg / n_het),
           beagle_flip = (n_flip_cpg_beagle / n_flip_beagle)/(n_het_cpg / n_het),
           eagle_flip = (n_flip_cpg_eagle / n_flip_eagle)/(n_het_cpg / n_het),
           shapeit_flip = (n_flip_cpg_shapeit / n_flip_shapeit)/(n_het_cpg / n_het)) %>%
    select(chrom, beagle_switch:shapeit_flip) %>%
    pivot_longer(cols = -chrom, names_to = "Metric", values_to = "Enrichment") %>%
    separate(Metric, into = c("Method", "Error_Type"), sep = "_") %>%
    mutate(Error_Type = factor(Error_Type, levels = c("switch", "flip"))) %>%
    # recode switch to Switch  and flip to Flip 
    mutate(Error_Type = recode(Error_Type, "switch" = "Switch", "flip" = "Flip")) %>%
    mutate(Method = factor(Method, levels = c("beagle", "eagle", "shapeit"))) %>%
    mutate(Method = recode(Method, "beagle" = "Beagle", "eagle" = "Eagle", "shapeit" = "SHAPEIT")) %>%
    mutate(chrom = factor(chrom, levels = c(1:22, 23), labels = c(1:22, "X"))) %>%
    ggplot(aes(x = chrom, y = Enrichment, color = Method)) +
    geom_boxplot() +
    facet_grid(rows = vars(Error_Type)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("Chromosome") +
    ylab("Enrichment at CpG Sites") +
    scale_color_manual(values = friendly_pal(pal_choice)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
    theme(strip.text = element_text(size = 14)) +
    theme(legend.position = "top", legend.box = "horizontal") + # remove legend title
    theme(legend.text = element_text(size = 14)) +
    theme(legend.title = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 5, 
        shape = 22, 
        linetype = "blank",
        fill = friendly_pal(pal_choice)[1:3]),
        title = ""))

# Get fraction of hets at CpG in each chromosome
df_all %>%
    mutate(prop_het_cpg = n_het_cpg / n_het) %>%
    select(chrom, prop_het_cpg) %>%
    mutate(chrom = factor(chrom, levels = c(1:22, 23), labels = c(1:22, "X"))) %>%
    ggplot(aes(x = chrom, y = prop_het_cpg)) +
    geom_violin(fill = friendly_pal(pal_choice)[1], alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.5) +
    labs(title = "",
         x = "Chromosome",
         y = "Proportion of Hets at CpG Sites") +
    theme_minimal() +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16))

# Get average proportion of hets at CpG by chromosome, along with 25th and 75th percentiles
df_all %>%
    mutate(prop_het_cpg = n_het_cpg / n_het) %>%
    select(chrom, prop_het_cpg) %>%
    group_by(chrom) %>%
    summarize(mean_prop_het_cpg = mean(prop_het_cpg),
              p25_prop_het_cpg = quantile(prop_het_cpg, 0.25),
              p75_prop_het_cpg = quantile(prop_het_cpg, 0.75)) %>%
    arrange(chrom) %>%
    knitr::kable(digits = 4, col.names = c("Chromosome", "Mean Proportion of Hets at CpG Sites", "25th Percentile", "75th Percentile"))

# Get proportion of hets in each maf category
# data are in output/autosome_summaries/het_maf_props_chr{CHR}.csv
output_dir <- "output"
df_maf_props <- vector("list", 23)
for (i in 1:23) {
    if(i == 23){
        chr_str <- "X"
    } else {
        chr_str <- as.character(i)
    }

    df_maf_props[[i]] <- read.csv(file.path(output_dir, "autosome_summaries", paste0("het_maf_props_chr", chr_str, ".csv")))
    
    if(i == 23){
        df_maf_props[[i]] <- df_maf_props[[i]] %>%
            rename(SUPER = pop)
    }

    df_maf_props[[i]]$chrom <- i
}

df_maf_props_all <- bind_rows(df_maf_props)

# Pivot from wide to long (prop_rare, prop_low, prop_common), and plot boxplots by chromosome
df_maf_props_long <- df_maf_props_all %>%
    pivot_longer(cols = starts_with("prop"), names_to = "maf_category", values_to = "proportion")

ggplot(df_maf_props_long, aes(x = as.factor(chrom), y = proportion, fill = maf_category)) +
    geom_boxplot() +
    labs(title = "",
         x = "Chromosome",
         y = "Proportion") +
    theme_minimal()

# Now just plot violin plots for prop_rare
df_maf_props_all %>%
    ggplot(aes(x = as.factor(chrom), y = prop_rare)) +
    geom_violin(fill = friendly_pal(pal_choice)[1], alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.5) +
    labs(title = "",
         x = "Chromosome",
         y = "Proportion of Rare Hets (MAF < 0.01)") +
    theme_minimal() +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16))

# Compute average proportion of heterozygous sites in rare category by chromosome
df_maf_props_all %>%
    group_by(chrom) %>%
    summarize(mean_prop_rare = mean(prop_rare)) %>%
    arrange(chrom) %>%
    knitr::kable(digits = 4, col.names = c("Chromosome", "Mean Proportion of Rare Hets (MAF < 0.01)"))

# Get range of values excluding 23 (X)
df_maf_props_all %>%
    filter(chrom != 23) %>%
    group_by(chrom) %>%
    summarize(mean_prop_rare = mean(prop_rare)) %>%
    summarize(min_prop_rare = min(mean_prop_rare),
              max_prop_rare = max(mean_prop_rare))

# Now lets look at the distributions of MAFs at errors
# These are in files output/autosome_summaries/error_maf_props_chr{CHR}.csv
df_error_maf_props <- vector("list", 23)
for (i in 1:23) {
    if(i == 23){
        chr_str <- "X"
    } else {
        chr_str <- as.character(i)
    }

    df_error_maf_props[[i]] <- read.csv(file.path(output_dir, "autosome_summaries", paste0("error_maf_props_chr", chr_str, ".csv")))

    if(i == 23){
        df_error_maf_props[[i]] <- df_error_maf_props[[i]] %>%
            rename(SUPER = pop)
    }

    df_error_maf_props[[i]]$chrom <- i
}

df_error_maf_props_all <- bind_rows(df_error_maf_props)

# plots violin plot of prop_rare by chrom for error_type == "switch", facet by method
df_error_maf_props_all %>%
    filter(error_type == "switch") %>%
    ggplot(aes(x = as.factor(chrom), y = prop_rare, fill = method)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.5) +
    labs(title = "",
         x = "Chromosome",
         y = "Proportion of Rare Hets (MAF < 0.01)") +
    theme_minimal() +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
    theme(legend.position = "top", legend.box = "horizontal") + # remove legend title
    theme(legend.text = element_text(size = 14)) +
    theme(legend.title = element_blank()) +
    guides(fill = guide_legend(override.aes = list(size = 5, 
        shape = 22, 
        linetype = "blank",
        fill = friendly_pal(pal_choice)[1:3]),
        title = "")) +
    facet_wrap(~method)

# Same for flips
df_error_maf_props_all %>%
    filter(error_type == "flip") %>%
    ggplot(aes(x = as.factor(chrom), y = prop_rare, fill = method)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.5) +
    labs(title = "",
         x = "Chromosome",
         y = "Proportion of Rare Hets (MAF < 0.01)") +
    theme_minimal() +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
    theme(legend.position = "top", legend.box = "horizontal") + # remove legend title
    theme(legend.text = element_text(size = 14)) +
    theme(legend.title = element_blank()) +
    guides(fill = guide_legend(override.aes = list(size = 5, 
        shape = 22, 
        linetype = "blank",
        fill = friendly_pal(pal_choice)[1:3]),
        title = "")) +
    facet_wrap(~method)



# Get table of mean number of hets per chromosome
df_proband %>%
    group_by(chrom) %>%
    summarize(mean_hets = mean(n_het)) %>%
    arrange(chrom) %>%
    knitr::kable()

# Get mean number of hets on chrX
df_sd_errors %>%
    summarize(mean_hets = mean(n_hets))


# Store mean number of hets in a data frame
df_mean_hets <- df_proband %>%
    group_by(chrom) %>%
    summarize(mean_hets = mean(n_het)) %>%
    arrange(chrom) %>% 
    bind_rows(data.frame(chrom = 23, mean_hets = mean(df_sd_errors$n_hets))) %>%
    arrange(chrom)

df_sizes <- data.frame(chrom = 1:23,
                       size = c(248956422, 242193529, 198295559,
                                190214555, 181538259, 170805979,
                                159345973, 145138636, 138394717,
                                133797422, 135086622, 133275309,
                                114364328, 107043718, 101991189,
                                90338345, 83257441, 80373285,
                                58617616, 64444167, 46709983,
                                50818468,156040895))
df_mean_hets <- df_mean_hets %>%
    left_join(df_sizes, by = "chrom") %>%
    mutate(hets_per_mb = mean_hets / (size / 1e6))

df_mean_hets %>%
    mutate(is_X = ifelse(chrom == 23, "X", "Autosome")) %>%
    ggplot(aes(x = size, y = mean_hets, color = is_X)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    theme_classic() +
    labs(x = "Chromosome Size (bp)",
         y = "Mean Number of Heterozygous Positions",
         color = "") +
    theme(text = element_text(size = 16)) + # larger axis labels
    theme(legend.text = element_text(size = 14), legend.title = element_text(size = 16)) +
    scale_color_manual(values = c(friendly_pal(pal_choice)[1], friendly_pal(pal_choice)[2])) + # add label next to each point
    geom_text(aes(label = ifelse(chrom == 23, "X", chrom)), vjust = -1, size = 5) + # remove legend
    theme(legend.position = "none")

# Mean number of hets by population
df_proband %>%
    group_by(SUPER, chrom) %>%
    summarize(mean_hets = mean(n_het)) %>%
    arrange(SUPER) %>%
    pivot_wider(names_from = SUPER, values_from = mean_hets) %>%
    knitr::kable()

df_sd_errors %>%
    group_by(pop) %>%
    summarize(mean_hets = mean(n_hets)) %>%
    arrange(pop) %>%
    pivot_wider(names_from = pop, values_from = mean_hets) %>%
    knitr::kable()


### 
df_all %>% head()

# I want to use the Friedman test to see if there are significant differences in switch rates between methods
# I want to do this for each chromosome separately
# We'll also limit our analysis to EUR samples

friedman_results <- vector("list", 23)
for (i in 1:23) {
    df_chr <- df_all %>%
        filter(chrom == i, pop == "EUR") %>%
        mutate(beagle_switch_rate = n_switch_beagle / n_het,
               eagle_switch_rate = n_switch_eagle / n_het,
               shapeit_switch_rate = n_switch_shapeit / n_het) %>%
        select(pair_id, beagle_switch_rate, eagle_switch_rate, shapeit_switch_rate) %>%
        pivot_longer(-pair_id, names_to = "Method", values_to = "Switch_Rate")
    friedman_results[[i]] <- df_chr %>%
        friedman_test(Switch_Rate ~ Method | pair_id) %>%
        mutate(chrom = i)
}

friedman_results_df <- bind_rows(friedman_results)
sum(friedman_results_df$p < 0.05)  # number of chromosomes with significant differences
friedman_results_df %>%
    arrange(chrom) %>%
    knitr::kable(digits = 4)

# As an alternative, let's try using a GLMM with a random effect for sample
# Let's also include a parameter for population -> would we also want to consider a population by method interaction?
library(lme4)
glmm_results <- vector("list", 23)
for (i in 1:23) {
    df_chr <- df_all %>%
        filter(chrom == i) %>% # get number of failures and successes for each method
        select(pair_id, pop, n_het, n_switch_beagle, n_switch_eagle, n_switch_shapeit) %>%
        pivot_longer(cols = -c(pair_id, pop, n_het), names_to = "Metric", values_to = "Count") %>% # recode Metric to Method
        separate(Metric, into = c("n", "type", "Method"), sep = "_") %>%
        select(-n, -type) %>%
        mutate(Success = n_het - Count) # compute number of successes
    glmm_obj <- glmer(cbind(Count, Success) ~ Method + pop + (1 | pair_id), data = df_chr, family = binomial(link = "logit"))
    null_glmm_obj <- glmer(cbind(Count, Success) ~ pop + (1 | pair_id), data = df_chr, family = binomial(link = "logit"))
    glmm_results[[i]] <- anova(null_glmm_obj, glmm_obj) %>% broom::tidy() %>%
        filter(term == "glmm_obj") %>%
        mutate(chrom = i)
}

glmm_results_df <- bind_rows(glmm_results)
sum(glmm_results_df$p.value < 0.05)  # number of chromosomes with
glmm_results_df %>%
    arrange(chrom) %>%
    knitr::kable(digits = 4)


# Enrichment of errors at CpG in X chromosome by population
df_enrich_X <- df_all %>%
    filter(chrom == 23) %>%
    mutate(beagle_switch = (n_switch_cpg_beagle / n_switch_beagle)/(n_het_cpg / n_het),
           eagle_switch = (n_switch_cpg_eagle / n_switch_eagle)/(n_het_cpg / n_het),
           shapeit_switch = (n_switch_cpg_shapeit / n_switch_shapeit)/(n_het_cpg / n_het),
           beagle_flip = (n_flip_cpg_beagle / n_flip_beagle)/(n_het_cpg / n_het),
           eagle_flip = (n_flip_cpg_eagle / n_flip_eagle)/(n_het_cpg / n_het),
           shapeit_flip = (n_flip_cpg_shapeit / n_flip_shapeit)/(n_het_cpg / n_het)) %>%
    select(pop, pair_id, beagle_switch:shapeit_flip) %>%
    pivot_longer(cols = -c(pop, pair_id), names_to = "Metric", values_to = "Enrichment") %>%
    separate(Metric, into = c("Method", "Error_Type"), sep = "_") %>%
    mutate(Error_Type = factor(Error_Type, levels = c("switch", "flip"))) %>%
    # recode switch to Switch  and flip to Flip 
    mutate(Error_Type = recode(Error_Type, "switch" = "Switch", "flip" = "Flip")) %>%
    mutate(Method = factor(Method, levels = c("beagle", "eagle", "shapeit"))) %>%
    mutate(Method = recode(Method, "beagle" = "Beagle", "eagle" = "Eagle", "shapeit" = "SHAPEIT"))

# Do we see a difference in enrichment across populations for switches and flips on X?
ggplot(df_enrich_X, aes(x = pop, y = Enrichment, color = Method)) +
    geom_boxplot() +
    facet_grid(rows = vars(Error_Type)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("Population") +
    ylab("Enrichment at CpG Sites") +
    scale_color_manual(values = friendly_pal(pal_choice)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
    theme(strip.text = element_text(size = 14)) +
    theme(legend.position = "top", legend.box = "horizontal") + # remove legend title
    theme(legend.text = element_text(size = 14)) +
    guides(color = guide_legend(override.aes = list(size = 5, 
        shape = 22, 
        linetype = "blank",
        fill = friendly_pal(pal_choice)[1:3]),
        title = "")) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") 
