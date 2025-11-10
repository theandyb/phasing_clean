source("code/figures/load_results.R")

# load mb error rates
df_sd_errors_mb <- load_sd_errors_mb()

# Range of switch rates across bins
df_sd_errors_mb %>%
    group_by(bin) %>%
    mutate(Beagle_Switch =  sum(beagle_switch),
            Eagle_Switch = sum(eagle_switch),
            Shapeit_Switch = sum(shapeit_switch),
            Beagle_Flip = sum(beagle_flip),
            Eagle_Flip = sum(eagle_flip),
            Shapeit_Flip = sum(shapeit_flip),
           n_het = sum(n_het)) %>%
    ungroup() %>%
    mutate(Beagle_Switch_Rate = Beagle_Switch/n_het,
           Eagle_Switch_Rate = Eagle_Switch/n_het,
           Shapeit_Switch_Rate = Shapeit_Switch/n_het,
           Beagle_Flip_Rate = Beagle_Flip/n_het,
           Eagle_Flip_Rate = Eagle_Flip/n_het,
           Shapeit_Flip_Rate = Shapeit_Flip/n_het) %>%
    summary()

# Same, but exclude bins 58 to 63
df_sd_errors_mb %>%
    filter(bin < 58 | bin > 63) %>% 
    group_by(bin) %>%
    mutate(Beagle_Switch =  sum(beagle_switch),
           Eagle_Switch = sum(eagle_switch),
           Shapeit_Switch = sum(shapeit_switch),
           Beagle_Flip = sum(beagle_flip),
           Eagle_Flip = sum(eagle_flip),
           Shapeit_Flip = sum(shapeit_flip),
           n_het = sum(n_het)) %>%
    ungroup() %>%
    mutate(Beagle_Switch_Rate = Beagle_Switch/n_het,
           Eagle_Switch_Rate = Eagle_Switch/n_het,
           Shapeit_Switch_Rate = Shapeit_Switch/n_het,
           Beagle_Flip_Rate = Beagle_Flip/n_het,
           Eagle_Flip_Rate = Eagle_Flip/n_het,
           Shapeit_Flip_Rate = Shapeit_Flip/n_het) %>%
    summary()
