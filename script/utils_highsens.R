# Helper functions to get high-sensitivity list across 3 methods
# date: 20240227

library(tidyverse)
# get combined scores & only filter taxa w/ p-val < 0.05 in >= 1 method
# input: depending on which fisher's test data
get_combined_scores <- function(all_glm_assoc, sg90_cpe_glm_assoc, ca_glm_assoc) {
  fisher_table_with_dir <- ca_glm_assoc %>% 
    # if OR < 1 -> negative, OR > 1 -> positive
    mutate(betaestimate = case_when(betaestimate > 1 ~ 1,
                                    betaestimate == Inf ~ 1,
                                    betaestimate == 0 ~ -1,
                                    betaestimate < 1 ~ -1))
  
  combined_scores <- rbind(
    all_glm_assoc %>% mutate(Group = 'All Cohorts\nEM'), 
    sg90_cpe_glm_assoc %>% mutate(Group = 'SG90 + CPE\nEM'),
    fisher_table_with_dir %>% mutate(Group = 'All Cohorts\nPresence\nAbsence')) 
  
  combined_scores <- combined_scores %>% 
    # 0. remove stderror column
    dplyr::select(-stderror, -`2.5%`, -`97.5%`) %>% 
    # 1. filter only significant taxa
    #filter(org_list %in% signif_taxa) %>% 
    # 2. filter pval < 0.05
    filter(pval < 0.05) %>% 
    # 3. make long to plot
    pivot_longer(-c('org_list', 'Group')) %>% 
    # 4. remove space, _ from org_list
    mutate(org_list = str_replace_all(org_list, '_', ' '))
  return(combined_scores)
  
}

# calculate consistency score
get_cons_score <- function(combined_scores) {
  
  ## Add scores
  signif_taxa_per_group <- combined_scores %>% 
    filter(name == 'padj', value < 0.05)
  
  # get betaestimate only & convert to sign for direction
  score_table <- combined_scores %>% 
    filter(name == 'betaestimate') %>% 
    right_join(signif_taxa_per_group %>% 
                 dplyr::select(org_list, Group, value), 
               by = c('org_list', 'Group')) %>% 
    dplyr::select(-name, -value.y) %>% 
    mutate(value.x = sign(value.x)) %>% 
    pivot_wider(names_from = 'Group', values_from = 'value.x') %>% 
    replace(is.na(.), 0) %>% 
    mutate(across(where(is.numeric), ~ . * 2)) %>% 
    mutate(Score = rowMeans(abs(across(c(where(is.numeric)))))/2)
  
  return(score_table)
}

# filter with only those with equal signs for the different methods
filter_same_dir <- function(score_table) {
  score_table_filt <- score_table %>% 
    filter(Score >= 0.5) %>% 
    dplyr::select(-Score) %>%
    # get sum
    mutate(check_sign = rowSums(across(c(where(is.numeric))))) %>% 
    # only sum = 6/-6/4/-4 are considered
    filter(check_sign == 6 | check_sign == -6 | check_sign == 4 | check_sign == -4)
  
  return(score_table_filt)
}

# plot combined heatmap
# input: depending on which fisher's test data
plot_heatmap <- function(all_glm_assoc, sg90_cpe_glm_assoc, ca_glm_assoc) {
  combined_scores <- get_combined_scores(all_glm_assoc, sg90_cpe_glm_assoc,
                                         ca_glm_assoc)
  
  # get no of signif taxa per Group
  combined_scores %>% 
    filter(name == 'padj', value < 0.05) %>% 
    count(Group)
  
  # get beta only
  combined_scores_beta <- combined_scores %>% filter(name == 'betaestimate') %>% 
    pivot_wider(names_from = 'Group', values_from = 'value') %>% 
    dplyr::select(-name) %>% 
    pivot_longer(-org_list, names_to = 'Group') %>% 
    mutate(value = as.factor(sign(value)),
           org_list = factor(org_list, levels = rev(sort(unique(combined_scores$org_list)))),
           Group = factor(Group, levels = c('All Cohorts\nEM',
                                            'SG90 + CPE\nEM',
                                            'All Cohorts\nPresence\nAbsence')))
  
  # get significant taxa only >= 1 method
  combined_scores_padj <- combined_scores %>% 
    filter(name == 'padj', value < 0.05) %>% 
    dplyr::select(-value)
  
  ## get consistency score
  score_table <- get_cons_score(combined_scores)
  
  # filter score >= 0.5
  score_filter <- score_table %>% 
    filter(Score >= 0.5)
  
  # and filter only those with same direction
  methods_ca_combined_assoc <- get_combined_scores(all_glm_assoc,
                                                   sg90_cpe_glm_assoc,
                                                   ca_glm_assoc)
  methods_score_ca_combined_assoc <- get_cons_score(methods_ca_combined_assoc)
  methods_score_ca_combined_same_dir_assoc <- filter_same_dir(methods_score_ca_combined_assoc)
  
  # list used to filter: either score >= 0.5 or score >= 0.5 + same direction
  filtered_list <- sort(methods_score_ca_combined_same_dir_assoc$org_list)
  
  # create heatmap: filter based on significant >= 1 method or >= 2 methods
  methods_heatmap <- combined_scores_beta %>% 
    # filter only taxa significant in >= 1 method
    #filter(org_list %in% combined_scores_padj$org_list) %>%
    
    # filter taxa w/ scores >= 0.5
    filter(org_list %in% filtered_list) %>%
    
    ggplot(aes(x = Group, y = org_list, fill = value)) + 
    geom_tile() + 
    geom_tile(colour = 'black', linewidth = 0.3, show.legend = F, alpha = 0.5) +
    geom_text(data = combined_scores_padj %>% 
                filter(org_list %in% filtered_list),
              aes(x = Group, y = org_list, label = 'X'), size = 4,
              color = 'white', inherit.aes = F) +
    scale_fill_manual(values = c('blue', 'red'), na.value = 'white',
                      labels = c('Negative', 'Positive', NA),
                      na.translate=F) +
    scale_x_discrete(expand = c(0, 0)) + 
    #scale_y_discrete(expand = c(0, 0)) +
    theme_bw() + 
    #coord_fixed(ratio=0.1) +
    labs(fill = 'Direction') + 
    # scale_x_discrete(labels = c('Cont.', 'Strat.',
    #                             'Cont.', 'Strat.',
    #                             'Cont.', 'Strat.',
    #                             'Cont.', 'Strat.',
    #                             'Strat.',
    #                             'Cont.', 'Strat.',
    #                             'Cont.', 'Strat.')) +
    theme(axis.text.y = element_text(face = 'italic', size = 12,
                                     margin = margin(r=0),
                                     hjust = 1),
          axis.text.x = element_text(size = 14),
          axis.title = element_blank(), 
          plot.background = element_blank(), 
          panel.border = element_blank(), 
          panel.grid = element_blank(), 
          text = element_text(family = 'Arial'), 
          legend.position = 'right', 
          legend.text = element_text(size = 14),
          axis.ticks = element_blank(),
          #plot.margin=grid::unit(c(0,0,0,0), "mm"))
          #plot.margin=margin(0, 0, 0, 0, 'cm'),
          panel.spacing = unit(0, 'cm'))
  
  return(methods_heatmap)
}

# get median p-adj values for combined results
get_median_padj <- function(all_glm_assoc, sg90_cpe_glm_assoc, ca_glm_assoc) {
  combined_padj <- rbind(all_glm_assoc %>% dplyr::select(org_list, padj) %>% 
                           mutate(Group = 'All Cohorts\nEM'),
                         sg90_cpe_glm_assoc %>% dplyr::select(org_list, padj) %>% 
                           mutate(Group = 'SG90 + CPE\nEM'),
                         ca_glm_assoc %>% dplyr::select(org_list, padj) %>% 
                           mutate(Group = 'All Cohorts\nPrevalence'))
  # calculate median
  combined_padj <- combined_padj %>% 
    group_by(org_list) %>% 
    summarize(n = n(),
              padj_median = if(n > 2) median(padj) else max(padj))
  return(combined_padj)
}

# plot distribution of values based on different thresholds
plot_distn <- function(tax_profiles, thresholds) {
  df_long <- reshape2::melt(tax_profiles)
  
  # Filter out rows with 0 values
  df_filtered <- df_long[df_long$value > 0, ]
  
  # Set up the ggplot with log-transformed x-axis and specific breaks
  p <- ggplot(df_filtered, aes(x = log10(value))) +
    geom_density() +
    facet_wrap(~org_list, scales = "free") +
    scale_x_continuous(
      breaks = log10(c(0.0005, 0.001, 0.005, 0.01)),
      labels = c("0.05%", "0.1%", "0.5%", "1%")
    ) +
    labs(x = "Log-transformed Relative Abundance", y = "Density") +
    theme_minimal()
  
  # Add smooth-fitted line
  p <- p + geom_density(aes(y = ..scaled..), color = "red", linetype = "dashed")
  
  print(p)
}
