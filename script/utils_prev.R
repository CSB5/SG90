# Helper functions to run Cochran-Armitage trend test
# date: 20240227

library(ggpubr)
library(rstatix)
library(DescTools)
library(tidyverse)


# remove low abundance taxa
remove_low_abund <- function(taxonomic_profiles, species_thr) {
  # convert colnames to rownames
  taxonomic_profiles <- taxonomic_profiles %>% 
    tibble::column_to_rownames('org_list')
  taxonomic_profiles[is.na(taxonomic_profiles)| taxonomic_profiles < species_thr] <- 0 
  
  ## remove taxa with 0s in all samples
  taxonomic_profiles <- taxonomic_profiles[rowSums(taxonomic_profiles) != 0, ] 
  
  ## renormalize to 100%
  taxonomic_profiles <- data.frame(apply(taxonomic_profiles, 2, function(x) x/sum(x)) * 100)
  
  return(taxonomic_profiles)
  
}

# filter only original list of taxa
get_profiles_for_test <- function(taxonomic_profiles_filt, org_list) {
  taxonomic_profiles_for_test <- 
    taxonomic_profiles_filt %>% 
    filter(row.names(.) %in% org_list) %>% 
    tibble::rownames_to_column('org_list')
  return(taxonomic_profiles_for_test)
}

# get count per taxa per age group
get_count_age_group <- function(tax_profiles, metadata) {
  tax_count_age_group <- tax_profiles %>% 
    pivot_longer(-org_list, names_to = 'LibraryID') %>% 
    left_join(metadata %>% dplyr::select(LibraryID, AgeG, Cohort)) %>% 
    group_by(org_list, AgeG) %>%
    # check for presence/absence
    filter(value > 0) %>%
    dplyr::count() %>% 
    pivot_wider(names_from = 'AgeG', values_from = 'n') %>% 
    replace(is.na(.), 0) %>% 
    mutate(org_list = str_remove(org_list, 's__'))
  return(tax_count_age_group)
}


# run Cochran-Armitage test
run_cohran_armitage <- function(tax_count_age_group, metadata_to_test) {
  org_list <- tax_count_age_group$org_list
  
  # get total count for each age group
  n21_40 <- metadata_to_test %>% dplyr::count(AgeG) %>% filter(AgeG == '21-40') %>% pull(n)
  n41_60 <- metadata_to_test %>% dplyr::count(AgeG) %>% filter(AgeG == '41-60') %>% pull(n)
  n61_70 <- metadata_to_test %>% dplyr::count(AgeG) %>% filter(AgeG == '61-70') %>% pull(n)
  n71_80 <- metadata_to_test %>% dplyr::count(AgeG) %>% filter(AgeG == '71-80') %>% pull(n)
  n81_90 <- metadata_to_test %>% dplyr::count(AgeG) %>% filter(AgeG == '81-90') %>% pull(n)
  n91_100 <- metadata_to_test %>% dplyr::count(AgeG) %>% filter(AgeG == '91-100') %>% pull(n)
  
  # generate empty table to be populated later
  pval_table <- data.frame(org_list=c(), pval = c())
  
  # calculate fisher's exact test for each species
  for (taxa in org_list) {
    #print(taxa)
    tax_count <- tax_count_age_group %>% filter(org_list == taxa)
    #print(tax_count)
    count_21_40 <- tax_count[['21-40']]
    count_41_60 <- tax_count[['41-60']]
    count_61_70 <- tax_count[['61-70']]
    count_71_80 <- tax_count[['71-80']]
    count_81_90 <- tax_count[['81-90']]
    count_91_100 <- tax_count[['91-100']]
    
    # cont_table <- data.frame('presence' = c(count_21_40, count_41_60, count_61_70,
    #                                         count_71_80, count_81_90, count_91_100),
    #                          'total' = c(n21_40, n41_60, n61_70, n71_80, n81_90,
    #                                      n91_100),
    #                          row.names = c('21-40', '41-60', '61-70', '71-80', '81-90',
    #                                        '91-100'))
    
    cont_table <- data.frame('presence' = c(count_21_40, count_41_60, count_61_70,
                                            count_71_80, count_81_90, count_91_100),
                             'absence' = c(n21_40 - count_21_40, n41_60 - count_41_60,
                                           n61_70 - count_61_70, n71_80 - count_71_80,
                                           n81_90 - count_81_90, n91_100 - count_91_100),
                             row.names = c('21-40', '41-60', '61-70', '71-80', '81-90',
                                           '91-100'))    
    # run Cochran-Armitage test
    # ca_pval <- prop.trend.test(t(cont_table)[1, ], 
    #                            t(cont_table)[2, ])$p.value
    ca_test <- CochranArmitageTest(cont_table, alternative = 'one.sided')
    
    pval_table <- rbind(pval_table, data.frame(org_list = taxa, pval = ca_test$p.value,
                                               betaestimate = ca_test$statistic))
  }
  
  pval_table$padj <- p.adjust(pval_table$pval, 'fdr')
  pval_table <- pval_table %>% 
    # dummy column to match glm format
    mutate(stderror = 1,
           `2.5%` = 1,
           `97.5%` = 1) %>% 
    dplyr::select(org_list, pval, betaestimate, stderror, padj, `2.5%`, `97.5%`)
  rownames(pval_table) <- NULL
  
  #print(pval_table)
  return(pval_table)
  
}

# run JoncheereTerpstraTest
run_jt <- function(tax_count_age_group, metadata_to_test) {
  org_list <- tax_count_age_group$org_list
  
  # get total count for each age group
  n21_40 <- metadata_to_test %>% dplyr::count(AgeG) %>% filter(AgeG == '21-40') %>% pull(n)
  n41_60 <- metadata_to_test %>% dplyr::count(AgeG) %>% filter(AgeG == '41-60') %>% pull(n)
  n61_70 <- metadata_to_test %>% dplyr::count(AgeG) %>% filter(AgeG == '61-70') %>% pull(n)
  n71_80 <- metadata_to_test %>% dplyr::count(AgeG) %>% filter(AgeG == '71-80') %>% pull(n)
  n81_90 <- metadata_to_test %>% dplyr::count(AgeG) %>% filter(AgeG == '81-90') %>% pull(n)
  n91_100 <- metadata_to_test %>% dplyr::count(AgeG) %>% filter(AgeG == '91-100') %>% pull(n)
  
  # generate empty table to be populated later
  pval_table <- data.frame(org_list=c(), pval = c())
  
  # calculate fisher's exact test for each species
  for (taxa in org_list) {
    #print(taxa)
    tax_count <- tax_count_age_group %>% filter(org_list == taxa)
    #print(tax_count)
    count_21_40 <- tax_count[['21-40']]
    count_41_60 <- tax_count[['41-60']]
    count_61_70 <- tax_count[['61-70']]
    count_71_80 <- tax_count[['71-80']]
    count_81_90 <- tax_count[['81-90']]
    count_91_100 <- tax_count[['91-100']]
    
    # cont_table <- data.frame('presence' = c(count_21_40, count_41_60, count_61_70,
    #                                         count_71_80, count_81_90, count_91_100),
    #                          'total' = c(n21_40, n41_60, n61_70, n71_80, n81_90,
    #                                      n91_100),
    #                          row.names = c('21-40', '41-60', '61-70', '71-80', '81-90',
    #                                        '91-100'))
    
    cont_table <- data.frame('presence' = c(count_21_40, count_41_60, count_61_70,
                                            count_71_80, count_81_90, count_91_100),
                             'absence' = c(n21_40 - count_21_40, n41_60 - count_41_60,
                                           n61_70 - count_61_70, n71_80 - count_71_80,
                                           n81_90 - count_81_90, n91_100 - count_91_100),
                             row.names = c('21-40', '41-60', '61-70', '71-80', '81-90',
                                           '91-100'))    
    # run JT test
    # update cont_table
    cont_table <- cont_table %>% tibble::rownames_to_column('Age') %>% 
      mutate(propn = presence/(presence + absence)) %>% 
      dplyr::select(Age, propn)
    # check decreasing
    jt_decreasing <- JonckheereTerpstraTest(propn ~ Age, cont_table, 
                                            alternative = 'decreasing',
                                            nperm = 1000)
    
    # check increasing
    jt_increasing <- JonckheereTerpstraTest(propn ~ Age, cont_table, 
                                            alternative = 'increasing',
                                            nperm = 1000)
    
    pval <- min(jt_increasing$p.value, jt_decreasing$p.value)
    sign <- ifelse(jt_increasing$p.value < jt_decreasing$p.value, 1, -1)
    pval_table <- rbind(pval_table, data.frame(org_list = taxa, pval = pval,
                                               betaestimate = sign))
  }
  
  pval_table$padj <- p.adjust(pval_table$pval, 'fdr')
  pval_table <- pval_table %>% 
    # dummy column to match glm format
    mutate(stderror = 1,
           `2.5%` = 1,
           `97.5%` = 1) %>% 
    dplyr::select(org_list, pval, betaestimate, stderror, padj, `2.5%`, `97.5%`)
  rownames(pval_table) <- NULL
  
  #print(pval_table)
  return(pval_table)
  
}

plot_desc <- function(tax_count_age_group, metadata_to_test) {
  org_list <- tax_count_age_group$org_list
  
  # get total count for each age group
  n21_40 <- metadata_to_test %>% dplyr::count(AgeG) %>% filter(AgeG == '21-40') %>% pull(n)
  n41_60 <- metadata_to_test %>% dplyr::count(AgeG) %>% filter(AgeG == '41-60') %>% pull(n)
  n61_70 <- metadata_to_test %>% dplyr::count(AgeG) %>% filter(AgeG == '61-70') %>% pull(n)
  n71_80 <- metadata_to_test %>% dplyr::count(AgeG) %>% filter(AgeG == '71-80') %>% pull(n)
  n81_90 <- metadata_to_test %>% dplyr::count(AgeG) %>% filter(AgeG == '81-90') %>% pull(n)
  n91_100 <- metadata_to_test %>% dplyr::count(AgeG) %>% filter(AgeG == '91-100') %>% pull(n)
  
  # generate empty table to be populated later
  pval_table <- data.frame(org_list=c(), pval = c())
  
  # calculate fisher's exact test for each species
  for (taxa in org_list) {
    print(taxa)
    #print(taxa)
    tax_count <- tax_count_age_group %>% filter(org_list == taxa)
    #print(tax_count)
    count_21_40 <- tax_count[['21-40']]
    count_41_60 <- tax_count[['41-60']]
    count_61_70 <- tax_count[['61-70']]
    count_71_80 <- tax_count[['71-80']]
    count_81_90 <- tax_count[['81-90']]
    count_91_100 <- tax_count[['91-100']]
    
    cont_table <- data.frame(age = rep(c(0, 1, 2, 3, 4, 5), 
                                       c(n21_40, n41_60, n61_70, 
                                         n71_80, n81_90, n91_100)),
                       presence = c(rep(c(0, 1), c(count_21_40, n21_40-count_21_40)), 
                                    rep(c(0, 1), c(count_41_60, n41_60-count_41_60)), 
                                    rep(c(0, 1), c(count_61_70, n61_70-count_61_70)), 
                                    rep(c(0, 1), c(count_71_80, n71_80-count_71_80)), 
                                    rep(c(0, 1), c(count_81_90, n81_90-count_81_90)),
                                    rep(c(0, 1), c(count_91_100, n91_100-count_91_100))))
    print(table(cont_table$age, cont_table$presence))
    print(DescTools::Desc(table(cont_table$age, cont_table$presence), main = taxa))
  }
}
