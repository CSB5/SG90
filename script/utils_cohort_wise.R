pacman::p_load(poolr)


# get padj from association test for species of interest
get_padj_from_signif <- function(cohort_assoc_res, species_of_interest) {
  cohort_assoc_res <- cohort_assoc_res %>% 
    #mutate(org_list = str_remove(org_list, 's__')) %>% 
    filter(org_list %in% species_of_interest) 
  cohort_assoc_res$padj <- p.adjust(cohort_assoc_res$pval, 'fdr')
  return(cohort_assoc_res)
}

# get profiles for specific cohort & filter according to either using org_list_to_use or different n_nonzero
get_cohort_profiles <- function(profiles, cohort, n_nonzero) {
  metadata_cohort <- metadata %>% filter(Cohort == cohort)#, Race == 'Chinese')
  print(dim(metadata_cohort))
  
  cohort_profiles <- all_profiles_with_cre %>% dplyr::select(Index, metadata_cohort %>% pull(LibraryID))
  # based on org_list_to_use
  # cohort_profiles <- cohort_profiles %>% 
  #   filter(Index %in% org_list_to_use) %>% column_to_rownames('Index')
  # 
  # cohort_profiles_association <- merge(cohort_profiles %>% t() %>% as.data.frame() 
  #                                      %>% rownames_to_column('LibraryID'),
  #                                      metadata_cohort, by = 'LibraryID')
  # based on n_nonzero
  cohort_profiles_association <- transform_tax_profiles(cohort_profiles %>% dplyr::rename(org_list = Index),
                                                        metadata_cohort, n_nonzero)
  
  return(cohort_profiles_association)
}


# get profiles for pairwise cohorts without batch correction
get_paired_cohort_profiles <- function(profiles, cohort1, cohort2, n_nonzero) {
  metadata_cohort <- metadata %>% filter(Cohort == cohort1 | Cohort == cohort2)
                                                  # ,
                                                  # Race == 'Chinese')
  print(dim(metadata_cohort))
  
  cohort_profiles <- profiles %>% dplyr::select(Index, metadata_cohort %>% pull(LibraryID))
  
  
  # based on org_list_to_use
  # cohort_profiles <- cohort_profiles %>% 
  #   filter(Index %in% org_list_to_use) %>% column_to_rownames('Index')
  # 
  # cohort_profiles_association <- merge(cohort_profiles %>% t() %>% as.data.frame() 
  #                                      %>% rownames_to_column('LibraryID'),
  #                                      metadata_cohort, by = 'LibraryID')
  # based on n_nonzero
  cohort_profiles_association <- transform_tax_profiles(cohort_profiles %>% dplyr::rename(org_list = Index),
                                                        metadata_cohort, n_nonzero)
  
  return(cohort_profiles_association)
}


# label significant results for combined tables
get_signif_res <- function(combined_padjs, combined_betaestimates) {
  combined_padj_res <- combined_padjs %>% pivot_longer(-org_list, values_to = 'padj') %>% 
    mutate(signif = case_when(padj < 0.1 ~ 'Significant',
                              TRUE ~ 'Not Significant')) %>% 
    inner_join(combined_betaestimates %>% pivot_longer(-org_list, values_to = 'beta'),
               by = c('org_list', 'name'))
  return(combined_padj_res)
}


# get sign from significant taxa of glm
get_sign_from_signif_species <- function(glm_res) {
  return(glm_res %>% 
           mutate(sign = sign(betaestimate)) %>% 
           dplyr::select(org_list, sign))
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

