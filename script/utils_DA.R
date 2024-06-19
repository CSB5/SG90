library(ALDEx2)
library(ANCOMBC)
library(Maaslin2)
library(VennDiagram)
library(ggrepel)
library(tidyverse)

# read cre profiles - replace x with org_list & '.' in CON.XXX
read_cre <- function(filepath) {
  profiles <- read.csv(filepath) %>% 
    rename(org_list = 'X')
  names(profiles) <- gsub("\\.", "", names(profiles))
  return(profiles)
}

# remove NA & low abundance; then renormalize & reorder colnames
preprocess_tax_profiles <- function(taxonomic_profiles, species_thr) {
  
  #Since some organisms are not present in CRE, there will be NA introduced in such cases
  #We substitute NA for 0
  taxonomic_profiles[is.na(taxonomic_profiles)] <- 0
  
  # Filtering the profiles by setting abundance < 0.1 (1) as 0 and 
  #removing the organisms which are zero across all samples
  original_profiles_copy <- taxonomic_profiles
  taxonomic_profiles[taxonomic_profiles < species_thr] <- 0
  taxonomic_profiles_filtered <-
    taxonomic_profiles[rowSums(taxonomic_profiles) != 0, ]
  dim(taxonomic_profiles_filtered) 
  
  # Renormalising to 1 after removing organisms
  taxonomic_profiles_filtered_2 <- taxonomic_profiles_filtered %>% apply(2, function(x) x / sum(x))
  #write.csv(taxonomic_profiles_filtered_2, 'Original_profiles_aarthi_filter_362_orgs.csv')
  # taxonomic_profiles_filtered_1 <- taxonomic_profiles_filtered/colSums(taxonomic_profiles_filtered)
  # tax3 <- sweep(taxonomic_profiles_filtered, 2, colSums(taxonomic_profiles_filtered), FUN="/")
  
  #tax_profiles_prop <- taxonomic_profiles_incre_filtered/100
  
  tax_profiles_prop <- taxonomic_profiles_filtered_2
  tax_profiles_prop <-
    tax_profiles_prop[, order(colnames(tax_profiles_prop))]
  (dim(tax_profiles_prop))
  return(tax_profiles_prop)
}

# filter profiles & merge w/ metadata (not normalized to 100%)
transform_tax_profiles <- function(tax_profiles, metadata, n_nonzero) {
  tax_profiles_t <- tax_profiles %>% tibble::column_to_rownames('org_list') %>% 
    t() %>% 
    as.data.frame() 
  
  # Check number of non zero entries
  nonzero_in_samples <- data.frame(nnzorg = base::colSums(tax_profiles_t !=0))
  ggplot(nonzero_in_samples, aes(x=nnzorg)) + 
    geom_histogram(color="darkblue", fill="lightblue") + 
    scale_x_continuous(breaks = scales::pretty_breaks(n=7)) + 
    xlab('Number of samples') + ylab('Count of organisms with non-zero abundances')
  
  orgsinmorethanfiftysamples <- nonzero_in_samples %>% filter(nnzorg >= n_nonzero)
  
  # # get taxa means >= 0.1%
  # taxa_means <- colMeans(tax_profiles_t) %>% 
  #   as.data.frame() %>% 
  #   filter(. > 0.1)
  # 
  # taxa_filtered <- intersect(rownames(taxa_means),
  #                            rownames(orgsinmorethanfiftysamples))
  taxa_with_morethan50 <- tax_profiles_t %>% 
   dplyr::select(rownames(orgsinmorethanfiftysamples)) %>%
    # dplyr::select(taxa_filtered) %>% 
    rownames_to_column('LibraryID')
  
  data_for_association <- merge(taxa_with_morethan50, metadata, by = 'LibraryID')
  return(data_for_association)
}

get_median_of_taxa <- function(tax_profiles, n_nonzero) {
  tax_profiles_t <- tax_profiles %>% column_to_rownames('org_list') %>% 
    t() %>% 
    as.data.frame() 
  
  # Check number of non zero entries
  nonzero_in_samples <- data.frame(nnzorg = colSums(tax_profiles_t !=0))
  orgsinmorethanfiftysamples <- nonzero_in_samples %>% filter(nnzorg >= n_nonzero)
  
  median_of_taxa <- tax_profiles_t %>% 
    dplyr::select(rownames(orgsinmorethanfiftysamples)) %>%  
    as.matrix() %>% matrixStats::colMedians()
  median_of_taxa_with_name <- data.frame(rownames(orgsinmorethanfiftysamples), 
                                         median_of_taxa)
  names(median_of_taxa_with_name)[1] <- 'org_list'
  median_of_taxa_with_name <- median_of_taxa_with_name %>% 
    mutate(org_list = str_remove(org_list, 's__'))
  return(median_of_taxa_with_name)
}

# Function to parse GLM results
param_with_age_correction <- function(org_list_to_use, resultssum, fname){
  pval <- NULL
  org_list <- NULL
  estimate <- NULL
  stderror <- NULL
  betaestimate <- NULL
  stderror <- NULL
  for (i in seq_along(1:length(resultssum)))
  {
    #print(paste0(org_list_to_use[i], ' ', resultssum[[i]][35]))
    if ("x" %in% rownames(resultssum[[i]])){
      pval <-
        rbind(pval, resultssum[[i]]["x",][4])
      org_list <- rbind(org_list, org_list_to_use[i])
      betaestimate <-
        rbind(betaestimate, resultssum[[i]]["x",][1])
      stderror <-
        rbind(stderror, resultssum[[i]]["x",][2])
    }
    else{
      
      print(org_list_to_use[i])}
  }
  
  all_df <- data.frame(org_list, pval, betaestimate, stderror)
  colnames(all_df) <- c('org_list', 'pval', 'betaestimate', 'stderror')
  all_df$padj <- p.adjust(all_df$pval, 'fdr')
  all_df$org_list <- str_remove(all_df$org_list, 's__')
  
  #write.csv(all_df, here('results_0.1',fname), row.names = FALSE)
  return(all_df)
}

param_with_age_correction_flip <- function(org_list_to_use, resultssum, fname){
  pval <- NULL
  org_list <- NULL
  estimate <- NULL
  stderror <- NULL
  betaestimate <- NULL
  stderror <- NULL
  for (i in seq_along(1:length(resultssum)))
  {
    #print(paste0(org_list_to_use[i], ' ', resultssum[[i]][35]))
    if ("data_for_association$Age" %in% rownames(resultssum[[i]])){
      pval <-
        rbind(pval, resultssum[[i]]["data_for_association$Age",][4])
      org_list <- rbind(org_list, org_list_to_use[i])
      betaestimate <-
        rbind(betaestimate, resultssum[[i]]["data_for_association$Age",][1])
      stderror <-
        rbind(stderror, resultssum[[i]]["data_for_association$Age",][2])
    }
    else{
      
      print(org_list_to_use[i])}
  }
  
  all_df <- data.frame(org_list, pval, betaestimate, stderror)
  colnames(all_df) <- c('org_list', 'pval', 'betaestimate', 'stderror')
  all_df$padj <- p.adjust(all_df$pval, 'fdr')
  
  #write.csv(all_df, here('results',fname), row.names = FALSE)
  return(all_df)
}

param_with_ageG_correction_flip <- function(org_list_to_use, resultssum, age_group){
  pval <- NULL
  org_list <- NULL
  estimate <- NULL
  stderror <- NULL
  betaestimate <- NULL
  stderror <- NULL
  for (i in seq_along(1:length(resultssum)))
  {
    #print(paste0(org_list_to_use[i], ' ', resultssum[[i]][35]))
    age_group_data <- paste("data_for_association$AgeG", age_group, sep='')
    #print(age_group_data)
    if (age_group_data %in% rownames(resultssum[[i]])){
      pval <-
        rbind(pval, resultssum[[i]][age_group_data,][4])
      org_list <- rbind(org_list, org_list_to_use[i])
      betaestimate <-
        rbind(betaestimate, resultssum[[i]][age_group_data,][1])
      stderror <-
        rbind(stderror, resultssum[[i]][age_group_data,][2])
    }
    else{
      
      print(org_list_to_use[i])}
  }
  
  all_df <- data.frame(org_list, pval, betaestimate, stderror)
  colnames(all_df) <- c('org_list', 'pval', 'betaestimate', 'stderror')
  all_df$padj <- p.adjust(all_df$pval, 'fdr')
  
  #write.csv(all_df, here('results',fname), row.names = FALSE)
  return(all_df)
}

# run GLM - age as dependent variable
run_glm <- function(data_for_association) {
  print(dim(data_for_association))
  org_list_to_use <- colnames(data_for_association)[grepl('__', colnames(data_for_association))]
  print(length(org_list_to_use))
  species_with_age <-
    lapply(data_for_association[, org_list_to_use], function(x)
      coefficients(summary(
        glm(
          data_for_association$Age ~ x +  data_for_association$FBGlu +  data_for_association$Tgly +
            data_for_association$TCHOL + 
            data_for_association$HDL + 
            #data_for_association$LDL + 
            data_for_association$BMI + data_for_association$Gender,
          data = data_for_association
        )
      )))
  
  # adjust for multiple corrections
  species_with_age_results <- param_with_age_correction(org_list_to_use, species_with_age, fname) %>%
    mutate(org_list = str_remove(org_list, "s__"))
  
  # get statistically significant taxa
  # species_with_age_results <- species_with_age_results %>%
  #   filter(padj < 0.05)
  
  # add 95% confint
  species_with_age_results <- species_with_age_results %>% 
    mutate(`2.5%` = betaestimate - 1.96 * stderror, 
           `97.5%` = betaestimate + 1.96 * stderror)
  return(species_with_age_results)
}

run_glm_age_only <- function(data_for_association, fname) {
  print(dim(data_for_association))
  org_list_to_use <- colnames(data_for_association)[grepl('__', colnames(data_for_association))]
  print(length(org_list_to_use))
  species_with_age <-
    lapply(data_for_association[, org_list_to_use], function(x)
      coefficients(summary(
        glm(
          data_for_association$Age ~ x + data_for_association$Gender,
          data = data_for_association
        )
      )))
  
  # adjust for multiple corrections
  species_with_age_results <- param_with_age_correction(org_list_to_use, species_with_age, fname) %>%
    mutate(org_list = str_remove(org_list, "s__"))
  
  # get statistically significant taxa
  # species_with_age_results <- species_with_age_results %>%
  #   filter(padj < 0.05)
  return(species_with_age_results)
}

# Volcano plots
plot_volcano <- function(taxa_profiles, n_nonzero, species_results, title) {
  median_of_taxa_with_name <- get_median_of_taxa(taxa_profiles, n_nonzero)
  species_results$Significant <-
    ifelse(species_results$padj < 0.05, "FDR < 0.05", "Not Sig")
  colnames(median_of_taxa_with_name)[1] <- c('org_list')
  species_results_merge <-
    merge(species_results, median_of_taxa_with_name, by = 'org_list')
  #species_results_merge <- merge(species_results_merge, median_abun, by='org_list')
  # To plot mean > 1
  #species_results_merge$meantohigh <- ifelse(species_results_merge$meanabundance > 1, "Mean > 1", "NIL") 
  species_results_merge$highlight1 <-
    ifelse(species_results_merge$median_of_taxa > 0 &
             species_results$padj < 0.05, "2", 
           ifelse(species_results_merge$median_of_taxa == 0 &
                    species_results$padj < 0.05 |species_results$padj == 0.05 , "1", "0"))
  # species_results_merge$highlight1 <-
  #   ifelse(species_results_merge$median_of_taxa > 0 &
  #            species_results$padj < 0.05, "2", ifelse(species_results_merge$median_of_taxa == 0 &
  #            species_results$padj < 0.05 |species_results$padj == 0.05 , "1", species_results$padj < 0.05 |species_results$padj == 0.05 , "1", "0"))
  
  #species_to_mark <- c('s__Ruminococcus_obeum','s__Coprococcus_catus', 's__Alistipes_senegalensis', 's__Dorea_formicigenerans')    
  species_results_merge$highlight2 <-  ifelse(species_results_merge$betaestimate < -12 , "3", "0")
  #ifelse(species_results_merge$betaestimate > 20, "3",                                
  species_results_merge$highlight3 <- ifelse(species_results_merge$org_list == 'Alistipes_senegalensis', '-1', 0)            
  
  
  #species_results_merge$org_list <- str_replace_all(species_results_merge$org_list, 's__', '')
  species_results_merge$org_list <- str_replace_all(species_results_merge$org_list, '_', ' ')
  
  volcanoplot <- ggplot(species_results_merge, aes(x = betaestimate, y = -log10(padj))) +
    geom_point(aes(color = highlight1), size=4) +
    geom_vline(
      xintercept = c(0),
      col = "black",
      linetype = "solid",
      size = 0.5
    ) +
    theme_bw() + theme(legend.position = "bottom") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(face="bold", size=12),
          axis.text.x = element_text(face="bold",size=14),
          axis.text.y = element_text(face="bold",size=14),
          axis.text = element_text(face="bold",size=14),
          legend.position="none",
          plot.title = element_text(hjust = 0.5, size=14),
          legend.text = element_text(size=10)) +
    geom_text_repel(
      data = subset(species_results_merge, highlight1 ==2 ), #Changed from ==2
      aes(label = org_list),
      size = 6, xlim=c(-45, 20),ylim=c(0,10),
      box.padding = unit(0.75, "lines"),
      point.padding = unit(0.35, "lines"), max.overlaps = 30, fontface = 'italic',nudge_x = .15,
      nudge_y = .5
    ) +
    geom_text_repel(
      data = subset(species_results_merge, highlight2 >2),
      aes(label = org_list),
      size = 6, ylim=c(0,10), xlim=c(-45, 20),
      box.padding = unit(0.65, "lines"),
      point.padding = unit(0.65, "lines"), max.overlaps = 30, fontface = 'italic',nudge_x = .15,
      nudge_y = .5
    )+
    geom_text_repel(
      data = subset(species_results_merge, highlight3 <0),
      aes(label = org_list),
      size = 6, ylim=c(0,10), xlim=c( 20,0),
      box.padding = unit(0.65, "lines"),
      point.padding = unit(0.65, "lines"), max.overlaps = 30, fontface = 'italic', nudge_x = .15,
      nudge_y = .65
    )+ geom_hline(
      yintercept = -log10(0.05),
      col = "black",
      linetype = "dashed",
      size = 0.5
    ) + geom_hline(
      yintercept = -log10(0.1),
      col = "black",
      linetype = "dashed",
      size = 0.5
    ) + scale_y_continuous(breaks = scales::pretty_breaks(5)) +
    scale_x_continuous(breaks = scales::pretty_breaks(2))+ ggtitle(title) +
    labs(x= "\u03b2 coefficient", y= expression(paste(-log[10]," p-value")))+ scale_color_manual('', values = c("0"="grey","1"="green", "2"= "red"), labels = c('p-value>0.05', 'p-value \u2264 0.05 (med = 0)', 'p-value \u2264 0.05 (med > 0)' ))
  print(volcanoplot)
  }

# taxa as dependent variable
run_glm_flip <- function(data_for_association, fname) {
  org_list_to_use <- colnames(data_for_association)[grepl('__', colnames(data_for_association))]
  species_with_age_flip <-
    lapply(data_for_association[, org_list_to_use], function(x)
      coefficients(summary(
        glm(
          x ~ data_for_association$Age +  data_for_association$FBGlu + 
            data_for_association$Tgly + data_for_association$TCHOL + 
            data_for_association$HDL + data_for_association$LDL + 
            data_for_association$BMI + data_for_association$Gender,
          data = data_for_association
        )
      )))
  
  species_results_flip <- param_with_age_correction_flip(org_list_to_use, species_with_age_flip, 'species_age_associations_flip.csv')
  species_results_flip <- species_results_flip %>% 
    mutate(org_list = str_remove(org_list, "s__")) 
  # %>% 
  #   filter(padj < 0.05)
  return(species_results_flip)
  }

run_glm_flip_ageG <- function(data_for_association, age_group) {
  org_list_to_use <- colnames(data_for_association)[grepl('__', colnames(data_for_association))]
  species_with_age_flip <-
    lapply(data_for_association[, org_list_to_use], function(x)
      coefficients(summary(
        glm(
          x ~ data_for_association$AgeG +  data_for_association$FBGlu + 
            data_for_association$Tgly + data_for_association$TCHOL + 
            data_for_association$HDL + data_for_association$LDL + 
            data_for_association$BMI + data_for_association$Gender,
          data = data_for_association
        )
      )))
  
  species_results_flip <- param_with_ageG_correction_flip(org_list_to_use, 
                                                          species_with_age_flip, 
                                                          age_group)
  species_results_flip <- species_results_flip %>% 
    mutate(org_list = str_remove(org_list, "s__")) %>% 
    filter(padj < 0.05) %>% 
    mutate(age_group)
  return(species_results_flip)
}

run_glm_flip_ageG_only <- function(data_for_association, age_group) {
  org_list_to_use <- colnames(data_for_association)[grepl('__', colnames(data_for_association))]
  species_with_age_flip <-
    lapply(data_for_association[, org_list_to_use], function(x)
      coefficients(summary(
        glm(
          x ~ data_for_association$AgeG + data_for_association$Gender,
          data = data_for_association
        )
      )))
  
  species_results_flip <- param_with_ageG_correction_flip(org_list_to_use, 
                                                          species_with_age_flip, 
                                                          age_group)
  species_results_flip <- species_results_flip %>% 
    mutate(org_list = str_remove(org_list, "s__")) %>% 
    filter(padj < 0.05) %>% 
    mutate(age_group)
  return(species_results_flip)
}

# https://stackoverflow.com/questions/24237399/how-to-select-the-rows-with-maximum-values-in-each-group-with-dplyr
get_oldest_beta <- function(ageG_glm_list) {
  return(ageG_glm_list %>% 
    dplyr::select(org_list, betaestimate, age_group) %>%
    group_by(org_list) %>% 
    mutate(age_group = factor(age_group,
                              levels = c('91-100',
                                         '81-90', 
                                         '71-80',
                                         '61-70', 
                                         '41-60'))) %>% 
    mutate(rank = rank(age_group, ties.method='random')) %>% 
    filter(rank == 1) %>% dplyr::select(-rank, -age_group))
}

# Build ps object
build_ps <- function(data_for_association) {
  org_list_to_use <- colnames(data_for_association)[grepl('s__', colnames(data_for_association))]
  tax_profiles_t <- data_for_association %>% dplyr::select(LibraryID, org_list_to_use) %>% 
    column_to_rownames('LibraryID') %>% 
    t() %>% as.data.frame()

  metadata <- data_for_association %>% dplyr::select(-org_list_to_use) %>% 
    column_to_rownames('LibraryID')
  
  tax_ps <- metaphlanToPhyloseq(tax_profiles_t, metadata, simplenames = T, roundtointeger = F)
  
}

# CSS Normalization
# CSSnorm = function(data_for_association) {
#   org_list_to_use <- colnames(data_for_association)[grepl('__', colnames(data_for_association))]
#   tax_profiles <- data_for_association %>% dplyr::select(LibraryID, org_list_to_use) %>% 
#     column_to_rownames('LibraryID')
# Convert to Matrix from Data Frame
#   features_norm = as.matrix(t(tax_profiles))
CSSnorm <- function(tax_profiles) {  
  features_norm <- tax_profiles %>% 
    tibble::column_to_rownames('org_list')
  org_list_to_use <- rownames(features_norm)
  dd <- colnames(features_norm)
  
  # CSS Normalizing the Data
  # Create the metagenomeSeq object
  MGS = metagenomeSeq::newMRexperiment(
    t(features_norm),
    featureData = NULL,
    libSize = NULL,
    normFactors = NULL
  )
  # Trigger metagenomeSeq to calculate its Cumulative Sum scaling factor.
  MGS = metagenomeSeq::cumNorm(MGS, p = metagenomeSeq::cumNormStat(MGS))
  # Save the normalized data as data.frame
  features_CSS = as.data.frame(t(
    metagenomeSeq::MRcounts(MGS, norm = TRUE, log = FALSE)))
  
  # Rename the True Positive Features - Same Format as Before
  colnames(features_CSS) <- dd
  features_CSS <- cbind(org_list_to_use, features_CSS)
  
  rownames(features_CSS) <- NULL
  names(features_CSS)[1] <- 'org_list'
  return(features_CSS)
}

# ALDEX2
aldex2.fun <- function(otu.tab, meta, formula, alpha) {
  tt <- FALSE
  if(length(all.vars(as.formula(paste('~', formula)))) == 1) {
    if(!is.numeric(meta[, formula])) {
      tt <- TRUE
    }
  }
  if(tt) {
    res <- aldex(otu.tab, meta[, formula], test = "t", effect = FALSE, denom = "all")
    pval <- res$wi.ep
  } else {
    design <- model.matrix(as.formula(paste('~', formula)), meta)
    x <- aldex.clr(otu.tab, design, denom = "all")
    res <- aldex.glm(x, design)
    pval <- res[, 8]
  }
  qval <- p.adjust(pval, method = 'BH')
  rej <- which(qval <= alpha)
  if(length(pval) != nrow(otu.tab)) {
    rej <- as.numeric(gsub('taxon', '', rownames(res)))[rej]
  }
  return(rej)
}

# ANCOMII
get_ancom2_da <- function(tax_ps, ancom2_out) {
  ancom2_out_res = ancom2_out$res
  
  # Visualization for W stats
  q_val = ancom2_out$q_data
  beta_val = ancom2_out$beta_data
  # Only consider the effect sizes with the corresponding q-value less than alpha
  beta_val = beta_val * (q_val < 0.05) 
  # Choose the maximum of beta's as the effect size
  beta_pos = apply(abs(beta_val), 2, which.max) 
  beta_max = vapply(seq_along(beta_pos), function(i) beta_val[beta_pos[i], i],
                    FUN.VALUE = double(1))
  # Number of taxa except structural zeros
  n_taxa = ifelse(is.null(ancom2_out$zero_ind), 
                  nrow(otu_table(tax_ps)), 
                  sum(apply(ancom2_out$zero_ind, 1, sum) == 0))
  # Cutoff values for declaring differentially abundant taxa
  cut_off = 0.7 * (n_taxa - 1)
  
  ancom2_da = ancom2_out_res %>%
    dplyr::mutate(beta = beta_max,
                  direct = case_when(
                    detected_0.7 == TRUE & beta > 0 ~ "Positive",
                    detected_0.7 == TRUE & beta <= 0 ~ "Negative",
                    TRUE ~ "Not Significant"
                  )) %>%
    dplyr::arrange(W)
  ancom2_da$taxon_id = factor(ancom2_da$taxon_id, levels = ancom2_da$taxon_id)
  ancom2_da$W = replace(ancom2_da$W, is.infinite(ancom2_da$W), n_taxa - 1)
  ancom2_da$direct = factor(ancom2_da$direct, 
                           levels = c("Negative", "Positive", "Not Significant"))
  
  
  return(list(ancom2_da = ancom2_da,
         cut_off = cut_off))
  
}

plot_ancom2_da <- function(ancom2_da, cut_off) {
  p_w = ancom2_da %>%
    ggplot(aes(x = taxon_id, y = W, color = direct)) +
    geom_point(size = 2, alpha = 0.6) +
    labs(x = "Taxon", y = "W") +
    scale_color_discrete(name = NULL) + 
    geom_hline(yintercept = cut_off, linetype = "dotted", 
               color = "blue", size = 1.5) +
    geom_text(aes(x = 2, y = cut_off + 0.5, label = "W[0.7]"), 
              size = 5, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank())
  
  print(p_w)
}

# ANCOM-BC
get_ancombc_da <- function(ancombc_out) {
  ancombc_out_res <- ancombc_out$res
  ancombc_out_lfc <- ancombc_out_res$lfc
  
  # Visualization for age
  # adding this line for new ancombc
  ancombc_out_res$diff_abn <- ancombc_out_res$diff_abn %>%
    column_to_rownames('taxon') %>%
    mutate_all(~ as.numeric(.))
  ancombc_out_res$lfc <- ancombc_out_res$lfc %>% column_to_rownames('taxon')
  
  df_lfc = data.frame(ancombc_out_res$lfc * ancombc_out_res$diff_abn, 
                      check.names = FALSE) %>% 
    rownames_to_column("taxon_id")
  ancombc_out_res$se <- ancombc_out_res$se %>% column_to_rownames('taxon')
  df_se = data.frame(ancombc_out_res$se * ancombc_out_res$diff_abn, check.names = FALSE) %>% 
    rownames_to_column("taxon_id")
  colnames(df_se)[-1] = paste0(colnames(df_se)[-1], "SE")
  
  ancombc_da = df_lfc %>% 
    dplyr::left_join(df_se, by = "taxon_id") %>%
    dplyr::transmute(taxon_id, Age, AgeSE) %>%
    dplyr::filter(Age != 0) %>% 
    dplyr::arrange(desc(Age)) %>%
    dplyr::mutate(direct = ifelse(Age > 0, "Positive LFC", "Negative LFC"))
  ancombc_da$taxon_id = factor(ancombc_da$taxon_id, levels = ancombc_da$taxon_id)
  ancombc_da$direct = factor(ancombc_da$direct, 
                             levels = c("Positive LFC", "Negative LFC"))
  ancombc_da <- mutate(ancombc_da, taxon_id = str_remove(taxon_id, 's__'))
  return(ancombc_da)
  
  # Use global test on Age Group
  # global test: determine taxa that are DA btw >= 2 groups across 3/more diff groups
  # ancombc_age_global <- ancombc(phyloseq = tax_ps, formula = "AgeG + FBGlu + Tgly + TCHOL + HDL + LDL + BMI + gender", 
  #                               p_adj_method = "BH", prv_cut = 0, lib_cut = 0, 
  #                               group = "AgeG", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
  #                               max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
  # ancombc_age_global_res <- ancombc_age_global$res
  # ancombc_age_global_lfc <- ancombc_age_global_res$lfc
}

# convert ancombc output to glm output
convert_ancombc_to_glm <- function(ancombc_ageG_out) {
  
  ancombc_out_res <- ancombc_ageG_out$res
  ancombc_out_pval <- ancombc_out_res$p_val %>% 
    dplyr::select(Age, taxon) %>% 
    dplyr::rename(pval = 'Age', org_list='taxon') 
  
  ancombc_out_lfc <- ancombc_out_res$lfc %>% 
    dplyr::select(Age, taxon) %>% 
    dplyr::rename(betaestimate = 'Age', org_list='taxon') 
  
  ancombc_out_padj <- ancombc_out_res$q_val %>% 
    dplyr::select(Age, taxon) %>% 
    dplyr::rename(padj = 'Age', org_list='taxon') 
  return(ancombc_out_pval %>% 
           left_join(ancombc_out_lfc, by = 'org_list') %>% 
           left_join(ancombc_out_padj, by = 'org_list') %>% 
           mutate(org_list = str_remove(org_list, 's__')) %>% 
           mutate(stderror = 1))
  
}

plot_ancombc_da <- function(ancombc_da) {
  p_lfc = ggplot(data = ancombc_da %>% 
                   mutate(org_list = fct_reorder(org_list, desc(Age))), 
                 aes(x = org_list, y = Age, fill = direct, color = direct)) + 
    geom_bar(stat = "identity", width = 0.7, 
             position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(ymin = Age - AgeSE, ymax = Age + AgeSE), width = 0.2,
                  position = position_dodge(0.05), color = "black") + 
    labs(x = NULL, y = "Log fold change", 
         title = "Waterfall Plot of Age") + 
    scale_fill_discrete(name = NULL) +
    scale_color_discrete(name = NULL) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1, size=11))
  print(p_lfc)
  
}

get_oldest_lfc_ancombc <- function(ancombc_ageG_global_out) {
  oldest_da <- ancombc_ageG_global_out$res$diff_abn %>% 
    tibble::rownames_to_column('org_list') %>% 
    dplyr::select(org_list, `AgeG41-60`:`AgeG91-100`) %>% 
    pivot_longer(-org_list) %>% 
    mutate(org_list = str_remove(org_list, 's__'),
           name = str_remove(name, 'AgeG')) %>% 
    filter(value == TRUE) %>%
    mutate(name = factor(name,
                         levels = c('91-100',
                                    '81-90', 
                                    '71-80',
                                    '61-70', 
                                    '41-60'))) %>% 
    group_by(org_list) %>% 
    mutate(rank = rank(name, ties.method='random')) %>% 
    filter(rank == 1)
  
  oldest_lfc <- ancombc_ageG_global_out$res$lfc %>% 
    tibble::rownames_to_column('org_list') %>% 
    dplyr::select(org_list, `AgeG41-60`:`AgeG91-100`) %>% 
    pivot_longer(-org_list) %>% 
    mutate(org_list = str_remove(org_list, 's__'),
           name = str_remove(name, 'AgeG')) %>% 
    rename(betaestimate = value)
  
  return(oldest_da %>% left_join(oldest_lfc, by = c('org_list', 'name')) %>% 
           dplyr::select(org_list, betaestimate))
}

get_sign_from_ancombc <- function(ancombc_da) {
  return(ancombc_da  %>% 
           mutate(sign = sign(Age)) %>% 
           dplyr::select(org_list, sign))
}

# run Maaslin2
run_maaslin2 <- function(data_for_association, normalization, transform, output_dir) {
  org_list_to_use <- colnames(data_for_association)[grepl('s__', colnames(data_for_association))]
  tax_profiles <- data_for_association %>% dplyr::select(LibraryID, org_list_to_use) %>% 
    column_to_rownames('LibraryID')
  
  metadata <- data_for_association %>% dplyr::select(-org_list_to_use) %>% 
    column_to_rownames('LibraryID')
  maaslin2_out <- Maaslin2(input_data = tax_profiles,
                           normalization = normalization,
                           transform = transform,
                           analysis_method = "LM",
                           min_prevalence = 0,
                           input_metadata = metadata,
                           output = output_dir,
                           fixed_effects = c('Age', 'Gender', 'FBGlu', 'Tgly',
                                             #LDL,
                                             'TCHOL', 'HDL', 'BMI'),
                           correction = 'BH',
                           reference = NULL, 
                           standardize = FALSE)
  maaslin2_res <- maaslin2_out$results %>% filter(value == 'Age')
  maaslin2_res$padj <- p.adjust(maaslin2_res$pval, 'fdr')
  maaslin2_res_sig <- maaslin2_res %>% #filter(padj < 0.05) %>% 
    mutate(feature = str_remove(feature, 's__')) %>% 
    dplyr::rename(org_list = feature)
  return(maaslin2_res_sig)
}

run_maaslin2_ageG <- function(data_for_association, normalization, transform, output_dir) {
  org_list_to_use <- colnames(data_for_association)[grepl('s__', colnames(data_for_association))]
  tax_profiles <- data_for_association %>% dplyr::select(LibraryID, org_list_to_use) %>% 
    column_to_rownames('LibraryID')
  
  metadata <- data_for_association %>% dplyr::select(-org_list_to_use) %>% 
    column_to_rownames('LibraryID')
  maaslin2_out <- Maaslin2(input_data = tax_profiles,
                           normalization = normalization,
                           transform = transform,
                           analysis_method = "LM",
                           min_prevalence = 0,
                           input_metadata = metadata,
                           output = output_dir,
                           fixed_effects = c('AgeG', 'Gender'),
                           correction = 'BH',
                           reference = NULL, 
                           standardize = FALSE)
  maaslin2_res <- maaslin2_out$results %>% filter(metadata == 'AgeG')
  maaslin2_res$padj <- p.adjust(maaslin2_res$pval, 'fdr')
  maaslin2_res_sig <- maaslin2_res %>% #filter(padj < 0.05) %>% 
    mutate(feature = str_remove(feature, 's__')) %>% 
    rename(org_list = feature)
  return(maaslin2_res_sig)
}

# run Linear mixed-effect model
run_lme <- function(data_for_association) {
  ranef_function <- lme4::ranef
  org_list_to_use <- colnames(data_for_association)[grepl('__', colnames(data_for_association))]
  species_with_age <-
    lapply(data_for_association[, org_list_to_use], function(x)
      coefficients(summary(
        lmerTest::lmer(
          data_for_association$Age ~ x +  data_for_association$FBGlu +  data_for_association$Tgly +
            data_for_association$TCHOL + data_for_association$HDL + data_for_association$LDL + 
            data_for_association$BMI + data_for_association$Gender + (1 | AgeG),
          data = data_for_association
        )
      )))
  
  pval <- NULL
  org_list <- NULL
  estimate <- NULL
  stderror <- NULL
  betaestimate <- NULL
  stderror <- NULL
  for (i in seq_along(1:length(species_with_age)))
  {
    if ("x" %in% rownames(species_with_age[[i]])){
      pval <-
        rbind(pval, species_with_age[[i]]["x",][5])
      org_list <- rbind(org_list, org_list_to_use[i])
      betaestimate <-
        rbind(betaestimate, species_with_age[[i]]["x",][1])
      stderror <-
        rbind(stderror, species_with_age[[i]]["x",][2])
    }
    else{
      print(org_list_to_use[i])}
  }
  
  all_df <- data.frame(org_list, pval, betaestimate, stderror)
  colnames(all_df) <- c('org_list', 'pval', 'betaestimate', 'stderror')
  all_df$padj <- p.adjust(all_df$pval, 'fdr')
  
  # filter using p-val < 0.05
  all_df <- all_df %>% filter(pval < 0.05)
  return(all_df)
}

# Venn diagram
## for 5 lists
plot_venn <- function(list_1, list_2, list_3, list_4, list_5, Title) {
  venn.plot <- venn.diagram(
    x = list(list_1$org_list, 
             list_2$org_list,
             list_3$org_list,
             list_4$org_list,
             list_5$org_list),
    category.names = c("Original", 
                       "Age-adjusted No CRE",
                       "Age-adjusted With CRE",
                       "AgeG-adjusted No CRE", 
                       "AgeG-adjusted With Cre"),
    fill = c("#8DD3C7", "#FFD700", "#80B1D3", "#FB8072", "#BEBADA"),
    margin=0.2,
    filename = NULL, disable.logging = TRUE,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    compression = "lzw",
    lwd = 1,
    col=c("#8DD3C7", "#FFD700", "#80B1D3", "#FB8072", "#BEBADA"),
    cex = 1,
    fontfamily = "sans",
    cat.cex = 1,
    cat.default.pos = "outer",
    cat.pos = c(-25, 0, 200, 200, 0),
    cat.dist = c(0.25, 0.3, 0.25, 0.25, 0.25),
    cat.fontfamily = "sans",
    cat.col = c("#8DD3C7", "#FFD700", "#80B1D3", "#FB8072", "#BEBADA")#,
    #rotation = 1
  )
  
  
  #grid.newpage()
  #grid.draw(venn.plot)
  require(gridExtra)
  grid.arrange(gTree(children=venn.plot), top=Title)
}

get_beta_of_signif <- function(glm_res, colname) {
  return(glm_res %>% filter(padj < 0.05) %>% 
           dplyr::select(org_list, betaestimate) %>% 
           dplyr::rename(`colname` = betaestimate))
}

# merge multiple results
merge_multiple_results <- function(results, quantity) {
  combined_res <- results$original_ageG_glm %>% dplyr::select(org_list, quantity) %>% rename(ori_ageG = quantity) %>% 
    full_join(results$original_age_glm %>% dplyr::select(org_list, quantity) %>% rename(ori_age = quantity)) %>% 
    full_join(results$ageG_with_cre_glm %>% dplyr::select(org_list, quantity) %>% rename(with_cre_ageG = quantity)) %>%
    full_join(results$age_with_cre_glm %>% dplyr::select(org_list, quantity) %>% rename(with_cre_age = quantity)) %>%
    full_join(results$ageG_no_cre_glm %>% dplyr::select(org_list, quantity) %>% rename(no_cre_ageG = quantity)) %>%
    full_join(results$age_no_cre_glm %>% dplyr::select(org_list, quantity) %>% rename(no_cre_age = quantity)) %>%
    full_join(results$ageG_nosg90_glm %>% dplyr::select(org_list, quantity) %>% rename(no_sg90_ageG = quantity)) %>%
    full_join(results$ageG_no_cre_nosg90_glm %>% dplyr::select(org_list, quantity) %>% rename(no_cre_no_sg90_ageG = quantity)) %>%
    full_join(results$age_nosg90_glm %>% dplyr::select(org_list, quantity) %>% rename(no_sg90_age = quantity)) %>%
    full_join(results$age_no_cre_nosg90_glm %>% dplyr::select(org_list, quantity) %>% rename(no_cre_no_sg90_age = quantity)) %>% 
    mutate(org_list = str_replace(org_list, '_', ' '))
  
}
# plot heatmap of beta color-coded according to significance
plot_signif <- function(combined_res) {
  combined_res %>% filter(org_list %in% desired_species) %>% 
    mutate(beta = round(beta, 2),
           name = factor(name, levels = c('ori_ageG', 'with_cre_ageG', 'no_cre_ageG',
                                          'ori_age', 'with_cre_age', 'no_cre_age',
                                          'no_sg90_ageG', 'no_cre_no_sg90_ageG', 
                                          'no_sg90_age', 'no_cre_no_sg90_age'))) %>% 
    ggplot(aes(x = name, y = org_list, color = signif, label = beta)) + 
    geom_text(size=10) + labs(x = 'Data', y = 'Species') +
    theme_bw() +
    theme(axis.text = element_text(size=14), axis.text.x = element_text(angle = 10),
          axis.title = element_text(size=18))
}

library(Biobase)

metaphlanToPhyloseq <- function(
    metaphlandir,
    metadat=NULL,
    simplify=TRUE){
  ## tax is a matrix or data.frame with the table of taxonomic abundances, rows are taxa, columns are samples
  ## metadat is an optional data.frame of specimen metadata, rows are samples, columns are variables
  ## if simplify=TRUE, use only the most detailed level of taxa names in the final object
  ## metaphlanToPhyloseq("~/Downloads/metaphlan_bugs_list")
  .getMetaphlanTree <- function(removeGCF=TRUE, simplify=TRUE){
    if (!requireNamespace("ape")) {
      stop("Please install the ape package to read Newick trees")
    }
    nwkfile <- bzfile(system.file("extdata/metaphlan2_selected.tree.reroot.nwk.bz2",
                                  package="curatedMetagenomicData"))
    tree <- ape::read.tree(nwkfile)
    close(nwkfile)
    if(removeGCF)
      tree$tip.label <- sub("\\|GCF_[0-9]+$", "", tree$tip.label)
    if(simplify)
      tree$tip.label <- gsub(".+\\|", "", tree$tip.label)
    return(tree)
  }
  .joinListOfMatrices <- function(obj) {
    rnames <- Reduce(union, lapply(obj, rownames))
    cnames <- names(obj)
    if (!all(isUnique(cnames))) {
      stop("Column names are not unique.")
    }
    output <- matrix(0,
                     nrow = length(rnames),
                     ncol = length(cnames),
                     dimnames = list(rnames, cnames)
    )
    for (i in seq_along(obj)) {
      output[match(rownames(obj[[i]]), rownames(output)), i] <- obj[[i]][, 1]
    }
    return(output)
  }
  fnames <- list.files(metaphlandir)
  bug.l <- lapply(fnames, function(x){
    res <- read.delim(file.path(metaphlandir, x), stringsAsFactors = FALSE, row.names = 1)
    colnames(res) <- x
    return(res)
  })
  names(bug.l) <- fnames
  tax <- .joinListOfMatrices(bug.l)
  xnames = rownames(tax)
  shortnames = gsub(".+\\|", "", xnames)
  if(simplify){
    rownames(tax) = shortnames
  }
  x2 = strsplit(xnames, split="|", fixed=TRUE)
  taxmat = matrix(NA, ncol=max(sapply(x2, length)), nrow=length(x2))
  colnames(taxmat) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:ncol(taxmat)]
  rownames(taxmat) = rownames(tax)
  for (i in 1:nrow(taxmat)){
    taxmat[i, 1:length(x2[[i]])] <- x2[[i]]
  }
  taxmat = gsub("[a-z]__", "", taxmat)
  taxmat = phyloseq::tax_table(taxmat)
  otutab = phyloseq::otu_table(tax, taxa_are_rows=TRUE)
  tree <- .getMetaphlanTree(simplify=simplify)
  if(is.null(metadat)){
    metadat <- data.frame(file=fnames, row.names=fnames, stringsAsFactors = FALSE)
  }
  res = phyloseq::phyloseq(taxmat, otutab, phyloseq::sample_data(metadat), tree)
  return(res)
}

# from https://github.com/waldronlab/presentations/blob/b6c6cbc812b64c1bb371dd9ec261300562089cf2/Waldron_2016-06-07_EPIC/metaphlanToPhyloseq.R
metaphlanToPhyloseq <- function(
    tax,
    metadat=NULL,
    simplenames=TRUE,
    roundtointeger=FALSE,
    split="|"){
  ## tax is a matrix or data.frame with the table of taxonomic abundances, rows are taxa, columns are samples
  ## metadat is an optional data.frame of specimen metadata, rows are samples, columns are variables
  ## if simplenames=TRUE, use only the most detailed level of taxa names in the final object
  ## if roundtointeger=TRUE, values will be rounded to the nearest integer
  xnames = rownames(tax)
  shortnames = gsub(paste0(".+\\", split), "", xnames)
  if(simplenames){
    rownames(tax) = shortnames
  }
  if(roundtointeger){
    tax = round(tax * 1e4)
  }
  x2 = strsplit(xnames, split=split, fixed=TRUE)
  taxmat = matrix(NA, ncol=max(sapply(x2, length)), nrow=length(x2))
  colnames(taxmat) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:ncol(taxmat)]
  rownames(taxmat) = rownames(tax)
  for (i in 1:nrow(taxmat)){
    taxmat[i, 1:length(x2[[i]])] <- x2[[i]]
  }
  taxmat = gsub("[a-z]__", "", taxmat)
  taxmat = phyloseq::tax_table(taxmat)
  otutab = phyloseq::otu_table(tax, taxa_are_rows=TRUE)
  if(is.null(metadat)){
    res = phyloseq::phyloseq(taxmat, otutab)
  }else{
    res = phyloseq::phyloseq(taxmat, otutab, phyloseq::sample_data(metadat))
  }
  return(res)
}

krakenToPhyloseq <- function(
    tax,
    metadat=NULL,
    simplenames=TRUE,
    roundtointeger=FALSE,
    split="|"){
  ## tax is a matrix or data.frame with the table of taxonomic abundances, rows are taxa, columns are samples
  ## metadat is an optional data.frame of specimen metadata, rows are samples, columns are variables
  ## if simplenames=TRUE, use only the most detailed level of taxa names in the final object
  ## if roundtointeger=TRUE, values will be rounded to the nearest integer
  xnames = rownames(tax)
  shortnames = gsub(paste0(".+\\", split), "", xnames)
  if(simplenames){
    rownames(tax) = shortnames
  }
  if(roundtointeger){
    tax = round(tax * 1e4)
  }
  x2 = strsplit(xnames, split=split, fixed=TRUE)
  taxmat = matrix(NA, ncol=max(sapply(x2, length)), nrow=length(x2))
  colnames(taxmat) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:ncol(taxmat)]
  rownames(taxmat) = rownames(tax)
  for (i in 1:nrow(taxmat)){
    taxmat[i, 1:length(x2[[i]])] <- x2[[i]]
  }
  taxmat = gsub("[a-z]__", "", taxmat)
  taxmat = phyloseq::tax_table(taxmat)
  otutab = phyloseq::otu_table(tax, taxa_are_rows=TRUE)
  if(is.null(metadat)){
    res = phyloseq::phyloseq(taxmat, otutab)
  }else{
    res = phyloseq::phyloseq(taxmat, otutab, phyloseq::sample_data(metadat))
  }
  return(res)
}

# get median p-adj values for combined results
get_median_padj <- function(all_cohorts_res, sg90_cpe_res, prevalence_res) {
  combined_padj <- rbind(all_cohorts_res %>% dplyr::select(org_list, padj) %>% 
                           mutate(Group = 'All Cohorts\nEM'),
                         sg90_cpe_res %>% dplyr::select(org_list, padj) %>% 
                           mutate(Group = 'SG90 + CPE\nEM'),
                         prevalence_res %>% dplyr::select(org_list, padj) %>% 
                           mutate(Group = 'All Cohorts\nPrevalence'))
  # calculate median
  combined_padj <- combined_padj %>% 
    group_by(org_list) %>% 
    summarize(n = n(),
              padj_median = if(n > 2) median(padj) else max(padj))
  return(combined_padj)
}
