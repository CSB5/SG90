## Helper functions
create_ps <- function(tax.profile, metadata) {
  taxmat <- matrix(sample(letters, 7*nrow(tax.profile), replace = T), 
                   nrow = nrow(tax.profile), ncol = 7)
  rownames(taxmat) <- rownames(tax.profile)
  colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  otu.s <- otu_table(tax.profile, taxa_are_rows = T)
  tax.s <- tax_table(taxmat)
  gut.s.fil <- phyloseq(otu.s, tax.s, metadata)
  metadata.ps <- sample_data(metadata)
  # add metadata
  gut.s.fil <- merge_phyloseq(gut.s.fil, metadata.ps)
  return(gut.s.fil)
}

# run envfit to get direction of covariates
get_nmds <- function(tax_profiles, metadata, distance_metric) {
  ps_ts <- create_ps(tax_profiles, metadata %>% 
                       mutate(Cohort = case_when(Cohort == 'iOMICS' ~ 'SPMP',
                                                 TRUE ~ Cohort)) %>% 
                       tibble::column_to_rownames('LibraryID'))
  nmds_bray <- ordinate(ps_ts, method = "NMDS", distance = distance_metric,
                        try_max = 50)
  return(list(ps_ts = ps_ts, 
              nmds_bray = nmds_bray))
}

# get scores for continuous covariate
get_envfit_score <- function(ts_envfit) {
  ts_envfit_scores <- scores(ts_envfit, "vectors")[ts_envfit$vectors$pvals<0.05, ]
  ts_envfit_coord <- as.data.frame(ts_envfit_scores) * ordiArrowMul(ts_envfit) * 1.25 # only plot arrows that are significant
  ts_envfit_length <- sort(sqrt((ts_envfit_scores[,1])^2 + (ts_envfit_scores[,2])^2))
  return(list(scores = ts_envfit_scores,
              coord = ts_envfit_coord,
              length = ts_envfit_length))
}

get_envfit_factor <- function(ts_envfit) {
  # for cohort
  ts_envfit_factors <- scores(ts_envfit, "factors")[ts_envfit$factors$pvals<5, ]
  ts_envfit_factors_coord <- as.data.frame(ts_envfit_factors) * ordiArrowMul(ts_envfit) * 1.25 # only plot arrows that are significant
  ts_envfit_factors_length <- sort(sqrt((ts_envfit_factors[, 1])^2 + (ts_envfit_factors[, 2])^2))
  return(list(factors = ts_envfit_factors,
              coord = ts_envfit_factors_coord,
              length = ts_envfit_factors_length))
  
}

## add ellipse from betadisper
get_ellipse_data <- function(ord, metadata, scaling = 1, choices = c(1,2),
                             kind = c("sd", "se", "ehull"), conf=NULL, show.groups="all", 
                             ellipse = TRUE, label = FALSE, pt.size = 3) {
  # get groups
  rownames(metadata) <- NULL
  metadata_order <- metadata %>% 
    mutate(Cohort = case_when(Cohort == 'iOMICS' ~ 'SPMP',
                              TRUE ~ Cohort)) %>% 
    tibble::column_to_rownames('LibraryID')
  groups <- as.factor(metadata_order$Cohort)
  if (show.groups[1]=="all") {
    show.groups <- as.vector(levels(groups))
  }
  
  # Get site coordinates to plot.
  df_ord <- vegan::scores(ord, display = "sites", scaling=scaling, choices=choices)
  axis.labels <- ord_labels(ord)[choices]
  df_ord <- data.frame(x=df_ord[ , 1], y=df_ord[ , 2], Group=groups)
  # Get ellipse centers to annotate.
  df_mean.ord <- aggregate(df_ord[,1:2], by=list(df_ord$Group),mean)
  colnames(df_mean.ord) <- c("Group", "x", "y")
  df_mean.ord <- df_mean.ord[df_mean.ord$Group %in% show.groups, ]
  
  # Get parameters from the ordiellipse function.
  if (is.null(conf)) {
    rslt <- vegan::ordiellipse(ord, groups=groups, display = "sites", scaling=scaling, choices=choices, kind = kind, show.groups = show.groups, draw = "none", label = label)
  } else {
    rslt <- vegan::ordiellipse(ord, groups=groups, display = "sites", scaling=scaling, choices=choices, kind = kind, show.groups = show.groups, draw = "none", conf = conf, label = label)
  }
  
  # Get points to plot for the ellipses.
  df_ellipse <- data.frame()
  for(g in show.groups) {
    df_ellipse <- rbind(df_ellipse, cbind(as.data.frame(with(df_ord[df_ord$Group==g,],
                                                             vegan:::veganCovEllipse(rslt[[g]]$cov,rslt[[g]]$center, rslt[[g]]$scale))),Group=g))
  }
  colnames(df_ellipse) <- c("x", "y", "Group")
  df_ellipse <- df_ellipse[ , c(3,1,2)]
  return(list('df_ellipse' = df_ellipse,
              'df_center' = df_mean.ord))
}

# plot ordination with envfit arrows
plot_ord <- function(ps_obj, nmds_bray_mat, envfit_scores_coord, envfit_factors_coord,
                     ellipse) {
  # plot_ordination(physeq = ps_obj, ordination = nmds_bray_mat,
  #                 color = 'Cohort') +
 # browser()
  rownames(envfit_factors_coord) %<>% str_replace_all(c("Cohort"=" "))
  
  coords <- scores(nmds_bray_mat, display = 'sites', choices = 1:2, physeq = ps_obj) %>% 
    as.data.frame() %>% tibble::rownames_to_column('LibraryID') %>% 
    left_join(data.frame(sample_data(ps_obj)) %>% tibble::rownames_to_column('LibraryID'),
              by = 'LibraryID') 
  ggplot(coords, aes(x = NMDS1, y = NMDS2, color = Cohort)) + 
    geom_point(alpha = 0.8, size = 2.5) +
    #stat_ellipse(aes(color = Cohort), linetype = 'dashed', lwd = 1.2) +
    # geom_path(data = ellipse,
    #           aes(x = x, y = y, color = Group),
    #           size = 1.5, linetype = 2) +
    geom_segment(data = envfit_scores_coord, mapping = aes(x = 0, y = 0, 
                                                           xend = NMDS1, 
                                                           yend = NMDS2), 
                 size = 0.75, colour = "black", arrow = arrow(length = unit(0.05, "npc"))) +
    geom_text(data = envfit_scores_coord, mapping = aes(x = NMDS1,
                                                        y = NMDS2, 
                                                        label = rownames(envfit_scores_coord)), 
              colour = "black", fontface = "bold", size = 4) +
    geom_segment(data = envfit_factors_coord, mapping = aes(x = 0, y = 0, 
                                                            xend = NMDS1, 
                                                            yend = NMDS2), 
                 size = 0.75, colour = "black", arrow = arrow(length = unit(0.05, "npc"))) +
    geom_text(data = envfit_factors_coord, mapping = aes(x = NMDS1,
                                                         y = NMDS2), 
              label = row.names(envfit_factors_coord),
              colour = "black", fontface = "bold", size = 4) + 
    scale_colour_manual(values = c(
      '#e6194B', #red
      '#3cb44b', #green
      '#354a99', # blue
      '#ffe119', #yellow
      # original colors from aarthi
      "#a50023",
      "#f47d4a",
      "#fed98a",
      "#add39f",
      '#6da5cb',
      "#354a99"
    )) +
    theme_bw() + 
    theme(legend.background = element_blank(),
          legend.position = "right",
          legend.key = element_blank(),
          legend.text = element_text(size = 11),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title       = element_text(face = "bold", size = 14),
          plot.title = element_text(face = 'bold', size = 14), 
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          text = element_text(family = 'Arial'))
  
}

# PERMANOVA R2
#Checking if batch correction worked fine
get_adonis <- function(tax_profiles, metadata) {
  D <- vegdist(t(tax_profiles), method = 'bray')
  
  set.seed(1)
  fit_adonis <- adonis2(D ~ Age + Cohort + Gender + Reads + Host,
                        data = metadata)
  
  print(fit_adonis)
  
}


