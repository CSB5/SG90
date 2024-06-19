# Script to plot for figure 1
# date: 20230314

library(vegan)
library(ggrepel)
library(MASS)
library(reshape2)
library(RColorBrewer)
library(officer)
library(rvg)
library(here)
library(tidyverse)

# read cre profiles - replace x with org_list & '.' in CON.XXX
read_cre <- function(filepath) {
  profiles <- read.csv(filepath) %>% 
    rename(org_list = 'X')
  names(profiles) <- gsub("\\.", "", names(profiles))
  return(profiles)
}

# remove NA & low abundance; then renormalize & reorder colnames
preprocess_tax_profiles <- function(taxonomic_profiles) {
  
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

#Functions to plot PCoA

performpcoa <- function(tax_profiles, metadata_order) {
  distmetric <- vegdist(t(tax_profiles), method = "bray")
  cmds <- cmdscale(distmetric, k = 3, eig = TRUE)
  eigen <- cmds$eig / sum(cmds$eig) * 100
  cmdvals <- as.data.frame(cmds$points)
  row.names(metadata_order) <- metadata_order$LibraryID
  dat_merged_all <- data.frame(cmdvals, metadata_order)
  listofop <- list("dat_merged_all" = dat_merged_all, "eigen" = eigen, "distmetric"=distmetric)
  return(listofop)
  
}
plotpcoa <-
  function(dat_merged_all, eigen, title,legendneed, pointnum) {
    ggplot(dat_merged_all,
           aes(x = V1,
               y = V2)) +
      geom_point(
        data = subset(dat_merged_all, Twogroups == 'Elderly'),
        aes(color = AgeGNum),
        size = pointnum,
      ) + geom_point(
        data = subset(dat_merged_all,
                      Twogroups == 'Young'),
        aes(color = AgeGNum),
        size = pointnum,
      ) +
      scale_colour_manual(values = c(
        "#a50023",
        "#f47d4a",
        "#fed98a",
        "#add39f",
        '#6da5cb',
        "#354a99"
      )) +
      labs(
        x = paste0('PCoA1 (', round(eigen[1], 1), '%)'),
        y = paste0('PCoA2 (', round(eigen[2], 1), '%)'),
        colour = 'Age groups'
      ) +
      theme_bw() +
      theme(
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = legendneed,
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title       = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)
      ) + ggtitle('') +
      geom_abline(slope=1, intercept = 0, linetype = 'dashed') +
      geom_abline(slope=-1, intercept = 0, linetype = 'dashed') 
  }

violin_plot_pcoa <- function(dat_merged_decade_adj, pcoa_axis, pcoa_axis_label){
  ggplot(
    dat_merged_decade_adj$dat_merged_all,
    aes(y = dat_merged_decade_adj$dat_merged_all[, pcoa_axis], x = AgeG, fill = AgeG)
  ) + geom_violin() + geom_boxplot(width = 0.2, outlier.shape = NA) +
    scale_fill_manual(values = c(
      "#a50023",
      "#f47d4a",
      "#fed98a",
      "#add39f",
      '#6da5cb',
      "#354a99"
    )) + theme_bw() + theme(
      plot.title = element_text(
        color = "black",
        size = 14,
        hjust = 0.5
      ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(size = 12,face = 'bold'),
      axis.title.y = element_text(size = 12, face = 'bold'),
      axis.text.x = element_text(
        size = 12,
        angle = 40,
        hjust = 1
      ),
      axis.text.y = element_text(size = 12), legend.position = 'none'
    ) + xlab('Age group') + ylab(pcoa_axis_label)}

## Alpha diversity
# remove outliers for alpha diversity per age group
remove_outliers <- function(alpha, ageG) {
  alpha_ageG <- alpha %>% filter(AgeG == ageG)
  print(paste('##########', ageG, '###########', sep=' '))
  print(paste("Before removing outliers: ", dim(alpha_ageG)[[1]]))
  q1 <- quantile(alpha_ageG[, 2], .25)
  q3 <- quantile(alpha_ageG[, 2], .75)
  iqr <- IQR(alpha_ageG[, 2])
  no_outliers <- subset(alpha_ageG, alpha_ageG[, 2] > (q1 - 1.5*iqr) &
                          alpha_ageG[, 2] < (q3 + 1.5*iqr))
  print(paste("After removing outliers: ", dim(no_outliers)[[1]]))
  return(no_outliers)
}

# combine after removing for outliers
combine_no_outliers <- function(alpha_to_plot) {
  return(rbind(remove_outliers(alpha_to_plot, '21-40'),
               remove_outliers(alpha_to_plot, '41-60'),
               remove_outliers(alpha_to_plot, '61-70'),
               remove_outliers(alpha_to_plot, '71-80'),
               remove_outliers(alpha_to_plot, '81-90'),
               remove_outliers(alpha_to_plot, '91-100'))
  )
}

# Functions to plot alpha diversity violin
# plot_alpha <- function(alpha_to_plot, alpha_measure) {
#   ggplot(alpha_to_plot, aes_string(
#     x = 'AgeG',
#     y = alpha_measure,
#     outline = FALSE,
#     fill = 'AgeG'
#   )) +
#     geom_violin() + geom_boxplot(width = 0.2, outlier.shape = NA) +
#     scale_fill_manual(values = c(
#       "#a50023",
#       "#f47d4a",
#       "#fed98a",
#       "#add39f",
#       '#6da5cb',
#       "#354a99"
#     )) + theme_bw() + theme(
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       plot.title = element_text(
#         color = "black",
#         face = "bold",
#         size = 14,
#         hjust = 0.5
#       ),
#       axis.title.x = element_text(size = 14),
#       axis.title.y = element_text(size = 14, face = "bold"),
#       axis.text.x = element_text(
#         size = 14,
#         angle = 40,
#         hjust = 1
#       ),
#       axis.text.y = element_text(size = 14), legend.position = 'none',
#       text = element_text(family = 'Arial')
#     ) + xlab('') + ylab(paste(alpha_measure, ' diversity', sep=''))
#   
# }

# calculate ordinal logistic test
calc_olr <- function(alpha_to_plot) {
  richness_olr <- polr(as.factor(AgeG) ~ Richness, Hess = TRUE, data = richness_to_plot)
  ctable_richness <- coef(summary(richness_olr))
  p_richness <- pnorm(abs(ctable_richness[, "t value"]), lower.tail = FALSE) * 2
  ctable_richness <- cbind(ctable_richness, "p value" = p_richness)
  return(p_richness = p_richness, ctable_richness = ctable_richness)
}


## Function to plot boxplots relative abundance
# filter profiles & merge w/ metadata (not normalized to 100%)
transform_tax_profiles <- function(tax_profiles, metadata, n_nonzero) {
  tax_profiles_t <- tax_profiles %>% column_to_rownames('org_list') %>% 
    t() %>% 
    as.data.frame() 
  
  # Check number of non zero entries
  nonzero_in_samples <- data.frame(nnzorg = colSums(tax_profiles_t !=0))
  ggplot(nonzero_in_samples, aes(x=nnzorg)) + 
    geom_histogram(color="darkblue", fill="lightblue") + 
    scale_x_continuous(breaks = scales::pretty_breaks(n=7)) + 
    xlab('Number of samples') + ylab('Count of organisms with non-zero abundances')
  
  orgsinmorethanfiftysamples <- nonzero_in_samples %>% filter(nnzorg >= n_nonzero)
  
  taxa_with_morethan50 <- tax_profiles_t %>% 
    dplyr::select(rownames(orgsinmorethanfiftysamples)) %>% 
    rownames_to_column('LibraryID')
  
  data_for_association <- merge(taxa_with_morethan50, metadata, by = 'LibraryID')
  return(data_for_association)
}

# plotbox <- function(data_for_plot, xaxislabel, yaxislabel, title, foldername, fname, ylimitplot, species_results){ 
#   #ylimitplot <- boxplot.stats(data_for_plot[[yaxislabel]])$stats[c(1, 5)]
#   bac_plot <-
#     ggplot(data_for_plot,
#            aes_string(
#              x = xaxislabel,
#              y = (yaxislabel),
#              outline =  FALSE,
#              fill = xaxislabel
#            )) +
#     geom_boxplot(width = 1.0, outlier.shape=NA) +
#     ggtitle(title) + theme_bw() +
#     theme(
#       legend.position = 'none', 
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       axis.line = element_line(colour = "black"),
#       plot.title = element_text(
#         color = "black",
#         size = 14,
#         face = "bold.italic",
#         hjust = 0.5
#       ),
#       axis.title.x = element_text(size = 14),
#       axis.title.y = element_text(size = 14, face = 'bold'),
#       axis.text.x = element_text(
#         size = 14,
#         angle = 40,
#         hjust = 1
#       ),
#       axis.text.y = element_text(size = 14),
#       text = element_text(family = 'Arial')
#     ) + labs(x = "",
#              y = "Relative abundance",
#              face = "bold",
#              size = 14) + scale_fill_manual(values = c(
#                "#a50023",
#                "#f47d4a",
#                "#fed98a",
#                "#add39f",
#                '#6da5cb',
#                "#354a99"
#              )) + coord_cartesian(ylim = c(0, ylimitplot))
#   ylimit_bac <- layer_scales(bac_plot)$y$range$range
#   xlimit_bac <- layer_scales(bac_plot)$x$range$range
#   p_bac <- subset(species_results, org_list == yaxislabel)[, 'padj']
#   print(p_bac)
#   
#   if(p_bac <  0.001){
#     bac_plot_1<- bac_plot + annotate("text", x = 0.5+length(xlimit_bac)/2, y = ylimitplot+ 0.05, ymax =ylimitplot+ 1, label = "***", size=8)
#   } else if(p_bac <= 0.01){
#     bac_plot_1<- bac_plot + annotate("text", x = 0.5+length(xlimit_bac)/2, y = ylimitplot+ 0.05, ymax =ylimitplot+ 1, label = "**", size=6)
#   }
#   else if(p_bac <= 0.05){
#     bac_plot_1<- bac_plot + annotate("text", x = 0.5+length(xlimit_bac)/2, y = ylimitplot+ 0.05, label = "*", size=6)
#   } else {
#     bac_plot_1 <- bac_plot + annotate("text", x = 0.5+length(xlimit_bac)/2, y = ylimitplot+ 0.05, label = "paste(italic(n.s.))", size=6, parse=TRUE)
#   }
#   
#   #+  coord_cartesian(ylim = c(0, 8)),, outlier.shape = FALSE
#   bacplotfname <- paste0(fname,"_speciesthr", species_thr, "_", "genusthr", genus_thr,"_", hist_filter, ".png")
#   ggsave(here(foldername,bacplotfname), plot=bac_plot_1, width=10, height=10, units=c("cm"),limitsize = FALSE)
#   return(bac_plot_1)
# }
# 
plotbox_nonzero <- function(data_for_plot, yaxislabel, title, foldername, fname, ylimitplot, species_results){ 
  #ylimitplot <- boxplot.stats(data_for_plot[[yaxislabel]])$stats[c(1, 5)]
  # get nonzero values for the species of interest
  bac_for_plot <- data_for_plot %>% dplyr::select(LibraryID:AgeG) %>% 
    dplyr::select(-c('Patient.ID', 'Group')) %>% 
    pivot_longer(-c('LibraryID', 'AgeG')) %>% 
    filter(name == yaxislabel, value > 0)
  
  bac_plot <- bac_for_plot %>% 
    ggplot(aes(x = AgeG, y = value, outline =  FALSE, fill = AgeG)) +
    geom_boxplot(width = 1.0, outlier.shape=NA) +
    ggtitle(title) + theme_bw() +
    theme(
      legend.position = 'none', 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      plot.title = element_text(
        color = "black",
        size = 14,
        face = "bold.italic",
        hjust = 0.5
      ),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14, face = 'bold'),
      axis.text.x = element_text(
        size = 14,
        angle = 40,
        hjust = 1
      ),
      axis.text.y = element_text(size = 14),
      text = element_text(family = 'Arial')
    ) + labs(x = "",
             y = "Relative abundance",
             face = "bold",
             size = 14) + scale_fill_manual(values = c(
               "#a50023",
               "#f47d4a",
               "#fed98a",
               "#add39f",
               '#6da5cb',
               "#354a99"
             )) + coord_cartesian(ylim = c(0, ylimitplot))
  ylimit_bac <- layer_scales(bac_plot)$y$range$range
  xlimit_bac <- layer_scales(bac_plot)$x$range$range
  p_bac <- subset(species_results, org_list == yaxislabel)[, 'padj']
  print(p_bac)
  
  if(p_bac <  0.001){
    bac_plot_1<- bac_plot + annotate("text", x = 0.5+length(xlimit_bac)/2, y = ylimitplot+ 0.05, ymax =ylimitplot+ 1, label = "***", size=8)
  } else if(p_bac <= 0.01){
    bac_plot_1<- bac_plot + annotate("text", x = 0.5+length(xlimit_bac)/2, y = ylimitplot+ 0.05, ymax =ylimitplot+ 1, label = "**", size=6)
  }
  else if(p_bac <= 0.05){
    bac_plot_1<- bac_plot + annotate("text", x = 0.5+length(xlimit_bac)/2, y = ylimitplot+ 0.05, label = "*", size=6)
  } else {
    bac_plot_1 <- bac_plot + annotate("text", x = 0.5+length(xlimit_bac)/2, y = ylimitplot+ 0.05, label = "paste(italic(n.s.))", size=6, parse=TRUE)
  }
  
  #+  coord_cartesian(ylim = c(0, 8)),, outlier.shape = FALSE
  bacplotfname <- paste0(fname,"_speciesthr", species_thr, "_", "genusthr", genus_thr,"_", hist_filter, ".png")
  ggsave(here(foldername,bacplotfname), plot=bac_plot_1, width=10, height=10, units=c("cm"),limitsize = FALSE)
  return(bac_plot_1)
}

remove_outliers <- function(x) {
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  iqr <- q3 - q1
  upper_bound <- q3 + 1.5 * iqr
  lower_bound <- q1 - 1.5 * iqr
  x[x < lower_bound | x > upper_bound] <- NA
  return(x)
}

plotbox_nooutlier <- function(data_for_plot, yaxislabel, title, foldername, fname, ylimitplot, species_results){ 
  #ylimitplot <- boxplot.stats(data_for_plot[[yaxislabel]])$stats[c(1, 5)]
  # get nonzero values for the species of interest
  bac_for_plot <- data_for_plot %>% dplyr::select(LibraryID:AgeG) %>% 
    dplyr::select(-c('Patient.ID', 'Group')) %>% 
    pivot_longer(-c('LibraryID', 'AgeG')) %>% 
    filter(name == yaxislabel) %>% 
    group_by(AgeG) %>% 
    mutate(value = remove_outliers(value)) %>% 
    filter(!is.na(value))
  
  bac_plot <- bac_for_plot %>% 
    ggplot(aes(x = AgeG, y = value, outline =  FALSE, fill = AgeG)) +
    geom_boxplot(width = 1.0, outlier.shape=NA) +
    ggtitle(title) + theme_bw() +
    theme(
      legend.position = 'none', 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      plot.title = element_text(
        color = "black",
        size = 14,
        face = "bold.italic",
        hjust = 0.5
      ),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14, face = 'bold'),
      axis.text.x = element_text(
        size = 14,
        angle = 40,
        hjust = 1
      ),
      axis.text.y = element_text(size = 14),
      text = element_text(family = 'Arial')
    ) + labs(x = "",
             y = "Relative abundance",
             face = "bold",
             size = 14) + scale_fill_manual(values = c(
               "#a50023",
               "#f47d4a",
               "#fed98a",
               "#add39f",
               '#6da5cb',
               "#354a99"
             )) + coord_cartesian(ylim = c(0, ylimitplot))
  ylimit_bac <- layer_scales(bac_plot)$y$range$range
  xlimit_bac <- layer_scales(bac_plot)$x$range$range
  p_bac <- subset(species_results, org_list == yaxislabel)[, 'padj']
  print(p_bac)
  
  if(p_bac <  0.001){
    bac_plot_1<- bac_plot + annotate("text", x = 0.5+length(xlimit_bac)/2, y = ylimitplot+ 0.05, ymax =ylimitplot+ 1, label = "***", size=8)
  } else if(p_bac <= 0.01){
    bac_plot_1<- bac_plot + annotate("text", x = 0.5+length(xlimit_bac)/2, y = ylimitplot+ 0.05, ymax =ylimitplot+ 1, label = "**", size=6)
  }
  else if(p_bac <= 0.05){
    bac_plot_1<- bac_plot + annotate("text", x = 0.5+length(xlimit_bac)/2, y = ylimitplot+ 0.05, label = "*", size=6)
  } else {
    bac_plot_1 <- bac_plot + annotate("text", x = 0.5+length(xlimit_bac)/2, y = ylimitplot+ 0.05, label = "paste(italic(n.s.))", size=6, parse=TRUE)
  }
  
  #+  coord_cartesian(ylim = c(0, 8)),, outlier.shape = FALSE
  bacplotfname <- paste0(fname,"_speciesthr", species_thr, "_", "genusthr", genus_thr,"_", hist_filter, ".png")
  ggsave(here(foldername,bacplotfname), plot=bac_plot_1, width=10, height=10, units=c("cm"),limitsize = FALSE)
  return(bac_plot_1)
}
## Plot volcano

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

plot_volcano <- function(taxa_profiles, n_nonzero, species_results, title) {
  median_of_taxa_with_name <- get_median_of_taxa(taxa_profiles, n_nonzero)
  species_results$Significant <-
    ifelse(species_results$padj < 0.05, "FDR < 0.05", "Not Sig")
  colnames(median_of_taxa_with_name)[1] <- c('org_list')
  
  species_results_merge <-
    merge(species_results, median_of_taxa_with_name, by = 'org_list')
  
  species_results_merge$highlight1 <-
    ifelse(species_results_merge$median_of_taxa > 0 &
             species_results$padj < 0.05, "2", 
           ifelse(species_results_merge$median_of_taxa == 0 &
                    species_results$padj < 0.05 |species_results$padj == 0.05 , "1", "0"))
  
  #species_to_mark <- c('s__Ruminococcus_obeum','s__Coprococcus_catus', 's__Alistipes_senegalensis', 's__Dorea_formicigenerans')    
  species_results_merge$highlight2 <-  ifelse(species_results_merge$betaestimate < -12 , "3", "0")
  #ifelse(species_results_merge$betaestimate > 20, "3",                                
  species_results_merge$highlight3 <- ifelse(species_results_merge$org_list == 'Alistipes_senegalensis', '-1', 0)            
  
  
  #species_results_merge$org_list <- str_replace_all(species_results_merge$org_list, 's__', '')
  species_results_merge$org_list <- str_replace_all(species_results_merge$org_list, '_', ' ')
  
  volcanoplot <- ggplot(species_results_merge, aes(x = betaestimate, y = -log10(padj))) +
    geom_point(aes(color = highlight1), size=3) +
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
      size = 3.5, xlim=c(-45, 20),ylim=c(0,10),
      box.padding = unit(0.75, "lines"),
      point.padding = unit(0.35, "lines"), max.overlaps = 30, fontface = 'italic',nudge_x = .15,
      nudge_y = .5
    ) +
    geom_text_repel(
      data = subset(species_results_merge, highlight2 >2),
      aes(label = org_list),
      size = 3.5, ylim=c(0,10), xlim=c(-45, 20),
      box.padding = unit(0.65, "lines"),
      point.padding = unit(0.65, "lines"), max.overlaps = 30, fontface = 'italic',nudge_x = .15,
      nudge_y = .5
    )+
    geom_text_repel(
      data = subset(species_results_merge, highlight3 <0),
      aes(label = org_list),
      size = 3.5, ylim=c(0,10), xlim=c( 20,0),
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

# Function to plot editable ppt
create_dml <- function(plot){
  rvg::dml(ggobj = plot)
}

plot_editable <- function(obj, h, w, fname) {
  obj_dml <- create_dml(obj)
  officer::read_pptx() %>%
    # add slide ----
  officer::add_slide() %>%
    # specify object and location of object ----
  officer::ph_with(obj_dml, ph_location(width = w, height = h)) %>%
    # export slide -----
  base::print(
    target = here::here(
      editablepptfname,
      fname
    )
  )
}

# to plot many plots at the same time
create_pptx <- function(plot, path, h, w){
  
  # if file does not yet exist, create new PowerPoint ----
  if (!file.exists(path)) {
    out <- officer::read_pptx()
  }
  # if file exist, append slides to exisiting file ----
  else {
    out <- officer::read_pptx(path)
  }
  
  out %>% 
    officer::add_slide() %>% 
    officer::ph_with(plot, location = officer::ph_location(
      width = w, height = h)) %>% 
    base::print(target = path)
}


plotpcoa_bycohort <-
  function(dat_merged_all, eigen, title) {
    ggplot(dat_merged_all,
           aes(x = V1,
               y = V2)) +
      geom_point(
        data = subset(dat_merged_all, Twogroups == 'Elderly'),
        aes(color = Cohort),
        size = 4,
      ) + geom_point(
        data = subset(dat_merged_all,
                      Twogroups == 'Young'),
        aes(color = Cohort),
        size = 4,
      ) +
      scale_colour_manual(values = c(
        "#a50023",
        "#f47d4a",
        "#fed98a",
        "#add39f",
        '#6da5cb',
        "#354a99"
      )) +
      labs(
        x = paste0('PCoA1 (', round(eigen[1], 1), '%)'),
        y = paste0('PCoA2 (', round(eigen[2], 1), '%)'),
        shape = 'Age in decades',
        colour = 'Age groups'
      ) +
      theme_bw() +
      theme(
        legend.background = element_blank(),
        legend.position = "right",
        legend.key = element_blank(),
        legend.text = element_text(size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title       = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)
      ) + ggtitle(title)
  }

plotpcoa_colorbyagegp <- function(dat_merged_all, eigen, title, filtercriteria = NULL)
{
  #geom_density_2d()+
  #cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
  #   "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999")
  if (!is.null(filtercriteria))
  {
    dat_merged_all <- dat_merged_all %>% 
      dplyr::filter(Cohort %in% filtercriteria)
    #Leaving eigen for now
  }
  
  cbp2 = c(
    '#e6194B',
    '#3cb44b',
    '#ffe119',
    '#4363d8',
    '#f58231',
    '#42d4f4',
    '#f032e6',
    '#fabed4', 
    '#469990'
  )
  ggplot(dat_merged_all,
         aes(x = V1,
             y = V2, colour=as.factor(Age))) +
    geom_point(
      data = subset(dat_merged_all, AgeG == '0-20'),
      aes(color = AgeG, shape = Cohort),
      size = 4) +
    geom_point(
      data = subset(dat_merged_all, AgeG == '21-40'),
      aes(color = AgeG, shape = Cohort),
      size = 4) +
    geom_point(
      data = subset(dat_merged_all, AgeG == '41-60'),
      aes(color = AgeG, shape = Cohort),
      size = 4) +
    geom_point(
      data = subset(dat_merged_all, AgeG == '41-50'),
      aes(color = AgeG, shape = Cohort),
      size = 4) +
    geom_point(
      data = subset(dat_merged_all, AgeG == '51-60'),
      aes(color = AgeG, shape = Cohort),
      size = 4) +
    geom_point(
      data = subset(dat_merged_all, AgeG == '61-70'),
      aes(color = AgeG, shape = Cohort),
      size = 4) +
    geom_point(
      data = subset(dat_merged_all, AgeG == '71-80'),
      aes(color = AgeG, shape = Cohort),
      size = 4) +
    geom_point(
      data = subset(dat_merged_all, AgeG == '81-90'),
      aes(color = AgeG, shape = Cohort),
      size = 4) +
    geom_point(
      data = subset(dat_merged_all, AgeG == '91-100'),
      aes(color = AgeG, shape = Cohort),
      size = 4) +
    labs(
      x = paste0('PCoA1 (', round(eigen[1], 1), '%)'),
      y = paste0('PCoA2 (', round(eigen[2], 1), '%)'),
      colour = 'Age Group'
    ) +
    theme_bw() +
    theme(
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(size = 14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title       = element_text(face = "bold", size = 14),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14)
    ) + ggtitle(title)+   scale_color_manual(values=cbp2)
}
#'#cc79A7',
#'#geom_boxplot(width = 0.5, lwd = 1.2, outlier.shape = NA) +
plot_alpha <- function(alphadiv, title, y){
alpha_plot <-
  ggplot(alphadiv, aes(
    x = AgeG,
    y = .data[[y]],
    outline = FALSE,
    fill = AgeG
  )) +
  geom_violin() + geom_boxplot(width = 0.2, outlier.shape = NA) +
  scale_fill_manual(values = c(
    "#a50023",
    "#f47d4a",
    "#fed98a",
    "#add39f",
    '#6da5cb',
    "#354a99"
  )) + theme_bw() + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(
      color = "black",
      size = 14,
      hjust = 0.5
    ),
    axis.title.x = element_text(size = 14, face='bold'),
    axis.title.y = element_text(size = 14, face='bold'),
    axis.text.x = element_text(
      size = 14,
      angle = 40,
      hjust = 1
    ),
    axis.text.y = element_text(size = 14), legend.position = 'none'
  ) + xlab('Age group') + ylab(title)
return(alpha_plot)
}

plot_uniqueness <- function(uniqueness, title, y){
  beta_plot <-
    ggplot(uniqueness, aes(
      x = AgeG,
      y = .data[[y]],
      outline = FALSE,
      fill = AgeG
    )) +
   geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = c(
      "#a50023",
      "#f47d4a",
      "#fed98a",
      "#add39f",
      '#6da5cb',
      "#354a99"
    )) + coord_cartesian(ylim = c(0.1, 0.7))+
    theme_bw() + theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(
        color = "black",
        size = 14,
        hjust = 0.5
      ),
      axis.title.x = element_text(size = 14, face='bold'),
      axis.title.y = element_text(size = 14, face='bold'),
      axis.text.x = element_text(
        size = 14,
        angle = 40,
        hjust = 1
      ),
      axis.text.y = element_text(size = 14), legend.position = 'none'
    ) + xlab('Age group') + ylab(title) +
    scale_y_continuous(breaks=scales::pretty_breaks(n=3))
  return(beta_plot)}


finddist <- function(tax_profiles, metadata_order)
{
  Dist_mat <- vegdist(t(tax_profiles), method = 'bray', upper = TRUE)
  df_D <- melt(as.matrix(Dist_mat), varnames = c("row", "col")) %>% 
    distinct(row, col, .keep_all = TRUE)
  
  #df[as.numeric(df_D$row) > as.numeric(df_D$col), ]
  df_D <- unique(t(apply(df_D, 1, sort))) %>% data.frame()
  df_D <- subset(df_D, df_D$X2 != df_D$X3) %>%
    rename('LibraryID'='X2', 'LibraryID2'='X3')
  df_D <- left_join(df_D, metadata_order, by='LibraryID')
  df_D %<>% dplyr::select('X1','LibraryID',  'LibraryID2','AgeG') %>% 
    rename('LibraryID1'='LibraryID', 'LibraryID'='LibraryID2') 
  df_D <- left_join(df_D, metadata_order, by='LibraryID') %>%
    dplyr::select('X1','LibraryID1',  'LibraryID','AgeG.x', 'AgeG.y')
  toplot <- subset(df_D, df_D$AgeG.x ==df_D$AgeG.y)
  return(toplot)
}

plotbox <- function(data_for_plot, xaxislabel, yaxislabel, title, foldername, fname, ylimitplot, species_results){ 
  #ylimitplot <- boxplot.stats(data_for_plot[[yaxislabel]])$stats[c(1, 5)]
  # Initialize an empty list to store filtered data for each age group
  filtered_data_list <- list()
  
  # Iterate over unique age groups
  for (age_group in unique(data_for_plot$AgeG)) {
    # Subset data for the current age group
    data_age_group <- subset(data_for_plot, AgeG == age_group)
    
    # Calculate quartiles for the current age group
    Q1 <- quantile(data_age_group[[yaxislabel]], 0.25)
    Q3 <- quantile(data_age_group[[yaxislabel]], 0.75)
    IQR <- Q3 - Q1
    
    # Define outlier boundaries for the current age group
    lower_bound <- Q1 - 1.5 * IQR
    upper_bound <- Q3 + 1.5 * IQR
    
    # Filter data for the current age group
    filtered_data <- data_age_group[data_age_group[[yaxislabel]] >= lower_bound & data_age_group[[yaxislabel]] <= upper_bound, ]
    
    # Store filtered data in the list
    filtered_data_list[[age_group]] <- filtered_data
  }
  
  # Combine filtered data from all age groups
  filtered_data <- do.call(rbind, filtered_data_list)
  
  bac_plot <-
    ggplot(data_for_plot,
           aes_string(
             x = xaxislabel,
             y = (yaxislabel),
             outline =  FALSE,
             fill = xaxislabel
           )) +
    geom_boxplot(width = 1, outlier.shape = NA) +
    geom_point(data = filtered_data, 
               position = position_jitterdodge(dodge.width = 0.4),
                size = 4, alpha = 0.5, pch = 21,
                color = 'black') +
    ggtitle(title) + theme_bw() +
    theme(
      legend.position = 'none', 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      plot.title = element_text(
        color = "black",
        size = 14,
        face = "italic",
        hjust = 0.5
      ),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14, face = 'bold'),
      axis.text.x = element_text(
        size = 12,
        angle = 40,
        hjust = 1
      ),
      axis.text.y = element_text(size = 12)
    ) + labs(x = "",
             y = "Relative Abundance",
             face = "bold",
             size = 14) + scale_fill_manual(values = c(
               "#a50023",
               "#f47d4a",
               "#fed98a",
               "#add39f",
               '#6da5cb',
               "#354a99"
             )) +
    # scale_color_manual(values = c(
  #              "#a50023",
  #              "#f47d4a",
  #              "#fed98a",
  #              "#add39f",
  #              '#6da5cb',
  #              "#354a99"
  #            )) + 
  coord_cartesian(ylim = c(0, ylimitplot))
  ylimit_bac <- layer_scales(bac_plot)$y$range$range
  xlimit_bac <- layer_scales(bac_plot)$x$range$range
  p_bac <- subset(species_results, org_list == yaxislabel)[, 'padj']
  print(p_bac)
  
  # if(p_bac <  0.001){
  #   bac_plot_1<- bac_plot + annotate("text", x = 0.5+length(xlimit_bac)/2, y = ylimitplot+ 0.03, ymax =ylimitplot+ 1, label = "***", size=12)
  # } else if(p_bac <= 0.01){
  #   bac_plot_1<- bac_plot + annotate("text", x = 0.5+length(xlimit_bac)/2, y = ylimitplot+ 0.03, ymax =ylimitplot+ 1, label = "**", size=12)
  # }
  # else if(p_bac <= 0.05){
  #   bac_plot_1<- bac_plot + annotate("text", x = 0.5+length(xlimit_bac)/2, y = ylimitplot+ 0.03, label = "*", size=12)
  # }else{
  #   bac_plot_1 <- bac_plot + annotate("text", x = 0.5+length(xlimit_bac)/2, y = ylimitplot+ 0.03, label = "paste(italic(n.s.))", size=10, parse=TRUE)
  # }
  # 
  #+  coord_cartesian(ylim = c(0, 8)),, outlier.shape = FALSE
  bacplotfname <- paste0(fname, ".png")
  ggsave(here(foldername,bacplotfname), plot=bac_plot, width=10, height=10, units=c("cm"),limitsize = FALSE)
  return(bac_plot)
  }

# Function to parse the results of associations
param_with_age_correction <-
  function(org_list_to_use,
           resultssum,
           foldername,
           fname, data_for_association) {
    pval <- NULL
    org_list <- NULL
    estimate <- NULL
    stderror <- NULL
    betaestimate <- NULL
    stderror <- NULL
    for (i in seq_along(1:length(resultssum)))
    {
      #print(paste0(org_list_to_use[i], ' ', resultssum[[i]][35]))
      if ("x" %in% rownames(resultssum[[i]])) {
        pval <-
          rbind(pval, resultssum[[i]]["x",][4])
        org_list <- rbind(org_list, org_list_to_use[i])
        betaestimate <-
          rbind(betaestimate, resultssum[[i]]["x",][1])
        stderror <-
          rbind(stderror, resultssum[[i]]["x",][2])
      }
      else{
        print(org_list_to_use[i])
      }
    }
    
    #Confidence intervals
    confint_val <- data.frame()
    for (x in seq_along(org_list_to_use)) {
      print(x)
      formula_glm <-
        as.formula(sprintf(
          "Age ~ %s + FBGlu +Tgly + TCHOL + HDL + LDL + BMI + Gender",
          org_list_to_use[x]
        ))
      model_glm <- glm(formula_glm, data = data_for_association)
      confidence_int <-
        confint(model_glm, param = org_list_to_use[x], level = 0.95)
      confint_needed <-
        cbind(org_list_to_use[x], confidence_int[org_list_to_use[x],][1], confidence_int[org_list_to_use[x],][2])
      confint_val <- rbind(confint_val, confint_needed)
    }
    #browser()
    colnames(confint_val) <- c('org_list', '2.5%', '97.5%')
    all_df <- data.frame(org_list, pval, betaestimate, stderror)
    colnames(all_df) <-
      c('org_list', 'pval', 'betaestimate', 'stderror')
    all_df$padj <- p.adjust(all_df$pval, 'fdr')
    all_df_withconfint <- merge(all_df, confint_val)
    all_df_withconfint <- all_df_withconfint %>% 
      mutate(org_list = str_remove(org_list, 's__'))
    write.csv(all_df_withconfint, here(foldername, fname), row.names = FALSE)
    return(all_df_withconfint)
  }


plot_volcano_asscn <- function(species_results_merge){
  volcanoplot <- ggplot(species_results_merge, aes(x = betaestimate, y = -log10(padj))) +
    geom_point(aes(color = highlight1), size=2) +
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
      size = 3.1,xlim=c(-45, 20),ylim=c(0,10),
      box.padding = unit(0.75, "lines"),
      point.padding = unit(0.35, "lines"), max.overlaps = 30, fontface = 'italic',nudge_x = .15,
      nudge_y = .5
    ) +
    geom_text_repel(
      data = subset(species_results_merge, highlight2 >2),
      aes(label = org_list),
      size = 3.1, ylim=c(0,10), xlim=c(-45, 20),
      box.padding = unit(0.65, "lines"),
      point.padding = unit(0.65, "lines"), max.overlaps = 30, fontface = 'italic',nudge_x = .15,
      nudge_y = .5
    )+
    geom_text_repel(
      data = subset(species_results_merge, highlight3 <0),
      aes(label = org_list),
      size = 3.1, ylim=c(0,10), xlim=c( 20,0),
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
    scale_x_continuous(breaks = scales::pretty_breaks(5))+
    labs(x= "\u03b2 coefficient", y= expression(paste(-log[10]," p-value")))+ scale_color_manual('', values = c("0"="grey","1"="green", "2"= "red"), labels = c('p-value>0.05', 'p-value \u2264 0.05 (med = 0)', 'p-value \u2264 0.05 (med > 0)' ))
  
  #print(volcanoplot)
  return(volcanoplot)
}

plot_bray_curtis <- function(combined_dist){
 bc_plot <- combined_dist %>% 
    ggplot(aes(x = Group, y = value, fill = Group)) + 
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = c(
      "#a50023",
      "#f47d4a",
      "#fed98a",
      "#add39f",
      '#6da5cb',
      "#354a99")) +
    #ylim(c(0.4, 1)) +
    coord_cartesian(ylim = c(0.4, 1)) + 
    labs(x = 'Age Group', y = 'Bray-Curtis Dissimilarity Index') +
    theme_bw() +
    theme(
      legend.background = element_blank(),
      legend.position = "none",
      legend.key = element_blank(),
      legend.text = element_text(size = 11),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title       = element_text(face = "bold", size = 14),
      plot.title = element_text(face = 'bold', size = 14), 
      axis.text.x = element_text(
        size = 14,
        angle = 40,
        hjust = 1),
      axis.text.y = element_text(size = 14),
      text = element_text(family = 'Arial')
    )
  
  return(bc_plot)
}


plot_volcano_genus <- function(genus_results_merge){
  genusvolcanoplot <- ggplot(genus_results_merge, aes(x = betaestimate, y = -log10(padj))) +
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
          axis.title = element_text(face="bold", size=16),
          axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),
          axis.text = element_text(size=16),
          legend.position="bottom",
          legend.title = element_text(size=12, face="bold"),
          plot.title = element_text(hjust = 0.5, size=16),
          legend.text = element_text(size=8)) +
    geom_text_repel(
      data = subset(genus_results_merge, highlight1 >1),
      aes(label = org_list),
      size = 4,
      box.padding = unit(0.65, "lines"),
      point.padding = unit(0.65, "lines"), max.overlaps = 30, fontface = 'italic'
    ) +
    labs(x= "\u03b2 coefficient", y= expression(paste(-log[10]," p-value")))+ scale_color_manual('', values = c("0"="grey","1"="green", "2"= "red"), labels = c('p-value>0.05', 'p-value \u2264 0.05 (med = 0)', 'p-value \u2264 0.05 (med > 0)' )) 
  return(genusvolcanoplot)
  
}