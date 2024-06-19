library(tidyverse)
library(ggrepel)
pacman::p_load(ggridges)

plotpcoa <- function(dat_merged_all, eigen, title) {ggplot(dat_merged_all,
                                                           aes(x = V1,
                                                               y = V2)) +
    geom_point(
      data = subset(dat_merged_all, metadata_order.twogroups == 'Elderly'),
      aes(color = metadata_order.decade, shape= metadata_order.cohortdetails),
      size = 4,
    ) + geom_point(
      data = subset(dat_merged_all,
                    metadata_order.twogroups == 'Young'),
      aes(color = metadata_order.decade,shape= metadata_order.cohortdetails),
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
      legend.key = element_blank(),
      legend.text = element_text(size = 14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title       = element_text(face = "bold", size = 14),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14)
    ) + ggtitle(title)}
plotpcoa_colorbyage <- function(dat_merged_all, eigen, title) {
  #geom_density_2d()+
  cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
            "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999")
  ggplot(dat_merged_all,
         aes(x = V1,
             y = V2, colour=as.factor(metadata_order.age))) +
    
    geom_point(
      data = subset(dat_merged_all, metadata_order.twogroups == 'Elderly'),
      aes( shape= metadata_order.cohortdetails),
      size = 4,
    ) + geom_point(
      data = subset(dat_merged_all,
                    metadata_order.twogroups == 'Young'),
      aes(shape= metadata_order.cohortdetails),
      size = 4,
    ) +
    labs(
      x = paste0('PCoA1 (', round(eigen[1], 1), '%)'),
      y = paste0('PCoA2 (', round(eigen[2], 1), '%)'),
      shape = 'Cohort',
      colour = 'Age'
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
    ) + ggtitle(title)+   scale_color_brewer(palette = "Paired")
  #scale_fill_manual(values=cbp2)
}
pcoa2var <- function(dat_merged_all, title) {ggplot(dat_merged_all, aes(x= metadata_order.decade, y=V2)) + 
    geom_violin()  + geom_boxplot(width = 0.2, outlier.shape = NA) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    ) + labs(x = "Group",
             y = "Value",
             face = "bold",
             size = 14) +
    theme(
      plot.title = element_text(
        color = "black",
        size = 14,
        hjust = 0.5
      ),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(
        size = 14,
        angle = 40,
        hjust = 1
      ),
      axis.text.y = element_text(size = 14)
    ) + scale_fill_manual(values = c(
      "#a50023",
      "#f47d4a",
      "#fed98a",
      "#add39f",
      '#6da5cb',
      "#354a99"
    )) + ggtitle(title)}
performpcoa <- function(tax_profiles, metadata_order){
  distmetric <- vegdist(t(tax_profiles), method = "bray")
  cmds <- cmdscale(distmetric, k = 3, eig = TRUE)
  eigen <- cmds$eig / sum(cmds$eig) * 100
  cmdvals <- as.data.frame(cmds$points) %>% rownames_to_column('LibraryID')
  #browser()
  #meta_merged <-
   # data.frame(metadata_order$LibraryID,
  #             metadata_order$agegrp, metadata_order$decade, metadata_order$cohortdetails,metadata_order$twogroups, metadata_order$age) #Position= decade, site=cohortdetails
  #row.names(meta_merged) <- metadata_order$LibraryID
  #dat_merged_all <- data.frame(cmdvals, meta_merged)
  dat_merged_all <- merge(cmdvals, metadata_order, by='LibraryID')
  listofop <- list("dat_merged_all" = dat_merged_all, "eigen"=eigen, 
                   "distmetric"=distmetric)
  return(listofop)
  
}
plotpcoa_bycohort <- function(dat_merged_all, eigen, title)
  {ggplot(dat_merged_all, aes(x = V1, y = V2)) +
    geom_point(
      data = subset(dat_merged_all, twogroups == 'Elderly'),
      aes(color = cohortdetails),
      size = 4,
    ) + geom_point(
      data = subset(dat_merged_all,
                    twogroups == 'Young'),
      aes(color = cohortdetails),
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
    ) + ggtitle(title)}
plot_density <- function(dat_merged_age_adj) {
  pTop <-
    ggplot(
      data = dat_merged_age_adj$dat_merged_all,
      aes(
        x = dat_merged_age_adj$dat_merged_all$V1,
        fill = dat_merged_age_adj$dat_merged_all$metadata_order.cohortdetails
      )
    ) + geom_density(size = 0.2, alpha = 0.7) +
    ylab("Density") + theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = rel(0.8)),
      plot.title = element_text(hjust = 0.5,
                                size = rel(1.5)),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_blank(),
      legend.position = "none"
    ) +
    scale_fill_manual(values = c("#a50023",
                                 "#f47d4a",
                                 "#fed98a"))
  pRight <- ggplot(
    data = dat_merged_age_adj$dat_merged_all,
    aes(
      x = dat_merged_age_adj$dat_merged_all$V2,
      fill = dat_merged_age_adj$dat_merged_all$metadata_order.cohortdetails
    )
  ) + geom_density(size = 0.2, alpha = 0.7) +
    coord_flip() + ylab("Density") + theme(
      axis.title.x = element_text(size = rel(0.8)),
      axis.title.y = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_blank(),
      legend.position = "none"
    ) +
    scale_fill_manual(values = c("#a50023",
                                 "#f47d4a",
                                 "#fed98a"))
  listofden <-list("top"=pTop, "right"=pRight)
  return(listofden)
}

plot_pcoa_lab_study <- function(dat_merged_all, eigen, title)
{
  ggplot(dat_merged_all,
         aes(x = V1,
             y = V2)) +
    geom_point(
      data = subset(dat_merged_all, Study == 'SPMP'),
      aes(shape = Study, color = AgeG),
      size = 4,
    ) +geom_point(
      data = subset(dat_merged_all, Study == 'SG90'),
      aes(shape = Study, color = AgeG),
      size = 4,
    )+
    geom_point(
      data = subset(dat_merged_all, Study == 'CRE'),
      aes(shape = Study, color = AgeG),
      size = 4,
    )+
    geom_point(
      data = subset(dat_merged_all, Study == 'Amili'),
      aes(shape = Study, color = AgeG),
      size = 4,
    )+
    geom_point(
      data = subset(dat_merged_all, Study == 'ST131'),
      aes(shape = Study, color = AgeG),
      size = 4,
    )+
    scale_colour_manual(values = c(
      # "#a50023",
      # "#f47d4a",
      # "#fed98a",
      # "#add39f",
      # '#6da5cb',
      # "#354a99",
      '#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990'
    )) + scale_shape_manual(values=c(15,16,17,18,19)) +
    labs(
      x = paste0('PCoA1 (', round(eigen[1], 1), '%)'),
      y = paste0('PCoA2 (', round(eigen[2], 1), '%)'),
      
      colour = 'Cohorts'
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
    ) + ggtitle(title)}

plot_pcoa_lab_study_batcheffectadj <- function(dat_merged_all, eigen, title)
{
  ggplot(dat_merged_all,
         aes(x = V1,
             y = V2, label=Row.names)) +
    geom_point(
      data = subset(dat_merged_all, metadata_order.Study == 'SPMP'),
      aes(color = metadata_order.Study),
      size = 4,
    ) +geom_point(
      data = subset(dat_merged_all, metadata_order.Study == 'SG90'),
      aes(color = metadata_order.Study),
      size = 4,
    )+
    geom_point(
      data = subset(dat_merged_all, metadata_order.Study == 'CRE'),
      aes(color = metadata_order.Study),
      size = 4,
    )+
    geom_point(
      data = subset(dat_merged_all, metadata_order.Study == 'Amili'),
      aes(color = metadata_order.Study),
      size = 4,
    )+
    geom_point(
      data = subset(dat_merged_all, metadata_order.Study == 'T2D'),
      aes(color = metadata_order.Study),
      size = 4,
    )+
    geom_point(
      data = subset(dat_merged_all, metadata_order.Study == 'ST131'),
      aes(color = metadata_order.Study),
      size = 4,
    )+ geom_text_repel(
      min.segment.length = 0,
      seed = 42,
      box.padding = 0.5, max.overlaps = 20,
      aes(label = ifelse(
        metadata_order.Study == 'ST131', Row.names, ''
      )) 
    ) + geom_text_repel(
    min.segment.length = 0,
    seed = 42,
    box.padding = 0.5, max.overlaps = 20,
    aes(label = ifelse(
      metadata_order.Study == 'CRE', Row.names, ''
    )))+
    scale_colour_manual(values = c(
      # "#a50023",
      # "#f47d4a",
      # "#fed98a",
      # "#add39f",
      # '#6da5cb',
      # "#354a99",
      '#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990'
    )) + scale_shape_manual(values=c(15,16,17,18,19)) +
    labs(
      x = paste0('PCoA1 (', round(eigen[1], 1), '%)'),
      y = paste0('PCoA2 (', round(eigen[2], 1), '%)'),
      
      colour = 'Cohorts'
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
    ) + ggtitle(title)}

plot_density_filter_ridge <- function(dat_merged_age_adj) {
  pTop <-
    ggplot(
      data = dat_merged_age_adj,
      aes(
        x = dat_merged_age_adj$V1,
        fill = dat_merged_age_adj$metadata_order.Study,
        y= dat_merged_age_adj$metadata_order.Study
      )
    ) + geom_density_ridges(size = 0.2, alpha = 0.7) +
    ylab("Density") + theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = rel(0.8)),
      plot.title = element_text(hjust = 0.5,
                                size = rel(1.5)),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_blank(),
      legend.position = "none"
    ) +
    scale_fill_manual(values = c(
      # "#a50023",
      # "#f47d4a",
      # "#fed98a",
      # "#add39f",
      # '#6da5cb',
      # "#354a99",
      '#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990'
    ))
  pRight <- ggplot(
    data = dat_merged_age_adj,
    aes(
      x = dat_merged_age_adj$V2,
      fill = dat_merged_age_adj$metadata_order.Study,
      y = dat_merged_age_adj$metadata_order.Study
    )
  ) + geom_density_ridges(size = 0.2, alpha = 0.7) +
    coord_flip() + ylab("Density") + theme(
      axis.title.x = element_text(size = rel(0.8)),
      axis.title.y = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_blank(),
      legend.position = "none"
    ) +
    scale_fill_manual(values = c(
      # "#a50023",
      # "#f47d4a",
      # "#fed98a",
      # "#add39f",
      # '#6da5cb',
      # "#354a99",
      '#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990'
    ))
  listofden <-list("top"=pTop, "right"=pRight)
  return(listofden)
}

plot_density_filter <- function(dat_merged_age_adj) {
  pTop <-
    ggplot(
      data = dat_merged_age_adj,
      aes(
        x = V1,
        fill = metadata_order.Study
      )
    ) + geom_density(size = 0.2, alpha = 0.7) +
    ylab("Density") + theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = rel(0.8)),
      plot.title = element_text(hjust = 0.5,
                                size = rel(1.5)),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_blank(),
      legend.position = "none"
    ) +
    scale_fill_manual(values = c(
      # "#a50023",
      # "#f47d4a",
      # "#fed98a",
      # "#add39f",
      # '#6da5cb',
      # "#354a99",
      '#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990'
    ))
  pRight <- ggplot(
    data = dat_merged_age_adj,
    aes(
      x = V2,
      fill = metadata_order.Study
    )
  ) + geom_density(size = 0.2, alpha = 0.7) +
    coord_flip() + ylab("Density") + theme(
      axis.title.x = element_text(size = rel(0.8)),
      axis.title.y = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_blank(),
      legend.position = "none"
    ) +
    scale_fill_manual(values = c(
      # "#a50023",
      # "#f47d4a",
      # "#fed98a",
      # "#add39f",
      # '#6da5cb',
      # "#354a99",
      '#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990'
    ))
  listofden <-list("top"=pTop, "right"=pRight)
  return(listofden)
}



plot_pcoa_lab_study_batcheffectadj_wolabels <- function(dat_merged_all, eigen, title)
{
  ggplot(dat_merged_all,
         aes(x = V1,
             y = V2, label=Row.names)) +
    geom_point(
      data = subset(dat_merged_all, metadata_order.Study == 'SPMP'),
      aes(color = metadata_order.Study),
      size = 4,
    ) +geom_point(
      data = subset(dat_merged_all, metadata_order.Study == 'SG90'),
      aes(color = metadata_order.Study),
      size = 4,
    )+
    geom_point(
      data = subset(dat_merged_all, metadata_order.Study == 'CRE'),
      aes(color = metadata_order.Study),
      size = 4,
    )+
    geom_point(
      data = subset(dat_merged_all, metadata_order.Study == 'Amili'),
      aes(color = metadata_order.Study),
      size = 4,
    )+
    geom_point(
      data = subset(dat_merged_all, metadata_order.Study == 'T2D'),
      aes(color = metadata_order.Study),
      size = 4,
    )+
    geom_point(
      data = subset(dat_merged_all, metadata_order.Study == 'ST131'),
      aes(color = metadata_order.Study),
      size = 4,
    ) + 
    scale_colour_manual(values = c(
      # "#a50023",
      # "#f47d4a",
      # "#fed98a",
      # "#add39f",
      # '#6da5cb',
      # "#354a99",
      '#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990'
    )) + scale_shape_manual(values=c(15,16,17,18,19)) +
    labs(
      x = paste0('PCoA1 (', round(eigen[1], 1), '%)'),
      y = paste0('PCoA2 (', round(eigen[2], 1), '%)'),
      
      colour = 'Cohorts'
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
    ) + ggtitle(title)}


plot_box_pcoa2 <-   function(dat_merged_all, title)
  {
    ggplot(dat_merged_all,
           aes(y = V2,  x = metadata_order.Study, fill=metadata_order.Study)) +
    geom_boxplot()+
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
    ) + ggtitle(title)+ 
    scale_fill_manual(values = c(
      '#e6194B',
      '#3cb44b',
      '#ffe119',
      '#4363d8',
      '#f58231',
      '#42d4f4',
      '#f032e6',
      '#fabed4', 
      '#469990'
    )) + xlab('')}


plot_box_pcoa2_facet <-   function(dat_merged_all, title)
{
  ggplot(dat_merged_all,
         aes(y = V2,  x = metadata_order.Study, fill=metadata_order.Study)) +
    geom_boxplot()+ facet_wrap(~metadata_order.AgeG) +
    theme_bw() +
    theme(
      legend.background = element_blank(),
      legend.position = "none",
      legend.key = element_blank(),
      legend.text = element_text(size = 11),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title       = element_text(face = "bold", size = 14),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    ) + ggtitle(title)+ 
    scale_fill_manual(values = c(
      '#e6194B',
      '#3cb44b',
      '#ffe119',
      '#4363d8',
      '#f58231',
      '#42d4f4',
      '#f032e6',
      '#fabed4', 
      '#469990'
    )) + xlab('')}

plot_box_pcoa2_withage <-   function(dat_merged_all, title)
{
  ggplot(dat_merged_all,
         aes(y = V2,  x = metadata_order.AgeG,fill = metadata_order.AgeG )) +
    geom_boxplot()+ 
    theme_bw() +
    theme(
      legend.background = element_blank(),
      legend.position = "none",
      legend.key = element_blank(),
      legend.text = element_text(size = 11),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title       = element_text(face = "bold", size = 14),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    ) + ggtitle(title)+ 
    scale_fill_manual(values = c(
      '#e6194B',
      '#3cb44b',
      '#ffe119',
      '#4363d8',
      '#f58231',
      '#42d4f4',
      '#f032e6',
      '#fabed4', 
      '#469990'
    )) + xlab('')}



plotpcoa_colorbyage <- function(dat_merged_all, eigen, title) {
  #geom_density_2d()+
  #cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
         #   "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999")
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
             y = V2, colour=as.factor(metadata_order.age))) +
    geom_point(
      data = subset(dat_merged_all, metadata_order.AgeG == '0-20'),
      aes(color = metadata_order.AgeG),
      size = 4) +
  geom_point(
    data = subset(dat_merged_all, metadata_order.AgeG == '21-30'),
    aes(color = metadata_order.AgeG),
    size = 4) +
  geom_point(
    data = subset(dat_merged_all, metadata_order.AgeG == '31-40'),
    aes(color = metadata_order.AgeG),
    size = 4) +
  geom_point(
    data = subset(dat_merged_all, metadata_order.AgeG == '41-50'),
    aes(color = metadata_order.AgeG),
    size = 4) +
  geom_point(
    data = subset(dat_merged_all, metadata_order.AgeG == '51-60'),
    aes(color = metadata_order.AgeG),
    size = 4) +
  geom_point(
    data = subset(dat_merged_all, metadata_order.AgeG == '61-70'),
    aes(color = metadata_order.AgeG),
    size = 4) +
  geom_point(
    data = subset(dat_merged_all, metadata_order.AgeG == '71-80'),
    aes(color = metadata_order.AgeG),
    size = 4) +
  geom_point(
    data = subset(dat_merged_all, metadata_order.AgeG == '81-90'),
    aes(color = metadata_order.AgeG),
    size = 4) +
  geom_point(
    data = subset(dat_merged_all, metadata_order.AgeG == '91-100'),
    aes(color = metadata_order.AgeG),
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
