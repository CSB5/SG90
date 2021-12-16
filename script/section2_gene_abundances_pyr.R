pacman::p_load(tidyverse, plyr, ggrepel)

rm(list = ls())

pyr_data <-
  read.table(
    file = "../data/pyr_genes.csv",
    sep = ",",
    header = TRUE,row.names='Index')
pyr_data <- as.matrix((pyr_data))
pyr_data  <- pyr_data * 10^6 #IN CPM
sg_metadata <- read.table("../data/sg90_metadata.txt",sep = "\t",header = TRUE)
pyr_mapping_table <- read.table("../data/pyr_uniref90_ko_individual.txt", sep='\t', col.names=c('UniRef90','KIDS','GeneName'))
pyr_mapping_table <- pyr_mapping_table[-1,]
mdf <- reshape2::melt(pyr_data)
colnames(mdf) <- c('Uniref', 'libID', 'value')
mdf$sampletype <- sg_metadata$group[match(mdf$libID, sg_metadata$Sample.ID)]
mdf$kids <- pyr_mapping_table$KIDS[match(mdf$Uniref, pyr_mapping_table$UniRef90)]
mdf$geneids <- pyr_mapping_table$GeneName[match(mdf$Uniref, pyr_mapping_table$UniRef90)]
mdf$Group <- sg_metadata$Position[match(mdf$libID, sg_metadata$Sample.ID)]
mdf$Age <- sg_metadata$age[match(mdf$libID, sg_metadata$Sample.ID)]
req_mdf <- mdf[mdf$value !=0,]
req_mdf1 <- req_mdf[complete.cases(req_mdf),]


# Gene plots
genes <-list(unique(req_mdf$geneids) %>% as.character())
makegeneboxplot<- function(reqmdf, valgiven, geneid, byval){
  p <- ggplot(data=reqmdf, aes(x=Group, y=value, fill=Group)) + geom_boxplot(width=1,outlier.shape = NA) + scale_fill_manual(values=c("#a50023", "#f47d4a", "#fed98a", "#add39f", '#6da5cb', "#354a99"))+
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold", size = 35),
      strip.text = element_text(face = "bold", size = 13),
      plot.title = element_text(face = "italic", hjust=0.5,size=45),
      strip.background = element_rect(fill = "grey90"),
      legend.position = "none",
      legend.text = element_text(size=13),
      legend.title = element_text(size=15),
      axis.text.x = element_text(size = 35, color = "black"),
      axis.text.y = element_text(size = 35, color = "black"),
      panel.border = element_rect(colour = "black", fill=NA, size=1.5)
    ) + ylab('Normalized counts') + ggtitle(geneid) + xlab("") +  labs(fill="Group") 
  #Compute lower and upper whiskers
  ylim1 <- boxplot.stats(req_mdf1$value)$stats[c(1, 5)]
  limonaxis <- ylim1 * valgiven
  p1 <- p +  coord_cartesian(ylim = limonaxis) + scale_y_continuous(breaks = seq(0, limonaxis[2], by = byval))
  print(p1)}

wantedthl <- filter(req_mdf, geneids %in% c('thl') )
makegeneboxplot(wantedthl, 1.28, 'thl',5)
summary (glm(wantedthl$value~ wantedthl$Age, data=wantedthl))
wantedhbd <- filter(req_mdf, geneids %in% c('hbd') )
makegeneboxplot(wantedhbd, 1.0, 'hbd',5)
hbdglm <- summary (glm(wantedhbd$value~ wantedhbd$Age, data=wantedhbd))
wantedcro<- filter(req_mdf, geneids %in% c('croR') )
wantedcro<-rbind(wantedcro, data.frame(Uniref = 'dummy1', libID='DUMMY1', value= 0, sampletype='Elderly', kids='K00998', geneids='croR', Group='71-80', Age=75))
wantedcro<-rbind(wantedcro, data.frame(Uniref = 'dummy1', libID='DUMMY1', value= 0, sampletype='Elderly', kids='K00998', geneids='croR', Group='91-100', Age=95))
makegeneboxplot(wantedcro, 0.13, 'croR', 0.5)
croglm <- summary (glm(wantedcro$value~ wantedcro$Age, data=wantedcro))
fortt_thl <- data.frame(wantedthl$Age, wantedthl$value)
thlglm <- summary (glm(fortt_thl$wantedthl.Age~ fortt_thl$wantedthl.value, data=fortt_thl))
allpval <- data.frame(geneids=character(), pvalues=numeric()) %>% add_row(geneids='croR', pvalues=croglm$coefficients[8])%>% add_row(geneids='thl', pvalues=thlglm$coefficients[8])%>% add_row(geneids='hbd', pvalues=hbdglm$coefficients[8])
write.csv(allpval, file='../results/pval_pyr_pathway.csv')

# Cumulative of all genes
p_meds <- ddply(req_mdf1, .(Group), summarise, med = round(median(value),3))
agegroupplot <- ggplot(data=req_mdf1, aes(x=Group, y=value, fill=Group)) + geom_boxplot(outlier.shape = NA) + scale_fill_manual(values=c("#a50023", "#f47d4a", "#fed98a", "#add39f", '#6da5cb', "#354a99")) +
  geom_text_repel(data = p_meds, aes(x = Group, y = med, label = med), 
                  size = 6,vjust = -0.5, fontface='bold') +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold", size = 15),
    strip.text = element_text(face = "bold", size = 13),
    plot.title = element_text(face = "bold", hjust=0.5,size=15),
    strip.background = element_rect(fill = "grey90"),
    legend.position = "none",
    legend.text = element_text(size=13),
    legend.title = element_text(size=13),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black")
  ) + ylab('Normalised counts') + ggtitle("Pyruvate pathway") + xlab("") +  labs(fill="Group") 
print(agegroupplot)
ylim1 <-boxplot.stats(req_mdf1$value)$stats[c(1, 5)]
limonaxis = ylim1 *1.3
agegroupplot1 <- agegroupplot +  coord_cartesian(ylim = limonaxis)
print(agegroupplot1)
 

allglm <-summary (glm(req_mdf1$value~ req_mdf1$Age, data=req_mdf1))
allglm$coefficients[8]