pacman::p_load(tidyverse, plyr, ggrepel)

rm(list = ls())

glu_data <-
  read.table(
    file = "../data/glu_genes.csv",
    sep = ",",
    header = TRUE,row.names='Index')
glu_data <- as.matrix((glu_data))
glu_data  <- glu_data * 10^6 #IN CPM
sg_metadata <- read.table("../data/sg90_metadata.txt",sep = "\t",header = TRUE)
glu_mapping_table <- read.table("../data/glutarate_uniref90_ko_individual.txt", sep='\t', col.names=c('UniRef90','KIDS','GeneName'))
glu_mapping_table <- glu_mapping_table[-1,]
mdf <- reshape2::melt(glu_data)
colnames(mdf) <- c('Uniref', 'libID', 'value')
mdf$sampletype <- sg_metadata$group[match(mdf$libID, sg_metadata$Sample.ID)]
mdf$kids <- glu_mapping_table$KIDS[match(mdf$Uniref, glu_mapping_table$UniRef90)]
mdf$geneids <- glu_mapping_table$GeneName[match(mdf$Uniref, glu_mapping_table$UniRef90)]
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
      plot.title = element_text(face = "italic", hjust=0.5,size=32),
      strip.background = element_rect(fill = "grey90"),
      legend.position = "none",
      legend.text = element_text(size=13),
      legend.title = element_text(size=15),
      axis.text.x = element_text(size = 35, color = "black"),
      axis.text.y = element_text(size = 35, color = "black"),
      panel.border = element_rect(colour = "black", fill=NA, size=1.5)
    ) + ylab('Normalized counts') + ggtitle(geneid) + xlab("") +  labs(fill="Group") 
  #Compute lower and upper whiskers
  ylim1 = boxplot.stats(req_mdf1$value)$stats[c(1, 5)]
  limonaxis = ylim1 * valgiven
  p1 <- p +  coord_cartesian(ylim = limonaxis) + scale_y_continuous(breaks = seq(0, limonaxis[2], by = byval))
  print(p1)}

wantedgctAB = filter(req_mdf, geneids %in% c('gctA', 'gctB') )
makegeneboxplot(wantedgctAB, 0.23, 'gctAB',1)
gctabglm <- summary (glm(wantedgctAB$value~ wantedgctAB$Age, data=wantedgctAB))
wantedgcdAB = filter(req_mdf, geneids %in% c('gcdA', 'gcdB') )
makegeneboxplot(wantedgcdAB, 1.6, 'gcdAB',5)
gcdabglm <- summary (glm(wantedgcdAB$value~ wantedgcdAB$Age, data=wantedgcdAB))
wantedhgdABC = filter(req_mdf, geneids %in% c('hgdA', 'hgdB', 'hgdC') )
makegeneboxplot(wantedhgdABC, 1.2, 'hgcoAdABC', 5)
hgdabcglm <- summary (glm(wantedhgdABC$value~ wantedhgdABC$Age, data=wantedhgdABC))

allpval <- data.frame(geneids=character(), pvalues=numeric()) %>% add_row(geneids='gctAB', pvalues=gctabglm$coefficients[8])%>% 
          add_row(geneids='gcdAB', pvalues=gcdabglm$coefficients[8])%>% 
          add_row(geneids='hgdABC', pvalues=hgdabcglm$coefficients[8])
write.table(allpval, file='../results/pval_glu_pathway.csv', append=TRUE, sep=',', col.names = FALSE)


# Cumulative of all genes
p_meds <- ddply(req_mdf1, .(Group), summarise, med = round(median(value),3))
agegroupplot <- ggplot(data=req_mdf1, aes(x=Group, y=value, fill=Group)) + geom_boxplot(outlier.shape = NA) + scale_fill_manual(values=c("#a50023", "#f47d4a", "#fed98a", "#add39f", '#6da5cb', "#354a99")) +
  geom_text_repel(data = p_meds, aes(x = Group, y = med, label = med), 
                  size = 6,vjust = -0.5, fontface='bold') + 
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold", size = 18),
    strip.text = element_text(face = "bold", size = 13),
    plot.title = element_text(face = "bold", hjust=0.5,size=20),
    strip.background = element_rect(fill = "grey90"),
    legend.position = "none",
    legend.text = element_text(size=13),
    legend.title = element_text(size=15),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black")
  ) + ylab('Normalised counts') + ggtitle("Glutarate pathway") + xlab("") +  labs(fill="Group") 

print(agegroupplot)

ylim1 <-boxplot.stats(req_mdf1$value)$stats[c(1, 5)]
limonaxis = ylim1 * 1.4
agegroupplot1 <- agegroupplot +  coord_cartesian(ylim = limonaxis)
print(agegroupplot1)

allglm <-summary (glm(req_mdf1$value~ req_mdf1$Age, data=req_mdf1))
allglm$coefficients[8]
