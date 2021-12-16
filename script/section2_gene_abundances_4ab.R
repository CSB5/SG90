pacman::p_load(tidyverse, plyr, ggrepel)

rm(list = ls())

fourab_data <-
  read.table(
    file = "../data/4ab_genes.csv",
    sep = ",",
    header = TRUE,row.names='Index')
fourab_data <- as.matrix((fourab_data))
fourab_data  <- fourab_data * 10^6 #IN CPM
#fourab_data[is.infinite(fourab_data)] <- 0

sg_metadata <- read.table("../data/sg90_metadata.txt",sep = "\t",header = TRUE)
fourab_mapping_table <- read.table("../data/4ab_uniref90_ko_individual.txt", sep='\t', col.names=c('UniRef90','KIDS','GeneName'))
fourab_mapping_table <- fourab_mapping_table[-1,]
mdf <- reshape2::melt(fourab_data)
colnames(mdf) <- c('Uniref', 'libID', 'value')
mdf$sampletype <- sg_metadata$group[match(mdf$libID, sg_metadata$Sample.ID)]
mdf$kids <- fourab_mapping_table$KIDS[match(mdf$Uniref, fourab_mapping_table$UniRef90)]
mdf$geneids <- fourab_mapping_table$GeneName[match(mdf$Uniref, fourab_mapping_table$UniRef90)]
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
  ) + ylab('Normalised counts') + ggtitle(geneid) + xlab("") +  labs(fill="Group") 
#Compute lower and upper whiskers
ylim1 = boxplot.stats(req_mdf1$value)$stats[c(1, 5)]
limonaxis = ylim1 * valgiven
p1 <- p +  coord_cartesian(ylim = limonaxis) + scale_y_continuous(breaks = seq(0, limonaxis[2], by = byval))
print(p1)}
wantedabfH <- filter(req_mdf, geneids=='abfH')
makegeneboxplot(wantedabfH, 4.15, 'abfH',2)
abfhglm <- summary (glm(wantedabfH$value~ wantedabfH$Age, data=wantedabfH))
wantedabfT <- filter(req_mdf, geneids=='abfT')
makegeneboxplot(wantedabfT, 1.35, 'abfT',1)
abfTglm <- summary (glm(wantedabfT$value~ wantedabfT$Age, data=wantedabfT))
wantedabfD <- filter(req_mdf, geneids=='abfD')
makegeneboxplot(wantedabfD, 1.35, 'abfD',1)
abfdglm <- summary (glm(wantedabfD$value~ wantedabfD$Age, data=wantedabfD))
allpval <- data.frame(geneids=character(), pvalues=numeric()) %>% add_row(geneids='abfH', pvalues=abfhglm$coefficients[8])%>% add_row(geneids='abfT', pvalues=abfTglm$coefficients[8])%>% add_row(geneids='abfD', pvalues=abfdglm$coefficients[8])
write.table(allpval, file='../results/pval_4ab_pathway.csv', append=TRUE, sep=',', col.names = FALSE)

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
  ) + ylab('Normalised counts') + ggtitle("4 aminobutyrate pathway") + xlab("") +  labs(fill="Group") 

print(agegroupplot)

ylim1 <-boxplot.stats(req_mdf1$value)$stats[c(1, 5)]
limonaxis = ylim1 * 1.8
agegroupplot1 <- agegroupplot +  coord_cartesian(ylim = limonaxis)
print(agegroupplot1)

allglm <-summary (glm(req_mdf1$value~ req_mdf1$Age, data=req_mdf1))
allglm$coefficients[8]

