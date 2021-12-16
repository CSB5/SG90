library(ggplot2)
library(vegan)
metadata = read.table(
  "../data/PcoA_metadata.txt",
  header = T,
  stringsAsFactors = F,
  sep = '\t',
  row.names = 1
)
mydata = read.table(
  "../data/filtered_taxa.txt",
  header = T,
  stringsAsFactors = F,
  sep = '\t',
  row.names = 1
)
meta_merged <-
  data.frame(metadata$LibraryID, metadata$group, metadata$Position)
rownames(meta_merged) <- meta_merged$metadata.LibraryID
distmetric <- vegdist(t(mydata), method = "bray")
cmds <- cmdscale(distmetric, k = 3, eig = TRUE)
eigen <- cmds$eig / sum(cmds$eig) * 100
cmdvals <- as.data.frame(cmds$points)
dat_merged_all <- data.frame(cmdvals, meta_merged)

s2 <- ggplot(dat_merged_all,
             aes(x = V1,
                 y = V2)) +
  geom_point(
    data = subset(dat_merged_all, metadata.group == '2Elderly'),
    aes(color = metadata.Position),
    size = 4,
  ) + geom_point(
    data = subset(dat_merged_all,
                  metadata.group == '1Young'),
    aes(color = metadata.Position),
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
    colour = 'Age group'
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
  )
print(s2)