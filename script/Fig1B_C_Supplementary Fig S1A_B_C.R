require(foreign)
require(ggplot2)
require(MASS)
require(Hmisc)
require(reshape2)

#Dataset for Richness
Dataset = read.csv (
  "../data/Species richness_S.csv"
)
mydata = read.csv (
  "../data/Species richness_S.csv",
  header = TRUE,
  sep = ","
)


#Plot Richness by decade
#violin
richnessplot <-
  ggplot(Dataset, aes(
    x = decade,
    y = value,
    outline = FALSE,
    fill = decade
  )) + geom_violin()  + geom_boxplot(width = 0.2, outlier.shape = NA) + coord_cartesian(ylim = c(0.0, 70.0)) + ggtitle("Richness") + theme_bw() + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  ) + labs(x = "Group",
           y = "Value",
           face = "bold",
           size = 14)
p + theme(
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
))
#Ordinal logistic test by decade
Richness <- polr(decade ~ value, Hess = TRUE, data = mydata)
(ctable <- coef(summary(Richness)))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
(ctable <- cbind(ctable, "p value" = p))


#Dataset for J_evenness
Dataset = read.csv (
  "../data/Pielou evenness_J.csv"
)
mydata = read.csv (
  "../data/Pielou evenness_J.csv",
  header = TRUE,
  sep = ","
)
#Plot J_evenness by decade
#violin
evennessplot <-
  ggplot(Dataset, aes(
    x = decade,
    y = value,
    outline = FALSE,
    fill = decade
  )) + geom_violin()  + geom_boxplot(width = 0.2, outlier.shape = NA) + coord_cartesian(ylim = c(0.0, 1.0)) + ggtitle("J_Evenness") + theme_bw() + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  ) + labs(x = "Group",
           y = "Value",
           face = "bold",
           size = 14) + theme(
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
           ))
print(evennessplot)
#Ordinal logistic test by decade
J_evenness <- polr(decade ~ value, Hess = TRUE, data = mydata)
(ctable <- coef(summary(J_evenness)))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
(ctable <- cbind(ctable, "p value" = p))


#Dataset for MDS2
mydata_MDS = read.csv (
  "../data/MDS.csv",
  header = TRUE,
  sep = ","
)
Dataset = read.csv ("../data/MDS.csv")
#Plot MDS2 by decade
#violin
mdsplot <-
  ggplot(Dataset, aes(
    x = decade10,
    y = MDS2,
    outline = FALSE,
    fill = decade10
  )) + geom_violin() + geom_boxplot(width = 0.2, outlier.shape = NA) + coord_cartesian(ylim = c(-0.12, 0.11)) + ggtitle("MDS2") + theme_bw() + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  ) + labs(x = "Group",
           y = "Value",
           face = "bold",
           size = 14) + theme(
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
           ))
print(mdsplot)
#Ordinal logistic test by decade
MDS <- polr(decade10 ~ MDS2, Hess = TRUE, data = mydata_MDS)
(ctable <- coef(summary(MDS)))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
(ctable <- cbind(ctable, "p value" = p))


#Dataset for simpson
Dataset = read.csv (
  "../data/diversity_values_simpson.csv"
)
mydata = read.csv (
  "../data/diversity_values_simpson.csv",
  header = TRUE,
  sep = ","
)
#Plot simpson by decade
#violin
simpsonplot <-
  ggplot(Dataset, aes(
    x = decade,
    y = simpson,
    outline = FALSE,
    fill = decade
  )) + geom_violin() + geom_boxplot(width = 0.2, outlier.shape = NA) + coord_cartesian(ylim = c(0.25, 1.0)) + ggtitle("Simpson diversity") + theme_bw() + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  ) + labs(x = "Group",
           y = "Value",
           face = "bold",
           size = 14) + theme(
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
           ))
print(simpsonplot)
#Ordinal logistic test by decade
simpson <- polr(decade ~ simpson, Hess = TRUE, data = mydata)
(ctable <- coef(summary(simpson)))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
(ctable <- cbind(ctable, "p value" = p))


#Dataset for shannon
Dataset = read.csv (
  "../data/diversity_values_shannon.csv"
)
mydata = read.csv (
  "../data/diversity_values_shannon.csv",
  header = TRUE,
  sep = ","
)
#Plot shannon by decade
#violin
shannonplot <-
  ggplot(Dataset, aes(
    x = decade,
    y = value,
    outline = FALSE,
    fill = decade
  )) + geom_violin()  + geom_boxplot(width = 0.2, outlier.shape = NA) + coord_cartesian(ylim = c(0.2, 3.5)) + ggtitle("Shannon diversity") + theme_bw() + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  ) + labs(x = "Group",
           y = "Value",
           face = "bold",
           size = 14) + theme(
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
           ))
print(shannonplot)
shannon <- polr(decade ~ value, Hess = TRUE, data = mydata)
(ctable <- coef(summary(shannon)))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
(ctable <- cbind(ctable, "p value" = p))
