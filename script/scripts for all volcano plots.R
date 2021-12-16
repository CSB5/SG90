#Scripts for all volcano plots

library(ggplot2)
library(scales)
library(ggrepel)

#Figure 1E
Species <- read.csv("../data/Age_associated_results.csv", header = TRUE)
#Species$Significant <- ifelse(Species$padj < 0.05, "FDR < 0.05", "Not Sig")
ggplot (Species, aes(x = Estimate, y = -log10(padj))) +
  geom_point(aes(color = Significant, size = 14)) +
  scale_color_manual(values = c("red", "green", "grey")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text_repel(
    data = subset(Species, highlight > 0),
    aes(label = Species),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(1.0, "lines")
  )

#Figure 4A
Species <- read.csv("../data/FBGlu_species_results_estimate_filter4.csv", header = TRUE)
Species$Significant <- ifelse(Species$padj < 0.05, "FDR < 0.05", "Not Sig")
ggplot (Species, aes(x = Estimate, y = -log10(padj))) +
  geom_point(aes(color = Significant, size = 10)) +
  geom_vline(
    xintercept = c(0),
    col = "black",
    linetype = "solid",
    size = 0.5) +
  geom_hline(
    yintercept = c(-log10(0.05),-log10(0.1)),
    col = c("black"),
    linetype = ("dotted"),
    size = 1)+
  scale_color_manual(values = c("red", "grey")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text_repel(
    data = subset(Species, highlight > 3),
    aes(label = Species),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(1.0, "lines")
  )

#Figure S5A
Species <- read.csv("../data/FBGlu_no_age_species_results_estimate_filter4.csv", header = TRUE)
Species$Significant <- ifelse(Species$padj < 0.05, "FDR < 0.05", "Not Sig")
ggplot (Species, aes(x = Estimate, y = -log10(padj))) +
  geom_point(aes(color = Significant, size = 10)) +
  geom_vline(
    xintercept = c(0),
    col = "black",
    linetype = "solid",
    size = 0.5) +
  geom_hline(
    yintercept = c(-log10(0.05),-log10(0.1)),
    col = c("black"),
    linetype = ("dotted"),
    size = 1)+
  scale_color_manual(values = c("red", "grey")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text_repel(
    data = subset(Species, highlight > 3),
    aes(label = Species),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(1.0, "lines")
  )

#Figure 4B
Species <- read.csv("../data/TCHOL_with_age_species_results.csv", header = TRUE)
Species$Significant <- ifelse(Species$padj < 0.05, "FDR < 0.05", "Not Sig")
ggplot (Species, aes(x = Estimate, y = -log10(padj))) +
  geom_point(aes(color = Significant, size = 10)) +
  geom_vline(
    xintercept = c(0),
    col = "black",
    linetype = "solid",
    size = 0.5) +
  geom_hline(
    yintercept = c(-log10(0.05),-log10(0.1)),
    col = c("black"),
    linetype = ("dotted"),
    size = 1)+
  scale_color_manual(values = c("red", "grey")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text_repel(
    data = subset(Species, highlight > 3),
    aes(label = Species),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(1.0, "lines")
  )

#Figure S5B
Species <- read.csv("../data/LDL_with_age_species_results.csv", header = TRUE)
Species$Significant <- ifelse(Species$padj < 0.05, "FDR < 0.05", "Not Sig")
ggplot (Species, aes(x = Estimate, y = -log10(padj))) +
  geom_point(aes(color = Significant, size = 10)) +
  geom_vline(
    xintercept = c(0),
    col = "black",
    linetype = "solid",
    size = 0.5) +
  geom_hline(
    yintercept = c(-log10(0.05),-log10(0.1)),
    col = c("black"),
    linetype = ("dotted"),
    size = 1)+
  scale_color_manual(values = c("red", "grey")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text_repel(
    data = subset(Species, highlight > 3),
    aes(label = Species),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(1.0, "lines")
  )

Species <- read.csv("../data/HDL_with_age_species_results.csv", header = TRUE)
Species$Significant <- ifelse(Species$padj < 0.05, "FDR < 0.05", "Not Sig")
ggplot (Species, aes(x = Estimate, y = -log10(padj))) +
  geom_point(aes(color = Significant, size = 10)) +
  geom_vline(
    xintercept = c(0),
    col = "black",
    linetype = "solid",
    size = 0.5) +
  geom_hline(
    yintercept = c(-log10(0.05),-log10(0.1)),
    col = c("black"),
    linetype = ("dotted"),
    size = 1)+
  scale_color_manual(values = c("grey")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text_repel(
    data = subset(Species, highlight > 3),
    aes(label = Species),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(1.0, "lines")
  )

Species <- read.csv("../data/Tgly_with_age_species_results.csv", header = TRUE)
Species$Significant <- ifelse(Species$padj < 0.05, "FDR < 0.05", "Not Sig")
ggplot (Species, aes(x = Estimate, y = -log10(padj))) +
  geom_point(aes(color = Significant, size = 10)) +
  geom_vline(
    xintercept = c(0),
    col = "black",
    linetype = "solid",
    size = 0.5) +
  geom_hline(
    yintercept = c(-log10(0.05),-log10(0.1)),
    col = c("black"),
    linetype = ("dotted"),
    size = 1)+
  scale_color_manual(values = c("grey")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text_repel(
    data = subset(Species, highlight > 3),
    aes(label = Species),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(1.0, "lines")
  )

#Figure 4C
Species <- read.csv("../data/HSCRP_Age_gender_ALT_AST_others.csv", header = TRUE)
Species$Significant <- ifelse(Species$padj < 0.05, "FDR < 0.05", "Not Sig")
ggplot (Species, aes(x = Estimate, y = -log10(padj))) +
  geom_point(aes(color = Significant, size = 10)) +
  geom_vline(
    xintercept = c(0),
    col = "black",
    linetype = "solid",
    size = 0.5) +
  geom_hline(
    yintercept = c(-log10(0.05),-log10(0.1)),
    col = c("black"),
    linetype = ("dotted"),
    size = 1)+
  scale_color_manual(values = c("red", "grey")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text_repel(
    data = subset(Species, highlight > 3),
    aes(label = Species),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(1.0, "lines")
  )

#Figure 4D
Species <- read.csv("../data/AST_Age_gender_HSCRP_ALT_others.csv", header = TRUE)
Species$Significant <- ifelse(Species$padj < 0.05, "FDR < 0.05", "Not Sig")
ggplot (Species, aes(x = Estimate, y = -log10(padj))) +
  geom_point(aes(color = Significant, size = 10)) +
  geom_vline(
    xintercept = c(0),
    col = "black",
    linetype = "solid",
    size = 0.5) +
  geom_hline(
    yintercept = c(-log10(0.05),-log10(0.1)),
    col = c("black"),
    linetype = ("dotted"),
    size = 1)+
  scale_color_manual(values = c("red", "grey")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text_repel(
    data = subset(Species, highlight > 3),
    aes(label = Species),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(1.0, "lines")
  )

#Figure S5C
Species <- read.csv("../data/ALT_Age_gender_HSCRP_AST_others.csv", header = TRUE)
Species$Significant <- ifelse(Species$padj < 0.05, "FDR < 0.05", "Not Sig")
ggplot (Species, aes(x = Estimate, y = -log10(padj))) +
  geom_point(aes(color = Significant, size = 10)) +
  geom_vline(
    xintercept = c(0),
    col = "black",
    linetype = "solid",
    size = 0.5) +
  geom_hline(
    yintercept = c(-log10(0.05),-log10(0.1)),
    col = c("black"),
    linetype = ("dotted"),
    size = 1)+
  scale_color_manual(values = c("grey")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text_repel(
    data = subset(Species, highlight > 3),
    aes(label = Species),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(1.0, "lines")
  )

#Figure S6A
Species <- read.csv("../data/GaitT1_Age_gender_others_estimate30.csv", header = TRUE)
Species$Significant <- ifelse(Species$padj < 0.05, "FDR < 0.05", "Not Sig")
ggplot (Species, aes(x = Estimate, y = -log10(padj))) +
  geom_point(aes(color = Significant, size = 10)) +
  geom_vline(
    xintercept = c(0),
    col = "black",
    linetype = "solid",
    size = 0.5) +
  geom_hline(
    yintercept = c(-log10(0.05),-log10(0.1)),
    col = c("black"),
    linetype = ("dotted"),
    size = 1)+
  scale_color_manual(values = c("red", "grey")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text_repel(
    data = subset(Species, highlight > 3),
    aes(label = Species),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(1.0, "lines")
  )

#Figure S6B
Species <- read.csv("../data/HandgripR1_Age_gender_others.csv", header = TRUE)
Species$Significant <- ifelse(Species$padj < 0.05, "FDR < 0.05", "Not Sig")
ggplot (Species, aes(x = Estimate, y = -log10(padj))) +
  geom_point(aes(color = Significant, size = 10)) +
  geom_vline(
    xintercept = c(0),
    col = "black",
    linetype = "solid",
    size = 0.5) +
  geom_hline(
    yintercept = c(-log10(0.05),-log10(0.1)),
    col = c("black"),
    linetype = ("dotted"),
    size = 1)+
  scale_color_manual(values = c("red", "grey")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text_repel(
    data = subset(Species, highlight > 3),
    aes(label = Species),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(1.0, "lines")
  )


#Figure S6C
Species <- read.csv("../data/HandgripL1_Age_gender_others.csv", header = TRUE)
Species$Significant <- ifelse(Species$padj < 0.05, "FDR < 0.05", "Not Sig")
ggplot (Species, aes(x = Estimate, y = -log10(padj))) +
  geom_point(aes(color = Significant, size = 10)) +
  geom_vline(
    xintercept = c(0),
    col = "black",
    linetype = "solid",
    size = 0.5) +
  geom_hline(
    yintercept = c(-log10(0.05),-log10(0.1)),
    col = c("black"),
    linetype = ("dotted"),
    size = 1)+
  scale_color_manual(values = c("grey")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text_repel(
    data = subset(Species, highlight > 3),
    aes(label = Species),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(1.0, "lines")
  )


