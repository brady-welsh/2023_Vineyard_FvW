library(phyloseq)
library(qiime2R)
library(ggplot2)
library(scales)
library(microbiome)
library(tidyverse)
library(ggmap)
library(ggrepel)
library(ggpubr)
library(VennDiagram)
library(osmdata)
library(RColorBrewer)
library(ggh4x)
library(ggtext)
library(ggVennDiagram)


#ps.tree <- qza_to_phyloseq(tree = "Vineyard-ITS-rooted-tree.qza")

#ps.phyF.tree <- merge_phyloseq(ps.phyF.prevF.rar, ps.tree)

# Seperate vineyards
ps.rar.wi <- subset_samples(ps.phyF.tree, sample.type != "Inoculated")
ps.rar.in <- subset_samples(ps.phyF.tree, sample.type != "Wild")

#Calculate beta diversities
#Unweighted unifrac
ord.bray.wi <- ordinate(ps.rar.wi, "PCoA", "bray", weighted=F)
ord.bray.in <- ordinate(ps.rar.in, "PCoA", "bray", weighted=F)


#Axes 1/2
bray.1.2.wi <- plot_ordination(ps.rar.wi, 
                                       ord.bray.wi, color="vineyard", 
                                       axes = c(1, 2))
bray.1.2.in <- plot_ordination(ps.rar.in, 
                                       ord.bray.in, color="vineyard", 
                                       axes = c(1, 2))


#Axes 1/3
bray.1.3.wi <- plot_ordination(ps.rar.wi, 
                                        ord.bray.wi, color="vineyard", 
                                        axes = c(1, 3))
bray.1.3.in <- plot_ordination(ps.rar.in, 
                                        ord.bray.in, color="vineyard", 
                                        axes = c(1, 3))


#set a theme
theme_pcoa <- theme(axis.text.x = element_text(face="bold", size=16), 
                    axis.text.y = element_text(face="bold", size=16),
                    axis.title.x = element_text(size=20, face="bold"),
                    axis.title.y = element_text(size=20, face="bold"),
                    axis.line = element_line(colour = "black"),
                    #Background panel
                    panel.background = element_rect(fill = "White"),
                    panel.grid.major = element_line(colour = "white"), 
                    panel.grid.minor = element_line(colour = "white"),
                    #Legend
                    legend.title = element_blank(),
                    legend.text = element_text(size=16),
                    legend.key = element_rect(fill = "white", color = NA),
                    legend.key.size = unit(2.5, "line"))


#Plot!
#Wild
wi.br.1.2 <- bray.1.2.wi +
  scale_x_reverse() +
  geom_point(size=7.5) +
  geom_path() +
  geom_text(aes(label = stage), color = "black", fontface = "bold") +
  scale_colour_manual(values=colours) +
  theme_pcoa +
  theme(legend.position = "none")

wi.br.1.3 <- bray.1.3.wi +
  scale_x_reverse() +
  geom_point(size=7.5) +
  geom_path() +
  geom_text(aes(label = stage), color = "black", fontface = "bold") +
  scale_colour_manual(values=colours) +
  theme_pcoa +
  theme(legend.position = "none")

#Inoculated
in.br.1.2 <- bray.1.2.in +
  scale_x_reverse() +
  geom_point(size=7.5) +
  geom_path() +
  geom_text(aes(label = stage), color = "black", fontface = "bold") +
  scale_colour_manual(values=colours) +
  theme_pcoa +
  theme(legend.position = "none")

in.br.1.3 <- bray.1.3.in +
  scale_x_reverse() +
  geom_point(size=7.5) +
  geom_path() +
  geom_text(aes(label = stage), color = "black", fontface = "bold") +
  scale_colour_manual(values=colours) +
  theme_pcoa +
  theme(legend.position = "none")


#merge!
ggarrange(wi.br.1.2, wi.br.1.3, in.br.1.2, in.br.1.3,
          nrow = 2, ncol = 2, common.legend = TRUE, labels = c("A)", "B)", "C)", "D)"), 
          font.label = list(size=18, face="bold", color="black"))

ggsave(filename = "R_output/PCoA.png", 
       width = 15, height = 15, dpi = 300)
