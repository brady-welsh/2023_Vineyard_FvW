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


##################################
## TAXANOMIC BAR PLOTS (FAMILY) ##
##################################

#Figure 1

#Collapse to family level
ps.phyllo.family <- tax_glom(ps.phyF.prevF.rar, taxrank = "Family", NArm = FALSE)

#Extract top 20 most abundant family names, bind to ps sampledata
top20families = names(sort(taxa_sums(ps.phyllo.family), TRUE)[1:22])
taxtab20 = cbind(tax_table(ps.phyllo.family), family_20 = NA)
taxtab20[top20families, "family_20"] <- as(tax_table(ps.phyllo.family)
                                           [top20families, "Family"], "character")
tax_table(ps.phyllo.family) <- tax_table(taxtab20)

ps.phyllo.family.ra <- transform_sample_counts(ps.phyllo.family, function(x) 100 * x/sum(x))

#Melt into a dataframe
pd.phyllo.family <- psmelt(ps.phyllo.family.ra)

#Replace NA with 'other', for plotting purposes
pd.phyllo.family <- arrange(pd.phyllo.family, stage)
pd.phyllo.family$family_20[is.na(pd.phyllo.family$family_20)] <- c("Other")

#Relative abundance of top 20 families?
mean(sample_sums(prune_taxa(top20families, ps.phyllo.family.ra)))


library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 23
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(nb.cols)

mycolors2 <- colorRampPalette(Blue2DarkRed18Steps)(nb.cols)


#Plot em

ggplot(pd.phyllo.family, aes(x = (stage), y = (Abundance),
                                  fill = fct_reorder(family_20, -Abundance))) +
  geom_bar(width = 0.9, stat = "identity", position="fill") +
  facet_nested(~vineyard + sample.type, scales = "free", space = "free", ) +
  labs(x = "", y = "Relative abundance") +
  theme(
    axis.text.y = element_text(size=14, face = 'bold'),
    axis.title.y = element_text(size=14, face = 'bold'),
    axis.ticks.y = element_line(size = 1),
    axis.ticks.x = element_line(size = 1),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_textbox_simple(
      padding = margin(5, 0, 5, 0),
      margin = margin(5, 5, 5, 5),
      size = 10,
      face = "bold",
      halign = 0.5,
      fill = "white",
      box.color = "grey",
      linewidth = 1.5,
      linetype = "solid",),
    panel.background = element_blank()
  ) +
  scale_fill_manual(values = mycolors)+
  scale_y_continuous(labels = scales::percent, expand = c(0, 0))

#Save!
ggsave(filename = "R_output/Figure4.png", 
       width = 15, height = 10, dpi = 300)


###########################
## Alpha Diversity Plots ##
###########################

# Generate Alpha statistics:
ps.phyF.prevF.rar.div <- alpha(ps.phyF.prevF.rar, index = c('observed', 'shannon', 'simpson'))

#Append metadata to alpha diversity values
ps.phyF.prevF.rar.meta <- meta(ps.phyF.prevF.rar)
ps.phyF.prevF.rar.meta$name <- rownames(ps.phyF.prevF.rar.meta)
ps.phyF.prevF.rar.div$name <- rownames(ps.phyF.prevF.rar.div)
ps.phyF.prevF.rar.div.df <- merge(ps.phyF.prevF.rar.div, ps.phyF.prevF.rar.meta, by = "name")


#Split dataframe by vineyard
Org.alpha.div.df <- filter(ps.phyF.prevF.rar.div.df, vineyard == "Organic")
Bio.alpha.div.df <- filter(ps.phyF.prevF.rar.div.df, vineyard == "Biodynamic")
Con.alpha.div.df <- filter(ps.phyF.prevF.rar.div.df, vineyard == "Conventional")

##GGPLOT!
colours.con <- c("#FF7C3A", "#FF9D6B", "#ffbE9D")
colours.bio <- c("#94D852", "#AFE27D", "#CAECA9")
colours.org <- c("#5293D8", "#7DAEE2", "#A9C9EC")

#Reorder to match GI tract order
Con.alpha.div.df.ordered <- Con.alpha.div.df
Con.alpha.div.df.ordered$sampleshort <- factor(Con.alpha.div.df.ordered$stage,
                                               levels = c('T0', 'T1', 'T2', 'T3', 'T4'))

Bio.alpha.div.df.ordered <- Bio.alpha.div.df
Bio.alpha.div.df.ordered$sampleshort <- factor(Bio.alpha.div.df.ordered$stage,
                                               levels = c('T0', 'T1', 'T2', 'T3', 'T4'))

Org.alpha.div.df.ordered <- Org.alpha.div.df
Org.alpha.div.df.ordered$sampleshort <- factor(Org.alpha.div.df.ordered$stage,
                                               levels = c('T0', 'T1', 'T2', 'T3', 'T4'))



Con.shannon <- ggplot(Con.alpha.div.df.ordered,
                       aes(x=stage, y=diversity_shannon, fill=sample.type))
Bio.shannon <- ggplot(Bio.alpha.div.df.ordered,
                        aes(x=stage, y=diversity_shannon, fill=sample.type))
Org.shannon <- ggplot(Org.alpha.div.df.ordered,
                       aes(x=stage, y=diversity_shannon, fill=sample.type))

Con.richness <- ggplot(Con.alpha.div.df.ordered,
                      aes(x=stage, y=observed, fill=sample.type))
Bio.richness <- ggplot(Bio.alpha.div.df.ordered,
                      aes(x=stage, y=observed, fill=sample.type))
Org.richness <- ggplot(Org.alpha.div.df.ordered,
                      aes(x=stage, y=observed, fill=sample.type))



gg.Con.shannon <- Con.shannon +
  #Boxplot
  facet_wrap(~stage, scales = "free_x", ncol = 5) +
  geom_boxplot() +
  #Jitter, size, colour
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, aes(colour=sample.type)) +
  #Custom manual colours
  scale_colour_manual(values=colours.con) +
  scale_fill_manual(values=colours.con) +
  #Tick labels
  theme(axis.text.y = element_text(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size=15, face="bold"),
        axis.title.y = element_text(size=15, face="bold"),
        axis.line = element_line(colour = "black"),
        #Background panel
        panel.background = element_rect(fill = "White"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        #Legend
        legend.position = 'none',
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white")
        ) +
  #Axis labels
  ylim(1,2.5) +
  labs(x = "") +
  labs(y = "Shannon Diversity")

gg.Bio.shannon <- Bio.shannon +
  #Boxplot
  facet_wrap(~stage, scales = "free_x", ncol = 5) +
  geom_boxplot() +
  #Jitter, size, colour
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, aes(colour=sample.type)) +
  #Custom manual colours
  scale_colour_manual(values=colours.bio) +
  scale_fill_manual(values=colours.bio) +
  #Tick labels
  theme(axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size=15, face="bold"),
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    axis.line = element_line(colour = "black"),
    #Background panel
    panel.background = element_rect(fill = "White"),
    panel.grid.major = element_line(colour = "white"),
    panel.grid.minor = element_line(colour = "white"),
    strip.background = element_blank(),
    strip.text = element_blank(),
    #Legend
    legend.position = 'none',
    legend.title = element_blank(),
    legend.key = element_rect(fill = "white")
  ) +
  #Axis labels
  ylim(1,2.5) +
  labs(x = "")

gg.Org.shannon <- Org.shannon +
  #Boxplot
  facet_wrap(~stage, scales = "free_x", ncol = 5) +
  geom_boxplot() +
  #Jitter, size, colour
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, aes(colour=sample.type)) +
  #Custom manual colours
  scale_colour_manual(values=colours.org) +
  scale_fill_manual(values=colours.org) +
  #Tick labels
  theme(axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size=15, face="bold"),
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    axis.line = element_line(colour = "black"),
    #Background panel
    panel.background = element_rect(fill = "White"),
    panel.grid.major = element_line(colour = "white"),
    panel.grid.minor = element_line(colour = "white"),
    strip.background = element_blank(),
    strip.text = element_blank(),
    #Legend
    legend.position = 'none',
    legend.title = element_blank(),
    legend.key = element_rect(fill = "white")
  ) +
  #Axis labels
  ylim(1,2.5) +
  labs(x = "Time Point")

gg.Bio.shannon <- gg.Bio.shannon + ggtitle("Biodynamic") + theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
gg.Con.shannon <- gg.Con.shannon + ggtitle("Conventional") + theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
gg.Org.shannon <- gg.Org.shannon + ggtitle("Organic") + theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))


######

gg.Con.richness <- Con.richness +
  #Boxplot
  facet_wrap(~stage, scales = "free_x", ncol = 5) +
  geom_boxplot() +
  #Jitter, size, colour
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, aes(colour=sample.type)) +
  #Custom manual colours
  scale_colour_manual(values=colours.con) +
  scale_fill_manual(values=colours.con) +
  #Tick labels
  theme(axis.text.y = element_text(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size=15, face="bold"),
        axis.title.y = element_text(size=15, face="bold"),
        axis.line = element_line(colour = "black"),
        #Background panel
        panel.background = element_rect(fill = "White"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        #Legend
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white")
  ) +
  #Axis labels
  ylim(5,30) +
  labs(x = "") +
  labs(y = "ASV Richness")

gg.Bio.richness <- Bio.richness +
  #Boxplot
  facet_wrap(~stage, scales = "free_x", ncol = 5) +
  geom_boxplot() +
  #Jitter, size, colour
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, aes(colour=sample.type)) +
  #Custom manual colours
  scale_colour_manual(values=colours.bio) +
  scale_fill_manual(values=colours.bio) +
  #Tick labels
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size=15, face="bold"),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.line = element_line(colour = "black"),
        #Background panel
        panel.background = element_rect(fill = "White"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        #Legend
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white")
  ) +
  #Axis labels
  ylim(5,30) +
  labs(x = "")

gg.Org.richness <- Org.richness +
  #Boxplot
  facet_wrap(~stage, scales = "free_x", ncol = 5) +
  geom_boxplot() +
  #Jitter, size, colour
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, aes(colour=sample.type)) +
  #Custom manual colours
  scale_colour_manual(values=colours.org) +
  scale_fill_manual(values=colours.org) +
  #Tick labels
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size=15, face="bold"),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.line = element_line(colour = "black"),
        #Background panel
        panel.background = element_rect(fill = "White"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        #Legend
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white")
  ) +
  #Axis labels
  ylim(5,30) +
  labs(x = "Time Point")


#Arrange plots together
shannon.plots <- ggarrange(gg.Con.shannon, gg.Org.shannon, gg.Bio.shannon,
          nrow = 1, ncol = 3, common.legend = FALSE)
richness.plots <- ggarrange(gg.Con.richness, gg.Org.richness, gg.Bio.richness,
                           nrow = 1, ncol = 3, common.legend = FALSE)
ggarrange(shannon.plots, richness.plots, nrow = 2, ncol = 1, common.legend = FALSE, labels = c("A)", "B)"))

#Save
ggsave(filename = "R_output/Figure5.png", 
       width = 15, height = 10, dpi = 300)



#########################
## Beta Diversity PCoA ##
#########################

# Seperate vineyards
ps.rar.wi <- subset_samples(ps.phyF.prevF.rar, sample.type != "Inoculated")
ps.rar.in <- subset_samples(ps.phyF.prevF.rar, sample.type != "Wild")

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

ggsave(filename = "R_output/Figure6.png", 
       width = 15, height = 15, dpi = 300)

# Beta diversity stats
library(vegan)

ps.rar.wi <- subset_samples(ps.phyF.prevF.rar, sample.type != "Inoculated")
ps.rar.in <- subset_samples(ps.phyF.prevF.rar, sample.type != "Wild")

metadata.wi <- as(sample_data(ps.rar.wi), "data.frame")
metadata.in <- as(sample_data(ps.rar.in), "data.frame")

adonis2(distance(ps.rar.wi, method="bray") ~ vineyard, data = metadata.wi)
adonis2(distance(ps.rar.in, method="bray") ~ vineyard, data = metadata.in)

