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

#Figure 2

#Collapse to family level
ps.family <- tax_glom(ps.prevF.rar, taxrank = "Family", NArm = FALSE)

#Extract top 20 most abundant family names, bind to ps sampledata
top20families = names(sort(taxa_sums(ps.family), TRUE)[1:24])
taxtab20 = cbind(tax_table(ps.family), family_20 = NA)
taxtab20[top20families, "family_20"] <- as(tax_table(ps.family)
                                           [top20families, "Family"], "character")
tax_table(ps.family) <- tax_table(taxtab20)

ps.family.ra <- transform_sample_counts(ps.family, function(x) 100 * x/sum(x))

#Melt into a dataframe
pd.family <- psmelt(ps.family.ra)

#Replace NA with 'other', for plotting purposes
pd.family <- arrange(pd.family, stage)
pd.family$family_20[is.na(pd.family$family_20)] <- c("Other")

#Relative abundance of top 20 families?
mean(sample_sums(prune_taxa(top20families, ps.family.ra)))


library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 24
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(nb.cols)


#Plot em

ggplot(pd.family, aes(x = (stage), y = (Abundance),
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
ggsave(filename = "R_output/Figure2.png", 
       width = 15, height = 10, dpi = 300)


###############################
## Venn Diagram Phyllo + All ##
###############################

# Figure 3

#Pull out the first berry samples
ps.prevF.rar.berry <- subset_samples(ps.prevF.rar, sample.type == "Berry")

ps.Bio.berry <- subset_samples(ps.prevF.rar.berry, vineyard == "Biodynamic")
ps.Bio.berry.table <- otu_table(ps.Bio.berry)

ps.Con.berry <- subset_samples(ps.prevF.rar.berry, vineyard == "Conventional")
ps.Con.berry.table <- otu_table(ps.Con.berry)

ps.Org.berry <- subset_samples(ps.prevF.rar.berry, vineyard == "Organic")
ps.Org.berry.table <- otu_table(ps.Org.berry)


#Merge tables
ps.berry.all <- merge_phyloseq(ps.Bio.berry, ps.Con.berry, ps.Org.berry)

#Remove ASVs with counts of 0.
ps.berry.all <- prune_taxa(taxa_sums(ps.berry.all) > 0, ps.berry.all)

#For each ASV (row), if abundance > 2, print ASV (rowname) to a vector
venn.Con.berry.0 <- rownames(ps.Con.berry.table[ apply(ps.Con.berry.table, MARGIN = 1,
                                                       function(x) any(x > 2))])

venn.Bio.berry.0 <- rownames(ps.Bio.berry.table[ apply(ps.Bio.berry.table, MARGIN = 1,
                                                       function(x) any(x > 2))])

venn.Org.berry.0 <- rownames(ps.Org.berry.table[ apply(ps.Org.berry.table, MARGIN = 1,
                                                       function(x) any(x > 2))])

berry <- list(Conventional = venn.Con.berry.0, Biodynamic = venn.Bio.berry.0, Organic = venn.Org.berry.0)

#Plot em using venn.diagram (figure 5)
colours <- c("#FF7C3A","#94D852", "#5293D8")

venn.diagram(berry, cex = 1.5, cat.cex = 2.5, print.mode = c("raw","percent"), fill = colours, inverted = TRUE,
             imagetype = "tiff", filename = "R_output/Figure3.png", cat.pos = c(340,20,180), lty="blank",
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.dist = c(0.075, 0.075, 0.055),
             cat.fontfamily = "sans",
             rotation = 1)
