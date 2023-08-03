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
library(lme4)
library(multcomp)
library(dplyr)

# Don't forget to set working directory. (2023_Vineyard_FvW/ITS_Fungi_ITS1-ITS2R/data/)

#################
## IMPORT DATA ##
#################

ps<-qza_to_phyloseq(
  features="QIIME2_output/Vineyard-ITS-table_filtered.qza",
  taxonomy="QIIME2_output/Vineyard-ITS-UNITE.qza",
  metadata="QIIME2_output/metadata/Vineyard-ITS-metadata_nocontrols.tsv"
)


###########################
## PREVELANCE FILTERTING ##
###########################

##Code modified from Callahan et al. 2016: https://doi.org/10.12688/f1000research.8986.2

#What phyla contain the most ASVs?
sort(table(tax_table(ps)[,"Phylum"]), exclude=NULL)

#Compute prevalence of each feature, store as data.frame
prevdf <- apply(X = otu_table(ps),
                MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
#Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(Prevalence = prevdf,
                     TotalAbundance = taxa_sums(ps),
                     tax_table(ps))

#OK, how many ASVs are only found in 1 replicate (we have replicates for each sample).
count(prevdf, vars = Prevalence <= 1)

#376/628 ASVs. We'll remove them. What % of samples is 1 replicate?
nsamples(ps)
#81 samples
1/81
#1/81 = 0.01234 x 100 = 1.234%

#Plot prevalance/abundance of each ASV, grouped by phylum, with the proposed prevalence 
# filter threshold line
prevdf1 <- subset(prevdf, Phylum %in% get_taxa_unique(ps,"Phylum"))

ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps), color = Phylum)) +
  geom_hline(yintercept=0.01234, alpha=0.5, linetype=2) + geom_point(size=2, alpha=0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

ggsave("R_output/AbundancePrevalence_Per_phylum.png", 
       width = 19, height = 10, dpi = 300)


#Execute prevalence filter
keepTaxa <- rownames(prevdf1)[(prevdf1$Prevalence >= 2)] 
ps.prevF <- prune_taxa(keepTaxa, ps)


#################
## RAREFACTION ##
#################

#Set seed, check the sample with the lowest read count
set.seed(1997)

sort(sample_sums(ps.prevF))

#Rarefy n number of reads where n = read count of samples with fewest reads
ps.prevF.rar <- rarefy_even_depth(ps.prevF, sample.size = 8220)


######################
## ALPHA STATISTICS ##
######################

# Generate Alpha statistics:
ps.rar.div <- alpha(ps.prevF.rar, index = c('observed', 'shannon', 'simpson'))

#Append metadata to alpha diversity values
ps.rar.meta <- meta(ps.prevF.rar)
ps.rar.meta$name <- rownames(ps.rar.meta)
ps.rar.div$name <- rownames(ps.rar.div)
ps.rar.div.df <- merge(ps.rar.div, ps.rar.meta, by = "name")


##############
## Figure 2 ##
##############

#Taxonomic Bar Plots (Family)

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


##############
## Figure 3 ##
##############

# Venn Diagram Phyllo + All

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

#Plot em using venn.diagram
colours <- c("#FF7C3A","#94D852", "#5293D8")

venn.diagram(berry, cex = 1.5, cat.cex = 2.5, print.mode = c("raw","percent"), fill = colours, inverted = TRUE,
             imagetype = "tiff", filename = "R_output/Figure3.png", cat.pos = c(340,20,180), lty="blank",
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.dist = c(0.075, 0.075, 0.055),
             cat.fontfamily = "sans",
             rotation = 1)


# Run PhyllosphereFilter.R

##################################
## PHYLLO PREVELANCE FILTERTING ##
##################################

##Code modified from Callahan et al. 2016: https://doi.org/10.12688/f1000research.8986.2

#What phyla contain the most ASVs?
sort(table(tax_table(ps.phyF)[,"Phylum"]), exclude=NULL)

#Compute prevalence of each feature, store as data.frame
prevdf.phyF <- apply(X = otu_table(ps.phyF),
                MARGIN = ifelse(taxa_are_rows(ps.phyF), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
#Add taxonomy and total read counts to this data.frame
prevdf.phyF <- data.frame(Prevalence = prevdf.phyF,
                     TotalAbundance = taxa_sums(ps.phyF),
                     tax_table(ps.phyF))

#OK, how many ASVs are only found in 1 replicate (we have replicates for each sample).
count(prevdf.phyF, vars = Prevalence <= 1)

#0/112 ASVs. We'll remove them. What % of samples is 1 replicate?
nsamples(ps.phyF)
#81 samples
1/81
#1/81 = 0.01234 x 100 = 1.234%

#Plot prevalance/abundance of each ASV, grouped by phylum, with the proposed prevalence 
# filter threshold line
prevdf1.phyF <- subset(prevdf.phyF, Phylum %in% get_taxa_unique(ps.phyF,"Phylum"))

ggplot(prevdf1.phyF, aes(TotalAbundance, Prevalence / nsamples(ps.phyF), color = Phylum)) +
  geom_hline(yintercept=0.01234, alpha=0.5, linetype=2) + geom_point(size=2, alpha=0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

ggsave("R_output/Phyllosphere_AbundancePrevalence_Per_phylum.png", 
       width = 19, height = 10, dpi = 300)


#Execute prevalence filter
keepTaxa <- rownames(prevdf1.phyF)[(prevdf1.phyF$Prevalence >= 2)] 
ps.phyF.prevF <- prune_taxa(keepTaxa, ps.phyF)


########################
## PHYLLO RAREFACTION ##
########################

#Set seed, check the sample with the lowest read count
set.seed(1997)

sort(sample_sums(ps.phyF.prevF))

#Rarefy n number of reads where n = read count of samples with fewest reads
ps.phyF.prevF.rar <- rarefy_even_depth(ps.phyF.prevF, sample.size = 653)


##############
## FIGURE 4 ##
##############

# Phyllo Taxonomic Bar Plot (Family)

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


##############
## FIGURE 5 ##
##############

# Alpha Diversity Plots

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


##################################
## LME MODEL FOR ALPHA DIVERSITY##
##################################

ps.inoc <- subset(ps.phyF.prevF.rar.div.df, sample.type != "Wild")
ps.wild <- subset(ps.phyF.prevF.rar.div.df, sample.type != "Inoculated")

## LME for ASV richness ##
# Inoculated Wines #
# Fit mixed-effects models for each vineyard separately
model_organic.in <- lm(observed ~ stage, data = subset(ps.inoc, vineyard == "Organic"))
model_biodynamic.in <- lm(observed ~ stage, data = subset(ps.inoc, vineyard == "Biodynamic"))
model_conventional.in <- lm(observed ~ stage, data = subset(ps.inoc, vineyard == "Conventional"))

# View the summaries of each model
summary(model_organic.in)
summary(model_biodynamic.in)
summary(model_conventional.in)

# Wild Wines
# Fit mixed-effects models for each vineyard separately
model_organic.wi <- lm(observed ~ stage, data = subset(ps.wild, vineyard == "Organic"))
model_biodynamic.wi <- lm(observed ~ stage, data = subset(ps.wild, vineyard == "Biodynamic"))
model_conventional.wi <- lm(observed ~ stage, data = subset(ps.wild, vineyard == "Conventional"))

# View the summaries of each model
summary(model_organic.wi)
summary(model_biodynamic.wi)
summary(model_conventional.wi)

model.list <- list(summary(model_organic.wi)$coefficients, summary(model_organic.in)$coefficients, 
                   summary(model_biodynamic.wi)$coefficients, summary(model_biodynamic.in)$coefficients, 
                   summary(model_conventional.wi)$coefficients, summary(model_conventional.in)$coefficients)


# Create a data frame from the extracted data
model.list.df <- data.frame(
  Coefficient = rownames(model.list),
  Estimate = coef_data[, "Estimate"],
  Std_Error = coef_data[, "Std. Error"],
  P_Value = coef_data[, "Pr(>|t|)"]
)

write.csv(model.list, file = "model-list_ASV-richness.csv", row.names = FALSE)


## LME for Shannon diversity ##
# Inoculated Wines #
# Fit mixed-effects models for each vineyard separately
model_organic.in <- lm(diversity_shannon ~ stage, data = subset(ps.inoc, vineyard == "Organic"))
model_biodynamic.in <- lm(diversity_shannon ~ stage, data = subset(ps.inoc, vineyard == "Biodynamic"))
model_conventional.in <- lm(diversity_shannon ~ stage, data = subset(ps.inoc, vineyard == "Conventional"))

# View the summaries of each model
summary(model_organic.in)
summary(model_biodynamic.in)
summary(model_conventional.in)

# Wild Wines #
# Fit mixed-effects models for each vineyard separately
model_organic.wi <- lm(diversity_shannon ~ stage, data = subset(ps.wild, vineyard == "Organic"))
model_biodynamic.wi <- lm(diversity_shannon ~ stage, data = subset(ps.wild, vineyard == "Biodynamic"))
model_conventional.wi <- lm(diversity_shannon ~ stage, data = subset(ps.wild, vineyard == "Conventional"))

# View the summaries of each model
summary(model_organic.wi)
summary(model_biodynamic.wi)
summary(model_conventional.wi)

model.list <- list(summary(model_organic.wi)$coefficients, summary(model_organic.in)$coefficients, 
                   summary(model_biodynamic.wi)$coefficients, summary(model_biodynamic.in)$coefficients, 
                   summary(model_conventional.wi)$coefficients, summary(model_conventional.in)$coefficients)


# Create a data frame from the extracted data
model.list.df <- data.frame(
  Coefficient = rownames(model.list),
  Estimate = coef_data[, "Estimate"],
  Std_Error = coef_data[, "Std. Error"],
  P_Value = coef_data[, "Pr(>|t|)"]
)

write.csv(model.list, file = "model-list_Shannon-diversity.csv", row.names = FALSE)


##############
## FIGURE 6 ##
##############

# Beta Diversity Plots

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


##############################################
## CORRELATION TEST (SENSORY AND MICROBIAL) ##
##############################################

## Berry PCoA

#Filter out berries
ps.rar.berry <- subset_samples(ps.phyF.tree, sample.type == "Berry")


#ordination
ord.bray.berry <- ordinate(ps.rar.berry, "PCoA", "bray", weighted=F)


#Axes
bray.1.2.berry <- plot_ordination(ps.rar.berry, 
                               ord.bray.berry, color="vineyard", 
                               axes = c(1, 2))

bray.1.3.berry <- plot_ordination(ps.rar.berry, 
                               ord.bray.berry, color="vineyard", 
                               axes = c(1, 3))

#Plot them
berry.br.1.2 <- bray.1.2.berry +
  scale_x_reverse() +
  geom_point(size=7.5) +
  geom_path() +
  geom_text(aes(label = stage), color = "black", fontface = "bold") +
  scale_colour_manual(values=colours) +
  theme_pcoa +
  theme(legend.position = "none")

berry.br.1.3 <- bray.1.3.berry +
  scale_x_reverse() +
  geom_point(size=7.5) +
  geom_path() +
  geom_text(aes(label = stage), color = "black", fontface = "bold") +
  scale_colour_manual(values=colours) +
  theme_pcoa +
  theme(legend.position = "none")

#PERMANOVA

library(vegan)

ps.rar.berry <- subset_samples(ps.phyF.prevF.rar, sample.type == "Berry")

metadata.berry <- as(sample_data(ps.rar.berry), "data.frame")

adonis2(distance(ps.rar.berry, method="bray") ~ vineyard, data = metadata.berry)


#merge!
ggarrange(berry.br.1.2, berry.br.1.3,
          nrow = 1, ncol = 2, common.legend = TRUE, labels = c("A)", "B)"), 
          font.label = list(size=18, face="bold", color="black"))

ggsave(filename = "R_output/SupplementaryFigure3.png", 
       width = 15, height = 7.5, dpi = 600)


#Make the matrices
pcoa_matrix <- ord.bray.berry$vectors
pcoa_matrix <- pcoa_matrix[, 1:5]
pcoa_matrix <- pcoa_matrix[2:9,]

pca_matrix <- read.csv("pca_matrix.csv", row.names = 1)


## Mantel Test ##

# Calculate the distance matrices from the PCA and PCoA matrices
pca_dist_matrix <- vegdist(pca_matrix)
pcoa_dist_matrix <- vegdist(pcoa_matrix)

# Add a constant value to make all elements of PCA and PCoA matrices non-negative
constant_value <- abs(min(c(min(pca_matrix), min(pcoa_matrix))))
pca_matrix_nonnegative <- pca_matrix + constant_value
pcoa_matrix_nonnegative <- pcoa_matrix + constant_value

# Calculate the distance matrices from the non-negative PCA and PCoA matrices
pca_dist_matrix <- vegdist(pca_matrix_nonnegative)
pcoa_dist_matrix <- vegdist(pcoa_matrix_nonnegative)

# Perform the Mantel test using the distance matrices
mantel_result <- mantel(pca_dist_matrix, pcoa_dist_matrix)

# Print the Mantel correlation coefficient and p-value
print(paste("Mantel correlation coefficient:", mantel_result$statistic))
print(paste("Mantel p-value:", mantel_result$signif))


## Procrustes Test ##

# Perform Procrustes analysis
procrustes_result <- procrustes(pca_matrix, pcoa_matrix)

# Get the Procrustes correlation coefficient
procrustes_corr <- 1 - procrustes_result$ss / sum(procrustes_result$X^2)

# Number of permutations (increase this for more accurate p-value estimation)
n_permutations <- 999

# Permutation testing to estimate p-value
permuted_correlations <- numeric(n_permutations)
for (i in 1:n_permutations) {
  shuffled_pcoa_matrix <- pcoa_matrix[sample(nrow(pcoa_matrix)), ]
  permuted_result <- procrustes(pca_matrix, shuffled_pcoa_matrix)
  permuted_corr <- 1 - permuted_result$ss / sum(permuted_result$X^2)
  permuted_correlations[i] <- permuted_corr
}

# Calculate the p-value (two-tailed test)
observed_correlation <- procrustes_corr
p_value <- sum(abs(permuted_correlations) >= abs(observed_correlation)) / n_permutations

# Print the Procrustes correlation coefficient and p-value
print(paste("Procrustes correlation coefficient:", observed_correlation))
print(paste("Procrustes p-value:", p_value))