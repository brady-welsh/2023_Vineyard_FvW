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

# Don't forget to set working directory. (2023_Vineyard_FvW/ITS_Fungi_ITS1-ITS2R/data/)

###########################
## PREVELANCE FILTERTING ##
###########################

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



#################
## RAREFACTION ##
#################

#Set seed, check the sample with the lowest read count
set.seed(1997)

sort(sample_sums(ps.phyF.prevF))

#Rarefy n number of reads where n = read count of samples with fewest reads
ps.phyF.prevF.rar <- rarefy_even_depth(ps.phyF.prevF, sample.size = 653)

