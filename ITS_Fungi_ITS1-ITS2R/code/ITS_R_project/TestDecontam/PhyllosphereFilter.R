library(phyloseq)
library(ggplot2)
library(decontam)
library(qiime2R)
library(devtools)
library(dplyr)
library(readr)

# Dont forget to set working directory to 2023_Vineyard_FvW/ITS_Fungi_ITS1-ITS2R/data/QIIME2_output/

# Put sample_data into a ggplot-friendly data.frame and plot
df <- as.data.frame(sample_data(ps))
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
plot1 <- ggplot(data=df, aes(x=Index, y=LibrarySize, color=sample.type)) + geom_point()


ps.bio <- subset_samples(ps, vineyard == "Biodynamic")
ps.con <- subset_samples(ps, vineyard == "Conventional")
ps.org <- subset_samples(ps, vineyard == "Organic")

#Run on Biodynamic samples (ps.bio)

sample_data(ps.bio)$is.neg <- sample_data(ps.bio)$sample.type == "Berry"
bio.phyllo <- isContaminant(ps.bio, method="prevalence", neg="is.neg")
table(bio.phyllo$contaminant)


bio.phyllo05 <- isContaminant(ps.bio, method="prevalence", neg="is.neg", threshold=0.6)
table(bio.phyllo05$contaminant)

ps.bio.pa <- transform_sample_counts(ps.bio, function(abund) 1*(abund>0))
ps.bio.pa.neg <- prune_samples(sample_data(ps.bio.pa)$sample.type == "Berry", ps.bio.pa)
ps.bio.pa.pos <- prune_samples(sample_data(ps.bio.pa)$sample.type != "Berry", ps.bio.pa)
# Make data.frame of prevalence in positive and negative samples
df.bio.phyllo.pa <- data.frame(pa.pos=taxa_sums(ps.bio.pa.pos), pa.neg=taxa_sums(ps.bio.pa.neg),
                    contaminant=bio.phyllo05$contaminant)
plot2 <- ggplot(data=df.bio.phyllo.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Berries)") + ylab("Prevalence (Ferments)")

#Create histogram of decontam scores
plot3 <- ggplot(bio.phyllo, aes(x=p)) +
  geom_histogram(binwidth = 0.01) +
  labs (x = "Decontam score", y="Number of species")

#Merge taxonomy to decontam output and create a tsv
tax <- as(tax_table(ps.bio.pa), "matrix")

tax.df <- data.frame(tax)

all.asvs <- merge(bio.phyllo05, tax.df, by=0, sort=FALSE)

bio.phyllo.asvs <- subset(all.asvs, contaminant == "TRUE")
bio.non.phyllo.asvs <- subset(all.asvs, contaminant == "FALSE")


#Run Conventional samples (ps.con)

sample_data(ps.con)$is.neg <- sample_data(ps.con)$sample.type == "Berry"
con.phyllo <- isContaminant(ps.con, method="prevalence", neg="is.neg")
table(con.phyllo$contaminant)


con.phyllo05 <- isContaminant(ps.con, method="prevalence", neg="is.neg", threshold=0.6)
table(con.phyllo05$contaminant)

ps.con.pa <- transform_sample_counts(ps.con, function(abund) 1*(abund>0))
ps.con.pa.neg <- prune_samples(sample_data(ps.con.pa)$sample.type == "Berry", ps.con.pa)
ps.con.pa.pos <- prune_samples(sample_data(ps.con.pa)$sample.type != "Berry", ps.con.pa)
# Make data.frame of prevalence in positive and negative samples
df.con.phyllo.pa <- data.frame(pa.pos=taxa_sums(ps.con.pa.pos), pa.neg=taxa_sums(ps.con.pa.neg),
                               contaminant=con.phyllo05$contaminant)
plot2 <- ggplot(data=df.con.phyllo.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Berries)") + ylab("Prevalence (Ferments)")

#Create histogram of decontam scores
plot3 <- ggplot(con.phyllo, aes(x=p)) +
  geom_histogram(binwidth = 0.01) +
  labs (x = "Decontam score", y="Number of species")

#Merge taxonomy to decontam output and create a tsv
tax <- as(tax_table(ps.con.pa), "matrix")

tax.df <- data.frame(tax)

all.asvs <- merge(con.phyllo05, tax.df, by=0, sort=FALSE)

con.phyllo.asvs <- subset(all.asvs, contaminant == "TRUE")
con.non.phyllo.asvs <- subset(all.asvs, contaminant == "FALSE")


#Run Organic samples (ps.org)

sample_data(ps.org)$is.neg <- sample_data(ps.org)$sample.type == "Berry"
org.phyllo <- isContaminant(ps.org, method="prevalence", neg="is.neg")
table(org.phyllo$contaminant)


org.phyllo05 <- isContaminant(ps.org, method="prevalence", neg="is.neg", threshold=0.6)
table(org.phyllo05$contaminant)

ps.org.pa <- transform_sample_counts(ps.org, function(abund) 1*(abund>0))
ps.org.pa.neg <- prune_samples(sample_data(ps.org.pa)$sample.type == "Berry", ps.org.pa)
ps.org.pa.pos <- prune_samples(sample_data(ps.org.pa)$sample.type != "Berry", ps.org.pa)
# Make data.frame of prevalence in positive and negative samples
df.org.phyllo.pa <- data.frame(pa.pos=taxa_sums(ps.org.pa.pos), pa.neg=taxa_sums(ps.org.pa.neg),
                               contaminant=org.phyllo05$contaminant)
plot2 <- ggplot(data=df.org.phyllo.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Berries)") + ylab("Prevalence (Ferments)")

#Create histogram of decontam scores
plot3 <- ggplot(org.phyllo, aes(x=p)) +
  geom_histogram(binwidth = 0.01) +
  labs (x = "Decontam score", y="Number of species")

#Merge taxonomy to decontam output and create a tsv
tax <- as(tax_table(ps.org.pa), "matrix")

tax.df <- data.frame(tax)

all.asvs <- merge(org.phyllo05, tax.df, by=0, sort=FALSE)

org.phyllo.asvs <- subset(all.asvs, contaminant == "TRUE")
org.non.phyllo.asvs <- subset(all.asvs, contaminant == "FALSE")

#save files
#write_tsv(org.non.phyllo.asvs, file="Decontam_output/org-non-phyllo.tsv")
#write_tsv(con.non.phyllo.asvs, file="Decontam_output/con-non-phyllo.tsv")
#write_tsv(bio.non.phyllo.asvs, file="Decontam_output/bio-non-phyllo.tsv")

#Execute phyllosphere filter
rownames(org.phyllo.asvs) <- org.phyllo.asvs[,1]
keepTaxa.org <- rownames(org.phyllo.asvs)
ps.org.phyllo <- prune_taxa(keepTaxa.org, ps.org)

rownames(bio.phyllo.asvs) <- bio.phyllo.asvs[,1]
keepTaxa.bio <- rownames(bio.phyllo.asvs)
ps.bio.phyllo <- prune_taxa(keepTaxa.bio, ps.bio)

rownames(con.phyllo.asvs) <- con.phyllo.asvs[,1]
keepTaxa.con <- rownames(con.phyllo.asvs)
ps.con.phyllo <- prune_taxa(keepTaxa.con, ps.con)

#Merge phyloseq objects
ps.phyF <- merge_phyloseq(ps.con.phyllo, ps.bio.phyllo, ps.org.phyllo)

