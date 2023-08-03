library(phyloseq)
library(ggplot2)
library(decontam)
library(qiime2R)
library(devtools)
library(dplyr)
library(readr)

# Dont forget to set working directory to 2023_Vineyard_FvW/ITS_Fungi_ITS1-ITS2R/data/QIIME2_output/

physeq<-qza_to_phyloseq(
  features="Vineyard-ITS-table_tidy.qza",
  taxonomy="Vineyard-ITS-UNITE.qza",
  metadata="metadata/Vineyard-ITS-metadata_tidy.tsv"
)

# Put sample_data into a ggplot-friendly data.frame and plot
df <- as.data.frame(sample_data(physeq))
df$LibrarySize <- sample_sums(physeq)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
plot1 <- ggplot(data=df, aes(x=Index, y=LibrarySize, color=sample.type)) + geom_point()

sample_data(physeq)$is.neg <- sample_data(physeq)$control == "TRUE"
contamdf.prev <- isContaminant(physeq, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)


contamdf.prev05 <- isContaminant(physeq, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)

physeq.pa <- transform_sample_counts(physeq, function(abund) 1*(abund>0))
physeq.pa.neg <- prune_samples(sample_data(physeq.pa)$control == "TRUE", physeq.pa)
physeq.pa.pos <- prune_samples(sample_data(physeq.pa)$control == "FALSE", physeq.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(physeq.pa.pos), pa.neg=taxa_sums(physeq.pa.neg),
                    contaminant=contamdf.prev05$contaminant)
plot2 <- ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

#Create histogram of decontam scores
plot3 <- ggplot(contamdf.prev, aes(x=p)) +
  geom_histogram(binwidth = 0.01) +
  labs (x = "Decontam score", y="Number of species")

#Merge taxonomy to decontam output and create a tsv
tax <- as(tax_table(physeq.pa), "matrix")

tax.df <- data.frame(tax)

merged <- merge(contamdf.prev05, tax.df, by=0, sort=FALSE)

contaminant.asvs <- subset(merged, contaminant == "TRUE")
biological.asvs <- subset(merged, contaminant == "FALSE")

write_tsv(merged, file="Decontam_output/all_decontam.tsv")
write_tsv(contaminant.asvs, file="Decontam_output/contaminant_asvs.tsv")
