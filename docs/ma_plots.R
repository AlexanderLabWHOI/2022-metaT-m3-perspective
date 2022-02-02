#
# Example log ratio (M) - mean average (A) plot (MA)
# and MANTA style plot with the same data but showing 
# pie charts representing the taxonomic breakdown
# for genes with significant differential expression
#
# Author: Rob Lampe, rlampe at ucsd dot edu
# 

# Load libraries 
library(ggplot2)
library(readr)
library(methods)
library(scatterpie)
library(RColorBrewer)
library(ggrastr)


# Static variables
#
# max_fc = max fold change value on plot
# alpha 
max_fc <- 5
min_fc <- max_fc * -1
alpha <- 0.05

# Read in deseq2 input
# Set padj = NA to 1
res <- read.csv('deseq2_means_results.csv')
res[is.na(res$padj),]$padj <- 1

# Calculate mean abundance and add pseudocount
res[,grep("mean", names(res))] <- res[,grep("mean", names(res))] + 1
res$mean <- log2(rowMeans(res[,c("treatment1_mean", "treatment2_mean")]))

#
# Plot settings
# adjust as needed
#
ylimits <- scale_y_continuous(breaks = seq(-4, 4, by = 2), limits = c(-5.5, 5.5))
xlimits <- scale_x_continuous(breaks = seq(0, 16, by=4), limits = c(0, max(res$mean))) 
xlabel <- xlab(bquote(~log[2]~"average normalized abundance"))
ylabel <- ylab(bquote(~log[2]~"fold change"))
hline <- geom_hline(yintercept = 0, linetype = 5)

theme_settings <- theme_bw() + theme(
  panel.border = element_rect(colour = "black"), 
  axis.line = element_line(colour = "black"), 
  # panel.grid.major = element_blank(), 
  # panel.grid.minor = element_blank(),
  axis.title=element_text(size=20), 
  axis.text.y = element_text(size=18, colour="black"),
  axis.text.x = element_text(size=18, colour="black"),
  plot.title = element_text(size=20, colour="black", face="bold"),
  plot.margin = unit(c(1,1,1,1), "cm"),
  legend.title = element_blank(),
  legend.position = "bottom",
  legend.text=element_text(size=16),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA)
)

# Set min p-value
res[res$padj < 0.01,]$padj <- 0.01

# Set max FC to max_fc
col <- "log2FoldChange"
if (nrow(res[res[,col] > max_fc,])) {
    res[res[,col] > max_fc, col] <- max_fc
}
if (nrow(res[res[,col] < min_fc,])) {
    res[res[,col] < min_fc, col] <- min_fc
}


#
# Regular MA plot
#
plot <- ggplot() + 
  rasterize(geom_point(data = res[res$padj > alpha & abs(res$log2FoldChange) < 5,],
                       aes(x = mean, y = log2FoldChange, size = -log10(padj)/15),
                       color = "grey")) +
  geom_point(data = res[res$padj > alpha & res$log2FoldChange == 5,],
             aes(x = mean, y = log2FoldChange, size = -log10(padj)/15),
             color = "grey", fill = "grey", shape = 24) +
  geom_point(data = res[res$padj > alpha & res$log2FoldChange == -5,],
             aes(x = mean, y = log2FoldChange, size = -log10(padj)/15),
             color = "grey", fill = "grey", shape = 25) +
  geom_point(data = res[res$padj <= alpha & abs(res$log2FoldChange) < 5,],
             aes(x = mean, y = log2FoldChange, size = -log10(padj)/15), color = "black") +
  geom_point(data = res[res$padj <= alpha & res$log2FoldChange == 5,],
             aes(x = mean, y = log2FoldChange, size = -log10(padj)/15), color = "black", fill = "black", shape = 24) +
  geom_point(data = res[res$padj <= alpha & res$log2FoldChange == -5,],
             aes(x = mean, y = log2FoldChange, size = -log10(padj)/15), color = "black", fill = "black", shape = 25) +
  hline + ylimits + xlimits + xlabel + ylabel + theme_settings +
  scale_size(range = c(1, 5)) +
  guides(size = "none")

plot
# uncomment to save as pdf
#ggsave("ma_plot.pdf", plot = plot, device = "pdf", width = 7, height = 8)


#
# Additional data for MANTA plot
#

# Select genes with padj < alpha
sig_genes <- res[res[grep("padj", names(res))] < alpha,]$gene


# Read in taxonomically annotated TPMs to generate pies. 
# Average TPMs for the treatments being compared were already calculated as meanTPM
tpms <- read_tsv('average_tpms.tsv')

# Select genus field in taxonomic string
tpms$taxonomy <- gsub("Eukaryota;Stramenopiles;Stramenopiles_X;Bacillariophyta;Bacillariophyta_X;[A-Za-z]*-centric-[A-Za-z-]*;", "", tpms$taxonomy)
tpms$taxonomy <- gsub(";[A-za-z0-9 ,-.()]*", "", tpms$taxonomy)

# Select genera of interest
# assign all others as "Other"

# This will select the top 4 most abundant genera
genera <- unique(tpms$taxonomy)
genera_totals <- c()
for (genus in genera) {
  genera_totals <- c(genera_totals, sum(tpms[grep(genus, tpms$taxonomy),]$meanTPM))
  names(genera_totals)[length(genera_totals)] <- genus
}
genera <- names(rev(sort(genera_totals))[1:4])

# Alternatively, genera can be selected manually
genera <- c("Chaetoceros", "Thalassiosira", "Eucampia", "Minutocellus")
tpms[!tpms$taxonomy %in% genera,]$taxonomy <- 'Other'

# Make dataframe with tpms for each genus of interest
# for all significantly DE genes and merge with DESeq2 results
df <- as.data.frame(matrix(0, nrow=nrow(res), ncol=(length(genera)+1)))
rownames(df) <- res$gene
names(df) <- c(genera, "Other")

for (gene in sig_genes) {
  for (genus in names(df)) {
    temp <- tpms[grep(gene, tpms$KO),]
    temp <- temp[grep(genus, temp$taxonomy),]
    df[gene, genus] <- sum(temp$meanTPM)
  }
}

df$gene <- rownames(df)
res <- merge(res, df, by="gene")

# Randomize order of genes
set.seed(42)
rows <- sample(nrow(res))
res <- res[rows,]
sig <- res[res$padj < alpha,]
# sig <- sig[sig$mean > 1,] # uncomment to set mean abundance cutoff

#
# Generate MANTA plot
#

# Color palette for pies
colors <- c(colorRampPalette(brewer.pal(4,"Set1"))(4), "#454545")
names(colors) <- c(genera, "Other")

# Base plot with rasterized grey circles for genes where padj > alpha
plot <- ggplot() + 
  rasterize(geom_point(data = res[res$padj > alpha,],
             aes(x = mean, y = log2FoldChange, size = -log10(padj)/15),
             color = "grey"))

# Iterate over significant genes to generate pies
for (i in seq(1, nrow(sig), by=5)) {
  plot <- plot + geom_scatterpie(data = sig[i:(i+4),],
                                   aes(x = mean, y = log2FoldChange, r = -log10(padj)/7.5),
                                   cols = c(genera, "Other"))
}

# Style plot
plot <- plot +
  scale_fill_manual(values = colors) +
  hline + ylimits + xlimits + xlabel + ylabel + theme_settings +
  scale_size(range = c(1, 5)) +
  guides(size = "none") +
  guides(fill = guide_legend(ncol=3))

plot
# uncomment to save as pdf
#ggsave("manta_plot.pdf", plot = plot, device = "pdf", width = 7, height = 8)
