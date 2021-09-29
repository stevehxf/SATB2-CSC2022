## Authour: Xiaofeng Steve Huang
## Please cite: https://doi.org/10.1016/j.stem.2021.09.004
## Paper: [SATB2 preserves colon stem cell identity and mediates ileum-colon conversion via enhancer remodeling]. Cell Stem Cell.
## Project: Satb2 knock out time course in vivo
## Date of update: June 26th 2021

# Load libraries
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(ggrepel)
library(data.table)
library(annotables)
library(org.Mm.eg.db)

#### Data cleaning ####
# load raw feature counts
data.raw <- read.table("data/satb2_time_course_cleancounts.txt", header=T, row.names=1, check.names = F)  # Put your unnormalized counts file in data directory

# extract protein_coding gene counts from the raw counts
## get protein coding genes symbol from grcm3
coding_genes <- grcm38 %>% dplyr::filter(biotype == "protein_coding") 
coding_genes <- coding_genes$symbol
## subset the count table by protein coding genes
idx <- row.names(data.raw) %in% c(coding_genes,"Inava") # Inava is missing in this package
data <- data.raw[idx,]

# load and organize meta
## load meta from the colleague who did the experiment (Wei Gu)
meta <- read.csv("meta/Time_course.csv", row.names=1, header=T) # Put your meta file in meta directory
## arrange the meta according to the column names of count data
meta <- meta[colnames(data),]
## set the time variable as factor
meta$time <- factor(meta$time)
## re-order the time facotr
ref <- c("colon_wt_tam", "colon_ko_0", "colon_ko_1", "colon_ko_2", "colon_ko_4", "colon_ko_6", "colon_ko_30", "ileum_wt_30")
meta$TissueGenotypeTime <- factor(meta$TissueGenotypeTime,ref)

# 3. Check the input files
##  (1) Check the data class to confirm they are "data.frame"
class(data)
class(meta)
##  (2) Check the colume names of data match the row names of meta: both should be TRUE
all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

#### Model QC ####

# create a directory for QC output
if(!dir.exists("./results/qc/")){
  dir.create("./results/qc/", recursive = TRUE)
}

# RNA distribution
## RNA distribution function
plotHistFunc <- function(dataframe, pre, na.rm = TRUE, ...) {     # make a function
  nm <- colnames(dataframe)              # put column names in the variable nm
  for (i in seq_along(nm)) {            # for loop, extract the index number from nm
    plots <- ggplot(dataframe,aes_string(x = nm[i])) +    # ggplot(file,aes_from_string)
      geom_histogram(stat ="bin", bins = 200, fill = "#446179", alpha = 0.8) + # plot histogram
      xlim(-5, 300)  + # x axis limit
      xlab("Raw expression counts") + # x lable
      ylab("Number of genes") + # y lable
      theme_light() + # theme
      theme(panel.border = element_rect(linetype = "solid", fill = NA, size = 1, colour = "black")) # border
    ggsave(plots,filename=paste("results/qc/",pre,nm[i],".png",sep="")) # export
  }}
## run the function
plotHistFunc(data,"") # plot counts distribution for RNA samples automatically. Call the function to check other dataset if you want.

# NB model (mean>variance)
## mean VS. variance # if var > mean, then model fits
mean_counts <- apply(data, MARGIN = 1, mean) # calculate mean for all genes
variance_counts <- apply(data, MARGIN = 1, var) # calculate variance for all genes
df <- data.frame(mean_counts, variance_counts) # make a data frame: x=mean, y=variance
meanVSvar <- ggplot(df) + # plot the dataframe
  geom_point(aes(x=mean_counts, y=variance_counts),alpha = 0.4) + #plot points with 60% transparency
  geom_line(aes(x=mean_counts, y=mean_counts), color = "red") +  # red line indicates mean=variance
  scale_y_log10() + scale_x_log10() +
  theme_light() +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, size = 1, colour = "black")) # border
ggsave(meanVSvar,filename="./results/qc/meanVSvar.jpeg") # export 

#### Overall exploratory analysis ####
# Genes that are barely detectable will contribute noise to the resulting clusters
# pre-filtering: only keep the genes that at least 3 (n.rep.min) samples have more than 5 counts
n.min <- 3
keep <- rowSums((data)>5) >= n.min
data <- data[keep,]
# create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = data,  
                              colData = meta,
                              design = ~tissue_genotype)  # Comparing factor 


# rlog Transformed counts for data visualization
rld <- rlog(dds, blind=TRUE) 

# PCA
plotPCA(rld, intgroup=c("tissue_genotype"), ntop=1000, returnData = F) +
  ggtitle ("PCA") + # main title  
  theme_light() +  # Theme
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(linetype = "solid", fill = NA, size = 1, colour = "black"),
        legend.title=element_blank()
  )+
  geom_text_repel(aes(label = meta$time)) +
  coord_fixed() +
  NULL

# create directory for figures output
if(!dir.exists("./figures/")){
  dir.create("./figures/")
}

# export the PCA plot
ggsave("figures/overall.pca.pdf")

#### Time course DE Analysis ####
# time course analysis for Satb2KO

# subset Satb2 knock out data and meta
data_satb2_ko <- data[,meta$tissue_genotype == "colon_ko"]
meta_satb2_ko <- meta[meta$tissue_genotype == "colon_ko",]

# Create DESeq2 object using time factor as design
dds_ko <- DESeqDataSetFromMatrix(countData = data_satb2_ko,  
                                 colData = meta_satb2_ko,
                                 design = ~time)  # Comparing factor

# DEseq analysis with LRT test to identify DEGs accorss different time points
dds_lrt_ko_time <- DESeq(dds_ko, test="LRT", reduced = ~1)

# extract the results
res_LRT_ko <- results(dds_lrt_ko_time)

# extract the coefficients
# betas_ko <- coef(dds_lrt_ko_time)
# normalization
# Extract normalized counts, check the size factors and total counts before and after normalization
dds_ko <- estimateSizeFactors(dds_ko)
sizeFactors(dds_ko) # check the size factors
normalized_counts_ko <- counts(dds_ko, normalized=TRUE)

# export normalized counts for downstream analysis
if(!dir.exists("./results/normalization/")){
  dir.create("./results/normalization/", recursive = TRUE)
}
write.table(normalized_counts_ko, file="results/normalization/normalized_counts_ko.txt", sep="\t", quote=F, col.names=NA)

# export the meta data for downstream analysis
write.csv(meta_satb2_ko, "meta/meta_ko.csv")

# Subset the LRT results to return genes with padj < 0.01
padj.cutoff <- 0.01
sig_res_LRT_ko <- res_LRT_ko %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff)

# Get significant gene list
sigLRT_genes.ko <- sig_res_LRT_ko %>% 
  pull(gene)

# how many DEGs?
length(sigLRT_genes.ko)

#### Time course clustering analysis ####

# VST Transform counts for data visualization
vsd.ko <- vst(dds_lrt_ko_time, blind=F)
# rename the column
if(all(colnames(vsd.ko) == rownames(meta_satb2_ko))){
  colnames(vsd.ko) <-  meta_satb2_ko$name
}

# extract the matrix and convert it into tibble format
vsd.ko.df <- vsd.ko %>% 
  assay() %>% 
  data.frame() %>% 
  rownames_to_column("gene") %>% 
  as_tibble()

# Summarise transformed counts to get the mean of scaled vsd for data visualization
vsd_mean <- vsd.ko.df %>% 
  ### convert to long format
  pivot_longer(cols = colnames(vsd.ko.df)[-1], names_to = "name", values_to = "vsd")  %>% 
  ### join with sample meta data by "name"
  full_join(meta_satb2_ko, by = ("name")) %>% 
  ### filter to retain only genes of interest (DEG across time points)
  filter(gene %in% sigLRT_genes.ko) %>% 
  ### for each gene
  group_by(gene) %>% 
  ### scale the vsd column
  mutate(vsd_scaled = (vsd - mean(vsd))/sd(vsd)) %>% 
  ### for each gene, strain and minute
  group_by(gene, time) %>%
  ### calculate the mean (scaled) cts
  summarise(mean_vsd_scaled = mean(vsd_scaled),
            nrep = n()) %>% 
  ungroup()

# hierarchical clustering ####

## scale the vst matrix of DEGs
hclust_matrix <- vsd.ko.df %>% 
  ### filter for DEG
  filter(gene %in% sigLRT_genes.ko) %>% 
  column_to_rownames("gene") %>% 
  ### convert to matrix
  as.matrix() %>% 
  ### transpose
  t() %>% 
  ### apply scaling to each column of the matrix (genes)
  scale() %>% 
  ### transpose back so genes are as rows again
  t()

## calculate distance of genes
gene_dist <- dist(hclust_matrix)

## clustering
gene_hclust <- hclust(gene_dist, method = "complete")
plot(gene_hclust, labels = FALSE)
#abline(h = 3.6, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram

## cut the tree
gene_cluster <- cutree(gene_hclust, k = 9) %>% 
  ### turn the named vector into a tibble
  enframe() %>% 
  ### rename the columns
  dplyr::rename(gene = name, cluster = value)

## join the summarized vsd with cluster infomration
vsd_cluster <- vsd_mean %>% 
  inner_join(gene_cluster, by = "gene")

# plot curves
vsd_cluster %>% 
  ggplot(aes(time, mean_vsd_scaled)) +
  geom_line(aes(group = gene), alpha = 0.1) +
  geom_line(stat = "summary", fun = "mean", colour = "firebrick", size = 1.5, 
            aes(group = 1)) +
  facet_grid(cols = vars(cluster))
# export the figure
ggsave("figures/timecourse.pdf", width = 10, height = 2)

# plot heatmap
## prepare wide matrix
vsd_mean_wide <- vsd_mean %>% 
  ### remove replicate numver
  dplyr::select(-nrep) %>% 
  ### convert to wide format
  pivot_wider(names_from = time, values_from = mean_vsd_scaled )%>% 
  ### join gene cluster information
  inner_join(gene_cluster, by = "gene") %>% 
  ### convert column to rownames
  column_to_rownames("gene") %>% 
  ### covert to dataframe format
  as.data.frame() %>% 
  ### order the cluster
  arrange(factor(cluster, levels = c(1,3,4,6,
                                     9,
                                     2,5, 8,
                                     7)))
## export the data 
vsd_mean_wide %>% write.csv("results/scaled_vsd_mean_wide.csv")

## plot
p1 <- pheatmap(vsd_mean_wide[,-7],
               cluster_cols = F,
               cluster_rows = F,
               fontsize_row = 1,
               color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(50)
)
### export
ggsave(p1, filename = "figures/timecourse_heatmap.pdf")
