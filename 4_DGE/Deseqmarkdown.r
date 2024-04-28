
---
title: "DESEQ2"
author: "Deekshya"
date: "2024-04-24"
output: pdf_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

```{r cars}
summary(cars)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2", force = TRUE)
library(DESeq2)
### You will need to set your working directory to the location you have your data.
# You can do this by using  the Session menu to set working directory To Source File Directory

#### Install the DESeq2 package if you have not already
## try http:// if https:// URLs are not supported


## Use the Session menu to set working directory To Source File Directory
```

````{r cars}
summary (cars)
countdata <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))
dim(countdata)
head(countdata)
```

```{r cars}
summary (cars)
coldata <-(read.table("Amphibian_PHENO_DATA.txt", header=TRUE, row.names=1))
dim(coldata)
head(coldata)
```

```{r cars}
summary  (cars)
all(rownames(coldata) %in% colnames(countdata))
countdata <- countdata[, rownames(coldata)]
all(rownames(coldata) == colnames(countdata))
```

``` {r cars}
summary(cars)

## Create the DESEQ dataset and define the statistical model (page 6 of the manual)
dds <- DESeqDataSetFromMatrix(countData = countdata, colData=coldata, design = ~treatment)


#look at it
dds
```

```{r cars}
summary(cars)
#Prefiltering manual
dds <- dds[ rowSums(counts(dds)) > 0, ]
dds
dds <- DESeq(dds)
res <- results(dds)
res
resOrdered <- res[order(res$padj),]
resOrdered
```

```{r cars}
summary(cars)
## set factors for statistical analyses

###### Note you need to change condition to treatment (to match our design above)
#  and levels to our treatment names in the PHENO_DATA: Ad_lib is the control, Caloric_Restriction is the treatment group
# example:
#dds$condition <- factor(dds$condition, levels=c("untreated","treated"))
dds$condition <- factor(dds$treatment, levels=c("control","heat_stress"))
```

#questions

# We can order our results table by the smallest adjusted p value:
resOrdered <- res[order(res$padj),]
resOrdered

# We can summarize some basic tallies using the summary function the default is p<0.1.
summary(res)

#How many adjusted p-values were less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)

#If the adjusted p value will be a value other than 0.1, alpha should be set to that value:
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

#MAplot
#chatgpt 
# Add a legend to the plot
plotMA(res05, main="DESeq2", ylim=c(-8,8))

legend("topright",  # Position of the legend
       legend=c("Significantly Up", "Significantly Down", "Not Significant"),  # Legend text
       col=c("red", "blue", "black"),  # Colors corresponding to your legend text
       pch=20,  # Symbol type (filled circle)
       cex=0.75)  # Adjust the size of the legend text
```

#After calling plotMA, one can use the function identify to interactively detect the row number of individual genes by clicking on the plot. 
# One can then recover the gene identiers by saving the resulting indices:
idx <- identify(res$baseMean, res$log2FoldChange)
# after selecting a gene. You need to press escape to move on
rownames(res)[idx]
```

##  Write your results to a file 
write.csv(as.data.frame(res05), file="DGESeq_results_greaterthan0-frfr.csv")

```{r cars}
summary(cars)
## 2.1.2 Extracting transformed values
rld <- rlog(dds)
vsd <- varianceStabilizingTransformation(dds)
head(assay(rld), 3)

### Heatmap of the count matrix
# Assuming 'vsd' is already defined and the necessary libraries are loaded

# Load the necessary libraries if not already loaded
library("genefilter")
library("pheatmap")

# Get the top variable genes
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)

# Extract the expression matrix for top variable genes
mat <- assay(vsd)[topVarGenes, ]
mat <- mat - rowMeans(mat)

# Check the column names in colData(vsd)
colnames(colData(vsd))

# Assuming the column names are "treatment", "type", "sizeFactor", and "condition"
# Create a data frame 'anno' containing columns 'treatment' and 'type'
anno <- as.data.frame(colData(vsd)[, c("treatment", "type")])
df <- as.data.frame(colData(dds)[,c("treatment","type")])
pheatmap(mat, annotation_col = anno)
````

## 2.2.1 Heatmap of the count matrix
#  library("pheatmap")
#  select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:50]
#  nt <- normTransform(dds) # defaults to log2(x+1)
#  log2.norm.counts <- assay(nt)[select,]
#   df <- as.data.frame(colData(dds)[,c("treatment","type")])
#  pheatmap(assay(vsd)[mat,], cluster_rows=TRUE, show_rownames=TRUE,
#          cluster_cols=TRUE, annotation_col=df)

# pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
#  pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,  cluster_cols=FALSE, annotation_col=df) 

```{r cars}
summary(cars)
#2.2.2 Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(rld)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$treatment)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
```

# 2.2.3 Principal component plot of the samples
plotPCA(rld, intgroup=c("treatment"))
```




Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

```{r cars}
summary(cars)
```

## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
