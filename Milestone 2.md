# Milestone 2
## 5. Generation of differential expression results
###5.1 Processing of dds data frames
####Pre-filtering
The aim is to remove low-count genes. Here, I keep genes with at least 10 counts.

```{r}
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
```
####Note on factor levels
The aim is to set a factor for comparing differences in gene expression, i.e. upperlobe lung cancer and lowerlobe lung cancer.

```{r}
dds$condition <- factor(dds$condition, levels = c("upperlobe","lowerlobe"))
```

```{r}
dds$condition <- droplevels(dds$condition)
```

###5.2 Obtain data frames for differential gene expression results
####Get Results
```{r}
dds <- DESeq(dds)
```
![5-1](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/5-1.png?raw=true)

Here I use Independent hypothesis weighting (IHW) to filter the p-values. I want to filter the differentially expressed genes by p-value < 0.05, so, I add alpha=0.05.

```{r}
library("IHW")
res <- results(dds, filterFun=ihw, contrast=c("condition","upperlobe","lowerlobe"), alpha=0.05)
res
```
![5-2](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/5-2.png?raw=true)
### 5.3 Change the Ensembl id in the result table to the gene name.
#### remove the version number of the gene Ensembl number
```{r}
enemble_id <- substr(row.names(res),1,15)
rownames(res) <- enemble_id
```
#### Add a colum to result table
```{r}
RawCounts <- res
Ensembl_ID <- data.frame(Ensembl_ID = row.names(RawCounts))
rownames(Ensembl_ID) <- Ensembl_ID[,1]
RawCounts <- cbind(Ensembl_ID,RawCounts)
```
#### Download gencode.v38.basic.annotation.gtf
```
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.basic.annotation.gtf.gz
```
```
gunzip gencode.v38.basic.annotation.gtf.gz
```
#### Create a file to associate the ensembl id and gene id
```{r}
get_map = function(input) {
  if (is.character(input)) {
    if(!file.exists(input)) stop("Bad input file.")
    message("Treat input as file")
    input = data.table::fread(input, header = FALSE)
  } else {
    data.table::setDT(input)
  }
  
  input = input[input[[3]] == "gene", ]
  
  pattern_id = ".*gene_id \"([^;]+)\";.*"
  pattern_name = ".*gene_name \"([^;]+)\";.*"
  
  gene_id = sub(pattern_id, "\\1",input[[9]])
  gene_name = sub(pattern_name, "\\1", input[[9]])
  
  Ensembl_ID_TO_Genename <- data.frame(gene_id = gene_id, gene_name = gene_name, stringsAsFactors = FALSE)
  return(Ensembl_ID_TO_Genename)
}
Ensembl_ID_TO_Genename <- get_map("gencode.v38.basic.annotation.gtf")

gtf_Ensembl_ID <- substr(Ensembl_ID_TO_Genename[,1],1,15)
Ensembl_ID_TO_Genename <- data.frame(gtf_Ensembl_ID, Ensembl_ID_TO_Genename[,2])
colnames(Ensembl_ID_TO_Genename) <- c("Ensembl_ID","gene_id")
write.csv(Ensembl_ID_TO_Genename,file = "Ensembl_ID_TO_Genename.csv")
```
[Ensembl ID TO Genename.csv](https://raw.githubusercontent.com/DZBohan/zhangboh_final_project/main/Files/Ensembl_ID_TO_Genename.csv)
####Replace ensembl id with gene id
```{r}
res_g <-merge(Ensembl_ID_TO_Genename,RawCounts,by="Ensembl_ID")
```
#### Remove unnecessary columns and duplicate gene ids
```{r}
res_g <- res_g[order(res_g[,"gene_id"]),]
index <- duplicated(res_g$gene_id)
res_g <- res_g[!index,]
rownames(res_g) <- res_g[,"gene_id"]
res_g <- res_g[,-c(1:2)]
```
#### Check the new table
```{r}
head(res_g)
```
![5-3](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/5-3.png?raw=true)
### 5.4 Optimization process
####Log fold change (LFC) shrinkage
I use the *apeglm* method for LFC shrinkage which is useful for gene visualization and ranking.

```{r}
resultsNames(dds)
```
![5-4](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/5-4.png?raw=true)

```{r}
resLFC <- lfcShrink(dds, coef="condition_lowerlobe_vs_upperlobe", type="apeglm")
resLFC
```
![5-5](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/5-5.png?raw=true)
####Parallelization
Splitting the work into 4 cores to speed it up.

```{r}
library("BiocParallel")
register(MulticoreParam(4))
```
###5.5 View the number of differentially expressed genes based on p-value
Sort the data in the result table according to the p-value.

```{r}
resOrdered <- res_g[order(res_g$pvalue),]
```

Check the number of differentially expressed genes at adjusted p-values less than 0.05.

```{r}
sum(res_g$padj < 0.05, na.rm=TRUE)
```
![5-6](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/5-6.png?raw=true)

Therefore, I finally screened 433 genes differentially expressed in upperlobe lung cancer and lowerlobe lung cancer.
## 6 Exploring and exporting results
### 6.1 MA-plot
```{r}
plotMA(res, ylim=c(-2,2))
```
![6-1](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/6-1.png?raw=true)

```{r}
plotMA(resLFC, ylim=c(-2,2))
```
![6-2](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/6-2.png?raw=true)

```{r}
resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")
```
![6-3](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/6-3.png?raw=true)

```{r}
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
```
![6-4](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/6-4.png?raw=true)
###6.2 Plot counts
```{r}
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
```
![6-5](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/6-5.png?raw=true)

```{r}
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))
```
![6-6](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/6-6.png?raw=true)
###6.3 More information on results colums
```{r}
mcols(res)$description
```
![6-7](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/6-7.png?raw=true)
###6.4 Write csv files
####Original table of all genes
```{r}
getwd()
write.csv(as.data.frame(res_g), 
          file="res_g.csv")
```
[Here is the file res_g.csv](https://raw.githubusercontent.com/DZBohan/zhangboh_final_project/main/Files/res_g.csv)
####Differential expression genes
Genes with p-values <0.05 were screened and determined as differentially expressed genes. And a new column was created to record the up- or down-regulation of genes.

```{r}
resSig_0.05 <- subset(resOrdered, padj < 0.05)
resSig_0.05[which(resSig_0.05$log2FoldChange > 0), "up_down"] <- "up"
resSig_0.05[which(resSig_0.05$log2FoldChange < 0), "up_down"] <- "down"
resSig_0.05
```

![6-8](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/6-8.png?raw=true)

This table contains the 433 differentially expressed genes I screened and their up- or down-regulation information.

```{r}
getwd()
write.csv(as.data.frame(resSig_p0.05), 
          file="differential_expression.csv")
```
[Here is the file differential_expression.csv](https://github.com/DZBohan/zhangboh_final_project/blob/main/Files/differential_expression.csv)
##7 Data transformations and visualization
###7.1 Extracting transformed values
```{r}
vsd <- vst(dds, blind=FALSE)
```
###7.2 Effects of transformations on the variance
```{r}
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))
```
![7-1](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/7-1.png?raw=true)

```{r}
meanSdPlot(assay(vsd))
```
![7-2](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/7-2.png?raw=true)
##8 Data quality assessment by sample clustering and visualization
###8.1 Heatmap of the count matrix
####Heatmap of ntd
```{r}
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition", "sizeFactor")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, show_colnames = FALSE)
```
![8-1](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/8-1.png?raw=true)
####Heatmap of vsd
```{r}
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, show_colnames = FALSE)
```
![8-2](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/8-2.png?raw=true)
###8.2 Heatmap of the sample-to-sample distances
```{r}
sampleDists <- dist(t(assay(vsd)))
```
```{r}
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, show_rownames = FALSE)
```
![8-3](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/8-3.png?raw=true)
###8.3 Principal component plot of the samples
```{r}
plotPCA(vsd, intgroup="condition")
```
![8-4](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/8-4.png?raw=true)