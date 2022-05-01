# <span id="jump"><center>Differential gene expression in TCGA within stage 1 lung cancers occurring in lower and upper lobe using DeSEQ2</center></span>

### <center>Bohan Zhang 12/03/2021</center>

## Index

- [1. Initial screening of data](#1)

	- [1.1 Screening of Files](#1.1)

	- [1.2 Screening of Upper Lobe Lung Cancer Cases](#1.2)

	- [1.3 Screening of Lower Lobe Lung Cancer Cases](#1.3)

- [2. Data Download and Collation](#2)

	- [2.1 Installation and Configuration of gdc-client](#2.1)

	- [2.2 Download of Count Files](#2.2)

	- [2.3 Organizing Count Files](#2.3)
	
	- [2.4 Metadata Json File Download](#2.4)

- [3. Secondary Screening of Data](#3)

	- [3.1 File Name & TCGA Number Mapping File](#3.1)

	- [3.2 Files Batch Rename](#3.2)

	- [3.3 Deletion of Unwanted Files](#3.3)

- [4. Load Count Files into Vignette](#4)

	- [4.1 Prefix Adding](#4.1)

	- [4.2 Loading Count Files](#4.2)

- [5. Generation of Differential Expression Results](#5)

	- [5.1 Processing of dds Dataframes](#5.1)

	- [5.2 Obtaining Dataframes](#5.2)

	- [5.3 Gene Names Changing](#5.3)

	- [5.4 Optimization](#5.4)

	- [5.5 The Number of Differentially Expressed Genes](#5.5)

- [6. Results Exploring & Exporting](#6)

	- [6.1 MA-plot](#6.1)

	- [6.2 Plot Counts](#6.2)

	- [6.3 More Information](#6.3)

	- [6.4 CSV Files Writing](#6.4)

- [7. Data Transformations & Visualization](#7)

	- [7.1 Extracting Transformed Values](#7.1)

	- [7.2 Effects of Transformations](#7.2)

- [8. Data Quality Assessment](#8)

	- [8.1 Heatmap of the Count Matrix](#8.1)

	- [8.2 Heatmap of the Sample-to-sample Distances](#8.2)

	- [8.3 Principal Component Plot](#8.3)

	- [6.4 CSV Files Writing](#6.4)

- [9. Evaluation of Differentially Expressed Genes](#9)

- [10. Conclusion](#10)

- [11. Known Issues](#11)

## <h2 id="1">1. Initial screening of data</h2>

### 1.1 Screening of Files
Data Category <- transcriptome profiling

Experimental Strategy <- RNA-Seq

Workflow Type <- HTSeq - Counts

Access <- open
### 1.2 Screening of Upper Lobe Lung Cancer Cases
Diagnoses Ajcc Pathologic Stage <- stage ia/stage ib/stage i

Diagnoses Tissue or Organ of Origin <- upper lobe, lung

Primary Site <- bronchus and lung

Program <- TCGA

Vital Status <- alive

Race <- white

Ethnicity <- not hispanic or latino


Now, I have filtered out 156 files and 136 cases. The number of files and cases is different because some cases have duplicate files; however, since I will be downloading files, the filtering in this step is incomplete. The second filtering will be done in the later steps to remove the duplicate files.

![Screen upper](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/Screen_upper.png?raw=true)

### 1.3 Screening of Lower Lobe Lung Cancer Cases

The filtering method is similar to Upper lobe lung cancer, except that the Diagnoses Tissue or Organ of Origin is changed to lower lobe, lung. here, I filtered 91 files and 81 cases. again, in the next steps I will do a secondary filter to remove duplicate files.

![Screen lower](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/Screen_lower.png?raw=true)

## 2.Data Download and Collation

### 2.1 Installation and Configuration of gdc-client

Gdc-client is a tool used to download files from the GDC website. I went to the [download page of gdc-client](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool) and selected GDC Data Transfer Tool Client's OSX version to download and install. 

Add the path of the software installation to `.zshrc` by adding a line to the .zshrc file

`export PATH="directory path:$PATH"`

Type g`dc-client -version` in the terminal to check if the software is installed successfully.

### 2.2 Download of Count Files

Here I click on the `Manifest` button to download a summary txt file with all the file names. Type the following command in the terminal to download all the files.

`gdc-client download -m gdc_manifest.2021-11-11.txt`

Here, I created two new directories `gdc_upper` and `gdc_lower` to download the files of upper lobe lung cancer and lower lobe lung cancer respectively.

### 2.3 Organizing Count Files

Here, I see that the downloaded files are not count.gz files but folders, so I use `R` to aggregate all the count.gz files into one folder. Taking upper lobe lung cancer as an example, I go to the `gdc_upper` directory in R, create a `gdc_upper_counts` directory, and run.

```
i <- list.dirs()

i
```

Now I can see all the folders in that directory and I find that the number of folders is greater than 156, this is because some of the downloaded folders contain subfolders. I'll ignore these subfolders to put all the files in my newly created `gdc_upper_counts` directory.

```{r}
m = i[2:183]
for(n in m){
  x.path=paste(n,list.files(n),sep='/')
  file.copy(x.path,'./GDC_upper',recursive = T)}
```

Now, organize the files in the `gdc_upper_counts` directory and keep only 156 Counts compressed files. Then do the same for lower lobe lung cancer and now I get two directories `gdc_upper_counts` and `gdc_lower_counts` which contain all the count.gz files I need.

### 2.4 Metadata Json File Download

Add the files of upper lobe lung cancer and lower lobe lung cancer to the `cart` separately, and click the `Metadata `button to download the `json` file. The role of this file is to convert the count file name to the data number of TCGA, which will be used in the secondary screening later.

[Metadata json file of lower](https://github.com/DZBohan/zhangboh_final_project/blob/main/Files/metadata.cart.2021-11-09_lower.json)

[Metadata json file of upper](https://github.com/DZBohan/zhangboh_final_project/blob/main/Files/metadata.cart.2021-11-09_upper.json)

## 3. Secondary Screening of Data

Now, I have obtained the count files for upper lobe lung cancer and lower lobe lung cancer. However, the duplicate files mentioned before still exist in them. I need to find them out and delete them. 

First, I need to know the meaning of TCGA data number.
As shown in the figure, the TCGA data number consists of 7 parts. The third part represents the patient number, so I first need to delete the extra files of the same patient according to the patient number. 

Second, if the number in the third part is greater than 10, it means the sample is a normal sample rather than a tumor sample; so I need to sun all samples with that number greater than 10. 

Third, the letter in the third part indicates the sample quality. I choose to keep only the samples with quality A.

![image](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/TCGA%20data%20number.png?raw=true)

I wrote a bash script to replace the filenames of the previously downloaded counts files with the TCGA data numbers, and compared these files in Finder and recorded the numbers of the files to be deleted. The following is an example of the operation with lower lobe lung cancer.

### 3.1 File Name & TCGA Number Mapping File

Write a file name and TCGA data number mapping file with R

**Lowerlobe mapping file building**

```{r}
meta <- jsonlite::fromJSON("metadata.cart.2021-11-09_lower.json")
View(meta)
```

```{r}
colnames(meta)
ids <- meta$associated_entities; class(ids)
ids[[1]]
ids[[1]][,1]
```

```{r}
ID = sapply(ids, function(x){x[,1]})
file2id_lower = data.frame(file_name = meta$file_name, ID= ID)
View(file2id_lower)
write.table(file2id_lower, file = "sample2id_lower.txt", sep = "\t", col.names = F, quote = F, row.names = F)
```

**Upper mapping file building**

```{r}
meta <- jsonlite::fromJSON("metadata.cart.2021-11-09_upper.json")
View(meta)
```

```{r}
colnames(meta)
ids <- meta$associated_entities; class(ids)
ids[[1]]
ids[[1]][,1]
```

```{r}
ID = sapply(ids, function(x){x[,1]})
file2id_upper = data.frame(file_name = meta$file_name, ID= ID)
View(file2id_upper)
write.table(file2id_lower, file = "sample2id_upper.txt", sep = "\t", col.names = F, quote = F, row.names = F)
```

[Here is the file of lower](https://github.com/DZBohan/zhangboh_final_project/blob/main/Files/3.1_sample2id_lower.txt)

[Here is the file of upper](https://github.com/DZBohan/zhangboh_final_project/blob/main/Files/3.1_sample2id_upper.txt)

### 3.2 Files Batch Rename

Use the mapping file to batch rename files via bash script

**Batch rename lowerlobe files**

```
#!/bin/bash

 cat $1 |while read line
 do
   arr=($line)
   filename=${arr[0]}
   submitterid=${arr[1]}
   gunzip -c ./${filename} > ./file_lower/${submitterid}.count
 done
```
**Batch rename upperlobe files**

```
#!/bin/bash

 cat $1 |while read line
 do
   arr=($line)
   filename=${arr[0]}
   submitterid=${arr[1]}
   gunzip -c ./${filename} > ./file_upper/${submitterid}.count
 done
```
**Bash script usage**

First, create a new `file2id_lower` directory in the `gdc_lower_counts `directory.

Terminal run `bash change_name.sh sample2id_lower.txt`

The files whose names are replaced are stored in the `file2id_lower `directory.

### 3.3 Deletion of Unwanted Files

The TCGA data numbers of the unwanted files were recorded and deleted according to the 3 selection principles described previously. 

22 files were selected from Upper lobe lung cancer that needed to be removed. 10 were selected from Lower lobe lung cancer. Therefore, the final number of Counts files for both is 134:81.

[The file name of the deleted files](https://github.com/DZBohan/zhangboh_final_project/blob/main/Files/3.3_deleted%20files.txt)

[215 count files](https://drive.google.com/drive/folders/1JpWD33umLiS46DZUPZbMImvW5mdi9KvM)

## 4. Load Count Files into Vignette

### 4.1 Prefix Adding

Adding prefix to the count files using a bash script

Because the files need to be loaded together when loading into the vignette. I use a bash script to add the prefix `upperlobe- `and `lowerlobe-` to the count files of upper lobe lung cancer and lower lobe lung cancer respectively. For example: `upperlobe-TCGA-NJ-A55R-01A-11R-A262-07.count`. This way I can put both sets of count files into the same directory `file_all` and use `regular expressions`  to take out "upperlobe" and "lowerlobe" from the file names as conditions.

**The script of lowerlobe**

```
#!/bin/sh
for files in $(ls *.count)
    do mv $files "lowerlobe-"$files
done
```
**The script of upperlobe**

```
#!/bin/sh
for files in $(ls *.count)
    do mv $files "upperlobe-"$files
done
```
### 4.2 Loading Count Files

Here, I am going to Load all count files into vignette

I used `R` to do the loading of the files. First I created a value and saved the path to the file_all directory in it. 

Import the filenames of the files in the directory into the `sampleFiles` value. 

Use a regular expression to get the upperlobe or lowerlobe of the file name into the `sampleCondition` value.

Create a dataframe named `sampleTable` with three columns for `sampleName`, `fileName` and `condition`.

Using `DESeq2` package, enter all count files into vignette and create `dds`.

```{r}
directory <- "~/myproject/file_all"
```

```{r}
sampleFiles <- grep("lobe",list.files(directory),value=TRUE)
```

```{r}
sampleCondition <- sub("(.*lobe).*","\\1",sampleFiles)
```

```{r}
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)
sampleTable$condition <- factor(sampleTable$condition)
```

```{r}
library("DESeq2")
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
dds
```
[Here is the sampleTable](https://github.com/DZBohan/zhangboh_final_project/blob/main/Files/4.2_sampleTable.csv)

**Information of dds**

![dds result](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/dds.png?raw=true)

![dds](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/dds_2.png?raw=true)

## 5. Generation of Differential Expression Results

### 5.1 Processing of dds Dataframes

#### Pre-filtering

The aim is to remove low-count genes. Here, I keep genes with at least 10 counts.

```{r}
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
```

#### Note on factor levels

The aim is to set a factor for comparing differences in gene expression, i.e. upperlobe lung cancer and lowerlobe lung cancer.

```{r}
dds$condition <- factor(dds$condition, levels = c("upperlobe","lowerlobe"))
```

```{r}
dds$condition <- droplevels(dds$condition)
```

### 5.2 Obtaining Dataframes

Here, I am going to obtain dataframes for differential gene expression results.

#### Get Results

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

I generated an [original res file](https://raw.githubusercontent.com/DZBohan/zhangboh_final_project/main/Files/original_res.csv) here.

```{r}
getwd()
write.csv(as.data.frame(res), 
          file="original_res.csv")
```

### 5.3 Gene Names Changing

Here, I am going to chane the ensembl id in the reult table to gene name.

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

#### Replace ensembl id with gene id

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

### 5.4 Optimization

#### Log fold change (LFC) shrinkage

I use the `apeglm` method for LFC shrinkage which is useful for gene visualization and ranking.

```{r}
resultsNames(dds)
```

![5-4](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/5-4.png?raw=true)

```{r}
resLFC <- lfcShrink(dds, coef="condition_lowerlobe_vs_upperlobe", type="apeglm")
resLFC
```

![5-5](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/5-5.png?raw=true)

#### Parallelization

Splitting the work into 4 cores to speed it up.

```{r}
library("BiocParallel")
register(MulticoreParam(4))
```

### 5.5 The Number of Differentially Expressed Genes

Here, I am going to view the number of differentially expressed genes based on p-values.

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

## 6. Results Exploring & Exporting

### 6.1 MA-plot

Use the code below to generate an MA-plot of the result, where the horizontal coordinates represent the mean value of the gene expressed in all samples and the vertical coordinates represent the fold difference in expression of the gene in the two samples. The blue points represent the differentially expressed genes. It can be seen that the mean values of differentially expressed genes in the samples are kept at relatively high levels so that the results are good.

```{r}
plotMA(res, ylim=c(-2,2))
```

![6-1](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/6-1.png?raw=true)

Now, using the shrunken results to generate the MA-plot again, theoretically, the differential gene expression ploidy should get shrunk, but in fact it does not. Therefore, I put this issue in Known issue. The good thing is that the M-plot I generated directly with the results is fine.

```{r}
plotMA(resLFC, ylim=c(-2,2))
```

![6-2](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/6-2.png?raw=true)

```{r}
resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")
```
I used the other three methods to try to shrink the resulting data and generate the MA-plot, and found that the MA-plot generated by shrinking the results in all three ways was not as good as the original results.

![6-3](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/6-3.png?raw=true)

```{r}
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
```

![6-4](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/6-4.png?raw=true)

### 6.2 Plot Counts

I found the genes associated with lung cancer on [NCBI](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3527990/), and they are *EGFR*, *KRAS*, *MET*, *LKB1*, *BRAF*, *PIK3CA*, *ALK*, *RET*, and *ROS1*. I found the Ensembl IDs corresponding to these genes from the file [Ensembl ID TO Genename.csv](https://raw.githubusercontent.com/DZBohan/zhangboh_final_project/main/Files/Ensembl_ID_TO_Genename.csv) I just generated. I selected 5 of these genes to generate plot counts. Then, I found their Ensembl ID with version number from the [original res table](https://raw.githubusercontent.com/DZBohan/zhangboh_final_project/main/Files/original_res.csv).

![lung cancer genes](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/lung_cancer_genes.png?raw=true)

```{r}
d <- plotCounts(dds,"ENSG00000146648.14", intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400)) +
ggtitle("EGFR") +
theme(plot.title = element_text(hjust = 0.5))

d <- plotCounts(dds,"ENSG00000133703.10", intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400)) +
ggtitle("KRAS") +
theme(plot.title = element_text(hjust = 0.5))

d <- plotCounts(dds,"ENSG00000157764.11", intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400)) +
ggtitle("BRAF") +
theme(plot.title = element_text(hjust = 0.5))

d <- plotCounts(dds,"ENSG00000121879.3", intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400)) +
ggtitle("PIK3CA") +
theme(plot.title = element_text(hjust = 0.5))

d <- plotCounts(dds,"ENSG00000047936.9", intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400)) +
ggtitle("ROS1") +
theme(plot.title = element_text(hjust = 0.5))
```

![EGFR](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/EGFR.png?raw=true)

![KRAS](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/KRAS.png?raw=true)

![BRAF](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/BRAF.png?raw=true)

![PIK3CA](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/PIK3CA.png?raw=true)

![ROS1](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/ROS1.png?raw=true)

The plots show that the expression of these lung cancer-associated genes do not differ significantly in lung cancers occurring in the upper and lower lobes of the lung.

### 6.3 More Information

```{r}
mcols(res)$description
```

![6-7](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/6-7.png?raw=true)

### 6.4 CSV Files Writing

#### Original table of all genes

```{r}
getwd()
write.csv(as.data.frame(res_g), 
          file="res_g.csv")
```

[Here is the file res_g.csv](https://raw.githubusercontent.com/DZBohan/zhangboh_final_project/main/Files/res_g.csv)

#### Differential expression genes

Genes with p-values <0.05 were screened and determined as differentially expressed genes. And a new column was created to record the up- or down-regulation of genes.

```{r}
resSig_0.05 <- subset(resOrdered, padj < 0.05)
resSig_0.05[which(resSig_0.05$log2FoldChange > 0), "up_down"] <- "up"
resSig_0.05[which(resSig_0.05$log2FoldChange < 0), "up_down"] <- "down"
resSig_up <- subset(resSig_0.05, log2FoldChange > 0)
resSig_down <- subset(resSig_0.05, log2FoldChange < 0)
```

![6-8](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/6-8.png?raw=true)

This table contains the 433 differentially expressed genes I screened and their up- or down-regulation information. Among the 433 differentially expressed genes, there are 345 up-regulated expressed genes and 88 down-regulated expressed genes.

```{r}
write.csv(as.data.frame(resSig_p0.05), 
          file="differential_expression_all.csv")
write.csv(as.data.frame(resSig_up), 
          file="differential_expression_up.csv")
write.csv(as.data.frame(resSig_down), 
          file="differential_expression_down.csv")
```

[Here is the file differential expression all.csv](https://github.com/DZBohan/zhangboh_final_project/blob/main/Files/differential_expression_all.csv)

[Here is the file differential expression down.csv](https://github.com/DZBohan/zhangboh_final_project/blob/main/Files/differential_expression_down.csv)

[Here is the file differential expression up.csv](https://github.com/DZBohan/zhangboh_final_project/blob/main/Files/differential_expression_down.csv)

## 7. Data Transformations & Visualization

### 7.1 Extracting Transformed Values

```{r}
vsd <- vst(dds, blind=FALSE)
```

### 7.2 Effects of Transformations

Here, I am going to see the effects of transformations on the variance.

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

As shown in the plots, the data transformation has little effect on the variance of the sample.

## 8 Data Quality Assessment

Here, I am going to assess the data quality by sample clustering and visualization.

### 8.1 Heatmap of the Count Matrix

#### Heatmap of ntd

```{r}
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:100]
df <- as.data.frame(colData(dds)[,c("condition", "sizeFactor")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, show_colnames = FALSE)
```

![8-1](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/8-1.png?raw=true)

I plotted heatmaps to show the differences in expression between the different subgroups, but no significant differences can be seen in the plots.

#### Heatmap of vsd

```{r}
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, show_colnames = FALSE)
```

![8-2](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/8-2.png?raw=true)

### 8.2 Heatmap of the Sample-to-sample Distances

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

The distances heatmap shows the similarity between samples.

### 8.3 Principal Component Plot

The PCA plot clusters the samples and can determine whether the clustering effect of the samples is consistent with the grouping set during the experimental design.

```{r}
plotPCA(vsd, intgroup="condition")
```

![8-4](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/8-4.png?raw=true)

From the plot, the samples are not clustered according to upperlobe and lowerlobe. It indicates that there are other factors that dominate the differential expression of genes in the sample.

Now, I need to find this factor that dominates the difference in gene expression of the samples. Assume that this factor is gender. I re-filter the count file from the TCGA database, this time I add another filter condition `male` and keep the other filters the same.

I filtered to get 80 upperlobe count files and 47 lowerlobe count files this time. since it was only used for validation, I did not filter the count files again. Repeating the procedure in the paper, these 127 files were entered into the vignette and a PCA plot was generated.

![pca_gender_check](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/PCA_gender_check.png?raw=true)

In the newl PCA plot, there is still another dominant factor in the gene expression differences of the samples at one location. So this factor is not gender. I recorded the issue in `Known issue`, and I can continue to explore and verify the dominant factor in the future.

## 9. Evaluation of Differentially Expressed Genes

I use the websit [Gene Set Enrichment Analysis(GSEA)](https://www.gsea-msigdb.org/gsea/msigdb/annotate.jsp) to find overlap of the differentially expressed genes in gene sets. Finally, I find 56 gene sets that have overlap of the differentially expressed genes. 52 have up-regulated genes, and 4 have down-regulated genes.

[Here is the table of gene sets that have overlap of the differentially expressed genes
](https://github.com/DZBohan/zhangboh_final_project/blob/main/Files/gene_sets_overlap.csv)

I choose the 9 gene sets that contain the most overlap to introduce them.

### `GOMF_SIGNALING_RECEPTOR_BINDING`

This gene sets have 35 overlap up-regulated genes. Interacting selectively and non-covalently with one or more specific sites on a receptor molecule, a macromolecule that undergoes combination with a hormone, neurotransmitter, drug or intracellular messenger to initiate a change in cell function.

### `GOBP_RESPONSE_TO_ENDOGENOUS_STIMULUS`

This gene sets have 31 overlap up-regulated genes. Any process that results in a change in state or activity of a cell or an organism (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a stimulus arising within the organism.

### `GOBP_REGULATION_OF_TRANSPORT`

This gene sets have 31 overlap up-regulated genes. Any process that modulates the frequency, rate or extent of the directed movement of substances (such as macromolecules, small molecules, ions) into, out of or within a cell, or between cells, by means of some agent such as a transporter or pore.

### `GOBP_CELL_CELL_SIGNALING`

This gene sets have 30 overlap up-regulated genes. Any process that mediates the transfer of information from one cell to another. This process includes signal transduction in the receiving cell and, where applicable, release of a ligand and any processes that actively facilitate its transport and presentation to the receiving cell.  Examples include signaling via soluble ligands, via cell adhesion molecules and via gap junctions.

### `BENPORATH_ES_WITH_H3K27ME3`

This gene sets have 29 overlap up-regulated genes. Set 'H3K27 bound': genes posessing the trimethylated H3K27 (H3K27me3) mark in their promoters in human embryonic stem cells, as identified by ChIP on chip.

### `GOBP_SECRETION`

This gene sets have 28 overlap up-regulated genes. The controlled release of a substance by a cell or a tissue.

### `YOSHIMURA_MAPK8_TARGETS_UP`

This gene sets have 25 overlap up-regulated genes. Genes up-regulated in vascular smooth muscle cells (VSMC) by MAPK8 (JNK1) [GeneID=5599].

### `GOBP_RESPONSE_TO_NITROGEN_COMPOUND`

This gene sets have 24 overlap up-regulated genes. Any process that results in a change in state or activity of a cell or an organism (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a nitrogen compound stimulus.

### `MIKKELSEN_MEF_HCP_WITH_H3K27ME3`

This gene sets have 22 overlap up-regulated genes. Genes with high-CpG-density promoters (HCP) bearing histone H3 trimethylation mark at K27 (H3K27me3) in MEF cells (embryonic fibroblast).

## 10. Conclusion

After screening, I found 215 cases of `Stage 1` lung cancer in the TCGA database, while I controlled for Race as `white` and Ethnicity as `not hispanic or latino`. 81 of these cases had cancer in the `lower lobe lung` and 134 patients had cancer in the `upper lobe lung`. The aim of this project was to look at the differences in gene expression between lung cancers occurring in the upper and lower lobes of the lung. By using Deseq2 and setting an adjusted p-value <0.05, I found 433 differentially expressed genes. Among these genes, there were 345 genes with up-regulated expression and 88 genes with down-regulated expression. Evaluating them in GSEA revealed that the up-regulated genes had overlaps in 52 gene sets and overlaps of more than 20 in 14 gene sets, suggesting that lung cancers occurring in the upper and lower lobes of the lung do have differential gene expression, especially up-regulated genes, and some of these genes affect the function and signaling of cells. In contrast, the overlaps of down-regulated genes in gene sets is not much.

## 11 Known Issues

### Resulting shrinkage and MA-plot

The MA plot plotted using the shrunken results is not satisfactory and further search for the cause is needed.


### The rlg method of data conversion

The rlg data transformation method is not applicable to the case of large sample size. Therefore, I choose to use vst's data conversion method only.### Heat mapping using differentially expressed genes

### Screening for differential genes using different criteria

In the screening of differentially expressed genes, my screening condition is p-value < 0.05, and if the differentially expressed gene pool needs to be further narrowed, the Log2change screening condition can be added.

### Another factor dominating differential gene expression in the samples

From the PCA plot, we know that there are other factors dominating differential gene expression in the samples. It was verified that this factor was not sex. Further exploration and validation is needed to find this dominant factor in the future.

### <center>[Back to Top](#jump)</center>
