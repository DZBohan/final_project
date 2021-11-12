# Milestone 1
## Section 1 Update
### 1.Initial screening of data
This step was completed in the previous stage, but I have documented the screening process to better reproduce the entire step.
#### 1.1 Screening of Files
Data Category <- transcriptome profiling

Experimental Strategy <- RNA-Seq

Workflow Type <- HTSeq - Counts

Access <- open
#### 1.2 Screening of Upper Lobe Lung Cancer Cases
Diagnoses Ajcc Pathologic Stage <- stage ia/stage ib/stage i

Diagnoses Tissue or Organ of Origin <- upper lobe, lung

Primary Site <- bronchus and lung

Program <- TCGA

Vital Status <- alive

Race <- white

Ethnicity <- not hispanic or latino

Now, I have filtered out 156 files and 136 cases. The number of files and cases is different because some cases have duplicate files; however, since I will be downloading files, the filtering in this step is incomplete. The second filtering will be done in the later steps to remove the duplicate files.
#### 1.2 Screening of Lower Lobe Lung Cancer Cases
The filtering method is similar to Upper lobe lung cancer, except that the Diagnoses Tissue or Organ of Origin is changed to lower lobe, lung. here, I filtered 91 files and 81 cases. again, in the next steps I will do a secondary filter to remove duplicate files.
### 2.Data download and collation
#### 2.1 Installation and Configuration of gdc-client
Gdc-client is a tool used to download files from the GDC website. I went to the [download page of gdc-client](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool) and selected GDC Data Transfer Tool Client's OSX version to download and install. 

Add the path of the software installation to `.zshrc` by adding a line to the .zshrc file

`export PATH="directory path:$PATH"`

Type g`dc-client -version` in the terminal to check if the software is installed successfully.

#### 2.1 Download of Count files
Here I click on the `Manifest` button to download a summary txt file with all the file names. Type the following command in the terminal to download all the files.

`gdc-client download -m gdc_manifest.2021-11-11.txt`

Here, I created two new directories `gdc_upper` and `gdc_lower` to download the files of upper lobe lung cancer and lower lobe lung cancer respectively.

#### 2.2 Organizing count files
Here, I see that the downloaded files are not count.gz files but folders, so I use `R` to aggregate all the count.gz files into one folder. Taking upper lobe lung cancer as an example, I go to the `gdc_upper` directory in R, create a `gdc_upper_counts` directory, and run.

`i <- list.dirs()`

`i`

Now I can see all the folders in that directory and I find that the number of folders is greater than 156, this is because some of the downloaded folders contain subfolders. I'll ignore these subfolders to put all the files in my newly created `gdc_upper_counts` directory.

[Here is the code](https://github.com/DZBohan/zhangboh_final_project/blob/main/scripts/2.2_organization.Rmd)

Now, organize the files in the `gdc_upper_counts` directory and keep only 156 Counts compressed files. Then do the same for lower lobe lung cancer and now I get two directories `gdc_upper_counts` and `gdc_lower_counts` which contain all the count.gz files I need.

#### 2.3 Metadata json file download
Add the files of upper lobe lung cancer and lower lobe lung cancer to the `cart` separately, and click the `Metadata `button to download the `json` file. The role of this file is to convert the count file name to the data number of TCGA, which will be used in the secondary screening later.

[Metadata json file of lower](https://github.com/DZBohan/zhangboh_final_project/blob/main/Files/metadata.cart.2021-11-09_lower.json)

[Metadata json file of upper](https://github.com/DZBohan/zhangboh_final_project/blob/main/Files/metadata.cart.2021-11-09_upper.json)
### 3. Secondary screening of data
Now, I have obtained the count files for upper lobe lung cancer and lower lobe lung cancer. However, the duplicate files mentioned before still exist in them. I need to find them out and delete them. 

First, I need to know the meaning of TCGA data number.
As shown in the figure, the TCGA data number consists of 7 parts. The third part represents the patient number, so I first need to delete the extra files of the same patient according to the patient number. 

Second, if the number in the third part is greater than 10, it means the sample is a normal sample rather than a tumor sample; so I need to sun all samples with that number greater than 10. 

Third, the letter in the third part indicates the sample quality. I choose to keep only the samples with quality A.

![image](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/TCGA%20data%20number.png)

I wrote a bash script to replace the filenames of the previously downloaded counts files with the TCGA data numbers, and compared these files in Finder and recorded the numbers of the files to be deleted. The following is an example of the operation with lower lobe lung cancer.

#### 3.1 Write a file name and TCGA data number mapping file with R
[Here is the code](https://github.com/DZBohan/zhangboh_final_project/blob/main/scripts/3.1_filename2TCGA.Rmd)

[Here is the file of lower](https://github.com/DZBohan/zhangboh_final_project/blob/main/Files/3.1_sample2id_lower.txt)

[Here is the file of upper](https://github.com/DZBohan/zhangboh_final_project/blob/main/Files/3.1_sample2id_upper.txt)
#### 3.2 Use this mapping file to batch rename files via bash script
[Here is the code of lower](https://github.com/DZBohan/zhangboh_final_project/blob/main/scripts/3.2_sample2id_rename_lower.sh)

[Here is the code of upper](https://github.com/DZBohan/zhangboh_final_project/blob/main/scripts/3.2_sample2id_rename_upper.sh)

Bash script usage.

First, create a new `file2id_lower` directory in the `gdc_lower_counts `directory.

Terminal run `bash change_name.sh sample2id_lower.txt`

The files whose names are replaced are stored in the `file2id_lower `directory.
#### 3.3 Deletion of unwanted files
The TCGA data numbers of the unwanted files were recorded and deleted according to the 3 selection principles described previously. 

22 files were selected from Upper lobe lung cancer that needed to be removed. 10 were selected from Lower lobe lung cancer. So, the final number of Counts files for both is 134:81.

[The file name of the deleted files](https://github.com/DZBohan/zhangboh_final_project/blob/main/Files/3.3_deleted%20files.txt)

### 4.Load the Counts file into the vignette
####4.1 Adding prefix to the count files using a bash script
Because the files need to be loaded together when loading into the vignette. I use a bash script to add the prefix `upperlobe- `and `lowerlobe-` to the count files of upper lobe lung cancer and lower lobe lung cancer respectively. For example: `upperlobe-TCGA-NJ-A55R-01A-11R-A262-07.count`. This way I can put both sets of count files into the same directory `file_all` and use `regular expressions`  to take out "upperlobe" and "lowerlobe" from the file names as conditions.

[Here is the code of lower](https://github.com/DZBohan/zhangboh_final_project/blob/main/scripts/4.1_add_prefix%20_lower.sh)

[Here is the code of upper](https://github.com/DZBohan/zhangboh_final_project/blob/main/scripts/4.1_add_prefix%20_upper.sh)
#### 4.2 Loading all count files into vignette
I used `R` to do the loading of the files. First I created a value and saved the path to the file_all directory in it. 

Import the filenames of the files in the directory into the `sampleFiles` value. 

Use a regular expression to get the upperlobe or lowerlobe of the file name into the `sampleCondition` value.

Create a dataframe named `sampleTable` with three columns for `sampleName`, `fileName` and `condition`.

Using `DESeq2` package, enter all count files into vignette and create `dds`.

[Here is the code](https://github.com/DZBohan/zhangboh_final_project/blob/main/scripts/4.2_htseq_input.Rmd)

[sampleTable](https://github.com/DZBohan/zhangboh_final_project/blob/main/Files/4.2_sampleTable.csv)

![dds](https://github.com/DZBohan/zhangboh_final_project/blob/main/Images/dds.png)

### 5.Read the counts file and create dataframe (Optional)
In the class on November 11, the professor said not to use the method; however, for the purpose of learning, I recorded the operation.

for upper lobe lung cancer and lower lobe lung cancer, respectively.
#### 5.1 Using R to build dataframe.

#### 5.2 Delete unnecessary rows and columns.

#### 5.3 Convert Counts file name to TCGA data number and remove the version number of the gene Ensembl number
Here I need to use the previously downloaded `json file` and install the `rjson` package.

#### 5.4 Change Ensembl number of gene in dataframe to gene id

[Ensembl to genename file](https://github.com/DZBohan/zhangboh_final_project/blob/main/Files/5_Ensembl_ID_TO_Genename.csv)
#### 5.5 Removal of id duplicate genes

#### 5.4 Deleting unwanted count files in dataframe
[Here is the code of 5 in lower](https://github.com/DZBohan/zhangboh_final_project/blob/main/scripts/5_creat_dataframe.Rmd)

[Here is the final datafram of lower lobe](https://github.com/DZBohan/zhangboh_final_project/blob/main/Files/5_datafram_lower.csv)
## Section 2
An initial completion of vignette. I will complete an entire first draft of analysis analyzed through the vignette. After that, I will send the results of the first draft to the professor for feedback.
