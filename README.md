# Subsetting secreted factors from a list of genes obtain from RNA Sequencing

## Load libraries
Needed libraries for wrangling data and plotting graphs.
```
library(readxl)
library(tidyverse)
library(ggplot2)
library(ggtext)
library(ggrepel)
```
## Read datasheets
The data sheets of all genes obrained from RNA seq data were queried with UniProt database and the information for their subcellular localization was obtained. 
```
y=read.table("NexCre_Foxg1cKO_resCutOff.txt")
up=read_excel("Upreg_NexCre_Foxg1_Subcellular_loc.xlsx")
down=read_excel("Downreg_NexCre_Foxg1_Subcellular_loc.xlsx")
```
## Subset genes with subcelluler localization as "Secreted"
Since we needed to identify genes which are secreted, we subset the data selecting only secreted molecules from the "subcellular localization" column. 
```
up_filtered = up %>% filter(grepl("Secreted", up$`Subcellular location [CC]`))
down_filtered = down %>% filter(grepl("Secreted", down$`Subcellular location [CC]`))
```
## Select only secreted factors from the RNA seq list 
The list was merged with the original RNA seq dataset to retain only the common elements.
```
secreted_p_up = y %>% inner_join(up_filtered, by=c("Gene"="From")) %>% 
  select(Gene:diffexpressed) %>% unique()
```
## Select only secreted factors from the RNA seq list 
```
secreted_p_down = y %>% inner_join(down_filtered, by=c("Gene"="From")) %>% 
  select(Gene:diffexpressed) %>% unique()
```  
## Merge both sheets into one

Since we queried the upregulated and downregulated genes speperately, we merged the final sheets into one.
```
secreted_p=rbind(secreted_p_up, secreted_p_down)
write.csv(secreted_p, file = "Secreted_protein_NexCre_Foxg1cKO.csv")
```     
