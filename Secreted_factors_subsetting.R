library(readxl)
library(tidyverse)
library(ggplot2)
library(ggtext)
library(ggrepel)
y=read.table("NexCre_Foxg1cKO_resCutOff.txt")
up=read_excel("Upreg_NexCre_Foxg1_Subcellular_loc.xlsx")
down=read_excel("Downreg_NexCre_Foxg1_Subcellular_loc.xlsx")

# Subset genes with subcelluler localization as "Secreted"
up_filtered = up %>% filter(grepl("Secreted", up$`Subcellular location [CC]`))
down_filtered = down %>% filter(grepl("Secreted", down$`Subcellular location [CC]`))

# Select only secreted factors from the RNA seq list 
secreted_p_up = y %>% inner_join(up_filtered, by=c("Gene"="From")) %>% 
  select(Gene:diffexpressed) %>% unique()

# Select only secreted factors from the RNA seq list 
secreted_p_down = y %>% inner_join(down_filtered, by=c("Gene"="From")) %>% 
  select(Gene:diffexpressed) %>% unique()
# Merge both sheets into one
secreted_p=rbind(secreted_p_up, secreted_p_down)
write.csv(secreted_p, file = "Secreted_protein_NexCre_Foxg1cKO.csv")

setwd("D:/ChIP seq and RNA seq analysis/P1_NexCre_Foxg1cKO_RNA_Seq/")
#resDF=read.table("NexCre_Foxg1cKO_RNASeq_results.txt")
resDF=read_excel("Secreted protein NexCre Foxg1cKO.xlsx")
resDF$diffexpressed <- "NO"
resDF$diffexpressed[resDF$log2FoldChange > 0.58 & resDF$padj < 0.05] <- "Upregulated"
resDF$diffexpressed[resDF$log2FoldChange < -0.58 & resDF$padj < 0.05] <- "Downregulated"
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("Downregulated", "Upregulated", "NO")
resDF$Label <- resDF[,1]
resDF$Label[resDF$diffexpressed=="NO"] <- "NA"
genes <- c("Fgf9", "Fgf18", "Bmp8b", "Bmp4", "Bmp6", "Wnt2", "Wnt5a")

p4 <- ggplot(resDF, aes(log2FoldChange, -log10(padj), col = diffexpressed, 
                        label=ifelse(Label %in% genes, as.character(Label), NA))) + 
  geom_point(alpha=ifelse(resDF$Label %in% genes, 1, 0.5)) + 
  theme_classic() + scale_colour_manual(values = mycolors) +
  geom_text_repel(color="black",
                  direction = "y", 
                  box.padding = 0.1, 
                  force = 10, 
                  nudge_y = 10, 
                  nudge_x = ifelse(resDF$log2FoldChange>0, 8, -8),
                  size=4) +  
  geom_vline(xintercept = 0, col="black") 
p4


