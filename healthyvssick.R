
## Specialized packages install code (if needed)

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")}

BiocManager::install("DESeq2")
BiocManager::install("fgsea")
BiocManager::install("clusterProfiler")

## libraries 

library(readxl)
library(DESeq2)
library(fgsea)
library(data.table)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(DT)

## load in the data, clean to proper format

setwd("C:\\Users\\Pat\\Downloads")
CountData2 <- read_excel("HealthyVsSick(c).xlsx")
CountData2 <- as.data.frame(CountData2)
CountData2 <- CountData2 %>% mutate_if(is.double, as.integer)


metaData2 <- read_excel("metaData(healthyvssick).xlsx")
metaData2 <- as.data.frame(metaData2)

## create the dds2 dataframw 
dds2 <- DESeqDataSetFromMatrix(countData=CountData2, 
                              colData=metaData2,
                              design=~dex, tidy = TRUE)

## Run the differential sequence model 
dds2 <- DESeq(dds2)

## See results
res2 <- results(dds2)
head(results(dds2, tidy=TRUE))
summary(res2)

SymbolSickvsNull <- read_excel("SymbolSickvsNull.xlsx")


resnew<-cbind(res2,SymbolSickvsNull)

resnew<-as.data.frame(resnew)
res2 <- resnew %>% 
  dplyr::select(Symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(Symbol) %>% 
  summarize(stat=mean(stat))
res2

library(fgsea)
library(tibble)

res2 <- as.data.frame(res2)
ranks <- res2[order(res2$stat,decreasing=TRUE),]
ranks <- deframe(ranks)
head(ranks, 20)


pathways.hallmark <- gmtPathways("C:\\Users\\Pat\\Documents\\custom(2.0).GMT")

# Look at them all if you want (uncomment)
# pathways.hallmark

# Show the first few pathways, and within those, show only the first few genes. 
pathways.hallmark %>% 
  head() %>% 
  lapply(head)

# Run FGSEA
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks,)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()







