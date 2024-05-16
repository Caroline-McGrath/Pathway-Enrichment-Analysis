
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

    
## load in the data, clean to proper format

setwd("C:\\Users\\Pat\\Downloads")
countData <- read_excel("ControlVsMenthol(c).xlsx")
countData <- as.data.frame(countData)
countData <- countData %>% mutate_if(is.double, as.integer)


metaData <- read_excel("metadata(controlvsmenthol).xlsx")
metaData <- as.data.frame(metaData)


#duplicate_rows <- countData[duplicated(countData[,1]) | duplicated(countData[,1], fromLast = TRUE), ]

# Remove duplicate rows
#countData <- countData[!duplicated(countData[,1]) & !duplicated(countData[,1], fromLast = TRUE), , drop = FALSE]




    ## create the dds dataframe
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData,
                              design=~dex, tidy = TRUE)
    ## Run the differential sequence model 
dds <- DESeq(dds)
    ## See results
res <- results(dds)
head(results(dds, tidy=TRUE))


##read in gene names
SymbolControlsvsMenthol <- read_excel("SymbolControlsvsMenthol.xlsx")


resnew<-cbind(res,SymbolControlsvsMenthol)

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

ranks <- deframe(res2)
head(ranks, 20)

pathways.hallmark <- gmtPathways("C:\\Users\\Pat\\Documents\\custom(2.0).GMT")

# Look at them all if you want (uncomment)
# pathways.hallmark

# Show the first few pathways, and within those, show only the first few genes. 
pathways.hallmark %>% 
  head() %>% 
  lapply(head)

# Run FGSEA
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>% 
  # dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip()
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

