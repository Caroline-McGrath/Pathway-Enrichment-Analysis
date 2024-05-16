


## install packages if needed 
BiocManager::install("KEGGREST")

## website with pathways id lists : https://www.genome.jp/kegg/pathway.html#disease

library("KEGGREST")

#Replace the number after 'hsa' in this line with the id of the pathway you want 
names <- keggGet("hsa05218")[[1]]$GENE

# strip it of excess info 
namesodd <-  names[seq(0,length(names),2)]
namestrue <- gsub("\\;.*","",namesodd)

# return list of genes in pathway
namestrue




