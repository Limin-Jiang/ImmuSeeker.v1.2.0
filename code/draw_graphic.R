#   ImmuSeeker
#
#   Authors:    Limin Jiang (lxj423@med.miami.edu)
#   Date:       June 1, 2024
#   Version:    0.1.0

rm(list = ls())
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggraph)
  library(igraph)
  library(docopt)
})
#filenames = "P001D00.bam"
#threhold_num = 0
doc <- "
Usage: myscript.R [-i <Imputname>] 
Options:
  -i, --Imputname   The file of list of cancers
"
options <- docopt(doc)


## parameters
filenames = options$Imputname 
source("./code/main_funs.R")

d1 = fread(paste0(filenames,".candidate.csv"),header = TRUE,sep = ",") %>%
  filter(Level != "Gene")
colnames(d1)[2] = "Num"
d1$ID = sapply(strsplit(d1$AlleleID, split=':', fixed=TRUE), function(x) (x[1]))


if (nrow(d1) > 1){
  genes = unique(d1$Gene)
  pdf(paste0(filenames,"_evolution_graphs.pdf"), width = 6, height = 4.5) 
  for (gene_ID in genes){
	  print(gene_ID)
    subd1 = d1%>%
      filter(Gene == gene_ID)
    
    data_all = getdata(subd1,gene_ID)
    
    mygraph <- graph_from_data_frame( data_all )
    
    plot <- ggraph(mygraph, layout = 'auto', circular = FALSE) + 
      geom_edge_diagonal(strength = 1) +
      geom_node_point(color = 'grey', size = 1)+
      geom_node_text(aes(label=name) ,colour = "blue", angle=90 ,size = 1.5, hjust=0.2,vjust=-0.7) +
      theme_void()+
      ggtitle(paste0("HLA/KIR ",gene_ID))
    
    print(plot)
  }
}
dev.off()

