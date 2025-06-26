#   ImmuSeeker
#
#   Authors:    Limin Jiang (lxj423@med.miami.edu)
#   Date:       June 1, 2024
#   Version:    0.1.0

rm(list = ls())
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(phyloseq)
  library(docopt)
})
doc <- "
Usage: myscript.R [-i <Imputname>]

Options:
  -i, --Imputname   The file of list of cancers
"
options <- docopt(doc)

## parameters filenames= "P001D00.bam"
filenames = options$Imputname 


#setwd("D:\\Projects\\HLA\\ImmuSeeker\\HLA_autoimmune\\results")

source("./code/main_funs.R")
d1 = fread(paste0(filenames,".candidate.csv"),header = TRUE,sep = ",") 
colnames(d1)[2] = "Num"
d1$ID = sapply(strsplit(d1$AlleleID, split=':', fixed=TRUE), function(x) (x[1]))



data_all = data.frame()
d1_temp = d1%>%
  filter(Level == "Gene")  
if (nrow(d1_temp)>=3){
  richness = getrich(d1_temp)
  data = data.frame(Level = "Gene", Diversity = colnames(richness), Value = t(richness[1,]) )
  rownames(data) <- NULL
  data_all = rbind(data_all,data)	  
}





d1_temp = d1%>%
  filter(Level == "One-field")
if (nrow(d1_temp)>=3){ 
  richness = getrich(d1_temp)
  data = data.frame(Level = "One-field", Diversity = colnames(richness), Value = t(richness[1,]) )
  rownames(data) <- NULL
  data_all = rbind(data_all,data)
}else{
  print("No HLA/KIR one-field allele was detected based on the provided conditions.")
}


d1_temp = d1%>%
  filter(Level == "Two-field")
if (nrow(d1_temp)>=3){ 
  richness = getrich(d1_temp)
  data = data.frame(Level = "Two-field", Diversity = colnames(richness), Value = t(richness[1,]) )
  rownames(data) <- NULL
  data_all = rbind(data_all,data)
}else{
  print("No HLA/KIR two-field allele was detected based on the provided conditions.")
}

if (nrow(data_all) > 0){
  colnames(data_all)[3] = "Value"
  write.csv(data_all,paste0(filenames,"_diversity_re.csv"),row.names = FALSE)
}else{
  print("No HLA/KIR alleles was detected based on the provided conditions.")
}

