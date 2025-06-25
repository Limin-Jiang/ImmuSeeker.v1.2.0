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

## parameters
filenames = options$Imputname 


#setwd("D:\\Projects\\HLA\\ImmuSeeker\\HLA_autoimmune\\results")

d1 = fread(paste0(filenames,".candidate.csv"),header = TRUE,sep = ",")
colnames(d1)[1] = "name"
subd1 = d1%>%
  filter(Level == "One-field")
aa = table(subd1$Gene) %>%
  as.data.frame()%>%
  filter(Freq == 1)
subd1 = subd1 %>%
  filter(Gene %in%aa$Var1 )
idx = match(subd1$HLAname,d1$HLAname)

d1$TPM[idx] = d1$TPM[idx]/2


data_all = data.frame()


getrich = function(d1_temp){
  d1_temp = d1_temp%>%
    select("name","TPM")%>%
    t()%>%
    as.data.frame()
  colnames(d1_temp) = d1_temp[1,]
  d1_temp = d1_temp[-1,]
  d1_temp <- data.frame(lapply(d1_temp, function(x) as.numeric(as.character(x))))
  d1_temp = round(d1_temp)
  
  d1_temp = t(d1_temp)
  OTU = otu_table(d1_temp, taxa_are_rows = TRUE)
  richness <- estimate_richness(OTU,measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"))
  
  return(richness)
  
}




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
}

d1_temp = d1%>%
  filter(Level == "Two-field")
if (nrow(d1_temp)>=3){ 
  richness = getrich(d1_temp)
  data = data.frame(Level = "Two-field", Diversity = colnames(richness), Value = t(richness[1,]) )
	rownames(data) <- NULL
	data_all = rbind(data_all,data)
	
}

if (nrow(data_all) > 0){
	colnames(data_all)[3] = "Value"
	write.csv(data_all,paste0(filenames,"_diversity_ex.csv"),row.names = FALSE)
}else{
  print("No HLA types were detected based on the provided conditions.")
}

