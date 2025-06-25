#   ImmuSeeker
#
#   Authors:    Limin Jiang (lxj423@med.miami.edu)
#   Date:       June 1, 2024
#   Version:    0.1.0
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(docopt)
})
doc <- "
Usage: myscript.R [-i <Imputname>] [-e <error>] [-p <probability>] [-s <probability1>]

Options:
  -i, --Imputname   The file of list of cancers
  -e, --error       Specify the sequencing error ratio (default: 0.02)
  -p, --probability Specify the probability of three genotypes with '-p (1/4,1/2,1/4)' (default: (1/3,1/3,1/3))
  -s, --probability1 Set the probability of observing an allele in genotype (default: 1/2).
"

options <- docopt(doc)
source("./code/main_funs.R")
## parameters
file_ID = options$Imputname #file_ID = "wt1203mo_sorted-003.bam_list.csv"
erp = as.numeric(options$error) #erp = 0.02
erp = 1 - erp
Pp = options$probability #Pp = "(1/3,1/3,1/3)"
Pp = eval(parse(text = paste0("c", Pp)))
pval = as.numeric(options$probability1)  #pval=0.9
#print(erp)
#print(threshold)
#print(Pp)

#setwd("C:\\Users\\lxj423\\OneDrive - University of Miami\\work_data\\HLA\\NUC_nEW\\R\\pipeline")





data_all = data.frame()
d1 = fread(paste0(file_ID,".candidate.csv"),header = T,sep = ",")

if (nrow(d1) > 0){
  
  geness = unique(d1$Gene)
  data_temp = data.frame()
  
  for (gene_ID in geness){
    subdata = subset(d1, Gene == gene_ID & Level == "One-field" & Num.Unique > 0)
    if (nrow(subdata) > 0){
      data_all = rbind(data_all,Getdata(subdata))
    }
  }
  

  
  d2 = data_all
  
  d2$ID = paste(d2$Sample,d2$gene,sep = "_")
  IDs = unique(d2$ID)
  
  
  data_all = data.frame()
  for (IDs_id in IDs){
    subd1 = subset(d2,ID == IDs_id)
    subd1 = subd1[which.max(subd1$Pro),]
    if (nrow(subd1) > 0){
      data_all = rbind(data_all,subd1)
    }
  }
  
  
  data_all = data_all[,-ncol(data_all)]
  write.csv(data_all, paste0(file_ID,"_genetype.csv") ,row.names = FALSE)
  

}else{
  print("No HLA genotype was detected under the specified conditions.")
}


