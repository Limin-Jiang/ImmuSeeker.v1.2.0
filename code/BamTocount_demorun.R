#   ImmuSeeker
#
#   Authors:    Limin Jiang (lxj423@med.miami.edu)
#   Date:       June 1, 2024
#   Version:    0.1.0

rm(list = ls())
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(dplyr)
  library(docopt)
})

doc <- "
Usage: myscript.R [-i <Imputname>] [-c <threshold>] [-x <expression>] [-q <coverage>] [-t <Type>] [-d <directory>]

Options:
  -i, --Imputname   The file of list of cancers
  -c, --threshold   Specify min number of supported reads (default: 0)
  -x, --expression  Specifies whether to output gene expression values (default: false)
  -q, --coverage  Specifies the minimum coverage threshold required for an allele to be retained (default: 98)
  -t, --Type   HLA or KIR
  -d, --directory   /app/data or /ImmuSeeker_data
"

options <- docopt(doc)

source("./code/main_funs.R")

## parameters
file_ID = options$Imputname #file_ID = "C:/Users/lxj423/Downloads/ImmuSeeker-main/ImmuSeeker-main/example.fq.output"
threshold = as.numeric(options$threshold) #threshold = 1
expression_lag = options$expression #expression_lag = 'true'
coverage_val = as.numeric(options$coverage) #coverage_val = 5
Type = options$Type #Type = 'HLA'
ref_dir = options$directory #ref_dir = "./data"
print(file_ID)

cov_id = paste0(file_ID, ".coverage.txt")
file_ID = paste0(file_ID, ".Nuc.txt")

d_cov = fread(cov_id)
d_cov$Percent = as.numeric(gsub("%","",d_cov$Percent))
d_cov = d_cov%>%
  filter(Percent >= coverage_val)

if (nrow(d_cov) == 0){
  print("No alleles meeting the specified coverage threshold were retained for analysis. You may consider lowering the read quality threshold using -q, reducing the coverage requirement, or trying a different input file.")
}else{
  print(paste0(nrow(d_cov)," alleles meeting the specified coverage threshold were retained for analysis."))
  
  Newbam = fread(file_ID,header = F)
  colnames(Newbam)[3] = "reads"
  
  Newbam = Newbam %>%
    filter(V1 %in% d_cov$Reference )%>%
    distinct()
  Newbam = Newbam[(!grepl("N",reads)),]
  
  colnames(Newbam) = c("V2","V1","reads")
  geneNameS = unique(sapply(strsplit(unique(Newbam$V2), split='*', fixed=TRUE), function(x) (x[1])))
  
  
  
  data_end = data.frame()
  data_all_out= data.frame()
  for (geneName in geneNameS){
    
    print(geneName)
    
    subd1_temp = Newbam[!grepl(paste0("^",geneName,"\\*"), Newbam$V2),]
    subd1_temp = unique(subd1_temp$reads)
    subd1 = Newbam[grepl(paste0("^",geneName,"\\*"), Newbam$V2),]
    subd1 = subset(subd1,!(reads %in% subd1_temp))
    
    if (nrow(subd1) > 0){
      data_all_uniq = getUnique(subd1)
      data_all_out = rbind(data_all_out,data_all_uniq)
      data_TEMP = data_all_uniq#%>% filter(Num > 0)
      data_TEMP$allele = sapply(strsplit(data_TEMP$HLAname, split=':', fixed=TRUE), function(x) (x[1]))
      data_class1 = data_TEMP %>%
        filter(Class != "Gene")
      all_allele = unique(data_class1$allele)
      data_class1 = data_class1 %>%
        filter(Class != "Class1")
      
      data_class1 = data_class1[order(data_class1$Total,decreasing = TRUE),]
      
      f1 <- function(x) { d <- !duplicated(x) ; data.frame(uniqueValue=x[d], firstIndex=which(d)) }
      
      dd = f1(data_class1$allele) 
      
      alleles_list = data_class1[dd$firstIndex,]
      aa = setdiff(all_allele,alleles_list$allele)
      alleles_list = alleles_list$HLAname
      if (length(aa)>0){
        alleles_list = c(alleles_list,aa)
      }
      
      if (length(alleles_list) > 1){
        alleles = getcandidate(subd1,alleles_list)
      }else{
        alleles = alleles_list
      }
      
      alleles = sapply(strsplit(alleles, split=':', fixed=TRUE), function(x) (x[1]))
      
      data_TEMP$allele = sapply(strsplit(data_TEMP$HLAname, split=':', fixed=TRUE), function(x) (x[1]))
      data_class1 = data_TEMP %>%
        filter(allele %in% alleles | Class == "Gene")
      data_end = rbind(data_end,data_class1)
    }
  }
  
  data_end = data_end[,-6]
  colnames(data_end) = c("HLAname","Num.Unique" , "Level","Gene", "Total"  )
  colnames(data_all_out)= c("HLAname","Num.Unique" , "Level","Gene", "Total"  )
}




if (nrow(data_end) > 0){
  data_end = subset(data_end,`Num.Unique`>=threshold & !is.na(`Num.Unique`))
  if (nrow(data_end) > 0){
    
    if (expression_lag == 'true'){
      
      if (Type == "HLA"){
        d_gene_length = fread(paste0(ref_dir, "/HLA-length.txt"),header = FALSE,sep = "\t")
        d_gene_length$V2 = paste0("HLA-",d_gene_length$V2)
      }
      
      if (Type == "KIR"){
        d_gene_length = fread(paste0(ref_dir, "/KIR-length.txt"),header = FALSE,sep = "\t")
      }
      
      
      colnames(d_gene_length) = c("ID","name","length")
      d_gene_length$length = d_gene_length$length/1000
      
      d_total = nrow(Newbam)/10^6
      rm(Newbam)
      data_all_out = merge(data_all_out,d_gene_length,by.x = "HLAname",by.y = "name",all.x = TRUE)
      
      for (ii in 1:nrow(data_all_out)){
        HLA_idx = data_all_out$HLAname[ii]
        d_gene_length_temp = d_gene_length[grepl(paste0("^",gsub("\\*","\\\\*",HLA_idx)),d_gene_length$name) ,]  
        data_all_out$length[ii] = median(d_gene_length_temp$length)
      }
      
      
      data_all_out$RPK = data_all_out$Total/data_all_out$length
      data_all_out$RPKM = (data_all_out$RPK*10^6)/rep(d_total,nrow(data_all_out))
      data_all_out$TPM = (data_all_out$RPKM*10^6)/(sum(data_all_out$RPKM,na.rm = TRUE))
      
      data_all_out$totalNum = rep(d_total*10^6,nrow(data_all_out))
    }
    
    idx = match(data_end$HLAname,data_all_out$HLAname)
    data_end$ID = data_all_out$ID[idx]
    data_end$length = data_all_out$length[idx]
    data_end$RPK = data_all_out$RPK[idx]
    data_end$RPKM = data_all_out$RPKM[idx]
    data_end$TPM = data_all_out$TPM[idx]
    data_end$totalNum = data_all_out$totalNum[idx]
    
    write.csv(data_end,paste0(options$Imputname, ".candidate.csv"),row.names = FALSE)
    write.csv(data_all_out,paste0(options$Imputname, ".all.csv"),row.names = FALSE)
  }
  
}










