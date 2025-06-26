getrich = function(d1_temp){
  d1_temp = d1_temp %>%
    select("AlleleID","Num")%>%
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

getdir = function(subdata_class3,subdata_class4){
  data = data.frame()
  for (ii in 1:nrow(subdata_class3) ){
    HLA_ID = subdata_class3$AlleleID[ii]
    subdata = subdata_class4[grepl( paste0("^",  gsub("\\*","\\\\*",HLA_ID),":" ),subdata_class4$AlleleID),]
    if (nrow(subdata)){
      data = rbind(data,data.frame(from = paste0(HLA_ID,"(",subdata_class3$Num[ii],")" ) ,to = paste0(subdata$AlleleID,"(",subdata$Num,")" )))
    }
    
  }
  return(data)
  
}


getdata = function(d1,gene_ID){
  
  
  subdata = subset(d1, Level == "One-field")
  if (nrow(subdata) > 0){
    subdata_class1 = subdata[order(subdata$Num,decreasing=TRUE),]
  }else{
    subdata_class1 = NULL
  }
  
  
  subdata = subset(d1, Level == "Two-field")
  if (nrow(subdata) > 0){
    subdata_class2 = subdata[order(subdata$Num,decreasing=TRUE),]
  }else{
    subdata_class2 = NULL
  }
  
  
  subdata = subset(d1, Level == "Three-field")
  if (nrow(subdata) > 0){
    subdata_class3 = subdata[order(subdata$Num,decreasing=TRUE),]
  }else{
    subdata_class3 = NULL
  }
  
  
  
  subdata = subset(d1, Level == "Four-field")
  if (nrow(subdata) > 0){
    subdata_class4 = subdata[order(subdata$Num,decreasing=TRUE),]
  }else{
    subdata_class4 = NULL
  }
  
  
  data_all = data.frame()
  
  if (!is.null(subdata_class1)){
    
    data = data.frame(from = gene_ID,to = paste0(subdata_class1$AlleleID,"(",subdata_class1$Num,")" ))
    data_all = rbind(data_all,data)
  }
  
  
  
  if (!is.null(subdata_class1) & !is.null(subdata_class2)){
    data_all = rbind(data_all,getdir(subdata_class1,subdata_class2))
    
  }
  
  if (!is.null(subdata_class2) & !is.null(subdata_class3)){
    data_all = rbind(data_all,getdir(subdata_class2,subdata_class3))
  }
  if (!is.null(subdata_class3) & !is.null(subdata_class4)){
    data_all = rbind(data_all,getdir(subdata_class3,subdata_class4))
  }
  return(data_all)
  
}


Getdata = function(subdata){
  
  colnames(subdata)[2] = "Num"
  if (max(subdata$Num) > 1000){
    subdata$Num =round(subdata$Num * (1000/max(subdata$Num)))
  }
  
  
  if (nrow(subdata) > 1){
    subdata = subdata[order(subdata$Num,decreasing=TRUE),]
    data_temp = subdata[c(1:2),]
    data_temp = getP(data_temp,Pp,erp,pval)
    data_temp$gene = gene_ID
    data_temp$Sample = basename(file_ID) 
  }
  if (nrow(subdata) == 1){
    subdata = subdata[order(subdata$Num,decreasing=TRUE),]
    data_temp = rbind(subdata,subdata)
    data_temp$AlleleID[2]= ""
    data_temp$Num[2] = 0
    data_temp = getP(data_temp,Pp,erp,pval)
    data_temp$gene = gene_ID
    data_temp$Sample = basename(file_ID) 
    data_temp = data_temp[1,]
  }
  return(data_temp)
  
}

getP = function(data_temp,Pp,erp,pval){
  numsum = sum(data_temp$Num )
  binomialcoff = data_temp$Num[1]
  Pd = c(dbinom(binomialcoff, size=numsum, prob=erp), dbinom(binomialcoff, size=numsum, prob=pval),dbinom(binomialcoff, size=numsum, prob=(1-erp) )  )
  
  PAAd = (Pd[1]* Pp[1])/(sum(Pd * Pp))
  PAGd = (Pd[2]* Pp[2])/(sum(Pd * Pp))
  PGGd = (Pd[3]* Pp[3])/(sum(Pd * Pp))
  
  data_end = data.frame(Genetype = c(paste0(data_temp$AlleleID[1],",",data_temp$AlleleID[1]),paste0(data_temp$AlleleID[1],",",data_temp$AlleleID[2]),paste0(data_temp$AlleleID[2],",",data_temp$AlleleID[2])),Pro = c(PAAd,PAGd,PGGd),
                        Number = c(paste0(data_temp$Num[1],"?",data_temp$Num[1]),paste0(data_temp$Num[1],"?",data_temp$Num[2]),paste0(data_temp$Num[2],"?",data_temp$Num[2])),stringsAsFactors = FALSE)
  data_end$PL = -10*log10(data_end$Pro + 10^-100)
  data_end$Normalized_PL = data_end$PL - min(data_end$PL)
  data_end$Type = c("homo","heter","homo")
  return(data_end)
}


getUnique = function(subd1){
  data_all = data.frame()
  aa = sapply(strsplit(subd1$V2, split=':', fixed=TRUE), function(x) (x[1:4]))
  
  subd1$Class1 = aa[1,]
  subd1$Class2 = paste(aa[1,],aa[2,],sep = ":")
  subd1$Class3 = paste(aa[1,],aa[2,],aa[3,],sep = ":")
  
  colnames(subd1)[2] = "Class4"
  
  
  subd1$Class1[grepl(":NA",subd1$Class1)] = NA
  subd1$Class2[grepl(":NA",subd1$Class2)] = NA
  subd1$Class3[grepl(":NA",subd1$Class3)] = NA
  subd1$Class4[grepl(":NA",subd1$Class4)] = NA
  
  data_temp4 = data.frame(table(subd1$Class4))
  data_temp4$Var1 =  as.character(data_temp4$Var1)
  data_temp3 = data.frame(table(subd1$Class3))
  data_temp3$Var1 =  as.character(data_temp3$Var1)
  data_temp2 = data.frame(table(subd1$Class2))
  data_temp2$Var1 =  as.character(data_temp2$Var1)
  data_temp1 = data.frame(table(subd1$Class1))
  data_temp1$Var1 =  as.character(data_temp1$Var1)
  
  Class4_list = data_temp4$Var1[str_count(data_temp4$Var1, ":") == 3]
  Class3_list = data_temp3$Var1[str_count(data_temp3$Var1, ":") == 2]
  Class2_list = data_temp2$Var1[str_count(data_temp2$Var1, ":") == 1]
  Class1_list = data_temp1$Var1[str_count(data_temp1$Var1, ":") == 0]
  
  
  data_temp = data.frame(AlleleID = geneName, Num = length(unique(subd1$reads)),Class = "Gene", Genes = geneName,  Total = length(unique(subd1$reads)), stringsAsFactors = FALSE)
  data_all = rbind(data_all,data_temp)
  
  
  for (Class1_list_ID in Class1_list){
    subsubset1 = subset(subd1, Class1 != Class1_list_ID)
    subsubset1 = unique(subsubset1$reads)
    subsubset2 = subset(subd1, Class1 == Class1_list_ID)
    subsubd1 = subset(subd1, Class1 == Class1_list_ID & !(reads %in% subsubset1) ) #
    
    Class4_list_temp = Class4_list[grepl(paste0("^",  gsub("\\*","\\\\*",Class1_list_ID),":" ),Class4_list)]
    Class3_list_temp = Class3_list[grepl(paste0("^",  gsub("\\*","\\\\*",Class1_list_ID),":" ),Class3_list)]
    Class2_list_temp = Class2_list[grepl(paste0("^",  gsub("\\*","\\\\*",Class1_list_ID),":" ),Class2_list)]
    
    data_temp = data.frame(AlleleID = Class1_list_ID, Num =  length(unique(subsubd1$reads )), Class = "One-field", Genes = geneName, Total = length(unique(subsubset2$reads)),stringsAsFactors = FALSE)
    data_all = rbind(data_all,data_temp)
    
    
    if (nrow(subsubd1) > 0){
      
      data_temp = data.frame()
      for (ii in Class2_list_temp){
        subsubset1 = subset(subd1, Class2 != ii )
        subsubset2 = subset(subd1, Class2 == ii )
        aa = setdiff(subsubset2$reads,subsubset1$reads)
        data_temp = rbind(data_temp,data.frame(AlleleID = ii, Num = length(aa),Total = length(unique(subsubset2$reads)),stringsAsFactors = FALSE))
      }
      if (nrow(data_temp) > 0){
        data_temp$Class = "Two-field"
        data_temp = data_temp[order(data_temp$Num,decreasing=TRUE),]
        data_temp$Genes = geneName
        data_temp = data_temp[,c(1,2,4,5,3)]
        data_all = rbind(data_all,data_temp)
      }
      
      data_temp = data.frame()
      for (ii in Class3_list_temp){
        subsubset1 = subset(subd1, Class3 != ii)
        subsubset2 = subset(subd1, Class3 == ii )
        aa = setdiff(subsubset2$reads,subsubset1$reads)
        data_temp = rbind(data_temp,data.frame(AlleleID = ii, Num = length(aa),Total = length(unique(subsubset2$reads)), stringsAsFactors = FALSE))
        
      }
      
      if (nrow(data_temp) > 0){
        data_temp$Class = "Three-field"
        data_temp = data_temp[order(data_temp$Num,decreasing=TRUE),]
        data_temp$Genes = geneName
        data_temp = data_temp[,c(1,2,4,5,3)]
        data_all = rbind(data_all,data_temp)
      }
      
      
      
      
      
      data_temp = data.frame()
      for (ii in Class4_list_temp){
        aa = paste(sapply(strsplit(ii, split=':', fixed=TRUE), function(x) (x[1:3])),collapse  = ":")
        subsubset1 = subset(subd1, Class4 != ii & Class3 == aa )
        subsubset2 = subset(subd1, Class4 == ii )
        aa = setdiff(subsubset2$reads,subsubset1$reads)
        data_temp = rbind(data_temp,data.frame(AlleleID = ii, Num = length(aa),Total = length(unique(subsubset2$reads)), stringsAsFactors = FALSE))
        
      }
      
      if (nrow(data_temp) > 0){
        data_temp$Class = "Four-field"
        data_temp = data_temp[order(data_temp$Num,decreasing=TRUE),]
        data_temp$Genes = geneName
        data_temp = data_temp[,c(1,2,4,5,3)]
        data_all = rbind(data_all,data_temp)
      }
      
    }
    
  }
  return(data_all)
}



getcandidate = function(subd1,alleles_list){
  
  aa = sapply(strsplit(subd1$V2, split=':', fixed=TRUE), function(x) (x[1:4]))
  
  subd1$Class1 = aa[1,]
  subd1$Class2 = paste(aa[1,],aa[2,],sep = ":")
  
  allele_reads <- list()
  for (allele_id in alleles_list) {
    allele_reads[[allele_id]] <- subd1$reads[subd1$Class2 == allele_id|subd1$Class1 == allele_id]
  }
  
  allele_pairs <- combn(names(allele_reads), 2, simplify = FALSE)
  
  # For each pair, compute the number of unique reads in their union
  pair_scores <- sapply(allele_pairs, function(pair) {
    reads_union <- union(allele_reads[[pair[1]]], allele_reads[[pair[2]]])
    length(reads_union)
  })
  
  # Find the pair with the maximum coverage
  best_pair <- allele_pairs[[which.max(pair_scores)]]
  
  return(best_pair)
}