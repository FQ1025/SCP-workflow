library(data.table)
library(tidyverse)
library(jsonlite)
library(languageserver)
library(GenomicRanges)
library(patchwork)
library(cowplot)
library(patchwork, help, pos = 2, lib.loc = NULL)
library(IRanges, help, pos = 2, lib.loc = NULL)
library(pheatmap, help, pos = 2, lib.loc = NULL)

library(magrittr)

library( scpdata )
library(QFeatures)
library(scp)
library(Seurat)

library(scales)
library(textshape)

library(impute)
library(scater)
library(sva)
library(RColorBrewer)

##### functions-------
runs.stats = function(tmt){
    all.runs = unique(tmt$Raw.file); length(all.runs)
    c10.runs = c( all.runs[grep("col19", all.runs)] , all.runs[grep("col20", all.runs)] ); length(unique(c10.runs))
    c100.runs = c( all.runs[grep("col21", all.runs)] , all.runs[grep("col22", all.runs)] ); length(unique(c100.runs))
    c1000.runs = c( all.runs[grep("col23", all.runs)] ); length(unique(c1000.runs))
    ladder.runs = c( all.runs[grep("col24", all.runs)] ); length(unique(ladder.runs))
    carrier.runs = c( all.runs[grep("arrier", all.runs)] ); length(unique(carrier.runs))
    other.runs = c( all.runs[c(grep("Ref", all.runs), grep("Master", all.runs), grep("SQC", all.runs), grep("blank", all.runs) )] ); length(unique(other.runs))
    sc.runs = all.runs[!all.runs%in%c(c10.runs,c1000.runs,c100.runs,ladder.runs, carrier.runs, other.runs)]; length(unique(sc.runs))

    colnames(anno)[1] = "Raw.file"
    colnames(batch)[1] = "Raw.file"
    batch.tmp = batch[batch$Raw.file %in% all.runs, ]
    anno.tmp = anno[anno$Raw.file %in% all.runs, ]
    sample.anno = merge( batch.tmp,anno.tmp,by = "Raw.file") %>% gather( key = "Channel",value = "SampleType", -lcbatch,-sortday,-digest,-Raw.file)
    ## count cells and proteins second time
    sc.left = length( unique(sample.anno$set[sample.anno$SampleType %in% c("sc_m0","sc_u")]) )
    prot.left = length(unique(tmt$Leading.razor.protein[tmt$Raw.file %in% sc.runs]))
    pep.left = length(unique(tmt$Modified.sequence[tmt$Raw.file %in% sc.runs]))
    print(paste0("all.runs = ",length(all.runs),", sc.runs = ",length(unique(sc.runs)),", peptides = ",pep.left,", protein = ", prot.left))
    return(list(sc = sc.runs,all = all.runs))
}

scp.stats = function(data){
  # sc.runs = unique(data$Set[data$SampleType %in% c("Macrophage","Monocyte")])
  # sc.runs = grep("^\\d+",names(data),value = T) ## cellline
  sc.runs = grep("^1",names(data),value = T)

  tmp = rbindRowData(data,i = sc.runs)
  
  all.runs = length(na.omit(unique(data$Raw.file)))
  n.peptide = length(na.omit(unique(tmp$Modified.sequence))) #tmp$Modified.sequence
  n.protein = length(na.omit(unique(tmp$Leading.razor.protein)))#tmp$Proteins
  n.cell = length( data$celltype[(data$celltype %in% c("Macrophage","Monocyte"))& data$Raw.file %in% sc.runs ] )
  sc.runs = length( sc.runs )
  print(paste0("all.runs = ",all.runs,", sc.runs = ",sc.runs,", cells = ",n.cell,
               ", peptides = ",n.peptide,", protein = ", n.protein))
  return(paste0("all.runs = ",all.runs,", sc.runs = ",sc.runs,", cells = ",n.cell,
                ", peptides = ",n.peptide,", protein = ", n.protein))
}
# scp.stats = function(data){
#   sc.runs = unique(data$Raw.file[data$celltype %in% c("HeLa","GM12878")])
#   tmp = rbindRowData(data,i = sc.runs)
  
#   all.runs = length(na.omit(unique(data$Raw.file)))
#   sc.runs = length( unique(data$Raw.file[data$celltype %in% c("HeLa","GM12878")]) )
  
#   n.peptide = length(na.omit(unique(tmp$Modified.sequence)))
#   n.protein = length(na.omit(unique(tmp$Proteins)))
#   n.cell = length( data$celltype[data$celltype %in% c("HeLa","GM12878")] )
#   print(paste0("all.runs = ",all.runs,", sc.runs = ",sc.runs,", cells = ",n.cell,
#                ", peptides = ",n.peptide,", protein = ", n.protein))
#   return(paste0("all.runs = ",all.runs,", sc.runs = ",sc.runs,", cells = ",n.cell,
#                 ", peptides = ",n.peptide,", protein = ", n.protein))
# }
# scp.stats = function(data){
#   sc.runs = unique(data$Raw.file[data$celltype %in% c("LPS","untreated")])
#   tmp = rbindRowData(data,i = sc.runs)
  
#   all.runs = length(na.omit(unique(data$Raw.file)))
#   sc.runs = length( unique(data$Raw.file[data$celltype %in% c("LPS","untreated")]) )
  
#   n.peptide = length(na.omit(unique(tmp$Modified.sequence)))
#   n.protein = length(na.omit(unique(tmp$Proteins)))
#   n.cell = length( data$celltype[data$celltype %in% c("LPS","untreated")] )
#   print(paste0("all.runs = ",all.runs,", sc.runs = ",sc.runs,", cells = ",n.cell,
#                ", peptides = ",n.peptide,", protein = ", n.protein))
#   return(paste0("all.runs = ",all.runs,", sc.runs = ",sc.runs,", cells = ",n.cell,
#                 ", peptides = ",n.peptide,", protein = ", n.protein))
# }

sc.run = function(data){
  # sc.runs = unique(data$Set[data$SampleType %in% c("Macrophage","Monocyte")])
  sc.runs =  grep("^\\d+",names(data),value = T)
  # sc.runs = unique(data$Raw.file[data$celltype %in% c("LPS","untreated")])
  return(sc.runs)
}
keep.run = function(data){
  sc.runs = unique(data$Raw.file[data$SampleType %in% c("Macrophage","Monocyte")])
  # sc.runs = unique(data$Raw.file[data$celltype %in% c("Blank")])
  return(sc.runs)
}

GSEA_Analysis = function( diff_data, gsea_pvalue = 0.05 ){
  
  gene_set = read.gmt("GSEA_imput_files/merge.gmt")
  
  gsea_data = diff_data[,which(colnames( diff_data )%in% colnames_needed)]
  gsea_data = rownames_to_column(gsea_data,"SYMBOL")
  # transfer_id = bitr( gsea_data$SYMBOL,
  #                      fromType = "SYMBOL",toType = "ENTREZID",
  #                      OrgDb = GO_database)
  # gsea_data = merge( gsea_data, transfer_id,by = 'SYMBOL')
  gsea_data = drop_na(gsea_data)
  
  gsea_input = gsea_data$log2FoldChange
  names(gsea_input) = gsea_data$SYMBOL
  gsea_input = sort( gsea_input,decreasing = T )
  
  gsea_output = GSEA( gsea_input,TERM2GENE = gene_set,verbose = F,pvalueCutoff = gsea_pvalue)
  head(gsea_output)
  return(gsea_output)
}

# sc.run = function(data){
#   # sc.runs = unique(data$Raw.file[data$SampleType %in% c("Macrophage","Monocyte")])
#   sc.runs = unique(data$Raw.file[data$celltype %in% c("HeLa","GM12878")])
#   return(sc.runs)
# }
# keep.run = function(data){
#   # sc.runs = unique(data$Raw.file[data$SampleType %in% c("Macrophage","Monocyte")])
#   sc.runs = unique(data$Raw.file[data$celltype %in% c("blank")])
#   return(sc.runs)
# }
