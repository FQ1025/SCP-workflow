setwd("/home/fangq/ffqq/proj/scp/cellLine/allfile")
source("/home/fangq/ffqq/proj/scp/pSCoPE/analysis/scp_packages.R")
#######-------------functions------------######
sc.run = function(data){
  sc.runs = unique(data$Raw.file[data$celltype %in% c("HeLa","U2OS")])
  # sc.runs = unique(data$Raw.file[data$celltype %in% c("LPS","untreated")])
  return(sc.runs)
}
scp.stats = function(data){
  sc.runs = unique(data$Raw.file[data$celltype %in% c("HeLa","U2OS")])

  tmp = rbindRowData(data,i = sc.runs)
  
  all.runs = length(na.omit(unique(data$Raw.file)))
  n.peptide = length(na.omit(unique(tmp$Modified.sequence))) #tmp$Modified.sequence
  n.protein = length(na.omit(unique(tmp$Leading.razor.protein)))#tmp$Proteins
  n.cell = length( data$celltype[(data$celltype %in% c("HeLa","U2OS"))& data$Raw.file %in% sc.runs ] )
  sc.runs = length( sc.runs )
  print(paste0("all.runs = ",all.runs,", sc.runs = ",sc.runs,", cells = ",n.cell,
               ", peptides = ",n.peptide,", protein = ", n.protein))
  return(paste0("all.runs = ",all.runs,", sc.runs = ",sc.runs,", cells = ",n.cell,
                ", peptides = ",n.peptide,", protein = ", n.protein))
}
#######-------------functions------------######

batch = read.csv("/home/fangq/ffqq/proj/scp/cellLine/allfile/metadata.csv")
batch = batch[1:22,1:6]

annotation = data.frame( Raw.file = rep(batch$Raw.file,each = 16),
                         tmt.label = rep(paste0("RI",1:16),time = 22))

annotation$celltype = ifelse(annotation$Raw.file%in%batch$Raw.file[batch$Order=="UH"],
                             c("Carrier","Reference","Unused","Unused",rep("U2OS",6),rep("HeLa",5),"Blank"),
                             ifelse(annotation$Raw.file=="20220712_SCP_HA_C_R",c("Carrier","Reference",rep("Unused",14)),
                             c("Carrier","Reference","Unused","Unused",rep("HeLa",6),rep("U2OS",5),"Blank")))

annotation$tmt.label = factor(annotation$tmt.label,levels = paste0("RI",1:16),ordered = T)
annotation = annotation %>% arrange(Raw.file,tmt.label)
metadata = merge(batch,annotation,by = "Raw.file") %>% arrange(Raw.file,tmt.label)

scp.cellline = readRDS("scp_peptides_aggregated.rds")

medianRI = colMeans( assay(scp.cellline[["peptides"]]),na.rm = T )
scp.cellline$medianRI = medianRI


scp.tmp = medianCVperCell( scp.cellline, 
                        i = sc.run(scp.cellline),
                        groupBy = "Leading.razor.protein",
                        nobs = 5,  # least peptide numbers per protein,2 and 5
                        norm = "SCoPE2",
                        na.rm = T, 
                        colDataName = "MedianCV" ) 
scp.cellline = scp.tmp[,!is.na(scp.tmp$MedianCV) & scp.tmp$MedianCV < 0.3 & scp.tmp$celltype %in% c("HeLa", "U2OS"), ]

runs = intersect(batch$Raw.file[batch$Carrier.x==200],sc.run(scp.cellline))
pep_runs = paste0("peptides_",runs)
scp.cellline = scp.cellline[,,c(runs,pep_runs)]
scp.cellline = joinAssays( scp.cellline,i = pep_runs,
                  name = "peptides" )

scp.cellline = sweep(scp.cellline, 
             i = "peptides", ## assay to be normalized
             MARGIN = 2,
             FUN = "/",
             STATS = colMedians(assay(scp.cellline[["peptides"]]), na.rm = TRUE), ## *****Can be changed by different methods for normalization (eg, quantile?)
             name = "peptides_norm_col")

scp.cellline = sweep(scp.cellline,
             i = "peptides_norm_col",  ## after col normalization do row normalization 
             MARGIN = 1,
             FUN = "/",
             STATS = rowMeans(assay(scp.cellline[["peptides_norm_col"]]),  na.rm = TRUE),
             name = "peptides_norm")

scp.cellline = filterNA(scp.cellline,
                i = "peptides_norm",
                pNA = 0.99)

scp.cellline = logTransform(scp.cellline,
                    base = 2,
                    i = "peptides_norm",
                    name = "peptides_log")

# tmp = rbindRowData(scp.cellline,"peptides_log")
# colnames(tmp)
scp.cellline = aggregateFeatures(scp.cellline,
                         i = "peptides_log",
                         name = "proteins",
                         fcol = "Leading.razor.protein",
                         fun = matrixStats::colMedians, 
                        na.rm = TRUE)


scp.cellline = sweep(scp.cellline, 
             i = "proteins", ## assay to be normalized
             MARGIN = 2,
             FUN = "-",
             STATS = colMedians(assay(scp.cellline[["proteins"]]), na.rm = TRUE), ## *****Can be changed by different methods for normalization (eg, quantile?)
             name = "proteins_norm_col")

scp.cellline = sweep(scp.cellline,
             i = "proteins_norm_col",  ## after col normalization do row normalization 
             MARGIN = 1,
             FUN = "-", ### initial : /
             STATS = rowMeans(assay(scp.cellline[["proteins_norm_col"]]),  na.rm = TRUE),
             name = "proteins_norm")
tmp = nNA(scp.cellline,i = "proteins_norm")
pros = tmp$nNArows[tmp$nNArows$pNA < 90.6,] %>% nrow()
scp.cellline = filterNA(scp.cellline,
                i = "proteins_norm",
                pNA = 0.906)

scp.cellline[["proteins_norm"]] %>%
    assay %>% 
    is.na %>%  mean() 
scp.cellline = impute(scp.cellline,
              i = "proteins_norm",
              method = "knn",
              k = 3, rowmax = 1, colmax= 1,
              maxp = Inf, rng.seed = 1025 )
scp.cellline[["proteins_norm"]] %>%
    assay %>%
    is.na %>%
    mean #### 0



sce = getWithColData(scp.cellline, "proteins_norm")

batch = colData(sce)$Raw.file
model = model.matrix(~ celltype, data = colData(sce))  # before combat, excluded the intersted variables
library(sva)
## Combat  bacth corrected
assay(sce) = ComBat(dat = assay(sce),
                     batch = batch,
                     mod = model)

scp.cellline = addAssay(scp.cellline,
                y = sce,
                name = "proteins_batchC")

scp.cellline = addAssayLinkOneToOne(scp.cellline, 
                            from = "proteins_norm",
                            to = "proteins_batchC")

scp.cellline = sweep(scp.cellline,
                 i = "proteins_norm_col",  ## after col normalization do row normalization 
                 MARGIN = 1,
                 FUN = "-",
                 STATS = rowMeans(assay(scp.cellline[["proteins_norm_col"]]),  na.rm = TRUE),
                 name = "proteins_norm_unimputed")
scp.cellline = filterNA(scp.cellline,
                i = "proteins_norm_unimputed",
                pNA = 0.906)

plot(scp.cellline)
sce = getWithColData(scp.cellline, "proteins_norm_unimputed")
library(limma)
batch = colData(sce)$Raw.file
assay(sce) = removeBatchEffect( assay(sce),
                                batch = batch)
scp.cellline = addAssay(scp.cellline,
                    y = sce,
                    name = "proteins_batchL")
scp.cellline = addAssayLinkOneToOne(scp.cellline, 
                                from = "proteins_norm_unimputed",
                                to = "proteins_batchL")
plot(scp.cellline,interactive = T)
scp.stats(scp.cellline)
save( scp.cellline,file = "for_seurat.RData",compress = T )
library(scater)
# library(loomR)
scp.tmp = getWithColData(scp.cellline,"proteins_batchC")

scp.tmp = runPCA(scp.tmp,
                ncomponents = 10, 
                ntop = Inf,
                scale = TRUE,
                exprs_values = 1,
                name = "PCA")
p = plotReducedDim(scp.tmp,
               dimred = "PCA",
               colour_by = "celltype",
               point_alpha = 1, text_size = 20,point_size = 4 ) + 
               theme_cowplot() 
p

mat = assay(scp.tmp)
pearson = cor(t(assay(scp.tmp)))
head(pearson[1:10,1:10])
rsum = rowSums(pearson^2)

# Calculate the weighted data matrix:
mat = diag(rsum) %*%  mat
pca.imp.cor = cor(mat)

# pca = prcomp(t(mat))
pca = prcomp(pca.imp.cor)
deviation = (pca$sdev)^2/sum((pca$sdev)^2)*100
sum(deviation[1:2]) ;length(deviation)

anno = scp.tmp@colData %>% as.data.frame()
pca_mat = cbind(pca$x,anno)

# pearson = cor(scp.tmp)
pca_mat$celltype = factor(pca_mat$celltype,levels = c("HeLa","U2OS"))
pca_mat$Raw.file = as.factor(pca_mat$Raw.file)

p8 = ggplot(pca_mat,aes(x = PC1,y = PC2,color = celltype)) + 
  geom_point(size =5) + 
  # scale_shape_manual(values = c(15:21,8,13)) + 
  scale_color_manual(values = c("#b08ed6","#81c2a2")) +
  xlab(paste0("PC1 (",round(deviation[1],digits = 1),"%)")) +
  ylab(paste0("PC2 (",round(deviation[2],digits = 1),"%)")) +
  theme(text = element_text(size = 26,family = "sans")) + 
  theme_cowplot() + 
  theme(text = element_text(size = 26),
        axis.text = element_text( size = 26),
        axis.title = element_text(size = 26),
        axis.line = element_line(size = 2),
         axis.ticks = element_line(size = 2)) 
p8

pca = prcomp( t(assay(scp.tmp) ))

anno = scp.tmp@colData %>% as.data.frame()
pca_mat = cbind(pca$x,anno)

# pearson = cor(scp.tmp)
pca_mat$celltype = factor(pca_mat$celltype,levels = c("HeLa","U2OS"))
pca_mat$Raw.file = as.factor(pca_mat$Raw.file)
deviation = (pca$sdev)^2/sum((pca$sdev)^2)*100
sum(deviation[1:2]) ;length(deviation)

p9 = ggplot(pca_mat,aes(x = PC1,y = PC2,color = celltype)) + 
  geom_point(size =5) + 
  # scale_shape_manual(values = c(15:21,8,13)) + 
  scale_color_manual(values = c("#b08ed6","#81c2a2")) +
  xlab(paste0("PC1 (",round(deviation[1],digits = 1),"%)")) +
  ylab(paste0("PC2 (",round(deviation[2],digits = 1),"%)")) +
  theme(text = element_text(size = 26,family = "sans")) + 
  theme_cowplot() + 
  theme(text = element_text(size = 26),
        axis.text = element_text( size = 26),
        axis.title = element_text(size = 26),
        line = element_line(size = 2)) 
p9


proteins = getWithColData( scp.cellline,i = "proteins_norm" )
pros = rownames(proteins) %>% as.data.frame()
write.table(pros,"uniprot_input.txt",sep = "\t",quote = F,col.names = F,row.names = F)
# uniprot_database = readxl::read_xlsx( "/lustre/user/liclab/fangq/proj/scp/pSCoPE/analysis/cellline_vector/uniprot-database_filtered-organism__Mus+musculus.xlsx" )
uniprot_search = readxl::read_xlsx( "uniprot_output.xlsx" )
colnames(uniprot_search)[1:2] = c("InputID","UniprotID")
colnames(uniprot_search)=gsub(" ",".",colnames(uniprot_search))
setdiff(uniprot_search$InputID,uniprot_search$UniprotID)

library(org.Hs.eg.db)
library(ensembldb)
library(clusterProfiler)
ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
wg.gene = data.frame(genes(ensdb)) %>%  dplyr::filter( gene_biotype == "protein_coding" )
wg.gene = wg.gene$symbol
###gene-protein match 
gene_pro_ref = uniprot_search[,c('UniprotID','Gene.Names.(primary)' )]
# colnames(gene_pro_ref) = c("UniprotID","Gene")
colnames(gene_pro_ref) = c("proteins","genes")
gene_pro_ref$gene = apply(gene_pro_ref,1,function(x){
  x = gsub(" ","",x[2])
  x = strsplit(x,";") %>% unlist()
  if(length(x)!=1){
  x = x[1]
  }
  return(x)
})

###GO use clusterProfiler
# all.genes = gsub(" ",",",unique(gene_pro_ref$gene))
# all = enrichGO( all.genes,
#                 OrgDb = org.Hs.eg.db,
#                 keyType = "SYMBOL",ont = "ALL",minGSSize = 5,
#                 pAdjustMethod = "BH",universe = wg.gene,
#                 pvalueCutoff = 0.05,qvalueCutoff = 0.05)
# all.mf = enrichGO( all.genes,
#                    OrgDb = org.Hs.eg.db,
#                    keyType = "SYMBOL",ont = "MF",minGSSize = 5,
#                    pAdjustMethod = "BH",universe = wg.gene,
#                    pvalueCutoff = 0.05,qvalueCutoff = 0.05)
# all.bp = enrichGO( all.genes,
#                    OrgDb = org.Hs.eg.db,
#                    keyType = "SYMBOL",ont = "BP",minGSSize = 5,
#                    pAdjustMethod = "BH",universe = wg.gene,
#                    pvalueCutoff = 0.05,qvalueCutoff = 0.05)
# library(clusterProfiler)
# all.bp.simplified = clusterProfiler::simplify(x = all.bp,cutoff=0.8,by="p.adjust",select_fun=min)  #去除冗余，可以调整cutoff值。
# all.mf.simplified = clusterProfiler::simplify(x = all.mf,cutoff=0.8,by="p.adjust",select_fun=min)  #去除冗余，可以调整cutoff值。
# #all.bp.simplified = clusterProfiler::simplify(x = all,cutoff=0.8,by="p.adjust",select_fun=min)  #去除冗余，可以调整cutoff值。
# save(all,all.mf,all.bp,all.bp.simplified,all.mf.simplified,file = "allgenes_GO.RData",compress = T)

load("allgenes_GO.RData",verbose = T)
go = rbind(all.bp.simplified@result %>% mutate( ONTOLOGY = "BP"),all.mf.simplified@result %>% mutate( ONTOLOGY = "MF"))

go$go_pros_match = apply( go,1,function(x){
  k = x[8]
  genes = strsplit("/",x = k) %>% unlist() %>% unique()
  pros = gene_pro_ref$proteins[gene_pro_ref$genes%in% genes] %>% unique()
  x = paste( c(pros),collapse = "/" )
  return(x)
} )
go$pros_num = apply(go,1,function(x){x=length( unlist(strsplit("/",x = x[11])))}) 
go$bgratio = apply(go,1,function(x){x=as.numeric(gsub("/\\d+","",x[3]))/as.numeric(gsub("/\\d+","",x[4]))}) 

go.filtered = go %>% dplyr::filter( pros_num >=5 & bgratio>=0.1  )
### proteins-GO term transformation
go_pro_ref = as.list( go.filtered$ID )
names(go_pro_ref) = go_pro_ref
go_pro_ref = lapply(go_pro_ref, function(x){
  k = go.filtered[go.filtered$ID==x,]
  pros  = k$go_pros_match
  pros = lapply( as.list(pros),function(x){ x = strsplit("/",x = x) %>% unlist } ) %>% unlist() %>% unique()
  return(pros)
})
go_pro_ref = plyr::ldply( go_pro_ref,cbind )
colnames(go_pro_ref) = c("GOterms","proteins")

scp.tmp = scp.cellline#scp.bmdm[["proteins_batchL"]]
dep_mat = assay(scp.tmp,"proteins_batchL") 
# dep_mat = dep_mat %>% t() %>% apply( .,2,function(x){x = (x-mean(x,na.rm = T))/sd(x,na.rm = T)} ) %>% t()
dep_mat = dep_mat %>% as.data.frame() %>% rownames_to_column(.,"proteins")

go_quant_mat = merge( go_pro_ref, dep_mat,by = "proteins",all.x = T )

merge = as.list(unique(go_quant_mat$GOterms))
names(merge) = unique(go_quant_mat$GOterms)

merge = lapply(merge,function(x){
  tmp = go_quant_mat %>% dplyr::filter( GOterms==x ) 
  median = matrixStats::colMedians(as.matrix(tmp[,-c(1:2)]),na.rm = T)
  names( median ) = colnames(tmp[,-c(1:2)])
  return(median)
})

merge = plyr::ldply(merge,rbind)
colnames(merge)[1] = "GOterms"
save( merge,file = "mergeGO_median_mat.RData",compress = T )

{
 normedassay1 = getWithColData(scp.cellline, "proteins_batchL")
  assay(normedassay1,"assay") = assay(normedassay1,"assay") %>% t() %>% apply( .,2,function(x){x = (x-mean(x,na.rm = T))/sd(x,na.rm = T)} ) %>% t()
  seur = as.Seurat( normedassay1,
                  counts = "aggcounts", data = "assay" )
  seur = RenameAssays(seur,originalexp = "batchLnorm")

  normedassay1 = getWithColData(scp.cellline, "proteins_batchL")
  # assay(normedassay1,"assay") = apply( assay(normedassay1,"assay") ,2,function(x){x = (x-mean(x,na.rm = T))/sd(x,na.rm = T)} )
  add = as.Seurat( normedassay1,
                    counts = "aggcounts", data = "assay" )
  ttt = GetAssayData( object = add,slot = "data" )
  seur[["batchLraw"]] = CreateAssayObject( data = ttt )


  normedassay2 = getWithColData(scp.cellline, "proteins_batchC")
  add = as.Seurat( normedassay2,
                   counts = "aggcounts", data = "assay" )
  ttt = GetAssayData( object = add,slot = "data" )
  seur[["imputedraw"]] = CreateAssayObject( data = ttt )
  assay(normedassay2,"assay") = assay(normedassay2,"assay")  %>% t() %>% apply( .,2,function(x){x = (x-mean(x,na.rm = T))/sd(x,na.rm = T)} ) %>% t()
  add = as.Seurat( normedassay2,
                    counts = "aggcounts", data = "assay" )
  ttt = GetAssayData( object = add,slot = "data" )
  seur[["imputedzscored"]] = CreateAssayObject( data = ttt )

  normedassay3 = scp.cellline[["proteins_batchC"]]
  # metadata = colData(scp.cellline)
  # idx = metadata %>% as.data.frame() %>% rownames_to_column("cell") %>% 
  #               select("cell","SampleType") %>% column_to_rownames("cell")
  pearson = cor(t(assay(normedassay3,"assay")))
  rsum = rowSums(pearson^2)
  assay(normedassay3,"assay",withDimnames = F) = diag(rsum) %*%  assay(normedassay3,"assay") 
  add = as.Seurat( normedassay3,
                    counts = "aggcounts", data = "assay" )
  ttt = GetAssayData( object = add,slot = "data" )
  seur[["weighted_imputed_batchC"]] = CreateAssayObject( data = ttt )

  tmp = assay(normedassay3,"assay") %>% cor()
  tmp = SingleCellExperiment(assays = list(assay = tmp,aggcounts = tmp))
  tmp@colData = normedassay3@colData

  add = as.Seurat( tmp,counts = "aggcounts", data = "assay" )
  ttt = GetAssayData( object = add,slot = "data" )
  seur[["weighted_cor_imputed_batchC"]] = CreateAssayObject( data = ttt )

  load("mergeGO_median_mat.RData",verbose = T)
  merge = column_to_rownames(merge,"GOterms")
  add = CreateSeuratObject( merge,assay = "GOterm",meta.data = seur@meta.data )
  ttt = GetAssayData( object = add,slot = "data" )
  seur[["GOterm"]] =  CreateAssayObject( data = ttt )
  Idents(seur) = "celltype"
  save(seur,file = "seurat.RData",compress = T)
}
load("seurat.RData",verbose = T)

library(Seurat)
# library(pacman)
# p_load(Seurat)
# p_unload(SeuratDisk)
# p_unload(Seurat)
# p_load(Seurat)
set.seed(1025)
cellline = seur

VlnPlot( cellline,features = c("nFeature_batchLnorm"),group.by = "Raw.file" )

ggplot(cellline@meta.data) + 
  geom_density(aes(x = nFeature_originalexp),fill = "grey",size = 1) + 
  # scale_fill_manual(values = alpha(brewer.pal(name = "Set2", n = 6),alpha = 0.6)) +
    geom_vline(xintercept = median(cellline@meta.data$nFeature_originalexp),size = 1.5,linetype = 2) +
    theme_cowplot() + xlab("Protein Number") +
    theme(text = element_text(size = 26),
          title = element_text(size = 16),
          axis.text = element_text(size = 32),
          axis.title = element_text(size = 32),
          axis.line = element_line(size = 2),
          axis.ticks = element_line(size = 2)) 

load("/lustre/user/liclab/fangq/proj/scp/scp/scope2_vector/seurat.RData",verbose = T)
scope2 = seur
ggplot(scope2@meta.data) + 
  geom_density(aes(x = nFeature_originalexp),fill = "grey",size = 1) + 
  # scale_fill_manual(values = alpha(brewer.pal(name = "Set2", n = 6),alpha = 0.6)) +
    geom_vline(xintercept = median(scope2@meta.data$nFeature_originalexp),size = 1.5,linetype = 2) +
    theme_cowplot() + xlab("Protein Number") +
    theme(text = element_text(size = 26),
          title = element_text(size = 16),
          axis.text = element_text(size = 32),
          axis.title = element_text(size = 32),
          axis.line = element_line(size = 2),
          axis.ticks = element_line(size = 2)) 
load("/lustre/user/liclab/fangq/proj/scp/pSCoPE/analysis/bmdm_vector/seurat.RData",verbose = T)
pscope = seur
ggplot(pscope@meta.data) + 
  geom_density(aes(x = nFeature_originalexp),fill = "grey",size = 1) + 
  # scale_fill_manual(values = alpha(brewer.pal(name = "Set2", n = 6),alpha = 0.6)) +
    geom_vline(xintercept = median(pscope@meta.data$nFeature_originalexp),size = 1.5,linetype = 2) +
    theme_cowplot() + xlab("Protein Number") +
    theme(text = element_text(size = 26),
          title = element_text(size = 16),
          axis.text = element_text(size = 32),
          axis.title = element_text(size = 32),
          axis.line = element_line(size = 2),
          axis.ticks = element_line(size = 2)) 



single_cell_stats = cellline@meta.data[,c("nFeature_originalexp","MedianCV")]
single_cell_stats$Data = "This study"
tmp = scope2@meta.data[,c("nFeature_originalexp","MedianCV")]
tmp$Data = "SCoPE2"

tmp2 = pscope@meta.data[,c("nFeature_originalexp","MedianCV")]
tmp2$Data = "pSCoPE"


single_cell_stats = rbind(rbind(single_cell_stats,tmp),tmp2)

data_summary <- function(x) { 
       m <- mean(x) 
       ymin <- m-sd(x) 
       ymax <- m+sd(x) 
       return(c(y=m,ymin=ymin,ymax=ymax)) 
}
protein_number = ggplot(single_cell_stats,aes(x = Data,y = nFeature_originalexp,fill = Data)) + 
  geom_violin(size = 1) + 
  scale_fill_manual(values = alpha(brewer.pal(name = "Set2", n = 6),alpha = 0.6)) +
  stat_summary(fun.data=data_summary,size = 1.2) + 
    theme_cowplot() + xlab(NULL) +ylab("Protein Number") +
    theme(text = element_text(size = 26),
          title = element_text(size = 16),
          legend.position = "none",
          axis.text = element_text(size = 32),
          axis.title = element_text(size = 32),
          axis.line = element_line(size = 2),
          axis.ticks = element_line(size = 2)) 

median_cv = ggplot(single_cell_stats,aes(x = Data,y = MedianCV,fill = Data)) + 
  geom_violin(size = 1) + 
  scale_fill_manual(values = alpha(brewer.pal(name = "Set2", n = 6),alpha = 0.6)) +
  stat_summary(fun.data=data_summary,size = 1.2) + 
    theme_cowplot() + xlab(NULL) +ylab("MedianCV") +
    theme(text = element_text(size = 26),
          title = element_text(size = 16),
          legend.position = "none",
          axis.text = element_text(size = 32),
          axis.title = element_text(size = 32),
          axis.line = element_line(size = 2),
          axis.ticks = element_line(size = 2)) 

top_tail= function(dataset,uniprot_path = character()){
  tmp = GetAssayData(dataset[["batchLraw"]])
  medianRI = apply(tmp,1,function(x){x = median(x,na.rm = T)}) %>% as.data.frame()
  colnames(medianRI) = "medianRI"
  medianRI$rank = rank((-1)*medianRI$medianRI)
  uniprot_search = readxl::read_xlsx( uniprot_path )
  colnames(uniprot_search)[1:2] = c("InputID","UniprotID")
  colnames(uniprot_search)=gsub(" ",".",colnames(uniprot_search))
  uniprot_search = uniprot_search%>% filter(InputID %in% rownames(tmp) )
  if(is_empty(grep("Gene.Names.(primary)",colnames(uniprot_search),value = T))){
    uniprot_search$`Gene.Names.(primary)`  = lapply( as.list(gsub(" ",",",uniprot_search$Gene.names)),function(x){ 
                                                x = strsplit(",",x = x) %>% unlist
                                                x = x[1] } ) %>% 
                                              unlist() 
  }
  medianRI = medianRI %>% rownames_to_column("InputID")
  proteins = intersect(medianRI$InputID,uniprot_search$InputID)
  medianRI= merge(uniprot_search[uniprot_search$InputID %in% proteins,],medianRI[medianRI$InputID %in% proteins,],by = "InputID")
  return(medianRI)
}

medianRI = top_tail(scope2,"/lustre/user/liclab/fangq/proj/scp/scp/scope2_vector/uniprot_search.xlsx")
medianRI = top_tail(pscope,"/lustre/user/liclab/fangq/proj/scp/pSCoPE/analysis/bmdm_vector/uniprot_bmdm_pros.xlsx")

Slope.cutoff = function(value = numeric(),slope.cutoff = 1,plot = logical()){
  ma = max(value,na.rm = T)
  mi = min(value,na.rm = T)
  n = length(value)
  value = apply(data.frame(value),1,function(x){(x-mi)/(ma-mi)})
  value = data.frame(y = sort(value),x = c(1:n)/n)  ## sort and rank
  value$intercept = value$y-slope.cutoff*value$x
  value$intercept2 = slope.cutoff*value$x-value$y

  cutoff.rank = rownames(value[which(value$intercept == min(value$intercept)),])
  cutoff.rank2 = rownames(value[which(value$intercept2 == min(value$intercept2)),])

  if(plot){
    value$type = ifelse(as.numeric(rownames(value)) >= as.numeric(cutoff.rank) ,"Significant","Non-signif")
    p = ggplot( value, aes(x = x,y = y,color = type) ) + geom_point(size = 1) + 
      geom_abline( slope = slope.cutoff,intercept = c(min(value$intercept),(-1)*min(value$intercept2)),
                   color = alpha("#e64047",alpha = 0.6),
                   show.legend = T ,linetype="dashed",size = 1) + 
      geom_hline( yintercept = c(value[rownames(value)==cutoff.rank,"y"],
                                 value[rownames(value)==cutoff.rank2,"y"]),
                  color = alpha("#42776a",alpha = 0.8),size = 1,
                  linetype="dashed" ) +
      scale_color_manual(values = alpha(c("grey","darkred"),alpha = 0.5),name = "Gene Type") + theme_cowplot()
    plot(p)
  }
  return(list(top = (nrow(value)-as.numeric(cutoff.rank)),
              tail = as.numeric(cutoff.rank2)))
}
Slope.cutoff(medianRI$medianRI,plot = T)


medianRI$cluster = ifelse(medianRI$rank<=Slope.cutoff(medianRI$medianRI,plot=F)[[1]],paste0("top",Slope.cutoff(medianRI$medianRI,plot=F)[[1]]),
                    ifelse(medianRI$rank<=(nrow(medianRI)-Slope.cutoff(medianRI$medianRI,plot=F)[[2]]),"middle",
                    paste0("tail",Slope.cutoff(medianRI$medianRI,plot=F)[[2]]))) %>% 
                    factor(.,levels = c(paste0("top",Slope.cutoff(medianRI$medianRI,plot=F)[[1]]),"middle",paste0("tail",Slope.cutoff(medianRI$medianRI,plot=F)[[2]])),ordered = T)

Protein_rank = ggplot(medianRI,aes(x=rank,y = medianRI,color = cluster)) + 
                geom_point(size= 7) +
                scale_color_manual(values = alpha(c("darkred","grey","navy"),0.6)) +
                    theme_cowplot() + xlab("Abndance Rank") +
                    theme(text = element_text(size = 26),
                          title = element_text(size = 16),
                          legend.title = element_blank(),
                          axis.text = element_text(size = 32),
                          axis.title = element_text(size = 32),
                          axis.line = element_line(size = 2),
                          axis.ticks = element_line(size = 2)) 
Protein_rank


species = "human"
species = "mouse"

database = function(species = character()){
  require(org.Hs.eg.db)
  require(org.Mm.eg.db)
  require(ensembldb)
  require(clusterProfiler)
  if(species=="human"){
    ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
    orgdb = org.Hs.eg.db
  }else if(species=="mouse"){
    ensdb = EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79
    orgdb = org.Mm.eg.db
  }
  wg.gene = data.frame(genes(ensdb)) %>%  dplyr::filter( gene_biotype == "protein_coding" )
  wg.gene = wg.gene$symbol
  return(list(wg.gene=wg.gene,orgdb=orgdb))
}


GO_Plot = function( GO_dataset, padj_threshold=0.05, show_num = 15, 
                    fill_color = "#e63f46",wid = 10, hgt = 8,
                    fontsize = 2,fontcolor = "black",
                    barwidth = 0.8,
                    GO_term = character(),
                    keywords = character(),
                    discard = character()){
                tmp = GO_dataset@result
                tmp = tmp[which(tmp$p.adjust<padj_threshold),] 
                tmp$logpval = (-log10( tmp$pvalue ))
                tmp  = tmp %>% arrange( dplyr::desc(logpval),dplyr::desc(Count)) 

                if(!is_empty( GO_term )){
                    tmp = tmp %>% filter( ONTOLOGY %in% GO_term )
                }
                if(!is_empty( keywords )){
                    seek = numeric()
                    for( i in keywords){
                    j = grep(i,tmp$Description,value = T)
                    seek = c(seek,j)
                    }
                    if(length(seek)==0){return("No aim GO terms was found!")}
                }
                if(show_num==0){
                    all_term = unique(seek)
                }else {
                    all_term = unique(c(tmp$Description[1:min(length(tmp$Description),show_num)],seek))
                }
                if(length(discard)!=0){
                    remove = numeric()
                    for( i in discard){
                    j = grep(i,tmp$Description,value = T)
                    remove = c(remove,j)
                    }
                    if(length(remove)==0){return("No discarded GO terms was found!")}
                    tmp = tmp[tmp$Description %in% setdiff(tmp$Description,remove),]
                }

                tmp = tmp[tmp$Description %in% all_term,] %>% arrange( dplyr::desc(logpval),dplyr::desc(Count)) 
                tmp$Description = factor(tmp$Description,levels = rev(tmp$Description))
                go = ggplot()+
                    geom_bar( data=tmp,
                            aes_string(x = "Description" ,y = "logpval"),
                            fill = fill_color,stat = "identity",width = barwidth) + ylab( "-log(p value)" ) +
                    geom_text( data=tmp,aes(x = Description ,y = 0.1, label = Description),
                                hjust = 0,size = fontsize,color = fontcolor,family = "sans" )+
                    xlab(NULL)+scale_y_continuous(expand = c(0,0)) + 
                    coord_flip( )+
                    # facet_grid( drop = T,rows = "ONTOLOGY",scales = "free",space = "free_y")+
                    theme_cowplot()+
                    theme( text = element_text( size = 30),
                        axis.title = element_text(size = 30),
                        axis.text.y = element_blank(),
                        axis.text.x= element_text(size = 30,family = "sans"),
                        axis.line = element_line(size = 1.5),
                        axis.ticks = element_line(size = 1.5),
                        axis.ticks.y = element_blank()

                    )
                return(go)
        }

top100 = medianRI$`Gene.Names.(primary)`[medianRI$rank<=Slope.cutoff(medianRI$medianRI,plot=F)[[1]]]
tail89 = medianRI$`Gene.Names.(primary)`[medianRI$rank>(nrow(medianRI)-Slope.cutoff(medianRI$medianRI,plot=F)[[2]])]

go_100 = enrichGO(top100,
                OrgDb = org.Mm.eg.db,
                keyType = "SYMBOL",ont = "BP",minGSSize = 5,
                pAdjustMethod = "BH",universe =database(species)[["wg.gene"]],
                pvalueCutoff = 0.05,qvalueCutoff = 0.05)
go_100_simplified = clusterProfiler::simplify(x = go_100 ,cutoff=0.8,by="p.adjust",select_fun=min)  #去除冗余，可以调整cutoff值。
go_89 = enrichGO(tail89,
                OrgDb = org.Mm.eg.db,
                keyType = "SYMBOL",ont = "BP",minGSSize = 5,
                pAdjustMethod = "BH",universe = database(species)[["wg.gene"]],
                pvalueCutoff = 0.05,qvalueCutoff = 0.05)
go_89_simplified = clusterProfiler::simplify(x = go_89 ,cutoff=0.7,by="p.adjust",select_fun=min)  #去除冗余，可以调整cutoff值。


p1 = GO_Plot(go_100,fontsize = 9,fill_color = "#f0cb5c",padj_threshold = 0.1,show_num = 10 )
p2 = GO_Plot(go_89_simplified,fontsize = 9,fill_color = "#f0cb5c",padj_threshold = 0.05,show_num = 10 )
p1
p2

p1 = GO_Plot(go_100,fontsize = 9,fill_color = "#f0cb5c",padj_threshold = 0.5,show_num = 10 )
p2 = GO_Plot(go_89,fontsize = 9,fill_color = "#f0cb5c",padj_threshold = 0.5,show_num = 10 )
p1
p2

pdf("top_100_go.pdf",width = 26,height = 12)
k='
AABBDDDD
#CC#EEEE
'
 protein_number + median_cv + Protein_rank + p1 + p2 + plot_layout(design = k)
dev.off()


pdf("scope2_go.pdf",width = 15,height = 14)
p1/p2
dev.off()

pdf("pscope_go.pdf",width = 13,height = 14)
p1/p2
dev.off()


DefaultAssay(seur)  = "batchLraw"
VariableFeatures( seur) = rownames(seur[['batchLraw']])
seur = ScaleData( seur )
# seur[["batchLnorm"]]@scale.data = seur[["batchLnorm"]]@data 
# seur = FindVariableFeatures( seur,selection.method = "vst" ) # default vst
# seur = ScaleData( seur,do.scale = F,do.center = F)
seur = Seurat::RunPCA( seur,assay = "batchLraw",npcs = 30)
ElbowPlot(seur,ndims = 30)




seur = Seurat::RunPCA( seur,assay = "batchLraw",npcs = 30)
# seur = Seurat::RunICA( seur,nics = 20)
seur = RunUMAP( seur,dims = 1:30,reduction = "pca",n.neighbors = 30,min.dist = 0.05 )#n.neighbors = 30,min.dist = 0.05 
# seur = RunTSNE( seur,dims = 1:20,reduction = "pca" )
DimPlot( seur,reduction = "pca",group.by = "celltype",pt.size = 4 )
DimPlot( seur,reduction = "umap",group.by = "celltype",pt.size = 4 )


library(clustree)

seur = FindNeighbors(object = seur,assay = "batchLraw")
cluster.test = FindClusters(object = seur,resolution = seq(0,2,by = 0.2))
clustree(cluster.test@meta.data, prefix = "batchLnorm_snn_res.") ### resolution = 0.4
seur = FindClusters(object = seur,resolution = 1)
# seur = FindNeighbors(seur, dims = 1:20)
# seur = FindClusters(seur, resolution = 0.8)
Idents(seur) = "seurat_clusters"
table(seur$seurat_clusters,seur$celltype)
signif_marker_pros = Seurat::FindAllMarkers( seur,
                                             test.use = "wilcox",
                                             assay = "batchLraw",
                                             min.cells.feature = 10,
                                             logfc.threshold = 0,only.pos = F,
                                             return.thresh = 0.05 )
display_seurat = signif_marker_pros %>%  dplyr::filter( p_val_adj < 0.05 ) %>% 
               group_by( cluster ) %>% top_n( wt = abs(avg_log2FC),n = 8 )

seur = ScaleData( seur,assay = "weighted_cor_imputed_batchC",do.scale = F,do.center = F)

VariableFeatures(seur[['weighted_cor_imputed_batchC']]) = rownames(seur[['weighted_cor_imputed_batchC']])

seur = Seurat::RunPCA( seur,assay = "weighted_cor_imputed_batchC",npcs = 50)
pca_weighted = DimPlot( seur,reduction = "pca",group.by = c("celltype"),
              pt.size = 3,
              cols = alpha(c("#8ac926","#6a4c93"),0.6)) + 
              theme(axis.text = element_text(size = 26),
              axis.line = element_line(size = 2),
              axis.ticks = element_line(size = 2),
              text = element_text(size = 26)) + ggtitle(NULL)
pca_weighted_cluster = DimPlot( seur,reduction = "pca",group.by = c("seurat_clusters"),
              pt.size = 3,
              cols = alpha(c( "0" = "#e97d72",
          "1" = "#bc7ff8",
          "2" = "#56bcc2"),0.6)) + 
              theme(axis.text = element_text(size = 26),
              axis.line = element_line(size = 2),
              axis.ticks = element_line(size = 2),
              text = element_text(size = 26)) + ggtitle(NULL)

cols =  c("HeLa"="#8ac926","U2OS"="#6a4c93")

cols2 = c( "0" = "#e97d72",
          "1" = "#bc7ff8",
          "2" = "#56bcc2")
  #   HeLa U2OS
  # 0   20   21
  # 1   36    3
  # 2    1   25
meta = seur@meta.data 
meta$seurat_clusters = factor(meta$seurat_clusters,levels = c(1,0,2),ordered = T)
meta = meta %>% arrange( seurat_clusters,celltype )
k = seur[["batchLraw"]]@scale.data[ unique(signif_marker_pros$gene),rownames(meta)]
d = dist(k,method = "euclidean")
#然后进行分层聚类(Hierarchical cluster)；
cluster2 <- hclust(d, method = "ward.D2")
#绘制聚类树；
plot(cluster2,hang = -1,cex=0.6,axes=FALSE,ann=FALSE)
cut = cutree(cluster2,3)

c1 = cut[cut==1] %>% names()
gene1 = gene_pro_ref$gene[gene_pro_ref$proteins%in%c1]
c2 = cut[cut==2] %>% names()
gene2 = gene_pro_ref$gene[gene_pro_ref$proteins%in%c2]
c3 = cut[cut==3] %>% names()
gene3 = gene_pro_ref$gene[gene_pro_ref$proteins%in%c3]

go1 = enrichGO( gene1,
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",ont = "ALL",minGSSize = 5,
                pAdjustMethod = "BH",universe = wg.gene,
                pvalueCutoff = 0.05,qvalueCutoff = 0.05)

go2 = enrichGO( c(gene1,gene2),
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",ont = "ALL",minGSSize = 5,
                pAdjustMethod = "BH",universe = wg.gene,
                pvalueCutoff = 0.05,qvalueCutoff = 0.05)

GO_Plot = function( GO_dataset, padj_threshold=0.05, show_num = 15, 
                    fill_color = "#e63f46",wid = 10, hgt = 8,
                    output_prefix,GO_term = character(),
                    keywords = character(),go_type = character()){
  tmp = GO_dataset
  if(!is_empty( GO_term )){
    tmp = tmp %>% filter( ONTOLOGY %in% GO_term )
  }
  if(length(keywords)!=0){
    seek = numeric()
    for( i in keywords){
      j = grep(i,tmp$Description)
      seek = c(seek,j)
    }
    if(length(seek)==0){return("No aim GO terms was found!")}
    tmp = tmp[seek,]
  }
  if(!is_empty(go_type)){
    tmp = tmp[tmp$ONTOLOGY %in% go_type,]
  }
  tmp = tmp[which(tmp$p.adjust<padj_threshold),] 
  tmp$logpval = (-log10( tmp$pvalue ))
  tmp  = tmp %>% arrange( dplyr::desc(logpval),dplyr::desc(Count)) 
  tmp = tmp[1:min(show_num,nrow(tmp)),] 
  tmp$Description = factor(tmp$Description,levels = rev(tmp$Description))
  
  go = ggplot()+
    geom_bar( data=tmp,
              aes_string(x = "Description" ,y = "logpval"),
              fill = fill_color,
              stat = "identity",width = 0.5) + ylab( "-log(p.value)" ) +
    coord_flip( )+
    facet_grid( drop = T,rows = "ONTOLOGY",scales = "free",space = "free_y")+
    theme_cowplot()+
    theme( text = element_text( size = 28),
           axis.title = element_text(size = 28),
           axis.text= element_text(size = 15,family = "sans"),
           line = element_line(size = 1),
    )
  go 
  #ggsave( go,filename = paste0( output_prefix,".pdf" ),width = wid, height = hgt, units = "in" )
}

go_plot = GO_Plot(go2,go_type = c("BP","MF"),fill_color = "#f19f39",
                  keywords = c("amide","pyruvate","oxidative","hypoxia","ubiquitination","posttranscriptional")) + 
              theme(axis.text.x = element_text(size = 26),
              axis.line = element_line(size = 2),
              axis.ticks = element_line(size = 2),
              text = element_text(size = 14)) + ggtitle(NULL)
#cytoplasmic


top = HeatmapAnnotation(bar = meta$celltype,
                        foo = meta$seurat_clusters,
                        col = list(bar = cols,foo = cols2),
                        name = "Type",
                        annotation_label = "Type",
                        show_legend = T,
                        show_annotation_name = F,
                        annotation_legend_param = list(labels_gp  = gpar(fontsize = 15),
                                                       title_gp = gpar(fontsize = 15))
)
summary(k)
k[k>quantile(k,0.99)] = quantile(k,0.99)
k[k<quantile(k,0.01)] = quantile(k,0.01)


a = Heatmap( k, 
         # heatmap_legend_param = list( title = "InsulationScore",
         #                              at = c(bot,0,(-1)*bot),
         #                              labels = c(round(bot,2),0,(-1)*round(bot,2)),
         #                              labels_gp  = gpar(fontsize = 15),
         #                              direction = "vertical",
         #                              title_gp = gpar(fontsize = 15)),
         show_heatmap_legend = T,,
         col = colorRampPalette(rev(brewer.pal("RdBu",n=10)))(50),
         cluster_columns = F,row_km = 3,
         cluster_rows = T ,
         show_column_names = F,
         clustering_method_rows = "ward.D2",#ward.D2
         clustering_distance_rows = "pearson",
         #("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"
         show_row_names = F,
         top_annotation = top,
         # left_annotation = left,
         column_names_rot = 45,
         column_names_gp = gpar(fontsize = 15) )
a
pdf("heatmap.pdf",width = 10,height = 4)
a
dev.off()

DefaultAssay(seur) = "GOterm" 
# seur = FindVariableFeatures( seur,selection.method = "vst",assay = "GOterm" ) # default vst,based on data not scale.data
# seur[["GOterm"]]@scale.data = seur[["GOterm"]]@data %>% t() %>% scale(.,center = T,scale = T) %>% t()
seur = ScaleData( seur,do.scale = T,do.center = T)
VariableFeatures(seur) = rownames(seur[['GOterm']])
# seur = Seurat::RunPCA( seur,assay = "GOterm",npcs = 20)
# ElbowPlot(seur)
# seur = FindNeighbors(seur, dims = 1:20,assay = "GOterm")
# seur = FindClusters(seur, resolution = 1)
# Idents(seur)
# seur = RunUMAP( seur,dims = 1:20,reduction = "pca",n.neighbors = 30,min.dist = 0.05 )
# 
# DimPlot( seur,reduction = "pca",group.by = c("celltype","seurat_clusters") )
# DimPlot( seur,reduction = "umap",pt.size = 2.5,group.by = c("seurat_clusters","celltype"))
# DimPlot( seur,reduction = "tsne",pt.size = 2) 

Idents(seur) = "seurat_clusters"
signif_marker_go = Seurat::FindAllMarkers( seur,slot = "data",
                                           features = VariableFeatures(seur[["GOterm"]]),
                                           test.use = "wilcox",
                                           assay = "GOterm",
                                           min.cells.feature = 10,
                                           logfc.threshold = 0,min.pct = 0.1 )
signif_marker_go = merge( signif_marker_go,go.filtered[,c("ID","Description","pros_num","bgratio","ONTOLOGY")], 
                          by.x = "gene",by.y = "ID",all.x = T,all.y = F) %>% 
                          arrange(cluster)
hela_go = c("GO:1901796","GO:1904667","GO:1990948","GO:0002181","GO:0044389","GO:0061684","GO:0034620","GO:0061077","GO:0046034")                      
u2os_go = c("GO:0005200","GO:0030506","GO:0098927","GO:0045022","GO:0006900")

display_go = signif_marker_go %>%  dplyr::filter( p_val_adj < 0.05 & gene %in% c(hela_go,u2os_go))
p1 = DimPlot( seur,reduction = "pca",group.by = "celltype",pt.size = 4 ) & 
  scale_colour_manual( values = alpha(cols,alpha = 0.7) )
p1

lim = min(abs(quantile(seur[["GOterm"]]@scale.data,0.99,na.rm = T)),abs(quantile(seur[["GOterm"]]@scale.data,0.01,na.rm = T)))
p_hela = FeaturePlot( seur,
             slot = "scale.data",
             order = T,
             # cols = alpha(rev(brewer.pal("RdBu",n = 3)), 0.6 ),
             reduction = "pca",
             features = hela_go,
             ncol = 5,pt.size = 2) & 
  scale_colour_gradientn(colours = alpha(rev(brewer.pal(n = 11, name = "RdBu")),alpha = 0.8),
                          limits = c((-1)*1.5,1.5)
                         ) &
  guides( colour = guide_colorbar(paste0("z-scored","\n","abundance")))
# titles = as.list( unique(display_go$gene) )
# names(titles) = unique(display_go$gene)
# titles = lapply(titles, function(x){x = go.filtered$Description[go.filtered$ID==x]})

p_hela$patches$plots = lapply(p_hela$patches$plots , function(x){
  name = go.filtered$Description[go.filtered$ID==x$labels$title]
  tmp = strsplit(name,split = " ")
  length = length(tmp[[1]])
  char = tmp %>% unlist
  if(length<=4){
    name.new = name
  }else{
    name.new = paste0(paste(char[1:4],collapse = " "),"\n",paste(char[5:length],collapse = " "))
  }
  x = x + ggtitle( name.new )
})
p_hela$labels$title = go.filtered$Description[go.filtered$ID==p_hela$labels$title]

p_hela = p_hela & theme(title = element_text(size = 18),
                        axis.title = element_text(size = 24),
                        axis.text = element_text(size = 26),
                        axis.line = element_line(size = 2),
                        axis.ticks = element_line(size = 2)
                        ) 
p_hela_key = p_hela + plot_layout( guides = 'collect' )




p_u2os = FeaturePlot( seur,
             slot = "scale.data",
             order = T,
             # cols = alpha(rev(brewer.pal("RdBu",n = 3)), 0.6 ),
             reduction = "pca",
             features = u2os_go,
             ncol = 5,pt.size = 2) & 
  scale_colour_gradientn(colours = alpha(rev(brewer.pal(n = 11, name = "RdBu")),alpha = 0.8),
                          limits = c((-1)*1.5,1.5)
                         ) &
  guides( colour = guide_colorbar(paste0("z-scored","\n","abundance")))

p_u2os$patches$plots = lapply(p_u2os$patches$plots , function(x){
  name = go.filtered$Description[go.filtered$ID==x$labels$title]
  tmp = strsplit(name,split = " ")
  length = length(tmp[[1]])
  char = tmp %>% unlist
  if(length<=4){
    name.new = name
  }else{
    name.new = paste0(paste(char[1:4],collapse = " "),"\n",paste(char[5:length],collapse = " "))
  }
  x = x + ggtitle( name.new )
})
p_u2os$labels$title = go.filtered$Description[go.filtered$ID==p_u2os$labels$title]

p_u2os = p_u2os & theme(title = element_text(size = 18),
                        axis.text = element_text(size = 26),
                        axis.title = element_text(size = 24),
                        axis.line = element_line(size = 2),
                        axis.ticks = element_line(size = 2)
                        ) 
p_u2os_key = p_u2os + plot_layout( guides = 'collect' )

pdf("GO_term.pdf",width = 28,height =20 )
k='
##AAABBB##
CCCCCCCCCC
CCCCCCCCCC
DDDDDDDDDD
'
pca_weighted + pca_weighted_cluster + p_hela_key + p_u2os_key + plot_layout(design = k,height=c(1.5,1,1,1))
dev.off()
pdf("GO.pdf",width = 15,height =6 )
go_plot
dev.off()

###### examples
seur$clusters  = ifelse(seur$seurat_clusters==1,"HeLa_sub",
                      ifelse(seur$seurat_clusters==2,"U2OS_sub","mixed"))
Idents(seur) = "clusters"
hela_2 = c("GO:1901796","GO:0061684")                 
u2os_2 = c("GO:0030506","GO:0098927")
lim = min(abs(quantile(seur[["GOterm"]]@scale.data,0.99,na.rm = T)),abs(quantile(seur[["GOterm"]]@scale.data,0.01,na.rm = T)))

double_plot = function(term = character()){
  hela_p = FeaturePlot( seur,
             slot = "scale.data",
             order = T,
             # cols = alpha(rev(brewer.pal("RdBu",n = 3)), 0.6 ),
             reduction = "pca",
             features = term,
             ncol = 1,pt.size = 4) & 
  scale_colour_gradientn(colours = alpha(rev(brewer.pal(n = 11, name = "RdBu")),alpha = 0.8),
                          limits = c((-1)*1.5,1.5)
                         ) &
  guides( colour = guide_colorbar(paste0("z-scored","\n","abundance")))
# titles = as.list( unique(display_go$gene) )
# names(titles) = unique(display_go$gene)
# titles = lapply(titles, function(x){x = go.filtered$Description[go.filtered$ID==x]})
hela_p$patches$plots = lapply(hela_p$patches$plots , function(x){
  name = go.filtered$Description[go.filtered$ID==x$labels$title]
  tmp = strsplit(name,split = " ")
  length = length(tmp[[1]])
  char = tmp %>% unlist
  if(length<=4){
    name.new = name
  }else{
    name.new = paste0(paste(char[1:4],collapse = " "),"\n",paste(char[5:length],collapse = " "))
  }
  x = x + ggtitle( name.new )
})
hela_p$labels$title = go.filtered$Description[go.filtered$ID==hela_p$labels$title]

hela_p = hela_p & theme(title = element_text(size = 18,hjust = 0.5,vjust = 0.5),
                        axis.title = element_text(size = 30),
                        legend.text = element_text(size = 22),
                        axis.text = element_text(size = 30),
                        axis.line = element_line(size = 2),
                        axis.ticks = element_line(size = 2)
                        ) 

data_summary <- function(x) {
   m <- mean(x)
   ymin <- m-sd(x)
   ymax <- m+sd(x)
   return(c(y=m,ymin=ymin,ymax=ymax))
}

hela_v = Seurat::VlnPlot(seur,slot = "scale.data",pt.size = 0,group.by = "clusters",
             features = term,idents = c("HeLa_sub","U2OS_sub"),adjust = 1,
             cols = alpha(c("#8ac926","#6a4c93"),0.7),
             ncol = 1) & stat_summary(fun.data=data_summary,size = 1.2) &
             theme(title = element_text(size = 18),
             legend.position = "none",
                         axis.title = element_text(size = 30),
                        legend.text = element_text(size = 22),
                        axis.text = element_text(size = 30),
                        axis.line = element_line(size = 2),
                        axis.ticks = element_line(size = 2)
                        ) 
k = '
AB
'
p_final = (hela_p) + (hela_v) + plot_layout(design = k,guides = 'collect')
return(p_final)
######
}

hela_p1 = double_plot(hela_2[1]) %>% ggplotify::as.ggplot()
hela_p2 = double_plot(hela_2[2])%>% ggplotify::as.ggplot()
u2os_p1 = double_plot(u2os_2[1])%>% ggplotify::as.ggplot()
u2os_p2 = double_plot(u2os_2[2])%>% ggplotify::as.ggplot()
pdf("double_plot.pdf",width = 28,height = 16)
k = '
AC
BD
'
hela_p1 + hela_p2 + u2os_p1 + u2os_p2 + plot_layout(design = k,guides = 'collect',widths = c(14,14))
dev.off()


Idents(seur) = "celltype"
signif_marker_pros = Seurat::FindAllMarkers( seur,
                                             slot = "data",
                                             test.use = "wilcox",
                                             assay = "batchLnorm",
                                             min.cells.feature = 20,
                                             logfc.threshold = 0.5,only.pos = T,
                                             return.thresh = 0.05 )
display_pros = signif_marker_pros %>%  dplyr::filter( p_val_adj < 0.05 ) %>% 
               group_by( cluster ) %>% top_n( wt = avg_log2FC,n = 20 )

go_pros = signif_marker_pros %>%  dplyr::filter( p_val < 0.05 ) %>% 
               group_by( cluster )


pdf("heat_test.pdf",width = 10,height = 6)                     
DoHeatmap( seur, features = signif_marker_pros$gene,group.by = "celltype",slot = "data",disp.min = -1.5,disp.max = 1.5)
# DoHeatmap( seur, features = display_seurat$gene,group.by = "seurat_clusters",slot = "data",
#            disp.min = -1.5,disp.max = 1.5)

k = seur[["batchLnorm"]]@scale.data[ signif_marker_pros$gene,]
top = HeatmapAnnotation(bar = seur@meta.data$celltype,
                        col = list(bar = cols),
                        name = "Type",
                        annotation_label = "Type",
                        show_legend = T,
                        show_annotation_name = F,
                        annotation_legend_param = list(labels_gp  = gpar(fontsize = 15),
                                                       title_gp = gpar(fontsize = 15))
)
summary(k)
k[k>quantile(k,0.95)] = quantile(k,0.95)
k[k<quantile(k,0.05)] = quantile(k,0.05)


Heatmap( k, 
         # heatmap_legend_param = list( title = "InsulationScore",
         #                              at = c(bot,0,(-1)*bot),
         #                              labels = c(round(bot,2),0,(-1)*round(bot,2)),
         #                              labels_gp  = gpar(fontsize = 15),
         #                              direction = "vertical",
         #                              title_gp = gpar(fontsize = 15)),
         show_heatmap_legend = T,
         col = c("navy","#3B4CC0", "#4F6AD9", "#6585EC", "#7B9FF9", "#93B5FF", 
                 "#AAC7FD", "#C0D4F5", "#D4DBE6", "#E5D8D1", "#F2CBB7",
                 "#F7B89C", "#F6A081", "#EE8568",
                 "#E0654F","#CC4039", "#B40426"),
         cluster_columns = T,row_km = 4,
         cluster_rows = T ,
         show_column_names = F,
         clustering_method_rows = "complete",#ward.D2
         clustering_distance_rows = "pearson",
         #("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"
         show_row_names = F,
         top_annotation = top,
         # left_annotation = left,
         column_names_rot = 45,
         column_names_gp = gpar(fontsize = 15) )


dev.off()

pca_batchLraw = DimPlot( seur,reduction = "pca",group.by = "celltype",pt.size = 3,
              cols = alpha(c("#8ac926","#6a4c93"),0.6))

DefaultAssay(seur)  = "imputedraw"
seur = FindVariableFeatures( seur,selection.method = "vst" ) 
seur = ScaleData( seur )
# seur[["batchLnorm"]]@scale.data = seur[["batchLnorm"]]@data 
# seur = FindVariableFeatures( seur,selection.method = "vst" ) # default vst
# seur = ScaleData( seur,do.scale = F,do.center = F)
seur = Seurat::RunPCA( seur,assay = "batchLnorm",npcs = 30)
ElbowPlot(seur,ndims = 30)
# VlnPlot( seur,features = c("nCount_batchLnorm","nFeature_batchLnorm"),group.by = "Set" )

seur = Seurat::RunPCA( seur,assay = "batchLnorm",npcs = 30)
pca_imputedraw = DimPlot( seur,reduction = "pca",group.by = "celltype",pt.size = 3,
              cols = alpha(c("#8ac926","#6a4c93"),0.6))

seur = ScaleData( seur,assay = "weighted_cor_imputed_batchC",do.scale = F,do.center = F)

VariableFeatures(seur[['weighted_cor_imputed_batchC']]) = rownames(seur[['weighted_cor_imputed_batchC']])

seur = Seurat::RunPCA( seur,assay = "weighted_cor_imputed_batchC",npcs = 50)
seur = Seurat::RunICA( seur,nics = 20)
seur = RunUMAP( seur,dims = 1:20,reduction = "pca",n.neighbors = 30,min.dist = 0.05 )#n.neighbors = 30,min.dist = 0.05 
seur = RunTSNE( seur,dims = 1:20,reduction = "pca" )


pca_weighted = DimPlot( seur,reduction = "pca",group.by = "celltype",pt.size = 3,
              cols = alpha(c("#8ac926","#6a4c93"),0.6))

pca_batchLnorm
pca_weighted
