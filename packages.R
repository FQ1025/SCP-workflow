if(pipe=="Seurat")
{
  # a = c('Seurat','survival','leiden','spatstat.sparse','BiocNeighbors','scater',
  #       'uwot','reticulate','scuttle','BiocSingular','mgcv','beachmat',
  #       'DelayedMatrixStats','rsvd','spatstat.explore','sctransform',
  #       'sparseMatrixStats','SummarizedExperiment','spatstat.data',
  #       'ScaledMatrix','DelayedArray','irlba','SeuratObject')
  pkgs =  c( 'scater', 'SingleCellExperiment', 'uwot', 'rtracklayer', # 'scran', 
             'GenomicFeatures', 'GenomicRanges', 'IRanges', 'clusterProfiler', 'rliger','SeuratWrappers',#'install',
             'ComplexHeatmap', 'ggsci', 'RColorBrewer', 'ggplotify', 'ggpubr', 'scales', 'cowplot', 'patchwork', 'tidyr', 'clustree',
             'VennDiagram', 'ggvenn', #'ggVennDiagram', 
             'tidyverse', 'plyr','data.table', 'dplyr', #'Seurat',
             'umap', 'reshape2' )
}else if( pipe=='scp' ){
  pkgs =  c( 'jsonlite', 'languageserver','IRanges','GenomicRanges','ComplexHeatmap',
             'pheatmap', 'data.table', 'tidyverse', 'magrittr', 'RColorBrewer', 'cowplot', 'patchwork', 'ggridges',
             'scales', 'textshape', 'stringr',"scuttle","SingleCellExperiment",
             'impute', 'scater', 'sva', 'GGally','limma',
             'scpdata', 'QFeatures', 'scp', 'Seurat')
}



Loadpkgs = function(i){
  # if(!require( i, quietly = TRUE )){
  #  BiocManager::install( i,update = F )
  # }
  library( i, character.only=TRUE )
}
sapply( pkgs,Loadpkgs )

