#!/usr/bin/env Rscript
library(optparse)
library(XML)
option_list = list(
  make_option(c("-i", "--filePath"), type="character", default=NULL,
              help="dataset file path , must be full path, must seperated by ',' if multiple paths are suppied", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output file path", metavar="character"),
  make_option(c("-t", "--numThreads"), type="integer", default=10, 
              help="number of Threads/Cores to run Maxquant. \n\t\t Default: %default", metavar="character"),
  # make_option(c("-t", "--taxonomy"), type="character", default="Human", 
  #             help="species, default Human", metavar="character"),
  make_option(c("-r", "--reference"), type="character", 
              default="/home/fangq/ffqq/ref/UniProt/uniprot-compressed_true_download_true_format_fasta_query__28Human_29-2023.03.27-10.36.15.54.fasta", 
              help="reference protein FASTA data. \n\t\t Default: %default", metavar="character"),
  make_option(c("-c", "--contaminant"), type="character", default="/lustre/user/liclab/fangq/tools/MaxQuant_v2.3.1.0/bin/conf/contaminants.fasta", 
              help="contaminant FASTA file, better to be concordant with your Maxquant version. \n\t\t Default: %default", metavar="character"),
  make_option(c("-x", "--xmlPath"), type="character", default="./", 
              help="mqpar.xml location", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$filePath)|is.null(opt$output)){
  print_help(opt_parser)
  stop("Inputs and Outputs must be supplied", call.=FALSE)
}


opt$xmlPath = ifelse(stringr::str_sub(opt$xmlPath,start = -1,end = -1) =="/",
                     stringr::str_sub(opt$xmlPath,start = 1,end = -2),opt$xmlPath)

if(purrr::is_empty(grep(opt$filePath,pattern = ","))){
  opt$filePath = ifelse(stringr::str_sub(opt$filePath,start = -1,end = -1) =="/",
                       stringr::str_sub(opt$filePath,start = 1,end = -2),opt$filePath)
  opt$filePath = list.files(opt$filePath,pattern = ".raw",full.names = T)
}else{
  opt$filePath = unlist(lapply(as.list(strsplit(opt$filePath,split = ",")[[1]]),
                               function(x){
                                 x = ifelse(stringr::str_sub(x,start = -1,end = -1) =="/",
                                            stringr::str_sub(x,start = 1,end = -2),x)
                                 x = list.files(path = x,pattern = ".raw",full.names = T)
                                 return(x)
                                 }))
}

opt$output = ifelse(stringr::str_sub(opt$output,start = -1,end = -1) =="/",
                    stringr::str_sub(opt$output,start = 1,end = -2),opt$output)

# print(opt)

file_num = length(opt$filePath)

mqpar = xmlParse("/lustre/user/liclab/fangq/tools/MaxQuant_v2.3.1.0/mqpar_2.3.xml",encoding = 'UTF-8')
# class(mqpar)
xmltop = xmlRoot(mqpar)
# class(xmltop)
# xmlName(xmltop[[1]])
# xmlSize(xmltop)
# getNodeSet(xmltop,path = "/MaxQuantParams/fastaFiles/FastaFileInfo/fastaFilePath")[[1]][[1]]


xmltop[['fastaFiles']][['FastaFileInfo']][['fastaFilePath']][[1]]=opt$reference
xmltop[['fastaFilesFirstSearch']][['FastaFileInfo']][['fastaFilePath']][[1]]=opt$contaminant
xmltop[['numThreads']][[1]]=opt$numThreads

for( i in 1:file_num){
  if(i==1){
    xmltop[['filePaths']][[1]][[1]] = opt$filePath[i]
    xmltop[['experiments']][[1]][[1]] = ""
    xmltop[['referenceChannel']][[1]][[1]] = ""
  }else{
    newXMLNode("string",opt$filePath[i], parent = xmltop[['filePaths']])
    newXMLNode("short",32767, parent = xmltop[['fractions']])
    newXMLNode("boolean","False", parent = xmltop[['ptms']])
    newXMLNode("int",0, parent = xmltop[['paramGroupIndices']])
    newXMLNode("string","",parent = xmltop[['referenceChannel']])
    newXMLNode("string","",parent = xmltop[['experiments']])
  }
}

xmltop[['fixedCombinedFolder']][[1]] = opt$output

saveXML(xmltop,file = paste0(opt$xmlPath,"/mqpar.xml"),prefix = '<?xml version="1.0" encoding="utf-8"?>\n',encoding = "UTF-8")
cat(paste0("mqpar.xml has been saved in ",paste0(opt$xmlPath,"/mqpar.xml")),fill = T)
