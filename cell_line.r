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

#####load evidence_updated data and generate matadata(****meta need to be arranged by RI1-16*****)
evidence = fread("/home/fangq/ffqq/proj/scp/cellLine/allfile/dart_out/ev_updated.txt",nThread = 10,header = T)
colnames(evidence)=gsub(" ",".",colnames(evidence))
colnames(evidence) = gsub("Reporter.intensity.corrected.","RI",colnames(evidence))
# evidence$Raw.file = gsub("_C_","_A2_",evidence$Raw.file)
# evidence$Raw.file = gsub("_B_","_B1_",evidence$Raw.file)

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

######scp analysis

scp.cellline = readSCP(featureData = evidence,
              colData = metadata,
              batchCol = "Raw.file", 
              channelCol = "tmt.label", 
              removeEmptyCols = TRUE)
scp.cellline = zeroIsNA(scp.cellline,i = names(scp.cellline))
#  filter low abundance runs
n.psm = dims(scp.cellline)[1, ] ;summary(n.psm)       # psm counts
n.runs = dims(scp.cellline)[2, ] ;summary(n.runs)      # run tmt counts

scp.cellline = pep2qvalue(scp.cellline,
                  i = names(scp.cellline), ## need to be presnt with single cell runs
                  PEP = "dart_PEP",
                  rowDataName = "qvalue_PSMs")

scp.cellline = pep2qvalue(scp.cellline,
                  i = names(scp.cellline), ## need to be presnt with single cell runs
                  PEP = "dart_PEP",
                  groupBy = "Leading.razor.protein",
                  rowDataName = "qvalue_proteins")

scp.stats(scp.cellline);dims(scp.cellline)
stats0 = dims(scp.cellline)[1,]
tab1 = dims(scp.cellline) %>% as.data.frame() %>% t 
sum1 = data.frame( cells = 231,peptides = 89173,proteins = 29058 ) %>% t
colnames(tab1) = c("PSM1","Cell1")
colnames(sum1) = c("Stats1")
# [1] "all.runs = 22, sc.runs = 21, cells = 231, peptides = 89173, protein = 29058"
#      20220629_SCP_B1_10 20220629_SCP_B1_13 20220704_SCP_B_7 20220704_SCP_B_8 20220704_SCP_B_9 20220704_SCP_C_10 20220704_SCP_C_7 20220704_SCP_C_9
# [1,]               8252               9396             7627             9141             9235              8207             9418             9152
# [2,]                 16                 16               16               16               16                16               16               16
#      20220704_SCP_HA_B_11 20220704_SCP_HA_B_6 20220704_SCP_HA_B_7 20220704_SCP_HA_S_12 20220704_SCP_LA_A1_11 20220705_SCP_HA_S_11 20220707_SCP_HA_B_12
# [1,]                 6059                9203                8646                10276                 10093                 7726                 9708
# [2,]                   16                  16                  16                   16                    16                   16                   16
#      20220707_SCP_HA_B_5 20220708_SCP_LA_B_10 20220708_SCP_LA_B2_6 20220712_SCP_HA_A3_5 20220712_SCP_HA_A3_6 20220712_SCP_HA_A3_8 20220712_SCP_HA_C_R
# [1,]                5364                 4863                10710                 8286                10074                 7732               12381
# [2,]                  16                   16                   16                   16                   16                   16                  16
scp.cellline = filterFeatures(scp.cellline, ~ qvalue_proteins < 0.01 &  qvalue_PSMs < 0.01  )
scp.stats(scp.cellline);dims(scp.cellline)
stats1 = dims(scp.cellline)[1,]
tab2 = dims(scp.cellline) %>% as.data.frame() %>% t 
sum2 = data.frame( cells = 231,peptides = 9471,proteins = 3423 ) %>% t
colnames(tab2) = c("PSM2","Cell2")
colnames(sum2) = c("Stats2")
# [1] "all.runs = 22, sc.runs = 21, cells = 231, peptides = 9471, protein = 3423"
#      20220629_SCP_B1_10 20220629_SCP_B1_13 20220704_SCP_B_7 20220704_SCP_B_8 20220704_SCP_B_9 20220704_SCP_C_10 20220704_SCP_C_7 20220704_SCP_C_9
# [1,]               1026               2688             1933             3286             3661              1230             3587             3063
# [2,]                 16                 16               16               16               16                16               16               16
#      20220704_SCP_HA_B_11 20220704_SCP_HA_B_6 20220704_SCP_HA_B_7 20220704_SCP_HA_S_12 20220704_SCP_LA_A1_11 20220705_SCP_HA_S_11 20220707_SCP_HA_B_12
# [1,]                  296                2340                 802                 3099                  4496                 1978                 2764
# [2,]                   16                  16                  16                   16                    16                   16                   16
#      20220707_SCP_HA_B_5 20220708_SCP_LA_B_10 20220708_SCP_LA_B2_6 20220712_SCP_HA_A3_5 20220712_SCP_HA_A3_6 20220712_SCP_HA_A3_8 20220712_SCP_HA_C_R
# [1,]                  22                   60                 4102                 2323                 3324                  906                5645
# [2,]                  16                   16                   16                   16                   16                   16                  16
scp.cellline = filterFeatures(scp.cellline, ~ Reverse != "+" & Potential.contaminant != "+" & !is.na(PIF) & PIF > 0.5 & !grepl("REV|CON", Leading.razor.protein))
scp.stats(scp.cellline);dims(scp.cellline)
stats2 = dims(scp.cellline)[1,]
tab3 = dims(scp.cellline) %>% as.data.frame() %>% t 
sum3 = data.frame( cells = 231,peptides = 5459,proteins = 1523 ) %>% t
colnames(tab3) = c("PSM3","Cell3")
colnames(sum3) = c("Stats3")
# [1] "all.runs = 22, sc.runs = 21, cells = 231, peptides = 5459, protein = 1523"
#      20220629_SCP_B1_10 20220629_SCP_B1_13 20220704_SCP_B_7 20220704_SCP_B_8 20220704_SCP_B_9 20220704_SCP_C_10 20220704_SCP_C_7 20220704_SCP_C_9 20220704_SCP_HA_B_11 20220704_SCP_HA_B_6
# [1,]                447               1440             1011             1931             2019               680             1942             1646                  110                1151
# [2,]                 16                 16               16               16               16                16               16               16                   16                  16
#      20220704_SCP_HA_B_7 20220704_SCP_HA_S_12 20220704_SCP_LA_A1_11 20220705_SCP_HA_S_11 20220707_SCP_HA_B_12 20220707_SCP_HA_B_5 20220708_SCP_LA_B_10 20220708_SCP_LA_B2_6 20220712_SCP_HA_A3_5
# [1,]                 317                 1504                  2748                  997                 1380                   3                   19                 2288                 1108
# [2,]                  16                   16                    16                   16                   16                  16                   16                   16                   16
#      20220712_SCP_HA_A3_6 20220712_SCP_HA_A3_8 20220712_SCP_HA_C_R
# [1,]                 1714                  345                3612
# [2,]                   16                   16                  16
save(scp.cellline,stats0,stats1,stats2,file = "raw_scp.RData",compress = T)
load("raw_scp.RData",verbose = T)
scp.cellline = computeSCR(scp.cellline,
                       i = sc.run(scp.cellline),
                  colvar = "celltype",
                  carrierPattern = "Carrier",
                  samplePattern = "HeLa|U2OS",
                  sampleFUN = "mean",
                  rowDataName = "MeanSCR")

tmp = rbindRowData(scp.cellline,sc.run(scp.cellline)) %>% as.data.frame() 
options(scipen = 200)
p1 = ggplot( tmp, aes( x = MeanSCR ) ) +
  geom_histogram() +
  scale_x_log10() +
  theme_cowplot(font_size = 26,font_family = "sans")  +
  # ggtitle( label = scp.stats(scp.cellline) ) +
  theme(text = element_text(size = 26),
        title = element_text(size = 13),
        axis.title = element_text(size = 26),
        axis.text.x = element_text(size = 32,hjust = 1),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        strip.background =element_blank(),strip.text = element_text(size = 20)) + 
        facet_wrap(~Raw.file, ncol = 5) 
p1 = p1 & geom_vline( xintercept = 0.1,size = 1,lty = 2) 
##### SCR comparison
scr = tmp[,c("assay","rowname","MeanSCR","Leading.razor.protein")]
colnames(scr)[1]="Raw.file"
head(scr)
scr = merge(scr,batch[1:21,],by = "Raw.file",all.x = T,all.y = F)
scr = scr[ scr$Raw.file!='20220707_SCP_HA_B_5' & scr$Raw.file !='20220712_SCP_HA_C_R' &scr$Raw.file !='20220708_SCP_LA_B_10',]
scr$Raw.file=gsub("\\d+_","",scr$Raw.file)
scr$Carrier.x = as.factor(scr$Carrier.x)
scr$Sonication = as.factor(scr$Sonication)
scr = scr %>% group_by(Raw.file) %>% mutate( median = median(MeanSCR,na.rm = T) )

mat1 = scr[scr$Sonication=="Normal frequency",]
mat1 = mat1 %>% arrange(Carrier.x,desc(median))
mat1$Raw.file = factor(mat1$Raw.file,levels = unique(mat1$Raw.file),ordered = T)

p_scr1 = ggplot( mat1, aes( x = Raw.file, y = MeanSCR, fill = Carrier.x ) ) +
  geom_boxplot(outlier.size = 1) + scale_fill_brewer( palette = "Set2" ) + 
  theme_cowplot(font_size = 26,font_family = "sans")   +
  scale_y_log10() +xlab(NULL) + ylab("MeanSCR") + 
  # ggtitle( label = scp.stats(scp.cellline) ) +
  theme(text = element_text(size = 26),
        axis.text.x = element_text( angle = 45,size = 18, hjust = 1, vjust = 1  ),
        title = element_text(size = 26),
        axis.ticks = element_line(size = 1.5),
        axis.title = element_text(size = 26),
        axis.line = element_line(size = 1.5)) 


mat2 = scr[scr$Sonication=="High frequency",]
mat2 = mat2 %>% arrange(Carrier.x,desc(median)) 
mat2$Raw.file = factor(mat2$Raw.file,levels = unique(mat2$Raw.file),ordered = T)
mat2$MeanSCR = log10(mat2$MeanSCR)

stats_p2 = na.omit(ungroup(mat2)) %>% 
          rstatix::anova_test(MeanSCR~Carrier.x) %>% 
          # rstatix::adjust_pvalue() %>% #rstatix::add_xy_position(dodge = 0.8) %>%  
          rstatix::add_significance(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns"))

p_scr2 = ggplot( mat2, aes( x = Raw.file, y =MeanSCR, fill = Carrier.x ) ) +
  geom_violin(aes(linetype=NA)) + scale_fill_brewer( palette = "Set2" ) + 
  geom_jitter(shape=21,aes(fill=Carrier.x),position = position_jitter(width = 0.1),
              color =alpha("#bebebe",0.6),size = 0.5)+
  theme_cowplot(font_size = 26,font_family = "sans") + ylab("log10(MeanSCR)") + xlab(NULL)+
  geom_text( aes(x = 5,y =0.4,label = "ANOVA test p < 10e-10"),family = "sans",size = 7) +
  # ggtitle( label = scp.stats(scp.cellline) ) +
  theme(text = element_text(size = 26),
        axis.text.x = element_text( angle = 45,size = 18, hjust = 1, vjust = 1  ),
        title = element_text(size = 26),
        axis.ticks = element_line(size = 1.5),
        axis.title = element_text(size = 26),
        axis.line = element_line(size = 1.5)) 

mat3 = scr[scr$Carrier.x=="200",]
mat3 = mat3 %>% arrange(Sonication,desc(median))
mat3$Raw.file = factor(mat3$Raw.file,levels = unique(mat3$Raw.file),ordered = T)
mat3$MeanSCR = log10(mat3$MeanSCR)
stats_p3 = ungroup(mat3) %>% 
          rstatix::wilcox_test(MeanSCR~Sonication) %>% 
          rstatix::adjust_pvalue() %>% rstatix::add_xy_position(x = "Sonication",dodge=0.8) %>%
          rstatix::add_significance(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns"))
          
p_scr3 = ggplot() +
  geom_violin( data = mat3 , aes( x = Sonication, y = MeanSCR,fill = Sonication,linetype=NA)) + 
  geom_jitter(data= mat3 , aes( x = Sonication, y = MeanSCR,fill = Sonication ) ,
              shape=21,position = position_jitter(width = 0.1),
              color = alpha("#bebebe",0.6),size = 0.3)+
  scale_fill_manual(values = alpha(brewer.pal("Accent",n=3),0.8)) +
  # geom_text( aes(x = 1.5,y =0.42,label = "wilcox test\n****"),family = "sans",size = 7) +
  stat_pvalue_manual(as.data.frame(stats_p3),label="p.adj.signif",size = 1,bracket.size=1,
                       tip.length=0.02,label.size = 10,family="sans") +
  theme_cowplot(font_size = 26,font_family = "sans")   +
  ylab("log10(MeanSCR)") + xlab(NULL)+ scale_y_continuous(limits = c(-3.5,0.8))+
  scale_x_discrete(label = c("High\nfrequency","Normal\nfrequency"))+
  # scale_y_log10() +
  # ggtitle( label = scp.stats(scp.cellline) ) +
  theme(text = element_text(size = 26),
        axis.text.x = element_text( size = 18  ),
        title = element_text(size = 26),
        axis.ticks = element_line(size = 1.5),
        axis.title = element_text(size = 26),
        axis.line = element_line(size = 1.5),
        legend.position = "none") 



scp.cellline = filterFeatures(scp.cellline, ~ !is.na(MeanSCR) & MeanSCR < 0.1 )
keepAssay =  names(dims(scp.cellline)[1, ][dims(scp.cellline)[1, ] > 0 ])
scp.cellline = scp.cellline[,,keepAssay]
scp.stats(scp.cellline);dims(scp.cellline)
tab4 = dims(scp.cellline) %>% as.data.frame() %>% t 
sum4 = data.frame( cells = 220,peptides = 5257,proteins = 1440 ) %>% t
colnames(tab4) = c("PSM4","Cell4")
colnames(sum4) = c("Stats4")
# [1] "all.runs = 20, sc.runs = 20, cells = 220, peptides = 5257, protein = 1440"
#      20220629_SCP_B1_10 20220629_SCP_B1_13 20220704_SCP_B_7 20220704_SCP_B_8 20220704_SCP_B_9 20220704_SCP_C_10 20220704_SCP_C_7 20220704_SCP_C_9 20220704_SCP_HA_B_11 20220704_SCP_HA_B_6
# [1,]                266               1348              946             1903             1991               662             1916             1611                  107                1144
# [2,]                 16                 16               16               16               16                16               16               16                   16                  16
#      20220704_SCP_HA_B_7 20220704_SCP_HA_S_12 20220704_SCP_LA_A1_11 20220705_SCP_HA_S_11 20220707_SCP_HA_B_12 20220708_SCP_LA_B_10 20220708_SCP_LA_B2_6 20220712_SCP_HA_A3_5 20220712_SCP_HA_A3_6
# [1,]                 302                 1365                  2709                  987                 1334                   15                 2271                 1018                 1635
# [2,]                  16                   16                    16                   16                   16                   16                   16                   16                   16
#      20220712_SCP_HA_A3_8
# [1,]                  308
# [2,]                   16

sc.runs = sc.run(scp.cellline)
n.psm = dims(scp.cellline)[1, ] ;summary(n.psm)       # psm counts
n.runs = dims(scp.cellline)[2, ] ;summary(n.runs)      # run tmt counts
keepAssay =  names(n.psm[n.psm > 200 ])
scp.cellline = scp.cellline[,,intersect(keepAssay,sc.runs)]
scp.stats(scp.cellline);dims(scp.cellline)
stats3 = dims(scp.cellline)[1,]
tab5 = dims(scp.cellline) %>% as.data.frame() %>% t 
sum5 = data.frame( cells = 198,peptides = 5252,proteins = 1439 ) %>% t
colnames(tab5) = c("PSM5","Cell5")
colnames(sum5) = c("Stats5")
# [1] "all.runs = 18, sc.runs = 18, cells = 198, peptides = 5252, protein = 1439"
#      20220629_SCP_B1_10 20220629_SCP_B1_13 20220704_SCP_B_7 20220704_SCP_B_8 20220704_SCP_B_9 20220704_SCP_C_10 20220704_SCP_C_7 20220704_SCP_C_9
# [1,]                266               1348              946             1903             1991               662             1916             1611
# [2,]                 16                 16               16               16               16                16               16               16
#      20220704_SCP_HA_B_6 20220704_SCP_HA_B_7 20220704_SCP_HA_S_12 20220704_SCP_LA_A1_11 20220705_SCP_HA_S_11 20220707_SCP_HA_B_12 20220708_SCP_LA_B2_6
# [1,]                1144                 302                 1365                  2709                  987                 1334                 2271
# [2,]                  16                  16                   16                    16                   16                   16                   16
#      20220712_SCP_HA_A3_5 20220712_SCP_HA_A3_6 20220712_SCP_HA_A3_8
# [1,]                 1018                 1635                  308
# [2,]                   16                   16                   16
stats = cbind(stats0,cbind(stats1,stats2)) %>% as.data.frame() %>%rownames_to_column("Raw.file")

stats3 = as.data.frame(stats3[grep("^\\d",names(stats3),value = T)])%>%rownames_to_column("Raw.file")
colnames(stats3)[2] = "stats3"
stats = merge(stats,stats3,all= T,by = "Raw.file")
stats = merge(stats,batch,by = "Raw.file") %>% filter(Raw.file!='20220707_SCP_HA_B_5' & Raw.file !='20220712_SCP_HA_C_R' &Raw.file !='20220708_SCP_LA_B_10')
stats$valid_ratio = stats$stats3/stats$stats0
stats$Raw.file=gsub("\\d+_","",stats$Raw.file)

stats$Carrier.x = as.factor(stats$Carrier.x)
stats$Sonication = as.factor(stats$Sonication)


mat0 = stats %>% arrange(Carrier.x,Sonication)
mat0$Raw.file = factor(mat0$Raw.file,levels = unique(mat0$Raw.file),ordered = T)

mat0 = mat0 %>% gather(value = Counts,key= "Filtration",-Raw.file,-Order,-Carrier.x,-Date,-Sonication,-Enzyme,-valid_ratio)
mat0$Filtration=ifelse(mat0$Filtration=="stats0","Raw",
                       ifelse(mat0$Filtration=="stats1","Filt1","Valid")) %>% factor(.,levels = c('Raw',"Filt1","Valid"),ordered=T)
                       
p_stats0 = ggplot( mat0 , aes( x = Raw.file, y = Counts ,fill = Filtration ) ) +
  geom_bar(position = "stack",stat = "identity") + 
  scale_fill_manual(values= alpha(brewer.pal("Accent",n=3),0.9)) + 
  theme_cowplot(font_size = 26,font_family = "sans") +
  ylab("Valid PSM ratio") + xlab(NULL)+ 
  # scale_y_log10() +
  # ggtitle( label = scp.stats(scp.cellline) ) +
  theme(text = element_text(size = 26),
        axis.text.x = element_text( angle = 45,size = 18, hjust = 1, vjust = 1  ),
        title = element_text(size = 26),
        axis.ticks = element_line(size = 1.5),
        axis.title = element_text(size = 26),
        axis.line = element_line(size = 1.5)) #+ facet_wrap(~Sonication+Carrier.x,drop = T)
p_stats0
mat1 = stats[stats$Sonication=="Norm
al frequency",]

p_stats1 = ggplot( mat1, aes( x = Carrier.x, y = valid_ratio , fill = Carrier.x ) ) +
  geom_boxplot(outlier.size = 2,size = 1) + 
  scale_fill_brewer( palette = "Set2" ) + 
  theme_cowplot(font_size = 26,font_family = "sans")   +
  scale_y_log10() +xlab(NULL) + ylab("Valid PSM ratio") +
   scale_x_discrete(label = c("<200\n(n=2)","200\n(n=8)")) +
  # ggtitle( label = scp.stats(scp.cellline) ) +
  theme(text = element_text(size = 26),
        axis.text.x = element_text(size = 22, hjust = 0.5, vjust = 0.5 ),
        title = element_text(size = 26),
        axis.ticks = element_line(size = 1.5),
        axis.title = element_text(size = 26),
        axis.line = element_line(size = 1.5),legend.position = "none") 
p_stats1

mat2 = stats[stats$Sonication=="High frequency",]
mat2 = mat2 %>% arrange(Carrier.x)
mat2$Raw.file = factor(mat2$Raw.file,levels = unique(mat2$Raw.file),ordered = T)
mat2$Raw.file = factor(mat2$Raw.file,levels = unique(mat2$Raw.file),ordered = T)

stats_p2 = na.omit(ungroup(mat2)) %>% 
          rstatix::anova_test(valid_ratio~Carrier.x) %>% 
        #   rstatix::adjust_pvalue() %>% 
          rstatix::add_significance(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns"))

p_stats2 = ggplot( mat2, aes( x = Raw.file, y = valid_ratio, fill = Carrier.x ) ) +
  geom_bar(position = "dodge",stat = "identity") + scale_fill_brewer( palette = "Set2" ) + 
  theme_cowplot(font_size = 26,font_family = "sans") + ylab("Valid PSM ratio") + xlab(NULL)+
  geom_text( aes(x = 5,y =0.3,label = "ANOVA test\nns"),family = "sans",size = 7) +
  scale_y_continuous(limits = c(0,0.35))+
  # ggtitle( label = scp.stats(scp.cellline) ) +
  theme(text = element_text(size = 26),
        axis.text.x = element_text( angle = 45,size = 18, hjust = 1, vjust = 1  ),
        title = element_text(size = 26),
        axis.ticks = element_line(size = 1.5),
        axis.title = element_text(size = 26),
        axis.line = element_line(size = 1.5)) 
p_stats2
mat3 = stats[stats$Carrier.x=="200",]
mat3 = mat3 %>% arrange(Sonication)
mat3$Raw.file = factor(mat3$Raw.file,levels = unique(mat3$Raw.file),ordered = T)

stats_p3 = ungroup(mat3) %>% 
          rstatix::t_test(valid_ratio ~Sonication) %>% 
          rstatix::adjust_pvalue() %>% rstatix::add_xy_position(x = "Sonication",dodge=0.8) %>%
          rstatix::add_significance(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns"))
        

p_stats3 = ggplot(  ) +
  geom_boxplot(data=mat3 , aes( x = Sonication, y = valid_ratio ,fill = Sonication ),size = 1) + 
  scale_fill_manual(values = alpha(brewer.pal("Accent",n=3),0.8)) +
  stat_pvalue_manual(as.data.frame(stats_p3),label="p.adj.signif",size = 1,bracket.size=1,
                       tip.length=0.02,label.size = 10,family="sans") +
  # geom_text( aes(x = 1.5,y =0.45),label = "wilcox test\np.adjust = 0.0299\n*",family = "sans",size = 7) +
  theme_cowplot(font_size = 26,font_family = "sans") + 
   scale_x_discrete(label = c("High\nfrequency\n(n=5)","Normal\nfrequency\n(n=8)")) +
  ylab("Valid PSM ratio") + xlab(NULL)+ 
  # scale_y_log10() +
  # ggtitle( label = scp.stats(scp.cellline) ) +
  theme(text = element_text(size = 26),
        axis.text.x = element_text( size = 20, hjust = 0.5, vjust = 1  ),
        title = element_text(size = 26),
        axis.ticks = element_line(size = 1.5),
        axis.title = element_text(size = 26),
        axis.line = element_line(size = 1.5),
        legend.position = "none") 
p_stats3

scp.cellline = joinAssays( scp.cellline,
                          i =  names(scp.cellline),
                          name = "psms_raw" )

scp.cellline = divideByReference( scp.cellline,
                         i = grep("^\\d+",names(scp.cellline),value = T), ### must have Reference or will be error
                         colvar = "celltype",
                         samplePattern = ".", 
                         refPattern = "Reference" ) 
scp.cellline = joinAssays( scp.cellline,
                          i =  grep("^\\d+",names(scp.cellline),value = T),
                          name = "psms" )
scp.cellline = aggregateFeaturesOverAssays(scp.cellline,
                                  i = grep("^\\d+",names(scp.cellline),value = T),
                                  fcol = "Modified.sequence",
                                  name = paste0("peptides_",grep("^\\d+",names(scp.cellline),value = T)),
                                  fun = matrixStats::colMedians, na.rm = TRUE)

peptide.runs = grep("^peptides_",names(scp.cellline),value = T)
scp.cellline = joinAssays( scp.cellline,
                  i = peptide.runs,
                  name = "peptides" )  ### may take long time
saveRDS(scp.cellline,"scp_peptides_aggregated.rds",compress = T)
scp.cellline = readRDS("scp_peptides_aggregated.rds")
# pep_to_protein = rbindRowData(scp.cellline, peptide.runs)[,c("Modified.sequence","Leading.razor.protein","Raw.file","razor_protein_fdr")] %>% 
#                  unique %>% as.data.frame() 
# # tmp = pep_to_protein[pep_to_protein$Modified.sequence %in% pep_to_protein$Modified.sequence[duplicated(pep_to_protein$Modified.sequence)],]

# unique(pep_to_protein$Modified.sequence) %>% length()
# rowdat = rbindRowData(scp.cellline, "peptides") %>% as.data.frame()
# unique(rowdat$Modified.sequence) %>% length()

# rowdat$protein = apply( rowdat,1,function(x){
#   id = pep_to_protein[pep_to_protein$Modified.sequence==x[6],]
#   if(length(unique(id$Leading.razor.protein))==1){
#     return(as.character(unique(id$Leading.razor.protein)))
#   }else if ( sum(id$razor_protein_fdr,na.rm = T)==0 ){
#      return(as.character(sample(id$Leading.razor.protein,size = 1)))
#   }else{
#     return(as.character(id$Leading.razor.protein[id$razor_protein_fdr==min(id$razor_protein_fdr,na.rm = T)]))
#   }
#   })#merge(rowdat,pep_to_protein,by = "Modified.sequence",all.x = T,all.y = F) %>% unique
# rowData(scp.cellline[["peptides"]]) = rowdat
psm = assay(scp.cellline[["psms_raw"]]) %>% as.data.frame() %>% 
      rownames_to_column( "psm" )
na_tbl = apply(psm[,-c(1)],2,function(x){
  na_num = length(x[is.na(x)])
  num = length(x[!is.na(x)])
  return(c(na_num,num))
}) %>% t %>% as.data.frame() %>% rownames_to_column("Channel")
colnames(na_tbl) = c("Channel","NA_number","Valid_number")


meta = metadata[,c("Raw.file","tmt.label","celltype","Carrier.x","Sonication")] %>% mutate(Channel = paste0(Raw.file,tmt.label)) %>%
        select( Channel,celltype,Carrier.x,Sonication ) 
meta = merge(meta,na_tbl,by = "Channel")
psm = psm %>% gather( key = "Channel",value = "RI", -psm )
psm = merge( psm, meta, by = c("Channel"),all.x = T )
psm = psm%>% arrange(Carrier.x,Sonication)
psm$batch = gsub("RI\\d+","",psm$Channel) %>% gsub("^\\d+_","",.) %>% factor(.)
psm$Channel = gsub("\\w+\\d+RI","RI",psm$Channel) %>% factor(.,levels = paste0("RI",1:16),order = T)
psm$TMTlabel = ifelse(psm$Channel=="RI1","126",
                      ifelse(psm$Channel=="RI2","127N",
                      ifelse(psm$Channel=="RI3","127C",
                      ifelse(psm$Channel=="RI4","128N",
                      ifelse(psm$Channel=="RI5","128C",
                      ifelse(psm$Channel=="RI6","129N",
                      ifelse(psm$Channel=="RI7","129C",
                      ifelse(psm$Channel=="RI8","130N",
                      ifelse(psm$Channel=="RI9","130C",
                      ifelse(psm$Channel=="RI10","131N",
                      ifelse(psm$Channel=="RI11","131C",
                      ifelse(psm$Channel=="RI12","132N",
                      ifelse(psm$Channel=="RI13","132C",
                      ifelse(psm$Channel=="RI14","133N",
                      ifelse(psm$Channel=="RI15","133C","134N"))))))))))))))) %>%
                      factor(.,levels =c("126","127N","127C","128N","128C","129N","129C","130N","130C","131N","131C","132N","132C","133N","133C","134N"),ordered = T)
psm$celltype = factor(psm$celltype,levels = c("Carrier","Reference","Unused","HeLa","U2OS","Blank"))

p2 = ggplot( psm, aes( x = TMTlabel, y = log10(RI), fill = celltype ) ) +
  geom_boxplot(outlier.size = 0.8) + scale_fill_brewer( palette = "Set2" ) + 
  theme_cowplot(font_size = 26,font_family = "sans")   +
  # ggtitle( label = scp.stats(scp.cellline) ) +
  theme(text = element_text(size = 26),
        axis.text.x = element_text( angle = 45,size = 14, hjust = 1, vjust = 1  ),
        title = element_text(size = 26),
        axis.ticks = element_line(size = 1.5),
        axis.title = element_text(size = 26),
        axis.line = element_line(size = 1.5),
        strip.background = element_blank(),strip.text = element_text(size = 18)) + 
        facet_wrap(~batch, ncol = 6)

p2.1 = ggplot( psm, aes( x = TMTlabel, y = log10(RI), fill = celltype ) ) +
  geom_boxplot(outlier.size = 1) + scale_fill_brewer( palette = "Set2" ) + 
  theme_cowplot(font_size = 26,font_family = "sans")   +
  # ggtitle( label = scp.stats(scp.cellline) ) +
  xlab(NULL)+ 
  theme(text = element_text(size = 26),
        axis.text.x = element_text( angle = 45,size = 30, hjust = 1, vjust = 1  ),
        title = element_text(size = 26),
        axis.ticks = element_line(size = 1.5),
        axis.title = element_text(size = 26),
        axis.line = element_line(size = 1.5),
        strip.background = element_blank(),strip.text = element_text(size = 18))

pep = assay(scp.cellline[["peptides"]]) %>% as.data.frame() %>% 
      rownames_to_column( "peptide" )
na_tbl = apply(pep[,-c(1)],2,function(x){
  na_num = length(x[is.na(x)])
  num = length(x[!is.na(x)])
  return(c(na_num,num))
}) %>% t %>% as.data.frame() %>% rownames_to_column("Channel")
colnames(na_tbl) = c("Channel","NA_number","Valid_number")


pep = pep %>% gather( key = "Channel",value = "RI", -peptide )
meta = metadata[,c("Raw.file","tmt.label","celltype","Carrier.x","Sonication")] %>% mutate(Channel = paste0(Raw.file,tmt.label)) %>%
        select( Channel,celltype,Carrier.x,Sonication ) 
meta = merge(meta,na_tbl,by = "Channel")
pep = merge( pep, meta, by = c("Channel"),all.x = T )
pep = pep %>% arrange(Carrier.x,Sonication)
pep$batch = gsub("RI\\d+","",pep$Channel) %>% gsub("^\\d+_","",.) %>% factor(.)
pep$Channel = gsub("\\w+\\d+RI","RI",pep$Channel) %>% 
              factor(.,levels = paste0("RI",1:16),order = T)
pep$TMTlabel = ifelse(pep$Channel=="RI1","126",
                      ifelse(pep$Channel=="RI2","127N",
                      ifelse(pep$Channel=="RI3","127C",
                      ifelse(pep$Channel=="RI4","128N",
                      ifelse(pep$Channel=="RI5","128C",
                      ifelse(pep$Channel=="RI6","129N",
                      ifelse(pep$Channel=="RI7","129C",
                      ifelse(pep$Channel=="RI8","130N",
                      ifelse(pep$Channel=="RI9","130C",
                      ifelse(pep$Channel=="RI10","131N",
                      ifelse(pep$Channel=="RI11","131C",
                      ifelse(pep$Channel=="RI12","132N",
                      ifelse(pep$Channel=="RI13","132C",
                      ifelse(pep$Channel=="RI14","133N",
                      ifelse(pep$Channel=="RI15","133C","134N"))))))))))))))) %>%
                      factor(.,levels =c("126","127N","127C","128N","128C","129N","129C","130N","130C","131N","131C","132N","132C","133N","133C","134N"),ordered = T)
                      
pep$celltype = factor(pep$celltype,levels = c("Carrier","Reference","Unused","HeLa","U2OS","Blank"))

p3 = ggplot( pep, aes( x = TMTlabel, y = log10(RI), fill = celltype ) ) +
  geom_boxplot(outlier.size = 1) + scale_fill_brewer( palette = "Set2" ) + 
  theme_cowplot(font_size = 26,font_family = "sans")   +
  xlab(NULL)+ 
  # ggtitle( label = scp.stats(scp.cellline) ) +
    theme(text = element_text(size = 26),
        axis.text.x = element_text( angle = 45,size = 14, hjust = 1, vjust = 1  ),
        title = element_text(size = 26),
        axis.ticks = element_line(size = 1.5),
        axis.title = element_text(size = 26),
        axis.line = element_line(size = 1.5),
        strip.background = element_blank(),strip.text = element_text(size = 18)) + facet_wrap(~batch, ncol = 6)

p3.1 = ggplot( pep, aes( x = TMTlabel, y = log10(RI), fill = celltype ) ) +
  geom_boxplot(outlier.size = 1) + scale_fill_brewer( palette = "Set2" ) + 
  theme_cowplot(font_size = 26,font_family = "sans")   +
  xlab(NULL)+
  # ggtitle( label = scp.stats(scp.cellline) ) +
  theme(text = element_text(size = 26),
        axis.text.x = element_text( angle = 45,size = 30, hjust = 1, vjust = 1  ),
        title = element_text(size = 26),
        axis.ticks = element_line(size = 1.5),
        axis.title = element_text(size = 26),
        axis.line = element_line(size = 1.5),
        strip.background = element_blank(),strip.text = element_text(size = 18))
p0 = (p2 |p3) + plot_layout(guides = 'collect')
p0.1 = (p2.1 /p3.1) + plot_layout(guides = 'collect')


meta = meta %>% gather( key = "Type",value = "Number",-celltype,-Channel,-Carrier.x,-Sonication )
meta$celltype = factor(meta$celltype,levels = c("Carrier","Reference","Unused","HeLa","U2OS","Blank"))
meta$batch = gsub("RI\\d+","",meta$Channel) %>% gsub("^\\d+_","",.) %>% factor(.,levels = levels(pep$batch))
meta$Channel = gsub("\\d+_\\w+_\\w+_\\d+RI","RI",meta$Channel) %>% 
              factor(.,levels = paste0("RI",1:16),order = T)

p4 = ggplot( meta,aes(x = celltype, y = Number/1000, fill = Type ) ) + 
  geom_bar(position = "dodge",stat = "identity") + 
  scale_fill_manual(values = c("#8f8e8e","#eb6059")) +
  theme_cowplot(font_size = 26,font_family = "sans")  + 
  xlab(NULL) + ylab(paste0("Pepetide counts (x1000)")) +
  # ggtitle( label = scp.stats(scp.cellline) ) +
  theme(text = element_text(size = 26),
        axis.text.x = element_text( angle = 45,size = 26, hjust = 1, vjust = 1  ),
        title = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.line = element_line(size = 1.5),
        strip.background = element_blank(),
        strip.text = element_text(size = 18)) + facet_wrap(facets = ~batch,ncol = 6)

pep_comp1 = meta %>% filter(Type =="Valid_number" & celltype %in% c("HeLa","U2OS") & Sonication=="Normal frequency" )

stats_p1 = pep_comp1 %>% 
          rstatix::t_test(Number~Carrier.x,paired = F) %>% 
          rstatix::adjust_pvalue() %>% rstatix::add_xy_position() %>% 
          rstatix::add_significance(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns"))

p_pep_num1 = ggplot( ) +
  geom_boxplot( data = pep_comp1,aes( x = Carrier.x, y = Number, fill = Carrier.x ) ,
                outlier.size = 1,size = 1) + 
  scale_fill_manual(values = ggsci::pal_npg(palette = c("nrc"), alpha = 0.8)(4)) + 
  theme_cowplot(font_size = 26,font_family = "sans")  +
  xlab(NULL)+ ylab("Peptides Number") + 
  # geom_text( aes(x = 3.5,y =2300,label = "ANOVA test\np = 0.00021\n***"),family = "sans",size = 7) +
    # ggtitle( label = scp.stats(scp.cellline) ) +
  stat_pvalue_manual(as.data.frame(stats_p1),label="p.adj.signif",size = 1,bracket.size=1,
                       tip.length=0.02,label.size = 10,family="sans") +
  theme(text = element_text(size = 26),
        axis.text.x = element_text( size = 22 ),
        title = element_text(size = 26),
        axis.ticks = element_line(size = 1.5),
        axis.title = element_text(size = 26),
        axis.line = element_line(size = 1.5)) 

pep_comp2 = meta %>% filter(Type =="Valid_number" & celltype %in% c("HeLa","U2OS") & Sonication=="High frequency" )

stats_p2 = pep_comp2 %>% 
          rstatix::t_test(Number~Carrier.x,paired = F) %>% 
          rstatix::adjust_pvalue() %>% rstatix::add_xy_position() %>% 
          rstatix::add_significance(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns"))

p_pep_num2 = ggplot() +
  geom_boxplot( data=pep_comp2, aes( x = Carrier.x, y = Number, fill = Carrier.x ) ,
              outlier.size = 1,size = 1) + 
  scale_fill_manual(values = ggsci::pal_npg(palette = c("nrc"), alpha = 0.8)(4)[2:4]) + 
  theme_cowplot(font_size = 26,font_family = "sans")  +
  xlab(NULL)+ ylab("Peptides Number") + 
  scale_y_continuous(limits = c(0,2000)) +
  # geom_text( aes(x = 3.5,y =2300,label = "ANOVA test\np = 0.00021\n***"),family = "sans",size = 7) +
    # ggtitle( label = scp.stats(scp.cellline) ) +
  stat_pvalue_manual(as.data.frame(stats_p2),label="p.adj.signif",
                        size = 1,bracket.size=1,step.increase = 0.05,
                       tip.length=0.02,label.size = 10,family="sans") +
  theme(text = element_text(size = 26),
        axis.text.x = element_text( size = 22 ),
        title = element_text(size = 26),
        axis.ticks = element_line(size = 1.5),
        axis.title = element_text(size = 26),
        axis.line = element_line(size = 1.5)) 

pep_comp3 = meta %>% filter(Type =="Valid_number" & celltype %in% c("HeLa","U2OS") & Carrier.x=="200" )

stats_p3 = pep_comp3 %>% 
          rstatix::t_test(Number~Sonication,paired = F) %>% 
          rstatix::adjust_pvalue() %>% rstatix::add_xy_position() %>% 
          rstatix::add_significance(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns"))

p_pep_num3 = ggplot() +
  geom_boxplot( data=pep_comp3, aes( x = Sonication, y = Number, fill = Sonication ) ,
                outlier.size = 1,size = 1) + 
  scale_fill_manual(values = brewer.pal("Accent",n=3)) + 
  theme_cowplot(font_size = 26,font_family = "sans")  +
  xlab(NULL)+ ylab("Peptides Number") +  
  scale_x_discrete(label = c("High\nfrequency","Normal\nfrequency")) + 
  stat_pvalue_manual(as.data.frame(stats_p3),label="p.adj.signif",size = 1,bracket.size=1,
                       tip.length=0.02,label.size = 10,family="sans") +
  theme(text = element_text(size = 26),
        axis.text.x = element_text( size = 22 ),
        title = element_text(size = 26),
        axis.ticks = element_line(size = 1.5),
        axis.title = element_text(size = 26),
        axis.line = element_line(size = 1.5),
        legend.position = "none") 

p_pep_num = (p_pep_num1 + p_pep_num2 + p_pep_num3 + plot_layout(guides  ='collect'))


plot(scp.cellline)
scp.stats(scp.cellline) 
# [1] "all.runs = 18, sc.runs = 18, cells = 198, peptides = 5252, protein = 1439"

medianRI = colMeans( assay(scp.cellline[["peptides"]]),na.rm = T )
scp.cellline$medianRI = medianRI

k = as.data.frame(colData( scp.cellline )) 
k$celltype = factor(k$celltype,levels = c("Carrier","Reference","Unused","HeLa","U2OS","Blank"))

p5 = ggplot(k) +
  geom_boxplot(aes( x = medianRI, y = celltype, fill =  celltype ),size = 1.5) + 
  scale_fill_brewer( palette = "Set2" ) + 
  geom_vline( xintercept = quantile(medianRI,0.95,na.rm = T),size = 2) + scale_x_log10( ) + 
  theme_cowplot() + ggtitle("MedianRI") +
  ggtitle(NULL) + xlab("log10(medianRI)") + ylab(NULL) + 
  theme(text = element_text(size = 26),
        title = element_text(size = 16),
        axis.text = element_text(size = 30),
        axis.title = element_text(size = 32),
        axis.line = element_line(size = 1.5),
        axis.ticks = element_line(size = 1.5),
        legend.title = element_blank()) 
p5

scp.tmp = medianCVperCell( scp.cellline, 
                        i = sc.run(scp.cellline),
                        groupBy = "Leading.razor.protein",
                        nobs = 5,  # least peptide numbers per protein,2 and 5
                        norm = "SCoPE2",
                        na.rm = T, 
                        colDataName = "MedianCV" ) #and MedianCV_nob5,
###Warning message:
# In medianCVperCell(scp.cellline, i = sc.runs, groupBy = "Leading.razor.protein",  :
#   The median CV could not be computed for one or more samples. You may want to try a smaller value for 'nobs'.
k2 = getWithColData(scp.tmp, "peptides") %>%
  colData %>%
  data.frame
k2$celltype = factor(k$celltype,levels = c("Carrier","Reference","Unused","HeLa","U2OS","Blank"))

p6 = k2 %>%
    ggplot(aes(x = MedianCV,
               fill = celltype)) +
    geom_density(size = 1) + scale_fill_manual(values = alpha(brewer.pal(name = "Set2", n = 6),alpha = 0.6)) +
    geom_vline(xintercept = 0.3,size = 1.5,linetype = 2) +
  theme_cowplot() + ggtitle("MedianCV") + xlab("MedianCV") +
  theme(text = element_text(size = 26),
        title = element_text(size = 16),
        axis.text = element_text(size = 32),
        axis.title = element_text(size = 32),
        axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 2)) 
p6  

cv_comp1 = k2 %>% filter( celltype %in% c("HeLa","U2OS") & Sonication=="Normal frequency" )

stats_p1 = cv_comp1 %>%
          rstatix::t_test(MedianCV~Carrier.x,paired = F) %>% 
          rstatix::adjust_pvalue() %>% rstatix::add_xy_position() %>% 
          rstatix::add_significance(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns"))

p_cv_num1 = ggplot( ) +
  geom_boxplot( data = cv_comp1,aes( x = Carrier.x, y = MedianCV, fill = celltype ) ,
                outlier.size = 1,size = 1) + 
  scale_fill_manual(values = alpha(c("#b08ed6","#81c2a2"),0.8)) + 
  theme_cowplot(font_size = 26,font_family = "sans")  +
  xlab(NULL)+ ylab(" MedianCV") +  
  # geom_text( aes(x = 3.5,y =2300,label = "ANOVA test\np = 0.00021\n***"),family = "sans",size = 7) +
    # ggtitle( label = scp.stats(scp.cellline) ) +
  stat_pvalue_manual(as.data.frame(stats_p1),label="p.adj.signif",size = 1,bracket.size=1,
                       tip.length=0.02,label.size = 10,family="sans") +
  theme(text = element_text(size = 26),
        axis.text.x = element_text( size = 22 ),
        title = element_text(size = 26),
        axis.ticks = element_line(size = 1.5),
        axis.title = element_text(size = 26),
        axis.line = element_line(size = 1.5)) 
p_cv_num1
cv_comp2 =  k2 %>% filter( celltype %in% c("HeLa","U2OS") &Sonication=="High frequency" )

stats_p2 = cv_comp2 %>% 
          rstatix::t_test(MedianCV~Carrier.x,paired = F) %>% 
          rstatix::adjust_pvalue() %>% rstatix::add_xy_position() %>% 
          rstatix::add_significance(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns"))

p_cv_num2 = ggplot() +
  geom_boxplot( data=cv_comp2, aes( x = Carrier.x, y = MedianCV, fill = celltype ) ,
              outlier.size = 1,size = 1) + 
  scale_fill_manual(values = alpha(c("#b08ed6","#81c2a2"),0.8)) + 
  theme_cowplot(font_size = 26,font_family = "sans")  +
  xlab(NULL)+ ylab(" MedianCV") +  
  scale_y_continuous(limits = c(0,0.75)) +
  # geom_text( aes(x = 3.5,y =2300,label = "ANOVA test\np = 0.00021\n***"),family = "sans",size = 7) +
    # ggtitle( label = scp.stats(scp.cellline) ) +
  stat_pvalue_manual(as.data.frame(stats_p2),label="p.adj.signif",
                        size = 1,bracket.size=1,step.increase = 0.05,
                       tip.length=0.02,label.size = 10,family="sans") +
  theme(text = element_text(size = 26),
        axis.text.x = element_text( size = 22 ),
        title = element_text(size = 26),
        axis.ticks = element_line(size = 1.5),
        axis.title = element_text(size = 26),
        axis.line = element_line(size = 1.5)) 
p_cv_num2
cv_comp3 = k2 %>% filter( celltype %in% c("HeLa","U2OS") & Carrier.x=="200" )

stats_p3 = cv_comp3 %>% group_by(celltype) %>%
          rstatix::t_test(MedianCV~Sonication,paired = F) %>% 
          rstatix::adjust_pvalue() %>% rstatix::add_xy_position() %>% 
          rstatix::add_significance(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns"))

p_cv_num3 = ggplot() +
  geom_boxplot( data=cv_comp3, aes( x = Sonication, y = MedianCV, fill = celltype),size = 1.2)+
  scale_fill_manual(values = alpha(c("#b08ed6","#81c2a2"),0.8)) + 
  theme_cowplot(font_size = 26,font_family = "sans")  +
  xlab(NULL)+ ylab(" MedianCV") +  
  scale_x_discrete(label = c("High\nfrequency","Normal\nfrequency")) + 
  stat_pvalue_manual(as.data.frame(stats_p3),label="p.adj.signif",
                      size = 1,bracket.size=1,step.increase = 0.4,
                       tip.length=0.02,label.size = 10,family="sans") +
  theme(text = element_text(size = 26),
        axis.text.x = element_text( size = 22 ),
        title = element_text(size = 26),
        axis.ticks = element_line(size = 1.5),
        axis.title = element_text(size = 26),
        axis.line = element_line(size = 1.5),legend.title = element_blank()) 
p_cv_num3

p_cv =  (p_cv_num1 + p_cv_num2 + p_cv_num3 + plot_layout(guides ='collect'))
pdf("carrier_plot.pdf",width = 20,height = 38)
k='
#AAAAA##
BBBBCCC#
#EEEEE##
FFGGGHHH
IIIIIIII
JJJJJJJJ
'
p_scr1 +p_scr2 + p_scr3 + p_stats0 + p_stats1 +p_stats2 + p_stats3 + p_pep_num + p_cv + plot_layout(design = k)
dev.off()


scp.cellline = scp.tmp[,!is.na(scp.tmp$MedianCV) & scp.tmp$MedianCV < 0.3 & scp.tmp$celltype %in% c("HeLa", "U2OS"), ]

scp.stats(scp.cellline);dims(scp.cellline)
tab6 = dims(scp.cellline[,,grep("peptides_",names(scp.cellline),value = T)]) %>% as.data.frame() %>% t 
sum6 = data.frame( cells = 153,peptides = 5252,proteins = 1439 ) %>% t
colnames(tab6) = c("Pep6","Cell6")
rownames(tab6) = gsub("peptides_","",rownames(tab6))
colnames(sum5) = c("Stats6")
# "all.runs = 18, sc.runs = 18, cells = 153, peptides = 5252, protein = 1439"
#      20220629_SCP_B1_10 20220629_SCP_B1_13 20220704_SCP_B_7 20220704_SCP_B_8 20220704_SCP_B_9 20220704_SCP_C_10 20220704_SCP_C_7 20220704_SCP_C_9 20220704_SCP_HA_B_6 20220704_SCP_HA_B_7
# [1,]                266               1348              946             1903             1991               662             1916             1611                1144                 302
# [2,]                  6                 10                6               10                7                 5                8                8                  11                   9
#      20220704_SCP_HA_S_12 20220704_SCP_LA_A1_11 20220705_SCP_HA_S_11 20220707_SCP_HA_B_12 20220708_SCP_LA_B2_6 20220712_SCP_HA_A3_5 20220712_SCP_HA_A3_6 20220712_SCP_HA_A3_8 psms_raw  psms
# [1,]                 1365                  2709                  987                 1334                 2271                 1018                 1635                  308    23716 23716
# [2,]                    9                     6                   11                   10                   11                    8                   11                    7      153   153
#      peptides_20220629_SCP_B1_10 peptides_20220629_SCP_B1_13 peptides_20220704_SCP_B_7 peptides_20220704_SCP_B_8 peptides_20220704_SCP_B_9 peptides_20220704_SCP_C_10 peptides_20220704_SCP_C_7
# [1,]                         260                        1316                       934                      1850                      1924                        646                      1829
# [2,]                           6                          10                         6                        10                         7                          5                         8
#      peptides_20220704_SCP_C_9 peptides_20220704_SCP_HA_B_6 peptides_20220704_SCP_HA_B_7 peptides_20220704_SCP_HA_S_12 peptides_20220704_SCP_LA_A1_11 peptides_20220705_SCP_HA_S_11
# [1,]                      1557                         1095                          293                          1220                           2627                           961
# [2,]                         8                           11                            9                             9                              6                            11
#      peptides_20220707_SCP_HA_B_12 peptides_20220708_SCP_LA_B2_6 peptides_20220712_SCP_HA_A3_5 peptides_20220712_SCP_HA_A3_6 peptides_20220712_SCP_HA_A3_8 peptides
# [1,]                          1292                          2175                           990                          1570                           301     5252
# [2,]                            10                            11                             8                            11                             7      153
plot(scp.cellline)

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
tab7 = tab6
sum7 = data.frame( cells = 153,peptides = 5243,proteins = 1439 ) %>% t
colnames(tab7) = c("Pep7","Cell7")
colnames(sum7) = c("Stats7")
scp.stats(scp.cellline)
# [1] "all.runs = 18, sc.runs = 18, cells = 153, peptides = 5252, protein = 1439"
scp.cellline
# An instance of class QFeatures containing 41 assays:
#  [1] 20220629_SCP_B1_10: SingleCellExperiment with 266 rows and 6 columns 
#  [2] 20220629_SCP_B1_13: SingleCellExperiment with 1348 rows and 10 columns 
#  [3] 20220704_SCP_B_7: SingleCellExperiment with 946 rows and 6 columns 
#  ...
#  [39] peptides: SingleCellExperiment with 5252 rows and 153 columns 
#  [40] peptides_norm_col: SingleCellExperiment with 5252 rows and 153 columns 
#  [41] peptides_norm: SingleCellExperiment with 5243 rows and 153 columns

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
  
plot(scp.cellline)

protein_norm = assay(scp.cellline[["proteins_norm"]]) %>% as.data.frame() %>% 
      rownames_to_column( "proteins" )
na_tbl = apply(protein_norm[,-c(1)],2,function(x){
  na_num = length(x[is.na(x)])
  num = length(x[!is.na(x)])
  return(c(na_num,num))
}) %>% t %>% as.data.frame() %>% rownames_to_column("Channel")
colnames(na_tbl) = c("Channel","NA_number","Valid_number")

meta = metadata[,c("Raw.file","tmt.label","celltype","Carrier.x","Sonication")] %>% 
        mutate(Channel = paste0(Raw.file,tmt.label)) %>%
        select( Channel,celltype,Sonication,Carrier.x ) 
meta = merge(meta,na_tbl,by = "Channel")
meta$Sonication = as.factor(meta$Sonication)
meta$Carrier.x = as.factor(meta$Carrier.x)
meta$batch = gsub("RI\\d+","",meta$Channel) %>% gsub('^\\d+_',"",.) 

meta  = meta %>% arrange(Carrier.x,Sonication)
meta$batch = factor(meta$batch ,levels = unique(meta$batch),ordered = T)
meta$Channel = gsub("\\d+_\\w+_\\w+_\\d+RI","RI",meta$Channel) %>% factor(., levels = paste0("RI",1:16),order = T)
meta = meta %>% gather( key = "Type",value = "Number",-batch,-Channel,-celltype,-Sonication,-Carrier.x )

p7 = ggplot( meta,aes(x = Channel, y = Number/1000, fill = Type ) ) + 
  geom_bar(position = "dodge",stat = "identity") +
  theme_cowplot(font_size = 26,font_family = "sans")  +
  # ggtitle( label = scp.stats(scp.cellline) ) +
  ylab( "Protein Nubmer(x1000)" )+ xlab(NULL)+
  scale_fill_manual(values = c("#8f8e8e","#eb6059")) +
  theme(text = element_text(size = 26),
        axis.text.x = element_text( angle = 45,size = 26, hjust = 1, vjust = 1  ),
        title = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 2),
        strip.background = element_blank(),
        strip.text = element_text(size = 18)
        ) + facet_wrap(~batch, ncol = 6)
p7 

scp.cellline[["proteins_norm"]] %>%
    assay %>% 
    is.na %>%  mean() ##[1] 0.6579384

library(impute)
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
save(scp.cellline,file = "imputed_scp.RData",compress = T)
load("imputed_scp.RData")
scp.stats(scp.cellline)
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
pca = p$data

mat = assay(scp.tmp)
pearson = cor(t(mat))
head(pearson[1:10,1:10])
rsum = rowSums(pearson^2)

# Calculate the weighted data matrix:
mat = diag(rsum) %*%  mat
pca.imp.cor <- cor(mat)

# pca = prcomp(t(mat))
pca = prcomp(pca.imp.cor)


pca = prcomp( t(assay(scp.tmp) %>% t %>% scale(.,center = T,scale = T) %>% t ))
deviation = (pca$sdev)^2/sum((pca$sdev)^2)*100
sum(deviation[1:2]) ;length(deviation)
# [1] 15.52804
# [1] 69

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
        line = element_line(size = 2)) 
p8
pdf("QCplot.pdf",width = 60,height = 46)
k = '
BBBBBBBBBB
AAAAAAACCC
DDDDDGGGGG
###EEFFF##
'
p1+p0+ p0.1+p4+p5+p6+p7 + plot_layout(design = k,height = c(12,14,12,8))
dev.off()
# batchL = assay(scp.cellline[["proteins_batchL"]])

# pears_cor = cor(batchL)
tab = list(tab1,tab2,tab3,tab4,tab5,tab6,tab7)
tab = lapply(tab,function(x){
  x = as.data.frame(x) %>% rownames_to_column(.,"Batch")
  return(x)
})
a = tab[[1]]
for(i in 2:7){
  tmp = tab[[i]]
  a = merge(a,tmp,by = "Batch",all.x = T)
}


sum = list(sum1,sum2,sum3,sum4,sum5,sum6,sum7)
sum = lapply(sum,function(x){
  x = as.data.frame(x) %>% rownames_to_column(.,"Stats")
  return(x)
})
b = sum[[1]]
for(i in 2:7){
  tmp = sum[[i]]
  b= merge(b,tmp,by = "Stats",all.x = T)
}

write.csv( a, file = "batch_state1.csv",quote = F,row.names = T,col.names = T)
write.csv( b, file = "batch_state2.csv",quote = F,row.names = T,col.names = T)
