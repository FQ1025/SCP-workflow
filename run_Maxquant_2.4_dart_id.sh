#!/bin/bash
#conda activate r4
filePath1=/lustre/user/liclab/fangq/proj/scp/cellLine/data_2023/raw_sc
filePath2=/lustre/user/liclab/fangq/proj/scp/scope_protocol_data/MSV000087041/raw/sc 
#filePath3=/lustre/user/liclab/fangq/proj/scp/cellLine/data_2023/raw_0329/
#filePath4=/lustre/user/liclab/fangq/proj/scp/cellLine/data_2023/raw_0330/
filePath=$filePath1','$filePath2
outputPath=/lustre/user/liclab/fangq/proj/scp/cellLine/data_2023/sc_maxquant2.4_test/maxquant/
para_xmlPath=/lustre/user/liclab/fangq/proj/scp/cellLine/data_2023/sc_maxquant2.4_test
#echo $filePath

Rscript /home/fangq/ffqq/tools/MaxQuant_2.4.0.0/Generate_xml.R -i ${filePath} \
	-o ${outputPath} \
	-t 32 \
	-x ${para_xmlPath}


dotnet /lustre/user/liclab/fangq/tools/MaxQuant_2.4.0.0/bin/MaxQuantCmd.exe ${para_xmlPath}/mqpar.xml 


dart_id -c config_cellline.yml 
