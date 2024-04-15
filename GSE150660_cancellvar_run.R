####需要的文件：D:/GRCh38_ENSG_Symbol.txt,D:/sample_type.txt
####需要更改Rdata aa_GSE123.Rdata->GSE123.Rdata
####需要更改以下两行
GSE_id <- "GSE150660"
input <- "H:/single_cell/LECA/"
#############################################################################################################################################################
#1.all_var_mat.R##################################################################################################################################################
#############################################################################################################################################################
dir_files<-list.dirs(path = paste(input,GSE_id,sep=""),recursive=F,full.names=F) 
for(dir_file in dir_files){
  #1-更换路径
  setwd(paste0(input,GSE_id,"/",dir_file))
  
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(tidyverse)
  
  #更换突变X46787,细胞X2356
  ##########################1-读取scsnv突变矩阵###############################
  barcodes <- read.table("final/barcodes.txt",header=F,sep="\t",as.is=T)  #2356
  snv_info <- read.table("final/snvs.csv",header=T,sep=",",as.is=T) #46787*35
  snv_info$pos <- snv_info$pos+1 #scsnv位置改为从1起始
  snv_info[,36] <- paste0(snv_info[,1],":",snv_info[,2],",",snv_info$strand,",",snv_info[,3],">",snv_info[,4])
  #第一列是snv，第二列细胞，第三列count
  ref_mat <- read.table("final/refs.mtx",sep=" ",header=T,comment.char = "%",as.is=T) 
  alt_mat <- read.table("final/alts.mtx",sep=" ",header=T,comment.char = "%",as.is=T) 
  colnames(ref_mat)<-c("snv_index","cell_index","read_count")
  colnames(alt_mat)<-c("snv_index","cell_index","read_count")
  
  #######################2-处理snvinfo用于annovar注释#####################
  options(scipen=999) #防止写出的时候数字是科学计数法
  snv_info$start <- snv_info$pos
  snv_info$end <- snv_info$pos
  #snv_info_new <- snv_info[,c(1,37,38,3,4)]
  #write.table(snv_info_new,"C:/Users/ttt/Desktop/data/5-MCC-gse117988/37/5_annovr/snv_info_5col.txt",sep="\t",row.names = F,col.names=F,quote=F)
  
  ########################3-加上annovar注释的结果#######################
  #setwd("C:/Users/ttt/Desktop/data/5-MCC-gse117988/37/5_annovr")
  annovar_anno <- read.csv("annovar/myanno.hg38_multianno.csv",header=T,as.is=T) #46787*32
  
  annovar_anno2 <- annovar_anno
  annovar_anno2$gnomad_af001 <- NA
  annovar_anno2$gnomad_af001[which(as.numeric(annovar_anno2$AF) >= 0.01 )] <- "True" #44990 
  
  annovar_anno2$ang1000_af001 <- NA
  annovar_anno2$ang1000_af001[which(as.numeric(annovar_anno2$ALL.sites.2015_08) >= 0.01 )] <- "True" #44000
  
  annovar_anno2 <- annovar_anno2[,-c(12:23)]
  names(annovar_anno2)[1:5] <- c("chrom", "start","end", "ref", "alt") 
  snv_info_ano <- full_join(snv_info,annovar_anno2) #55列
  write.table(snv_info_ano[snv_info_ano$g1000 %in% "False" & snv_info_ano$avsnp150 %in% ".",],paste0(dir_file,"_snv_info_ano.txt"),sep = "\t",quote = F,col.names = T,row.names=F)
  
  ###################统计突变注释情况
  anno_stat <- matrix(NA,1,13)
  colnames(anno_stat) <- c("s_edit","s_g1000_af001","a_g1000","a_g1000_af001","gnomod_genome","gnomod_af001","clivar","dbsnp","cosmic_coding","cosmic_noncoding","cosmic_isec","retain","retain_5cell") 
  #clivar不作为过滤，是看最终结果里的情况
  
  anno_stat[1,1:10] <- c(length(which(snv_info_ano$REDIportal=="True")),
                         length(which(snv_info_ano$g1000=="True")),
                         length(which(snv_info_ano$ALL.sites.2015_08!=".")),
                         length(which(snv_info_ano$ang1000_af001=="True")),
                         length(which(snv_info_ano$AF!=".")), #gnomad
                         length(which(snv_info_ano$gnomad_af001=="True")),
                         length(which(snv_info_ano$CLNDN!=".")),
                         length(which(snv_info_ano$avsnp150!=".")),
                         length(which(snv_info_ano$cosmic96_coding!=".")),
                         length(which(snv_info_ano$cosmic96_noncoding!=".")))
  
  ###############4-处理基因型######################
  #REF/ALT（相当于有突变或无突变）
  tmp_mat <- cbind(ref_mat[,3],alt_mat[,3])
  colnames(tmp_mat) <- c("ref_count","alt_count")
  tmp_mat <- as.data.frame(tmp_mat)
  tmp_mat$type <- 1 #ALT,有alt_count
  tmp_mat$type[which(tmp_mat$alt_count==0)] <- 0 #REF，没有alt_count , -1 no call
  tmp_mat$cell_AF<-ifelse(tmp_mat$alt_count!=0,tmp_mat$alt_count/(tmp_mat$ref_count+tmp_mat$alt_count),0)
  
  all_snp_mat <- ref_mat[,1:2]
  all_snp_mat <- cbind(all_snp_mat,tmp_mat)
  
  #把index还原为名字，因为后面会打乱
  all_snp_mat$mut_name <- snv_info_ano[all_snp_mat$snv_index,36]
  all_snp_mat$CB <- barcodes[all_snp_mat$cell_index,1]
  #all_snp_mat <- left_join(all_snp_mat,final_cell_type)
  less_snv_info_ano <- snv_info_ano[,c(1:4,29,35:45,46,50:53)]
  adc<-function(x){
    bridge_str<-paste0(x[20],x[21])
    cosmic_str<-strsplit(bridge_str,split="=|;")[[1]][2]
    return(cosmic_str)
  }
  less_snv_info_ano$cosmicID<-apply(less_snv_info_ano,1,function(x){adc(x)})
  names(less_snv_info_ano)[7] <- "mut_name"
  all_snp_mat <- left_join(all_snp_mat,less_snv_info_ano)
  all_snp_mat <- all_snp_mat[all_snp_mat$g1000 %in% "False" & all_snp_mat$avsnp150 %in% ".",]
  df2 <-all_snp_mat %>% as_tibble() %>% 
    separate_rows(Gene.ensGene, sep = ";")
  ID_symbol <- read.table("D:/GRCh38_ENSG_Symbol.txt",header=F,sep="\t",as.is=T) 
  bridge<-ID_symbol[,c(2,2)]
  colnames(bridge)<-colnames(ID_symbol)<-c("Gene.ensGene","gene_symbol")
  ID_symbol<-rbind(ID_symbol,bridge)
  df2<-merge(df2,ID_symbol)
  df2[is.na(df2)]<-"."
  write.table(df2,paste0(dir_file,"_all_var_mat.txt"),sep = "\t",quote = F,col.names = T,row.names=F)
}



#############################################################################################################################################################
#2.tsne_1.R##################################################################################################################################################
#############################################################################################################################################################
rm(list = ls()[-match(c("GSE_id","input"),ls())])
##使用的细胞为表达谱中的全部细胞colnames(tissue_exp)，全部细胞在突变|编辑谱的缺失值用NA补全，在突变|编辑谱中多余的细胞被删除
##使用的基因为突变|编辑谱和表达谱共有的基因exp_mut_G,仅有表达没有突变的基因或仅有突变没有表达的基因被删除

#########4.1 tsne
#Vento-Tormo Group
#五个表（横轴为基因，纵轴为组织，每一个格子是下面的内容）
#细胞类型
#细胞名
#表达值
#tsne坐标X
#tsne坐标y
get.paste <- function(x) {
  xx <- paste(x,collapse=",")
  return(xx)
}
get.paste2 <- function(x) {
  xx <- paste(x,collapse=";")
  return(xx)
}
library(Seurat)
library(dplyr)
setwd(paste(input,GSE_id,sep=""))
load(paste(GSE_id,".RData",sep=""))

tissue_seurat <- test.seu
tissue_exp <- tissue_seurat[["RNA"]]@data


samples <- list.dirs(path = ".",recursive=F,full.names=F)
pro<-samples[1]
file=read.table(paste0(pro,"/",pro,"_all_var_mat.txt"),header=T,sep="\t",as.is=T)
file$CB<-paste0(pro,"_",file$CB,"-1")
if(length(samples)>1){
  for(pro in samples[2:length(samples)]){
    file1=read.table(paste0(pro,"/",pro,"_all_var_mat.txt"),header=T,sep="\t",as.is=T)
    file1$CB<-paste0(pro,"_",file1$CB,"-1")
    file<-rbind(file,file1)
  }
}
library(reshape2)
file<-file[file$CB %in% colnames(tissue_exp),]##在突变|编辑谱中多余的细胞被删除
mut_freq<-as.data.frame(table(unique(file[,c("mut_name","CB")])$mut_name))##突变|编辑可在多个基因上，因此单独提突变|编辑会有重复,需对突变|编辑-细胞组合去重
file_50<-file[file$mut_name %in% mut_freq[which(mut_freq[,2]>50),1],]##保留大于50个细胞检测的突变;大概10000突变；涉及3000基因
exp_mut_G<-intersect(rownames(tissue_exp),file_50$gene_symbol)
file_50<-file_50[file_50$gene_symbol %in% exp_mut_G,]
file_50$cell_AF <- round(file_50$cell_AF,4)
#write.table(file_50,file=paste(GSE_id,"_50_var_mat.txt",sep=""),sep="\t",quote=F,row.names=F)

tissue_var<-dcast(data=unique(file_50[,c("mut_name","CB","cell_AF")]),mut_name ~CB, value.var ="cell_AF",fill="NONE")
rownames(tissue_var)<-tissue_var$mut_name
tissue_var<-tissue_var[,-1]
tissue_var[,setdiff(colnames(tissue_exp),colnames(tissue_var))]<-"NONE"#细胞突变|编辑谱var_cell_matrix
tissue_var_out<-as.data.frame(tissue_var)
tissue_var_out$gene_symbol<-rownames(tissue_var_out)
write.table(tissue_var_out,file=paste(GSE_id,"_filterCB_var_mat.txt",sep=""),sep="\t",quote=F,row.names=F)

tissue_exp <- tissue_exp[rownames(tissue_exp) %in% exp_mut_G,]#细胞表达谱G_cell_matrix
tissue_exp_out<-as.data.frame(tissue_exp)
tissue_exp_out$gene_symbol<-rownames(tissue_exp_out)
write.table(tissue_exp_out,file=paste(GSE_id,"_filterG_exp_mat.txt",sep=""),sep="\t",quote=F,row.names=F)

tsne <- as.data.frame(tissue_seurat@reductions$tsne@cell.embeddings)
tsne$cell <- rownames(tsne)
anno_value <- tissue_seurat@meta.data
inferCNV_celltype <- read.table("cell_type.inferCNV.txt",header=T,sep="\t",stringsAsFactors=F)
inferCNV_celltype_cancer <- inferCNV_celltype[which(inferCNV_celltype[,3] == "cancer"),]
anno_value[rownames(inferCNV_celltype_cancer),"celltype"] <- "malignant"
#anno_value$cli_type <- ifelse(anno_value$orig.ident == "SRR7722938","After T cell therapy","Before T cell therapy")
anno_value$cli_type <- anno_value$orig.ident

rm(list=c("tissue_seurat"))
anno_value$cell <- rownames(anno_value)
sample_type<-read.table("D:/sample_type.txt",header=F,sep="\t",as.is=T)
#anno_value$cli_type<-ifelse(anno_value$cli_type == "GSM4957683" | anno_value$cli_type =="GSM4957684", "primary", "metastasis") 
anno_value <- anno_value[,c("cell","celltype","cli_type")]
sample_type<-sample_type[sample_type[,1] %in% GSE_id,]
colnames(sample_type)[3]<-"cli_type"
anno_value<-merge(anno_value,sample_type)
anno_value<- anno_value[,c("cell","celltype","V4")]
colnames(anno_value)[3]<-"cli_type"
write.table(anno_value,file=paste(GSE_id,"_cell_type.txt",sep=""),sep = "\t",quote = F,col.names = T,row.names=F)

file_50_cell<-file_50
colnames(file_50_cell)[9]<-"cell"
file_50_cell<-merge(file_50_cell,anno_value)
write.table(file_50_cell,file=paste(GSE_id,"_cell_50_var_mat.txt",sep=""),sep="\t",quote=F,row.names=F)

tsne_result <- merge(tsne,anno_value,by="cell",all.x=T)
tsne_result <- tsne_result[order(tsne_result$celltype),]#细胞坐标和注释cell_x_y_type_Sample

tsne_result_1<-tsne_result
aa <- data.frame(table(tsne_result_1$celltype))
for (k in 1:nrow(aa)) {
  aaa <- subset(tsne_result_1,celltype == as.character(aa[k,1]))
  clitype_unique<-unique(aaa$cli_type)
  aa_cell_result <- c()
  #aa_cell_type_result <- c()
  aa_cli_type_result <- c()
  aa_tsne1_result <- c()
  aa_tsne2_result <- c()
  aa_exp_result <- c()
  aa_var_result <- c()
  aa_max_result <- c()
  for(l in 1:length(clitype_unique)){
    aaaa<-aaa[aaa$cli_type %in% clitype_unique[l],]
    
    aa_cell_result <- c(aa_cell_result,get.paste(aaaa$cell))
    #aa_cell_type_result <- paste(c(aa_cell_type_result,aa_cell_type),collapse=";")
    aa_cli_type_result <- c(aa_cli_type_result,get.paste(aaaa$cli_type))
    aa_tsne1_result <- c(aa_tsne1_result,get.paste(round(aaaa$tSNE_1,2)))
    aa_tsne2_result <- c(aa_tsne2_result,get.paste(round(aaaa$tSNE_2,2)))
    
    tissue_exp1 <- tissue_exp[,aaaa$cell]
    tissue_exp1 <- round(tissue_exp1,4)
    result1 <- apply(tissue_exp1,1,get.paste)
    aa_exp_result <- cbind(aa_exp_result,result1)
    #colnames(aa_exp_result)[k] <- as.character(aa[k,1])
    
    tissue_var1 <- tissue_var[,aaaa$cell]
    #tissue_var1 <- round(tissue_var1,4)
    result1 <- apply(tissue_var1,1,get.paste)
    aa_var_result <- cbind(aa_var_result,result1)
    #colnames(aa_var_result)[k] <- as.character(aa[k,1])
    
    tissue_exp1 <- tissue_exp[,aaaa$cell]
    tissue_exp1 <- round(tissue_exp1,4)
    result1 <- apply(tissue_exp1,1,max)
    aa_max_result <- cbind(aa_max_result,result1)
    #colnames(aa_max_result)[k] <- as.character(aa[k,1])
  }
  result3 <- apply(aa_exp_result,1,get.paste2)
  result_exp <- data.frame(Gene=rownames(aa_exp_result),result3)
  #result_exp <- merge(result_exp,result2,by="Gene",all=T)
  #colnames(result_exp)[2] <- Sample
  
  result3 <- apply(aa_var_result,1,get.paste2)
  result_var <- data.frame(Gene=rownames(aa_var_result),result3)
  #result_var <- merge(result_var,result2,by="Gene",all=T)
  #colnames(result_var)[2] <- Sample
  
  result3 <- apply(aa_max_result,1,max)
  result_max <- data.frame(Gene=rownames(aa_max_result),result3)
  #result_max <- merge(result_max,result2,by="Gene",all=T)
  #colnames(result_max)[2] <- Sample
  
  aa_cell_result <- get.paste2(aa_cell_result)
  aa_cli_type_result <- get.paste2(aa_cli_type_result)
  aa_tsne1_result <- get.paste2(aa_tsne1_result)
  aa_tsne2_result <- get.paste2(aa_tsne2_result)
  #result_cell_type <- c(result_cell_type,aa_cell_type_result)
  #result_cli_type <- c(result_cli_type,aa_cli_type_result)
  #result_tsne_x <- c(result_tsne_x,aa_tsne1_result)
  #result_tsne_y <- c(result_tsne_y,aa_tsne2_result)
  
  aa_cell_result <- c("cell",aa_cell_result)
  #result_cell_type <- c("cell_type",result_cell_type)
  aa_cli_type_result <- c("cli_type",aa_cli_type_result)
  aa_tsne1_result <- c("tsne_x",aa_tsne1_result)
  aa_tsne2_result <- c("tsne_y",aa_tsne2_result )
  
  result_exp <- rbind(aa_cell_result,rbind(aa_cli_type_result,rbind(aa_tsne1_result,rbind(aa_tsne2_result,rbind(result_exp,result_var)))))
  
  if(k ==1){
    final_result_exp<-result_exp
    final_result_max<-result_max
  }else{
    final_result_exp <- cbind(final_result_exp,result_exp[,2])
    final_result_max <- cbind(final_result_max,result_max[,2])
  }
  colnames(final_result_exp)[k+1]<-as.character(aa[k,1])
  colnames(final_result_max)[k+1]<-as.character(aa[k,1])
}
final_result_exp$Gene<-gsub(pattern = ":|>|,",replacement = "_",x=final_result_exp$Gene)
final_result_max$Gene<-gsub(pattern = ":|>|,",replacement = "_",x=final_result_max$Gene)
write.table(final_result_exp,file=paste(GSE_id,".cell.cell_type.tsne_x.tsne_y.exp.round_4y_1.txt",sep=""),sep="\t",quote=F,row.names=F)
write.table(final_result_max,file=paste(GSE_id,".cell.cell_type.tsne_x.tsne_y.exp.round_4y_max_1.txt",sep=""),sep="\t",quote=F,row.names=F)



#############################################################################################################################################################
#3.celltype_clitype_mut_cellcount.R##################################################################################################################################################
#############################################################################################################################################################
rm(list = ls()[-match(c("GSE_id","input"),ls())])
get.paste1 <- function(x) {
  xx <- paste(x,collapse=",")
  return(xx)
}
get.paste2 <- function(x) {
  xx <- paste(x,collapse=";")
  return(xx)
}

library(dplyr)

setwd(paste(input,GSE_id,sep=""))
file1=read.table(paste0(GSE_id,"_cell_50_var_mat.txt"),header=T,sep="\t",as.is=T)
##位点有检测的细胞
file1 <- unique(file1[,c("mut_name","cell_AF","cell","REDIportal","celltype","cli_type")])
##所有的细胞[表达谱、突变谱所用的细胞]
anno_value <-read.table(paste(GSE_id,"_cell_type.txt",sep=""),header=T,sep="\t",as.is=T)
celltype_count<-as.data.frame(table(anno_value$celltype))
bridge_merge<-data.frame(mut_type=c("Mut","Non_mut","Unknown"))
aa <- data.frame(table(file1$celltype))
for (k in 1:nrow(aa)) {
  aaa <- subset(file1,celltype == as.character(aa[k,1]))
  mut_unique<-unique(aaa$mut_name)
  aa_cli_type_result <- c()
  aa_mut_type_result <- c()
  aa_cell_count_result <- c()
  result_exp<-c()
  for(l in 1:length(mut_unique)){
    aaaa<-aaa[aaa$mut_name %in% mut_unique[l],]
    aaaa$mut_type<-"Unknown"
    aaaa$mut_type<-ifelse(aaaa$cell_AF>0,"Mut","Non_mut")
    aa_cli_type <- c()
    aa_mut_type<- c()
    aa_cell_count <- c()
    clitype_unique<-unique(aaaa$cli_type)
    for(m in 1:length(clitype_unique)){
      aaaaa<-aaaa[aaaa$cli_type %in% clitype_unique[m],]
      bbb<-aaaaa %>%
        group_by(mut_name,mut_type) %>%
        count()
      ccc<-aaaaa %>%
        group_by(mut_name) %>%
        count()
      ccc$n<-celltype_count[celltype_count$Var1 %in% as.character(aa[k,1]),2]-ccc$n
      ccc$mut_type<-"Unknown"
      ddd<-rbind(bbb,ccc)
      bridge<-merge(ddd,bridge_merge,all.y=T)
      bridge$mut_name<-mut_unique[l]
      bridge[is.na(bridge)]<-0
      bridge$cli_type<- clitype_unique[m]
      aa_cli_type <- c(aa_cli_type,get.paste1(bridge$cli_type))
      aa_mut_type<- c(aa_mut_type,get.paste1(bridge$mut_type))
      aa_cell_count <- c(aa_cell_count,get.paste1(bridge$n))
    }
    aa_cli_type_result <- c(mut_unique[l],get.paste2(aa_cli_type))
    aa_mut_type_result <- c(mut_unique[l],get.paste2(aa_mut_type))
    aa_cell_count_result <- c(mut_unique[l],get.paste2(aa_cell_count))
    result_1 <- rbind(aa_cli_type_result,rbind(aa_mut_type_result,aa_cell_count_result))
    result_exp<-rbind(result_exp,result_1)
  }
  result_exp<-as.data.frame(result_exp)
  colnames(result_exp)<-c("mut_name",as.character(aa[k,1]))
  if(k ==1){
    final_result_exp<-result_exp
  }else{
    ##相同的mut_name 用cbind,否则用merge
    same_M<-intersect(result_exp[,1],final_result_exp[,1])
    aa_result_1<-final_result_exp[final_result_exp$mut_name %in% same_M,]
    aaa_result_1<-result_exp[result_exp$mut_name %in% same_M,]
    bb<-cbind(aa_result_1[order(aa_result_1$mut_name),],aaa_result_1[order(aaa_result_1$mut_name),2])
    aa_result_2<-final_result_exp[!(final_result_exp$mut_name %in% same_M),]
    aaa_result_2<-result_exp[!(result_exp$mut_name %in% same_M),]
    cc<-merge(aa_result_2,aaa_result_2,by=c("mut_name"),all=T)
    colnames(bb)<-colnames(cc)
    final_result_exp<-rbind(bb,cc)
  }
}
final_result_exp$mut_name<-gsub(pattern = ":|>|,",replacement = "_",x=final_result_exp$mut_name)
write.table(final_result_exp,file=paste(GSE_id,"_celltype_clitype_mut_cellcount.txt",sep=""),sep="\t",quote=F,row.names=F)




#############################################################################################################################################################
#4.monocle.R##################################################################################################################################################
#############################################################################################################################################################
rm(list = ls()[-match(c("GSE_id","input"),ls())])
##########伪时间分析
##参考https://blog.csdn.net/m0_52069102/article/details/126051672
library(Seurat)
library(dplyr)
library(igraph)
setwd(paste(input,GSE_id,sep=""))
load(paste(GSE_id,".RData",sep=""))
##
# counts_1 <- Read10X(data.dir = "H:/single_cell/SKCM/GSE117988/SRR7722937/Monocle/")
# count_object <- CreateSeuratObject(counts = counts_1, project = "myproject", 
#                                    min.cells = 3, min.features = 200)
data <- as(as.matrix(test.seu@assays$RNA@counts), 'sparseMatrix')
library(monocle)
pd <- new('AnnotatedDataFrame', data = test.seu@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

HSMM <- newCellDataSet(data,
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size())
HSMM <- estimateSizeFactors(HSMM)   
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 3) ##cut 3

expressed_genes <- row.names(subset(fData(HSMM),num_cells_expressed >= 10))
##Step 1: choosing genes that define progress
diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],
                                      fullModelFormulaStr = "~percent.mt")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01)) ## 不要也写0.1 ，而是要写0.01。
HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)
#Step 2: reducing the dimensionality of the data 
HSMM <- reduceDimension(HSMM, max_components = 2,
                        method = 'DDRTree')
#Step 3: ordering the cells in pseudotime 
HSMM <- orderCells(HSMM)
plot_cell_trajectory(HSMM, color_by = "seurat_clusters")
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
##提取细胞坐标
bb_monocle<-data.frame(t(as.matrix(HSMM@reducedDimS)))
colnames(bb_monocle)<-c("tSNE_1","tSNE_2")
bb_monocle$Pseudotime<-HSMM@phenoData@data[,"Pseudotime"]
#plot(c(bb_monocle[,1]),c(bb_monocle[,2]))
##提取边坐标
##参考https://blog.csdn.net/wong2016/article/details/124847268
bb_location<-data.frame(t(as.matrix(HSMM@reducedDimK)))
bb_edges<-HSMM@minSpanningTree %>% ends(E(HSMM@minSpanningTree))
bb_edges<-data.frame(bb_edges)
#plot(c(bb[1,]),c(bb[2,]))

##########2.粘贴坐标数据，参考GSE117988_tsne.R
get.paste <- function(x) {
  xx <- paste(x,collapse=",")
  return(xx)
}
get.paste2 <- function(x) {
  xx <- paste(x,collapse=";")
  return(xx)
}
tsne_result<-bb_monocle
tsne_result$cell<-rownames(tsne_result)
anno_value<-read.table(file=paste(GSE_id,"_cell_type.txt",sep=""),sep="\t",header = T)
tsne_result <- merge(tsne_result,anno_value)

tissue_var <-read.table(file=paste(GSE_id,"_filterCB_var_mat.txt",sep=""),sep="\t",header = T)
rownames(tissue_var)<-tissue_var$gene_symbol
colnames(tissue_var)<-gsub(pattern = "\\.",replacement = "-",x=colnames(tissue_var))
tissue_exp <-read.table(file=paste(GSE_id,"_filterG_exp_mat.txt",sep=""),sep="\t",header = T)
rownames(tissue_exp)<-tissue_exp$gene_symbol
tissue_exp<-tissue_exp[,which(colnames(tissue_exp)!="gene_symbol")]
colnames(tissue_exp)<-gsub(pattern = "\\.",replacement = "-",x=colnames(tissue_exp))
tissue_pse<-data.frame(t(bb_monocle[colnames(tissue_exp),]))
colnames(tissue_pse)<-gsub(pattern = "\\.",replacement = "-",x=colnames(tissue_pse))
tissue_exp<-rbind(tissue_pse["Pseudotime",],tissue_exp)

Sample_id<-unique(tsne_result$cli_type)
for(ii in 0:length(Sample_id)){
  if(ii==0){
    Sample<-GSE_id
    tsne_result_1<-tsne_result
  }else{
    Sample<-Sample_id[ii]
    tsne_result_1<-tsne_result[tsne_result$cli_type %in% Sample,]
  }
  result_exp <- matrix(ncol=1,nrow=0)
  colnames(result_exp) <- "Gene"
  result_var <- matrix(ncol=1,nrow=0)
  colnames(result_var) <- "Gene"
  result_max <- matrix(ncol=1,nrow=0)
  colnames(result_max) <- "Gene"
  
  result_cell <- c()
  result_cell_type <- c()
  result_cli_type <- c()
  result_tsne_x <- c()
  result_tsne_y <- c()
  anno_result <- c()
  
  aa <- data.frame(table(tsne_result_1$celltype))
  aa_cell_result <- c()
  aa_cell_type_result <- c()
  aa_cli_type_result <- c()
  aa_tsne1_result <- c()
  aa_tsne2_result <- c()
  aa_exp_result <- c()
  aa_var_result <- c()
  aa_max_result <- c()
  for (k in 1:nrow(aa)) {
    aaa <- subset(tsne_result_1,celltype == as.character(aa[k,1]))
    aa_cell <- paste(aaa$cell,collapse=",")
    aa_cell_type <- paste(aaa$celltype,collapse=",")
    aa_cli_type <- paste(aaa$cli_type,collapse=",")
    aa_tsne1 <- paste(round(aaa$tSNE_1,2),collapse=",")
    aa_tsne2 <- paste(round(aaa$tSNE_2,2),collapse=",")
    
    aa_cell_result <- paste(c(aa_cell_result,aa_cell),collapse=";")
    aa_cell_type_result <- paste(c(aa_cell_type_result,aa_cell_type),collapse=";")
    aa_cli_type_result <- paste(c(aa_cli_type_result,aa_cli_type),collapse=";")
    aa_tsne1_result <- paste(c(aa_tsne1_result,aa_tsne1),collapse=";")
    aa_tsne2_result <- paste(c(aa_tsne2_result,aa_tsne2),collapse=";")
    
    tissue_exp1 <- tissue_exp[,aaa$cell]
    tissue_exp1 <- round(tissue_exp1,4)
    result1 <- apply(tissue_exp1,1,get.paste)
    aa_exp_result <- cbind(aa_exp_result,result1)
    colnames(aa_exp_result)[k] <- as.character(aa[k,1])
    
    tissue_var1 <- tissue_var[,aaa$cell]
    #tissue_var1 <- round(tissue_var1,4)
    result1 <- apply(tissue_var1,1,get.paste)
    aa_var_result <- cbind(aa_var_result,result1)
    colnames(aa_var_result)[k] <- as.character(aa[k,1])
    
    tissue_exp1 <- tissue_exp[,aaa$cell]
    tissue_exp1 <- round(tissue_exp1,4)
    result1 <- apply(tissue_exp1,1,max)
    aa_max_result <- cbind(aa_max_result,result1)
    colnames(aa_max_result)[k] <- as.character(aa[k,1])
  }
  result3 <- apply(aa_exp_result,1,get.paste2)
  result2 <- data.frame(Gene=rownames(aa_exp_result),result3)
  result_exp <- merge(result_exp,result2,by="Gene",all=T)
  colnames(result_exp)[2] <- Sample
  
  result3 <- apply(aa_var_result,1,get.paste2)
  result2 <- data.frame(Gene=rownames(aa_var_result),result3)
  result_var <- merge(result_var,result2,by="Gene",all=T)
  colnames(result_var)[2] <- Sample
  
  result3 <- apply(aa_max_result,1,max)
  result2 <- data.frame(Gene=rownames(aa_max_result),result3)
  result_max <- merge(result_max,result2,by="Gene",all=T)
  colnames(result_max)[2] <- Sample
  
  result_cell <- c(result_cell,aa_cell_result)
  
  result_cell_type <- c(result_cell_type,aa_cell_type_result)
  result_cli_type <- c(result_cli_type,aa_cli_type_result)
  result_tsne_x <- c(result_tsne_x,aa_tsne1_result)
  result_tsne_y <- c(result_tsne_y,aa_tsne2_result)
  
  result_cell <- c("cell",result_cell)
  result_cell_type <- c("cell_type",result_cell_type)
  result_cli_type <- c("cli_type",result_cli_type)
  result_tsne_x <- c("tsne_x",result_tsne_x)
  result_tsne_y <- c("tsne_y",result_tsne_y)
  
  result_exp <- rbind(result_cell,rbind(result_cell_type,rbind(result_cli_type,rbind(result_tsne_x,rbind(result_tsne_y,rbind(result_exp,result_var))))))
  #result_exp[is.na(result_exp)] <- "NONE"
  if(ii ==0){
    final_result_exp<-result_exp
    final_result_max<-result_max
  }else{
    final_result_exp <- cbind(final_result_exp,result_exp[,2])
    colnames(final_result_exp)[ii+2]<-Sample
    final_result_max <- cbind(final_result_max,result_max[,2])
    colnames(final_result_max)[ii+2]<-Sample
  }
}
write.table(final_result_max,file=paste(GSE_id,".cell.cell_type.Pse_x.tsne_y.Pse.round_4y_max.txt",sep=""),sep="\t",quote=F,row.names=F)

get.paste_pse <- function(x) {
  x_1 <- paste0("[",paste(bb_location[x[1],],collapse=","),"]")
  x_2 <- paste0("[",paste(bb_location[x[2],],collapse=","),"]")
  x_3 <- paste0(x_1,",",x_2)
  return(x_3)
}
bb_location<-apply(bb_location,2,FUN = function(x){round(x,2)})
edges_coordinate<-apply(bb_edges,1,get.paste_pse)
edges_coordinate<-paste(edges_coordinate,collapse=";")
final_result_exp<-data.frame(final_result_exp)
row_index<-nrow(final_result_exp)+1
final_result_exp[row_index,]<-c("edges_coordinate",rep(edges_coordinate,ncol(final_result_exp)-1))
final_result_exp[,1]<-gsub(pattern = ":|>|,",replacement = "_",x=final_result_exp[,1])
write.table(final_result_exp,file=paste(GSE_id,".cell.cell_type.Pse_x.tsne_y.Pse.round_4y.txt",sep=""),sep="\t",quote=F,row.names=F)




#############################################################################################################################################################
#5.cellchat_2.R##################################################################################################################################################
#############################################################################################################################################################
rm(list = ls()[-match(c("GSE_id","input"),ls())])
cellcell_LR<-function (object, sources.use = NULL, targets.use = NULL, signaling = NULL, 
                       pairLR.use = NULL, sort.by.source = FALSE, sort.by.target = FALSE, 
                       sort.by.source.priority = TRUE, color.heatmap = c("Spectral", 
                                                                         "viridis"), n.colors = 10, direction = -1, thresh = 0.05, 
                       comparison = NULL, group = NULL, remove.isolate = FALSE, 
                       max.dataset = NULL, min.dataset = NULL, min.quantile = 0, 
                       max.quantile = 1, line.on = TRUE, line.size = 0.2, color.text.use = TRUE, 
                       color.text = NULL, title.name = NULL, font.size = 10, font.size.title = 10, 
                       show.legend = TRUE, grid.on = TRUE, color.grid = "grey90", 
                       angle.x = 90, vjust.x = NULL, hjust.x = NULL, return.data = FALSE) 
{
  color.heatmap <- match.arg(color.heatmap)
  if (is.list(object@net[[1]])) {
    message("Comparing communications on a merged object \n")
  }
  else {
    message("Comparing communications on a single object \n")
  }
  if (is.null(vjust.x) | is.null(hjust.x)) {
    angle = c(0, 45, 90)
    hjust = c(0, 1, 1)
    vjust = c(0, 1, 0.5)
    vjust.x = vjust[angle == angle.x]
    hjust.x = hjust[angle == angle.x]
  }
  if (length(color.heatmap) == 1) {
    color.use <- tryCatch({
      RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
    }, error = function(e) {
      (scales::viridis_pal(option = color.heatmap, direction = -1))(n.colors)
    })
  }
  else {
    color.use <- color.heatmap
  }
  if (direction == -1) {
    color.use <- rev(color.use)
  }
  if (!is.null(pairLR.use)) {
    if (!is.data.frame(pairLR.use)) {
      stop("pairLR.use should be a data frame with a signle column named either 'interaction_name' or 'pathway_name' ")
    }
    else if ("pathway_name" %in% colnames(pairLR.use)) {
      pairLR.use$pathway_name <- as.character(pairLR.use$pathway_name)
    }
    else if ("interaction_name" %in% colnames(pairLR.use)) {
      pairLR.use$interaction_name <- as.character(pairLR.use$interaction_name)
    }
  }
  if (is.null(comparison)) {
    cells.level <- levels(object@idents)
    if (is.numeric(sources.use)) {
      sources.use <- cells.level[sources.use]
    }
    if (is.numeric(targets.use)) {
      targets.use <- cells.level[targets.use]
    }
    df.net <- subsetCommunication(object, slot.name = "net", 
                                  sources.use = sources.use, targets.use = targets.use, 
                                  signaling = signaling, pairLR.use = pairLR.use, 
                                  thresh = thresh)
    df.net$source.target <- paste(df.net$source, df.net$target, 
                                  sep = " -> ")
    source.target <- paste(rep(sources.use, each = length(targets.use)), 
                           targets.use, sep = " -> ")
    source.target.isolate <- setdiff(source.target, unique(df.net$source.target))
    if (length(source.target.isolate) > 0) {
      df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate), 
                                             ncol = ncol(df.net)))
      colnames(df.net.isolate) <- colnames(df.net)
      df.net.isolate$source.target <- source.target.isolate
      df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
      df.net.isolate$pval <- 1
      a <- stringr::str_split(df.net.isolate$source.target, 
                              " -> ", simplify = T)
      df.net.isolate$source <- as.character(a[, 1])
      df.net.isolate$target <- as.character(a[, 2])
      df.net <- rbind(df.net, df.net.isolate)
    }
    df.net$pval[df.net$pval > 0.05] = 1
    df.net$pval[df.net$pval > 0.01 & df.net$pval <= 0.05] = 2
    df.net$pval[df.net$pval <= 0.01] = 3
    df.net$prob[df.net$prob == 0] <- NA
    df.net$prob.original <- df.net$prob
    df.net$prob <- -1/log(df.net$prob)
    idx1 <- which(is.infinite(df.net$prob) | df.net$prob < 
                    0)
    if (sum(idx1) > 0) {
      values.assign <- seq(max(df.net$prob, na.rm = T) * 
                             1.1, max(df.net$prob, na.rm = T) * 1.5, length.out = length(idx1))
      position <- sort(prob.original[idx1], index.return = TRUE)$ix
      df.net$prob[idx1] <- values.assign[match(1:length(idx1), 
                                               position)]
    }
    df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% 
                                                                  unique(df.net$source)])
    df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% 
                                                                  unique(df.net$target)])
    group.names <- paste(rep(levels(df.net$source), each = length(levels(df.net$target))), 
                         levels(df.net$target), sep = " -> ")
    df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
    df.net <- with(df.net, df.net[order(interaction_name_2), 
    ])
    df.net$interaction_name_2 <- factor(df.net$interaction_name_2, 
                                        levels = unique(df.net$interaction_name_2))
    cells.order <- group.names
    df.net$source.target <- factor(df.net$source.target, 
                                   levels = cells.order)
    df <- df.net
  }
  else {
    dataset.name <- names(object@net)
    df.net.all <- subsetCommunication(object, slot.name = "net", 
                                      sources.use = sources.use, targets.use = targets.use, 
                                      signaling = signaling, pairLR.use = pairLR.use, 
                                      thresh = thresh)
    df.all <- data.frame()
    for (ii in 1:length(comparison)) {
      cells.level <- levels(object@idents[[comparison[ii]]])
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      df.net <- df.net.all[[comparison[ii]]]
      df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
      df.net$source.target <- paste(df.net$source, df.net$target, 
                                    sep = " -> ")
      source.target <- paste(rep(sources.use, each = length(targets.use)), 
                             targets.use, sep = " -> ")
      source.target.isolate <- setdiff(source.target, 
                                       unique(df.net$source.target))
      if (length(source.target.isolate) > 0) {
        df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate), 
                                               ncol = ncol(df.net)))
        colnames(df.net.isolate) <- colnames(df.net)
        df.net.isolate$source.target <- source.target.isolate
        df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
        df.net.isolate$pval <- 1
        a <- stringr::str_split(df.net.isolate$source.target, 
                                " -> ", simplify = T)
        df.net.isolate$source <- as.character(a[, 1])
        df.net.isolate$target <- as.character(a[, 2])
        df.net <- rbind(df.net, df.net.isolate)
      }
      df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% 
                                                                    unique(df.net$source)])
      df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% 
                                                                    unique(df.net$target)])
      group.names <- paste(rep(levels(df.net$source), 
                               each = length(levels(df.net$target))), levels(df.net$target), 
                           sep = " -> ")
      group.names0 <- group.names
      group.names <- paste0(group.names0, " (", dataset.name[comparison[ii]], 
                            ")")
      if (nrow(df.net) > 0) {
        df.net$pval[df.net$pval > 0.05] = 1
        df.net$pval[df.net$pval > 0.01 & df.net$pval <= 
                      0.05] = 2
        df.net$pval[df.net$pval <= 0.01] = 3
        df.net$prob[df.net$prob == 0] <- NA
        df.net$prob.original <- df.net$prob
        df.net$prob <- -1/log(df.net$prob)
      }
      else {
        df.net <- as.data.frame(matrix(NA, nrow = length(group.names), 
                                       ncol = 5))
        colnames(df.net) <- c("interaction_name_2", 
                              "source.target", "prob", "pval", "prob.original")
        df.net$source.target <- group.names0
      }
      df.net$group.names <- as.character(df.net$source.target)
      df.net$source.target <- paste0(df.net$source.target, 
                                     " (", dataset.name[comparison[ii]], ")")
      df.net$dataset <- dataset.name[comparison[ii]]
      df.all <- rbind(df.all, df.net)
    }
    if (nrow(df.all) == 0) {
      stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
    }
    idx1 <- which(is.infinite(df.all$prob) | df.all$prob < 
                    0)
    if (sum(idx1) > 0) {
      values.assign <- seq(max(df.all$prob, na.rm = T) * 
                             1.1, max(df.all$prob, na.rm = T) * 1.5, length.out = length(idx1))
      position <- sort(df.all$prob.original[idx1], index.return = TRUE)$ix
      df.all$prob[idx1] <- values.assign[match(1:length(idx1), 
                                               position)]
    }
    df.all$interaction_name_2[is.na(df.all$interaction_name_2)] <- df.all$interaction_name_2[!is.na(df.all$interaction_name_2)][1]
    df <- df.all
    df <- with(df, df[order(interaction_name_2), ])
    df$interaction_name_2 <- factor(df$interaction_name_2, 
                                    levels = unique(df$interaction_name_2))
    cells.order <- c()
    dataset.name.order <- c()
    for (i in 1:length(group.names0)) {
      for (j in 1:length(comparison)) {
        cells.order <- c(cells.order, paste0(group.names0[i], 
                                             " (", dataset.name[comparison[j]], ")"))
        dataset.name.order <- c(dataset.name.order, 
                                dataset.name[comparison[j]])
      }
    }
    df$source.target <- factor(df$source.target, levels = cells.order)
  }
  min.cutoff <- quantile(df$prob, min.quantile, na.rm = T)
  max.cutoff <- quantile(df$prob, max.quantile, na.rm = T)
  df$prob[df$prob < min.cutoff] <- min.cutoff
  df$prob[df$prob > max.cutoff] <- max.cutoff
  if (remove.isolate) {
    df <- df[!is.na(df$prob), ]
    line.on <- FALSE
  }
  if (!is.null(max.dataset)) {
    signaling <- as.character(unique(df$interaction_name_2))
    for (i in signaling) {
      df.i <- df[df$interaction_name_2 == i, , drop = FALSE]
      cell <- as.character(unique(df.i$group.names))
      for (j in cell) {
        df.i.j <- df.i[df.i$group.names == j, , drop = FALSE]
        values <- df.i.j$prob
        idx.max <- which(values == max(values, na.rm = T))
        idx.min <- which(values == min(values, na.rm = T))
        dataset.na <- c(df.i.j$dataset[is.na(values)], 
                        setdiff(dataset.name[comparison], df.i.j$dataset))
        if (length(idx.max) > 0) {
          if (!(df.i.j$dataset[idx.max] %in% dataset.name[max.dataset])) {
            df.i.j$prob <- NA
          }
          else if ((idx.max != idx.min) & !is.null(min.dataset)) {
            if (!(df.i.j$dataset[idx.min] %in% dataset.name[min.dataset])) {
              df.i.j$prob <- NA
            }
            else if (length(dataset.na) > 0 & sum(!(dataset.name[min.dataset] %in% 
                                                    dataset.na)) > 0) {
              df.i.j$prob <- NA
            }
          }
        }
        df.i[df.i$group.names == j, "prob"] <- df.i.j$prob
      }
      df[df$interaction_name_2 == i, "prob"] <- df.i$prob
    }
  }
  if (remove.isolate) {
    df <- df[!is.na(df$prob), ]
    line.on <- FALSE
  }
  if (nrow(df) == 0) {
    stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
  }
  if (!is.null(pairLR.use)) {
    interaction_name_2.order <- intersect(object@DB$interaction[pairLR.use$interaction_name, 
    ]$interaction_name_2, unique(df$interaction_name_2))
    df$interaction_name_2 <- factor(df$interaction_name_2, 
                                    levels = interaction_name_2.order)
  }
  df$source.target = droplevels(df$source.target, exclude = setdiff(levels(df$source.target), 
                                                                    unique(df$source.target)))
  if (sort.by.target & !sort.by.source) {
    if (!is.null(targets.use)) {
      df$target <- factor(df$target, levels = intersect(targets.use, 
                                                        df$target))
      df <- with(df, df[order(target, source), ])
      source.target.order <- unique(as.character(df$source.target))
      df$source.target <- factor(df$source.target, levels = source.target.order)
    }
  }
  if (sort.by.source & !sort.by.target) {
    if (!is.null(sources.use)) {
      df$source <- factor(df$source, levels = intersect(sources.use, 
                                                        df$source))
      df <- with(df, df[order(source, target), ])
      source.target.order <- unique(as.character(df$source.target))
      df$source.target <- factor(df$source.target, levels = source.target.order)
    }
  }
  if (sort.by.source & sort.by.target) {
    if (!is.null(sources.use)) {
      df$source <- factor(df$source, levels = intersect(sources.use, 
                                                        df$source))
      if (!is.null(targets.use)) {
        df$target <- factor(df$target, levels = intersect(targets.use, 
                                                          df$target))
      }
      if (sort.by.source.priority) {
        df <- with(df, df[order(source, target), ])
      }
      else {
        df <- with(df, df[order(target, source), ])
      }
      source.target.order <- unique(as.character(df$source.target))
      df$source.target <- factor(df$source.target, levels = source.target.order)
    }
  }
  return(df)
}

suppressMessages({
  library(Seurat)
  #library(SeuratWrappers)
  #library(SeuratDisk)
  library(CellChat)
  library(ggplot2)
  library(ggalluvial)
  library(NMF)
})


setwd(paste(input,GSE_id,sep=""))
load(paste(GSE_id,".RData",sep=""))
##添加cancer细胞注释
cell_type<-read.table(paste0(GSE_id,"_cell_type.txt"),header=T,sep="\t",as.is=T)
rownames(cell_type)<-cell_type$cell
Idents(test.seu)<-cell_type[names(Idents(test.seu)),"celltype"]

file1=read.table(paste0(GSE_id,"_cell_50_var_mat.txt"),header=T,sep="\t",as.is=T)
file1<-unique(file1[,c("cell","mut_name","type")])
##提取前100位点
mut_freq<-as.data.frame(table(unique(file1[,c("mut_name","cell")])$mut_name))
library(dplyr)
r2=mut_freq %>% top_n(n=100,wt=Freq)
file1<-file1[file1$mut_name %in% r2[,1],]
cell_mut_type<-merge(cell_type,file1)
#Sample<-unique(test.seu@meta.data$orig.ident)
sample_type<-read.table("D:/sample_type.txt",header=F,sep="\t",as.is=T)
sample_type<-sample_type[sample_type[,1] %in% GSE_id,]
Sample_id<-unique(cell_mut_type$cli_type)
chat_result<-data.frame(Mut_name=character(),cell_type=character(),cli_type=character(),From=character(),To=character(),edge_count=numeric(),From_size=numeric(),To_size=numeric())
LR_result<-data.frame(Mut_name=character(),cell_type=character(),cli_type=character(),source=character(),target=character(),ligand=character(),receptor=character(),prob=numeric(),pval=numeric(),interaction_name=character(),interaction_name_2=character(),pathway_name=character(),annotation=character(),evidence=character(),source.target=character(),prob.original=numeric())
for(sample_i in Sample_id){
  sample_type_1<-sample_type[sample_type[,4] %in% sample_i,]
  #获取seurat_obj的子集：基于metadata的值
  test.seu.1<-subset(x = test.seu, subset = orig.ident == sample_type_1[,3])
  cell_mut_type_1<-cell_mut_type[cell_mut_type$cli_type %in% sample_i,]
  cell_type_1<-cell_type[cell_type$cli_type %in% sample_i,]
  for(mut_i in unique(cell_mut_type_1$mut_name)){
    cell_mut_type_2<-cell_mut_type_1[cell_mut_type_1$mut_name %in% mut_i,]
    for(cell_type_i in unique(cell_mut_type_2$celltype)){
      cell_mut_type_3<-cell_mut_type_2[cell_mut_type_2$celltype %in% cell_type_i,]
      cell_mut_type_4<-merge(cell_type_1,cell_mut_type_3,all=T)
      cell_mut_type_4$mut_type<-apply(cell_mut_type_4,1,function(x){
        if(x[2]==cell_type_i & is.na(x[5])){
          mut_type<-"Unknown"
        }
        if(x[2]==cell_type_i & !is.na(x[5])){
          if(as.numeric(x[5])==0){
            mut_type<-"Non_mut"
          }
          if(as.numeric(x[5])==1){
            mut_type<-"Mut"
          }
        }
        if(x[2]!=cell_type_i){
          mut_type<-x[2]
        }
        return(mut_type)
      })
      if(length(which(cell_mut_type_4$mut_type=="Mut"))>50){
        ##添加mut_i定义的细胞注释
        rownames(cell_mut_type_4)<-cell_mut_type_4$cell
        Idents(test.seu.1)<-cell_mut_type_4[names(Idents(test.seu.1)),"mut_type"]
        ##1.创建cellchat分析对象，导入数据库
        cellchat <- createCellChat( object = test.seu.1) #创建cellchat对象
        #levels(cellchat@idents)#查看cellchat的分群情况
        groupSize <- as.numeric(table(cellchat@idents))
        groupSize#每个cluster中的细胞数
        CellChatDB <- CellChatDB.human #导入需要用的数据库
        #如果是鼠的数据 CellChatDB <- CellChatDB.mouse
        # colnames(CellChatDB$interaction)#可以看一下cellDB$interaction中都包含什么
        #Show the description of CellChatDB databse
        #showDatabaseCategory(CellChatDB) #可以看下数据库中包含什么
        
        ##2.初步分析
        #可以选择我们感兴趣的通路进行分析，在这里选择了Secreted Signaling，大家可以跟进自己的需求进行选择
        CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
        # 选择我们想使用的数据库
        cellchat@DB <- CellChatDB.use
        #提取signaling gene，结果保存在data.signaling
        cellchat <- subsetData(cellchat)
        #head(cellchat@data.signaling)
        future::plan("multisession", workers = 4) #实现并行，不着急的可以不用这个，一步一步跑
        #寻找每个细胞中高表达的供体和受体
        cellchat <- identifyOverExpressedGenes(cellchat)
        cellchat <- identifyOverExpressedInteractions(cellchat)#结果存在cellchat@LR$LRsig中 #上一步运行的结果存在cellchat@LR$LRsig中
        cellchat <- projectData(cellchat, PPI.mouse) #A diffusion process is used to smooth genes’ expression values based on their neighbors’ defined in a high-confidence experimentally validated protein-protein network.结果保存在@data.project
        # 推断细胞通讯网络
        cellchat <- computeCommunProb(cellchat, raw.use = TRUE) #Compute the communication probability/strength between any interacting cell groups
        # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
        cellchat <- filterCommunication(cellchat, min.cells = 10) #过滤掉10个一下细胞的群
        #推测细胞互作的概率
        cellchat <- computeCommunProbPathway(cellchat) #Compute the communication probability on signaling pathway level by summarizing all related ligands/receptors
        df.netp<-subsetCommunication(cellchat,slot.name='netP')
        #write.csv(df.netp,'net_pathway.csv')
        #Calculate the aggregated network by counting the number of links or summarizing the communication probability
        cellchat <- aggregateNet(cellchat) 
        
        ##3.细胞之间通讯的数量和强度
        #netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
        #netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
        aaa<-as.data.frame(cellchat@net$count)
        aaa$From<-rownames(aaa)
        library(reshape2)
        bbb<-melt(aaa,id.vars="From",variable.name="To",value.name="edge_count")
        From_size<-data.frame(table(cellchat@idents))
        colnames(From_size)<-c("From","From_size")
        To_size<-data.frame(table(cellchat@idents))
        colnames(To_size)<-c("To","To_size")
        bbb<-merge(merge(bbb,From_size),To_size)
        bbb$Mut_name<-mut_i
        bbb$cell_type<-cell_type_i
        bbb$cli_type<-sample_i
        bbb<-bbb[,c(6,7,8,2,1,3,4,5)]
        chat_result<-rbind(chat_result,bbb)
        ##受体配体互作可能性prob,pvalue
        ccc<-cellcell_LR(cellchat)
        ccc$Mut_name<-mut_i
        ccc$cell_type<-cell_type_i
        ccc$cli_type<-sample_i
        LR_result<-rbind(LR_result,ccc)
      }
    }
  }
}
chat_result$Mut_name<-gsub(pattern = ":|>|,",replacement = "_",x=chat_result$Mut_name)
write.table(chat_result,file=paste(GSE_id,"_chat_result.txt",sep=""),sep="\t",quote=F,row.names=F)
write.table(LR_result,file=paste(GSE_id,"_LR_result.txt",sep=""),sep="\t",quote=F,row.names=F)




#############################################################################################################################################################
#6.cellchat_2_1.R##################################################################################################################################################
#############################################################################################################################################################
rm(list = ls()[-match(c("GSE_id","input"),ls())])
get.paste1 <- function(x) {
  xx <- paste(x,collapse=",")
  return(xx)
}
get.paste2 <- function(x) {
  xx <- paste(x,collapse=";")
  return(xx)
}
get.paste3 <- function(x) {
  xx <- paste(x,collapse=":")
  return(xx)
}

setwd(paste(input,GSE_id,sep=""))
data<-read.table(paste(GSE_id,"_LR_result.txt",sep=""),header=T,sep="\t",as.is=T)
final_result<-data.frame(mut_name=character(),cell_type=character(),all_values=character())
Mut_unique<-unique(data$Mut_name)
for(i in 1:length(Mut_unique)){
  data_Mut<-data[data$Mut_name %in% Mut_unique[i],]
  celltype_unique<-unique(data_Mut$cell_type)
  for(j in 1:length(celltype_unique)){
    data_Mut_celltype<-data_Mut[data_Mut$cell_type %in% celltype_unique[j],]
    clitype_unique<-unique(data_Mut_celltype$cli_type)
    cell_re<-c()
    LR_re<-c()
    prob_re<-c()
    pval<-c()
    for(k in 1:length(clitype_unique)){
      data_Mut_celltype_clitype<-data_Mut_celltype[data_Mut_celltype$cli_type %in% clitype_unique[k],]
      re0 = dcast(data=data_Mut_celltype_clitype[,c(8,12,5)],source.target ~interaction_name_2 )
      re1 = melt(data = re0,id.vars=c("source.target"),variable.name="interaction_name_2",value.name="value")
      data_Mut_celltype_clitype<-merge(re1,data_Mut_celltype_clitype,all=T)
      data_Mut_celltype_clitype[is.na(data_Mut_celltype_clitype)]<-"NONE"
      st_unique<-unique(data_Mut_celltype_clitype$source.target)
      cell_re_1<-c()
      LR_re_1<-c()
      prob_re_1<-c()
      pval_1<-c()
      for(l in 1:length(st_unique)){
        bridge<-data_Mut_celltype_clitype[data_Mut_celltype_clitype$source.target %in% st_unique[l],]
        cell_re_1<-c(cell_re_1,get.paste1(bridge$source.target))
        LR_re_1<-c(LR_re_1,get.paste1(bridge$interaction_name_2))
        prob_re_1<-c(prob_re_1,get.paste1(bridge$prob))
        pval_1<-c(pval_1,get.paste1(bridge$pval))
      }
      cell_re<-c(cell_re,get.paste2(cell_re_1))
      LR_re<-c(LR_re,get.paste2(LR_re_1))
      prob_re<-c(prob_re,get.paste2(prob_re_1))
      pval<-c(pval,get.paste2(pval_1))
    }
    nrow_index<-nrow(final_result)+1
    final_result[nrow_index,1]<-Mut_unique[i]
    final_result[nrow_index,2]<-celltype_unique[j]
    final_result[nrow_index,3]<-get.paste3(clitype_unique)
    nrow_index<-nrow(final_result)+1
    final_result[nrow_index,1]<-Mut_unique[i]
    final_result[nrow_index,2]<-celltype_unique[j]
    final_result[nrow_index,3]<-get.paste3(cell_re)
    nrow_index<-nrow(final_result)+1
    final_result[nrow_index,1]<-Mut_unique[i]
    final_result[nrow_index,2]<-celltype_unique[j]
    final_result[nrow_index,3]<-get.paste3(LR_re)
    nrow_index<-nrow(final_result)+1
    final_result[nrow_index,1]<-Mut_unique[i]
    final_result[nrow_index,2]<-celltype_unique[j]
    final_result[nrow_index,3]<-get.paste3(prob_re)
    nrow_index<-nrow(final_result)+1
    final_result[nrow_index,1]<-Mut_unique[i]
    final_result[nrow_index,2]<-celltype_unique[j]
    final_result[nrow_index,3]<-get.paste3(pval)
  }
}
final_result$dataset<-"GSE117988"
final_result$cancer_type<-"NET"
final_result$mut_name<-gsub(pattern = ":|>|,",replacement = "_",x=final_result$mut_name)
write.table(final_result,file=paste(GSE_id,"_LR_result_1.txt",sep=""),sep="\t",quote=F,row.names=F)




#############################################################################################################################################################
#7.diffex_5.R##################################################################################################################################################
#############################################################################################################################################################
rm(list = ls()[-match(c("GSE_id","input"),ls())])
get.paste <- function(x) {
  xx <- paste(x,collapse=",")
  return(xx)
}
get.paste1 <- function(x) {
  xx <- paste(x,collapse=":")
  return(xx)
}
get.paste2 <- function(x) {
  xx <- paste(x,collapse=";")
  return(xx)
}

setwd(paste(input,GSE_id,sep=""))
##计算每个变异的平均水平，再把高于它的细胞作为高组，否则低组
data<-read.table(paste0(GSE_id,"_cell_50_var_mat.txt"),header=T,sep="\t",as.is=T)
dat_exp<-read.table(paste0(GSE_id,"_filterG_exp_mat.txt"),header=T,sep="\t",as.is=T)
library(reshape2)
dat_exp_mat<-melt(dat_exp,id.vars=c("gene_symbol"),variable.name="cell",value.name="exp_value")
dat_exp_mat$cell<-gsub(pattern = "\\.",replacement ="-",dat_exp_mat$cell)
dat_var_cell_G<-data[,c("cell_AF","mut_name","cell","REDIportal","gene_symbol")]#最小唯一组合：mut_name	cell  gene_symbol
##提取前1000位点
mut_freq<-as.data.frame(table(unique(dat_var_cell_G[,c("mut_name","cell")])$mut_name))
library(dplyr)
r2=mut_freq %>% top_n(n=1000,wt=Freq)
dat_var_cell_G<-dat_var_cell_G[dat_var_cell_G$mut_name %in% r2[,1],]

cell_type<-read.table(paste0(GSE_id,"_cell_type.txt"),header=T,sep="\t",as.is=T)
colnames(cell_type)[1]<-"cell"
dat_var_cell_type_G<-merge(dat_var_cell_G,cell_type)
dat_var_cell_type_G_exp<-merge(dat_var_cell_type_G,dat_exp_mat)
dat_exp_mat_clitype<-merge(dat_exp_mat,cell_type)

dat_var_cell_type_G_exp$exp_value<-round(dat_var_cell_type_G_exp$exp_value,4)
#dat_var_cell_type_G_exp<-dat_var_cell_type_G_exp[which(dat_var_cell_type_G_exp$cell_AF>0),]#使用编辑|突变水平大于0的数据
dat_var_cell_type_G_exp$dataset<-GSE_id
write.table(dat_var_cell_type_G_exp[which(dat_var_cell_type_G_exp$cell_AF>0),],file=paste(GSE_id,"dat_var_cell_type_G_exp_4.txt",sep=""),sep="\t",quote=F,row.names=F)

celltype_unique<-unique(dat_var_cell_type_G_exp$celltype)
for(ii in 1:length(celltype_unique)){
  exp_result<-dat_var_cell_type_G_exp[dat_var_cell_type_G_exp$celltype %in% celltype_unique[ii],]
  mut_G <- unique(exp_result[,c("mut_name","gene_symbol")])
  aa_result<-data.frame(mut_name=character(),gene_symbol=character(),value=character())
  for (k in 1:nrow(mut_G)) {
    aaa <- merge(exp_result,mut_G[k,])
    clitype_unique<-unique(aaa$cli_type)
    cell_mut_type_result<-c()
    cell_cli_type_result<-c()
    cell_G_exp_result<-c()
    cell_mut_AF_result<-c()
    #meanAF_result<-c()
    #meanGexp_result<-c()
    xy_logFC_result<-c()
    xy_p_value_result<-c()
    xz_logFC_result<-c()
    xz_p_value_result<-c()
    yz_logFC_result<-c()
    yz_p_value_result<-c()
    cor_result<-c()
    p_result<-c()
    for(j in 1:length(clitype_unique)){
      analysis_dat <- aaa[aaa$cli_type %in% clitype_unique[j],]
      analysis_dat$mut_type<-ifelse(analysis_dat$cell_AF>0,"Mut","Non_mut")
      x <- analysis_dat[analysis_dat$mut_type %in% "Mut",]$exp_value
      y <- analysis_dat[analysis_dat$mut_type %in% "Non_mut",]$exp_value
      bridge_3<-dat_exp_mat_clitype[(dat_exp_mat_clitype$gene_symbol %in% mut_G[k,2]) & !(dat_exp_mat_clitype$cell %in% analysis_dat$cell) & (dat_exp_mat_clitype$celltype %in% celltype_unique[ii]) & (dat_exp_mat_clitype$cli_type %in% clitype_unique[j]),]
      z <- bridge_3$exp_value
      #meanGexp <- round(mean(c(analysis_dat$exp_value,bridge_3$exp_value)),4)
      if(length(x)*length(y)>0){
        xy_logFC<-round(mean(x+1)/mean(y+1),4)
        xy_p_value<-wilcox.test(x,y)$p.value
        if(xy_p_value=="NaN"){
          xy_p_value<-"NONE"
        }else{
          xy_p_value<-round(xy_p_value,4)
        }
      }else{
        xy_logFC<-"NONE"
        xy_p_value<-"NONE"
      }
      if(length(x)*length(z)>0){
        xz_logFC<-round(mean(x+1)/mean(z+1),4)
        xz_p_value<-wilcox.test(x,z)$p.value
        if(xz_p_value=="NaN"){
          xz_p_value<-"NONE"
        }else{
          xz_p_value<-round(xz_p_value,4)
        }
      }else{
        xz_logFC<-"NONE"
        xz_p_value<-"NONE"
      }
      if(length(y)*length(z)>0){
        yz_logFC<-round(mean(y+1)/mean(z+1),4)
        yz_p_value<-wilcox.test(y,z)$p.value
        if(yz_p_value=="NaN"){
          yz_p_value<-"NONE"
        }else{
          yz_p_value<-round(yz_p_value,4)
        }
      }else{
        yz_logFC<-"NONE"
        yz_p_value<-"NONE"
      }
      analysis_dat_cor<-analysis_dat[analysis_dat$mut_type %in% "Mut",]
      # if(nrow(analysis_dat_cor)>=1){
      #   meanAF <- round(mean(analysis_dat_cor$cell_AF),4)
      # }else{
      #   meanAF <- "NONE"
      # }
      if(nrow(analysis_dat_cor)>=3){
        cor<-round(cor.test(analysis_dat_cor$cell_AF,analysis_dat_cor$exp_value)$estimate,4)
        p<-round(cor.test(analysis_dat_cor$cell_AF,analysis_dat_cor$exp_value)$p.value,4)
        if(is.na(cor)){
          cor<-"NONE"
          p<-"NONE"
        }
      }else{
        cor<-"NONE"
        p<-"NONE"
      }
      bridge_1<-analysis_dat[analysis_dat$mut_type %in% "Mut",]
      bridge_2<-analysis_dat[analysis_dat$mut_type %in% "Non_mut",]
      bridge_3<-bridge_3
      if(nrow(bridge_3)>0){
        bridge_3$mut_type<-"Unknown"
        if(nrow(bridge_1)>0 & nrow(bridge_2)>0){
          cell_mut_type<-get.paste1(c(get.paste(bridge_1$mut_type),get.paste(bridge_2$mut_type),get.paste(bridge_3$mut_type)))
          cell_cli_type<-get.paste1(c(get.paste(bridge_1$cli_type),get.paste(bridge_2$cli_type),get.paste(bridge_3$cli_type)))
          cell_G_exp<-get.paste1(c(get.paste(bridge_1$exp_value),get.paste(bridge_2$exp_value),get.paste(bridge_3$exp_value)))
          cell_mut_AF<-get.paste1(c(get.paste(bridge_1$cell_AF),get.paste(bridge_2$cell_AF),"NONE"))
        }
        if(nrow(bridge_1)==0 & nrow(bridge_2)>0){
          cell_mut_type<-get.paste1(c("Mut",get.paste(bridge_2$mut_type),get.paste(bridge_3$mut_type)))
          cell_cli_type<-get.paste1(c(bridge_2$cli_type[1],get.paste(bridge_2$cli_type),get.paste(bridge_3$cli_type)))
          cell_G_exp<-get.paste1(c("NONE",get.paste(bridge_2$exp_value),get.paste(bridge_3$exp_value)))
          cell_mut_AF<-get.paste1(c("NONE",get.paste(bridge_2$cell_AF),"NONE"))
        }
        if(nrow(bridge_1)>0 & nrow(bridge_2)==0){
          cell_mut_type<-get.paste1(c(get.paste(bridge_1$mut_type),"Non_mut",get.paste(bridge_3$mut_type)))
          cell_cli_type<-get.paste1(c(get.paste(bridge_1$cli_type),bridge_1$cli_type[1],get.paste(bridge_3$cli_type)))
          cell_G_exp<-get.paste1(c(get.paste(bridge_1$exp_value),"NONE",get.paste(bridge_3$exp_value)))
          cell_mut_AF<-get.paste1(c(get.paste(bridge_1$cell_AF),"NONE","NONE"))
        }
        if(nrow(bridge_1)==0 & nrow(bridge_2)==0){
          cell_mut_type<-get.paste1(c("Mut","Non_mut",get.paste(bridge_3$mut_type)))
          cell_cli_type<-get.paste1(c(bridge_3$cli_type[1],bridge_3$cli_type[1],get.paste(bridge_3$cli_type)))
          cell_G_exp<-get.paste1(c("NONE","NONE",get.paste(bridge_3$exp_value)))
          cell_mut_AF<-get.paste1(c("NONE","NONE","NONE"))
        }
      }else{
        if(nrow(bridge_1)*nrow(bridge_2)>0){
          cell_mut_type<-get.paste1(c(get.paste(bridge_1$mut_type),get.paste(bridge_2$mut_type),"Unknown"))
          cell_cli_type<-get.paste1(c(get.paste(bridge_1$cli_type),get.paste(bridge_2$cli_type),bridge_1$cli_type[1]))
          cell_G_exp<-get.paste1(c(get.paste(bridge_1$exp_value),get.paste(bridge_2$exp_value),"NONE"))
          cell_mut_AF<-get.paste1(c(get.paste(bridge_1$cell_AF),get.paste(bridge_2$cell_AF),"NONE"))
        }
        if(nrow(bridge_1)==0 & nrow(bridge_2)>0){
          cell_mut_type<-get.paste1(c("Mut",get.paste(bridge_2$mut_type),,"Unknown"))
          cell_cli_type<-get.paste1(c(bridge_2$cli_type[1],get.paste(bridge_2$cli_type),bridge_2$cli_type[1]))
          cell_G_exp<-get.paste1(c("NONE",get.paste(bridge_2$exp_value),"NONE"))
          cell_mut_AF<-get.paste1(c("NONE",get.paste(bridge_2$cell_AF),"NONE"))
        }
        if(nrow(bridge_1)>0 & nrow(bridge_2)==0){
          cell_mut_type<-get.paste1(c(get.paste(bridge_1$mut_type),"Non_mut","Unknown"))
          cell_cli_type<-get.paste1(c(get.paste(bridge_1$cli_type),bridge_1$cli_type[1],bridge_1$cli_type[1]))
          cell_G_exp<-get.paste1(c(get.paste(bridge_1$exp_value),"NONE","NONE"))
          cell_mut_AF<-get.paste1(c(get.paste(bridge_1$cell_AF),"NONE","NONE"))
        }
      }
      cell_mut_type_result<-get.paste2(c(cell_mut_type_result,cell_mut_type))
      cell_cli_type_result<-get.paste2(c(cell_cli_type_result,cell_cli_type))
      cell_G_exp_result<-get.paste2(c(cell_G_exp_result,cell_G_exp))
      cell_mut_AF_result<-get.paste2(c(cell_mut_AF_result,cell_mut_AF))
      #meanAF_result<-c()
      #meanGexp_result<-c()
      xy_logFC_result<-get.paste2(c(xy_logFC_result,xy_logFC))
      xy_p_value_result<-get.paste2(c(xy_p_value_result,xy_p_value))
      xz_logFC_result<-get.paste2(c(xz_logFC_result,xz_logFC))
      xz_p_value_result<-get.paste2(c(xz_p_value_result,xz_p_value))
      yz_logFC_result<-get.paste2(c(yz_logFC_result,yz_logFC))
      yz_p_value_result<-get.paste2(c(yz_p_value_result,yz_p_value))
      cor_result<-get.paste2(c(cor_result,cor))
      p_result<-get.paste2(c(p_result,p))
    }
    aaa_result<-data.frame(value=c(cell_mut_type_result,cell_cli_type_result, cell_G_exp_result,cell_mut_AF_result,xy_logFC_result,xy_p_value_result,xz_logFC_result,xz_p_value_result,yz_logFC_result,yz_p_value_result,cor_result,p_result))
    aaa_result$mut_name<-mut_G[k,1]
    aaa_result$gene_symbol<-mut_G[k,2]
    aa_result<-rbind(aa_result,aaa_result)
  }
  aa_result<-aa_result[,c("mut_name","gene_symbol","value")]
  colnames(aa_result)[3]<-celltype_unique[ii]
  #aa_result[is.na(aa_result)]<-"NONE"
  # aa_result$all_values<-apply(aa_result[,3:ncol(aa_result)],1,function(x){paste(x,collapse = ";")})
  # bridge_result<-rbind(bridge_result,aa_result[,c("mut_name","gene_symbol","all_values")])
  if(ii ==1){
    final_result<-aa_result
  }else{
    ##相同的mut_name  gene_symbol用cbind,否则用merge
    same_MG<-merge(final_result[,c("mut_name","gene_symbol")],aa_result[,c("mut_name","gene_symbol")])
    final_result_1<-final_result[final_result$mut_name %in% same_MG$mut_name & final_result$gene_symbol %in% same_MG$gene_symbol,]
    bridge_result_1<-aa_result[aa_result$mut_name %in% same_MG$mut_name & aa_result$gene_symbol %in% same_MG$gene_symbol,]
    bb<-cbind(final_result_1[order(final_result_1$mut_name,final_result_1$gene_symbol),],bridge_result_1[order(bridge_result_1$mut_name,bridge_result_1$gene_symbol),3])
    final_result_2<-final_result[!(final_result$mut_name %in% same_MG$mut_name & final_result$gene_symbol %in% same_MG$gene_symbol),]
    bridge_result_2<-aa_result[!(aa_result$mut_name %in% same_MG$mut_name & aa_result$gene_symbol %in% same_MG$gene_symbol),]
    cc<-merge(final_result_2,bridge_result_2,by=c("mut_name","gene_symbol"),all=T)
    colnames(bb)<-colnames(cc)
    final_result<-rbind(bb,cc)
  }
}
final_result[is.na(final_result)]<-"NONE"
final_result$mut_name<-gsub(pattern = ":|>|,",replacement = "_",x=final_result$mut_name)
write.table(final_result,file=paste(GSE_id,"var.var_type.cell_type.exp.p.round_4y_5.txt",sep=""),sep="\t",quote=F,row.names=F)
