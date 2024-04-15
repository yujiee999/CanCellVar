####需要的文件：H:/single_cell/GBM/GSE195682/aa_GSE123.Precily.DrugPredictions.IC50.txt
####需要更改IC50.txt aa_GSE123.Precily.DrugPredictions.IC50.txt->GSE123.Precily.DrugPredictions.IC50.txt
####需要更改以下两行
GSE_id <- "GSE150660"
input <- "H:/single_cell/LECA/"

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
dat_drug<-read.table(paste0(GSE_id,".Precily.DrugPredictions.IC50.txt"),header=T,sep="\t",as.is=T,check.names = F)
dat_drug$cell<-rownames(dat_drug)
library(reshape2)
dat_drug_mat<-melt(dat_drug,id.vars=c("cell"),variable.name="drug",value.name="IC50")
dat_var_cell<-unique(data[,c("cell_AF","mut_name","cell")])#最小唯一组合：mut_name	cell
##提取前100位点
mut_freq<-as.data.frame(table(unique(dat_var_cell[,c("mut_name","cell")])$mut_name))
library(dplyr)
r2=mut_freq %>% top_n(n=100,wt=Freq)
dat_var_cell<-dat_var_cell[dat_var_cell$mut_name %in% r2[,1],]

cell_type<-read.table(paste0(GSE_id,"_cell_type.txt"),header=T,sep="\t",as.is=T)
colnames(cell_type)[1]<-"cell"
dat_var_cell_type<-merge(dat_var_cell,cell_type)
#dat_var_cell_type_drug_IC50<-merge(dat_var_cell_type,dat_drug_mat,all.x=T)
#dat_var_cell_type_drug_IC50<-merge(dat_var_cell_type_G,dat_exp_mat)
dat_drug_mat_clitype<-merge(dat_drug_mat,cell_type)

#dat_var_cell_type_drug_IC50$IC50<-round(dat_var_cell_type_drug_IC50$IC50,4)
#dat_var_cell_type_drug_IC50<-dat_var_cell_type_drug_IC50[which(dat_var_cell_type_drug_IC50$cell_AF>0),]#使用编辑|突变水平大于0的数据
dat_var_cell_type$dataset<-GSE_id
#write.table(dat_var_cell_type_drug_IC50[which(dat_var_cell_type_drug_IC50$cell_AF>0),],file=paste(GSE_id,"dat_var_cell_type_drug_IC50_4.txt",sep=""),sep="\t",quote=F,row.names=F)

celltype_unique<-unique(dat_var_cell_type$celltype)
for(ii in 1:length(celltype_unique)){
  drug_result<-dat_var_cell_type[dat_var_cell_type$celltype %in% celltype_unique[ii],]
  mut_drug <- unique(drug_result[,c("mut_name")])
  aa_result<-data.frame(mut_name=character(),drug=character(),value=character())
  for(j in 1:length(mut_drug)){
    analysis_dat <- drug_result[drug_result$mut_name %in% mut_drug[j],]
    analysis_dat <- merge(analysis_dat,dat_drug_mat_clitype)
    if(nrow(analysis_dat)>0){
      drug_detail <- as.character(unique(analysis_dat$drug))
      for(l in 1:length(drug_detail)){
        analysis_dat_1<-analysis_dat[analysis_dat$drug %in% drug_detail[l],]
        clitype_unique<-unique(analysis_dat_1$cli_type)
        cell_mut_type_result<-c()
        cell_cli_type_result<-c()
        cell_drug_IC50_result<-c()
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
        for(m in 1:length(clitype_unique)){
          analysis_dat_2<-analysis_dat_1[analysis_dat_1$cli_type %in% clitype_unique[m],]
          analysis_dat_2$mut_type<-ifelse(analysis_dat_2$cell_AF>0,"Mut","Non_mut")
          x <- analysis_dat_2[analysis_dat_2$mut_type %in% "Mut",]$IC50
          y <- analysis_dat_2[analysis_dat_2$mut_type %in% "Non_mut",]$IC50
          bridge_3<-dat_drug_mat_clitype[(dat_drug_mat_clitype$drug %in% drug_detail[l]) & !(dat_drug_mat_clitype$cell %in% analysis_dat_2$cell) & (dat_drug_mat_clitype$celltype %in% celltype_unique[ii]) & (dat_drug_mat_clitype$cli_type %in% clitype_unique[m]),]
          z <- bridge_3$IC50
          #meanIC50 <- round(mean(c(analysis_dat_1$IC50,bridge_3$IC50)),4)
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
          analysis_dat_cor<-analysis_dat_2[analysis_dat_2$mut_type %in% "Mut",]
          if(nrow(analysis_dat_cor)>=3){
            cor<-round(cor.test(analysis_dat_cor$cell_AF,analysis_dat_cor$IC50)$estimate,4)
            p<-round(cor.test(analysis_dat_cor$cell_AF,analysis_dat_cor$IC50)$p.value,4)
          }else{
            cor<-"NONE"
            p<-"NONE"
          }
          bridge_1<-analysis_dat_2[analysis_dat_2$mut_type %in% "Mut",]
          bridge_2<-analysis_dat_2[analysis_dat_2$mut_type %in% "Non_mut",]
          bridge_3<-bridge_3
          if(nrow(bridge_3)>0){
            bridge_3$mut_type<-"Unknown"
            if(nrow(bridge_1)>0 & nrow(bridge_2)>0){
              cell_mut_type<-get.paste1(c(get.paste(bridge_1$mut_type),get.paste(bridge_2$mut_type),get.paste(bridge_3$mut_type)))
              cell_cli_type<-get.paste1(c(get.paste(bridge_1$cli_type),get.paste(bridge_2$cli_type),get.paste(bridge_3$cli_type)))
              cell_drug_IC50<-get.paste1(c(get.paste(bridge_1$IC50),get.paste(bridge_2$IC50),get.paste(bridge_3$IC50)))
              cell_mut_AF<-get.paste1(c(get.paste(bridge_1$cell_AF),get.paste(bridge_2$cell_AF),"NONE"))
            }
            if(nrow(bridge_1)==0 & nrow(bridge_2)>0){
              cell_mut_type<-get.paste1(c("Mut",get.paste(bridge_2$mut_type),get.paste(bridge_3$mut_type)))
              cell_cli_type<-get.paste1(c(bridge_2$cli_type[1],get.paste(bridge_2$cli_type),get.paste(bridge_3$cli_type)))
              cell_drug_IC50<-get.paste1(c("NONE",get.paste(bridge_2$IC50),get.paste(bridge_3$IC50)))
              cell_mut_AF<-get.paste1(c("NONE",get.paste(bridge_2$cell_AF),"NONE"))
            }
            if(nrow(bridge_1)>0 & nrow(bridge_2)==0){
              cell_mut_type<-get.paste1(c(get.paste(bridge_1$mut_type),"Non_mut",get.paste(bridge_3$mut_type)))
              cell_cli_type<-get.paste1(c(get.paste(bridge_1$cli_type),bridge_1$cli_type[1],get.paste(bridge_3$cli_type)))
              cell_drug_IC50<-get.paste1(c(get.paste(bridge_1$IC50),"NONE",get.paste(bridge_3$IC50)))
              cell_mut_AF<-get.paste1(c(get.paste(bridge_1$cell_AF),"NONE","NONE"))
            }
            if(nrow(bridge_1)==0 & nrow(bridge_2)==0){
              cell_mut_type<-get.paste1(c("Mut","Non_mut",get.paste(bridge_3$mut_type)))
              cell_cli_type<-get.paste1(c(bridge_3$cli_type[1],bridge_3$cli_type[1],get.paste(bridge_3$cli_type)))
              cell_drug_IC50<-get.paste1(c("NONE","NONE",get.paste(bridge_3$IC50)))
              cell_mut_AF<-get.paste1(c("NONE","NONE","NONE"))
            }
          }else{
            if(nrow(bridge_1)>0 & nrow(bridge_2)>0){
              cell_mut_type<-get.paste1(c(get.paste(bridge_1$mut_type),get.paste(bridge_2$mut_type),"Unknown"))
              cell_cli_type<-get.paste1(c(get.paste(bridge_1$cli_type),get.paste(bridge_2$cli_type),bridge_1$cli_type[1]))
              cell_drug_IC50<-get.paste1(c(get.paste(bridge_1$IC50),get.paste(bridge_2$IC50),"NONE"))
              cell_mut_AF<-get.paste1(c(get.paste(bridge_1$cell_AF),get.paste(bridge_2$cell_AF),"NONE"))
            }
            if(nrow(bridge_1)==0 & nrow(bridge_2)>0){
              cell_mut_type<-get.paste1(c("Mut",get.paste(bridge_2$mut_type),,"Unknown"))
              cell_cli_type<-get.paste1(c(bridge_2$cli_type[1],get.paste(bridge_2$cli_type),bridge_2$cli_type[1]))
              cell_drug_IC50<-get.paste1(c("NONE",get.paste(bridge_2$IC50),"NONE"))
              cell_mut_AF<-get.paste1(c("NONE",get.paste(bridge_2$cell_AF),"NONE"))
            }
            if(nrow(bridge_1)>0 & nrow(bridge_2)==0){
              cell_mut_type<-get.paste1(c(get.paste(bridge_1$mut_type),"Non_mut","Unknown"))
              cell_cli_type<-get.paste1(c(get.paste(bridge_1$cli_type),bridge_1$cli_type[1],bridge_1$cli_type[1]))
              cell_drug_IC50<-get.paste1(c(get.paste(bridge_1$IC50),"NONE","NONE"))
              cell_mut_AF<-get.paste1(c(get.paste(bridge_1$cell_AF),"NONE","NONE"))
            }
          }
          cell_mut_type_result<-get.paste2(c(cell_mut_type_result,cell_mut_type))
          cell_cli_type_result<-get.paste2(c(cell_cli_type_result,cell_cli_type))
          cell_drug_IC50_result<-get.paste2(c(cell_drug_IC50_result,cell_drug_IC50))
          cell_mut_AF_result<-get.paste2(c(cell_mut_AF_result,cell_mut_AF))
          xy_logFC_result<-get.paste2(c(xy_logFC_result,xy_logFC))
          xy_p_value_result<-get.paste2(c(xy_p_value_result,xy_p_value))
          xz_logFC_result<-get.paste2(c(xz_logFC_result,xz_logFC))
          xz_p_value_result<-get.paste2(c(xz_p_value_result,xz_p_value))
          yz_logFC_result<-get.paste2(c(yz_logFC_result,yz_logFC))
          yz_p_value_result<-get.paste2(c(yz_p_value_result,yz_p_value))
          cor_result<-get.paste2(c(cor_result,cor))
          p_result<-get.paste2(c(p_result,p))
        }
        aaa_result<-data.frame(value=c(cell_mut_type_result,cell_cli_type_result, cell_drug_IC50_result,cell_mut_AF_result,xy_logFC_result,xy_p_value_result,xz_logFC_result,xz_p_value_result,yz_logFC_result,yz_p_value_result,cor_result,p_result))
        aaa_result$mut_name<-mut_drug[j]
        aaa_result$drug<-drug_detail[l]
        aa_result<-rbind(aa_result,aaa_result)
      }
    }
  }
  aa_result<-aa_result[,c("mut_name","drug","value")]
  colnames(aa_result)[3]<-celltype_unique[ii]
  if(ii ==1){
    final_result<-aa_result
  }else{
    ##相同的mut_name  drug用cbind,否则用merge
    same_MG<-merge(final_result[,c("mut_name","drug")],aa_result[,c("mut_name","drug")])
    final_result_1<-final_result[final_result$mut_name %in% same_MG$mut_name & final_result$drug %in% same_MG$drug,]
    bridge_result_1<-aa_result[aa_result$mut_name %in% same_MG$mut_name & aa_result$drug %in% same_MG$drug,]
    bb<-cbind(final_result_1[order(final_result_1$mut_name,final_result_1$drug),],bridge_result_1[order(bridge_result_1$mut_name,bridge_result_1$drug),3])
    final_result_2<-final_result[!(final_result$mut_name %in% same_MG$mut_name & final_result$drug %in% same_MG$drug),]
    bridge_result_2<-aa_result[!(aa_result$mut_name %in% same_MG$mut_name & aa_result$drug %in% same_MG$drug),]
    cc<-merge(final_result_2,bridge_result_2,by=c("mut_name","drug"),all=T)
    colnames(bb)<-colnames(cc)
    final_result<-rbind(bb,cc)
  }
}
final_result[is.na(final_result)]<-"NONE"
final_result$mut_name<-gsub(pattern = ":|>|,",replacement = "_",x=final_result$mut_name)
final_result$drug_code<-final_result$drug
final_result$drug_code<-gsub(pattern = "-| |\\(|\\)",replacement = "_",x=final_result$drug_code)
write.table(final_result,file=paste(GSE_id,"var.var_type.cell_type.drug.p.round_4y_5.txt",sep=""),sep="\t",quote=F,row.names=F)

