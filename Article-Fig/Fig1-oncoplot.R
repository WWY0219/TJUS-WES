# ============================================== Loading Environment==========================================================
setwd("workspace")
library(tidyverse)
library(ggplot2)
library(maftools)
library(dplyr)
library(ComplexHeatmap)
library(reshape2)
library(paletteer)       ##连续型配色
library(RColorBrewer)
library(circlize)
library(ggsci)
rm(list=ls())

#============================================== DataReading==================================================================
TJUS <- read.maf(maf = './#TJUS_data/TJUS.maf',clinicalData = './#TJUS_data/clinicalall.tsv')
write.table(clinicaldata_TJUS, file='clinicalall.tsv', 
            quote=F, 
            row.names=F,        #不输出索引行 
            sep='\t')           #写出为maftools识别tsv文件

##-------------------- Complexmap_DATAinput -----------------------------
## oncoplot函数设置参数writeMatrix = TRUE会自动生成文件名为“onco_matrix.txt”的突变矩阵文件
oncoplot(maf = TJUS ,top= 40, writeMatrix = TRUE)                                             
matMut <- read.table("onco_matrix.txt", header = T, check.names = F, sep = "\t")     #读取突变文件
matMut$US058 <- ""
#matMut[matMut == "In-frame"] = "In_frame"#                                         
matMuttmp = matMut                         
matMuttmp$gene = row.names(matMuttmp)                                                #在matMuttmp增加gene
mat_long <- melt(matMuttmp, id.vars = "gene", value.name = "Variant_Classification") # mat_long数据框将包含三列：gene（标识符变量）Original_Column（原始数据框中非id.vars列的列名）以及Variant_Classification（这些列中的实际值）
levels(factor(mat_long$Variant_Classification))
##整理临床数据
pdata <- getClinicalData(TJUS)                                                       #通过maftools获取临床数据
pdata <- subset(pdata, pdata$Tumor_Sample_Barcode %in% colnames(matMut))     
pdata = as.data.frame(pdata)
##更改数据逻辑类型
pdata$Pathology      = factor(pdata$Pathology)
pdata$FIGO           = factor(pdata$FIGO)
pdata$Recurrent      = factor(pdata$Recurrent)
pdata$Radiotherapy   = factor(pdata$Radiotherapy)
pdata$Chemotherapy   = factor(pdata$Chemotherapy)
pdata$Prognosis      = factor(pdata$Prognosis)
pdata$Targetedtherapy=factor(pdata$Targetedtherapy)
pdata$TMB            =as.numeric(pdata$TMB)
pdata$Livetime   = as.numeric(pdata$Livetime)
pdata$Age        = as.numeric(pdata$Age)
pdata$MSI        = as.numeric(pdata$MSI)
pdata$Purity     = as.numeric(pdata$Purity)
pdata$Ploidy     = as.numeric(pdata$Ploidy)
pdata$Tumor_Size = as.numeric(pdata$Tumor_Size)
str(pdata)
matMut <- matMut[, pdata$Tumor_Sample_Barcode]

##确定图形参数##
alter_fun <- list(
  background = function(x, y, w, h) {                                 
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),             #grid.rect用来绘制矩形【x,y确定矩形位置；w,h确定矩形高度宽度
              gp = gpar(fill = "#dcddde", col = NA))                    #gp为图形参数，col=NA不设置矩形边框
  },
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Missense_Mutation"], col = NA))
  },
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),  
              gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
  },
  Multi_Hit = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Multi_Hit"], col = NA))
  },
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Frame_Shift_Del"], col = NA))
  },
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["In_Frame_Del"], col = NA))
  },
  In_Frame_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["In_Frame_Ins"], col = NA))
  }
  #  Splice_Site = function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(0.5, "mm"),h-unit(0.5, "mm"),
  #            gp = gpar(fill = col["Splice_Site"], col = NA))
  #}
)
heatmap_legend_param <- list(title = "Alternations", 
                             at = c("Missense_Mutation", "Nonsense_Mutation" , "Multi_Hit","Frame_Shift_Del","In_Frame_Del","Frame_Shift_Ins"), 
                             labels = c("Missense_Mutation", "Nonsense_Mutation" , "Multi_Hit","Frame_Shift_Del","In_Frame_Del","Frame_Shift_Ins"))   
###配色环节
display.brewer.all()   ##查看所有配色
brewer.pal(9,"YlOrRd")

col <- c(Missense_Mutation = "#ED7B61", 
         Nonsense_Mutation = "#FF8E10", 
         Multi_Hit = "#5679BA",
         Frame_Shift_Del="#FEC739",
         In_Frame_Del="#B47A93",
         Frame_Shift_Ins="#FFFFCC",
         In_Frame_Ins='#A89C87')        
#分类变量
color_group        = c(P = "#006093",R = "#CF1C29")

color_Pathology    = c(`ULMS`='#94B5D7',`HGESS` ='#9897BC',`LGESS` ='#C8B8D4',`AS` ='#D3E6EF')

color_FIGO         = c(I ='#29A15C', `II` ='#F49600', `III` ='#428DBF', `IV`  ='#BE0E23')


color_Radiotherapy = c(Yes = "#B8D4E6", No  = '#dcddde')
color_Chemotherapy = c(Yes = "#216FB0", No  = '#dcddde')
color_Targetedtherapy     = c(Yes = "#08519C", No  = '#dcddde')
color_status       = c(death = '#D3D3D3',live =  '#A5D6A7')
col_OS       = colorRamp2(c(0, 100), c("white", "red3"))
col_age      = colorRamp2(c(20,80),c("white","#CC4763"))
col_size     = colorRamp2(c(0,40),c("white","skyblue"))
color_MSI    = colorRamp2(c(0,20),c("white","#A16179"))
color_Purity = colorRamp2(c(0,1),c("white","purple3"))
color_ploidy = colorRamp2(c(0,5.5),c("white","#B80EC4"))
color_TMB    = colorRamp2(c(0,152),c("white","#9D4097"))

##底部注释
ha_bottom <- HeatmapAnnotation(Pathology = pdata$Pathology,
                               FIGO=pdata$FIGO,
                               Age = pdata$Age,
                               Tumor_size = pdata$Tumor_Size,
                               Radiotherapy=pdata$Radiotherapy,
                               Chemotherapy=pdata$Chemotherapy,
                               Immutherapy=pdata$Targetedtherapy,
                               Status=pdata$Prognosis,
                               OS = pdata$Livetime,
                               col = list(OS = col_OS,
                                          Age=col_age,
                                          Tumor_size=col_size,
                                          Pathology = color_Pathology,
                                          Radiotherapy=color_Radiotherapy,
                                          Chemotherapy=color_Chemotherapy,
                                          Immutherapy=color_Targetedtherapy,
                                          Status=color_status,
                                          FIGO=color_FIGO), 
                               show_annotation_name = TRUE,
                               annotation_name_side = "left",
                               annotation_name_gp = gpar(fontsize = 10)
)

##顶端注释
ha_top<- HeatmapAnnotation(Group = pdata$SeqSite,
                           TMB = pdata$TMB,
                           MSI_Score = pdata$MSI_Score,
                           Purity = pdata$Purity,
                           Ploidy= pdata$Ploidy,
                           col = list(MSI_Score=color_MSI,Purity=color_Purity,Ploidy=color_ploidy,Group = color_group,TMB=color_TMB), 
                           show_annotation_name = TRUE,
                           annotation_name_side = "left",
                           annotation_name_gp = gpar(fontsize = 10)
)
column_title <- "TJUS_ONCOPRINT "
##根据原发及复发进行分组
grouplist <- pdata[,1:2]
grouplist <- grouplist[order(grouplist$SeqSite), ]
sample_order <- grouplist %>%
  arrange(SeqSite) %>%  # 按 sequencing_site 升序排序（P 在前，R 在后）
  pull(Tumor_Sample_Barcode) 

##绘图
P_all <- oncoPrint(matMut,
          bottom_annotation = ha_bottom,
          top_annotation = ha_top, #注释信息在底部
          alter_fun = alter_fun, 
          col = col, 
          show_pct = T,                                                        #展示基因突变频率
          show_heatmap_legend = T,
          column_title = column_title, 
          heatmap_legend_param = heatmap_legend_param,
          row_names_side = "left",
          pct_side = "right",
          alter_fun_is_vectorized = FALSE,
          column_order = sample_order,
          column_split = factor(pdata$SeqSite, levels = c("P", "R"))
)
P_all
##save result
pdf("../03.Output/filtered_oncoprint_all.pdf", width = 10.33, height = 8.16)
P_all
dev.off()
png("../03.Output/filtered_oncoprint_all.png", width = 992, height = 783)
P_all
dev.off()


####复发图谱绘制####
gene <- cancerGeneList$`Hugo Symbol`
clinmaf <- subsetMaf(TJUS,genes= gene)
clin.P <- subset(pdata, sequencing_site=="P")$Tumor_Sample_Barcode
clin.R <- subset(pdata, sequencing_site=="R")$Tumor_Sample_Barcode
TJUS.P <- subsetMaf(maf=clinmaf, tsb=clin.P)
TJUS.R <- subsetMaf(maf=clinmaf, tsb=clin.R)
##复发图谱
oncoplot(maf = TJUS.R,top= 40, writeMatrix = TRUE)                                             
matMut <- read.table("onco_matrix_R.txt", header = T, check.names = F, sep = "\t")     #读取突变文件
matMuttmp = matMut                         
matMuttmp$gene = row.names(matMuttmp)                                                #在matMuttmp增加gene
mat_long <- melt(matMuttmp, id.vars = "gene", value.name = "Variant_Classification") # mat_long数据框将包含三列：gene（标识符变量）Original_Column（原始数据框中非id.vars列的列名）以及Variant_Classification（这些列中的实际值）
levels(factor(mat_long$Variant_Classification))
##整理临床数据
pdata <- getClinicalData(TJUS.R)                                               #通过maftools获取临床数据
pdata <- subset(pdata, pdata$Tumor_Sample_Barcode %in% colnames(matMut))     # %in%是一个二元操作符，用于测试左边的元素是否出现在右边的向量中
pdata = as.data.frame(pdata)
pdata$Pathol_Subtype = factor(pdata$Pathol_Subtype)
pdata$FIGO           = factor(pdata$FIGO)
pdata$Recurrent_site = factor(pdata$Recurrent_site)
pdata$Radiotherapy   = factor(pdata$Radiotherapy)
pdata$Chemotherapy   = factor(pdata$Chemotherapy)
pdata$Prognosis      = factor(pdata$Prognosis)
pdata$sequencing_site = factor(pdata$sequencing_site)
pdata$Immunity       =factor(pdata$Immunity)
pdata$TMB            =as.numeric(pdata$TMB)
pdata$Livetime  = as.numeric(pdata$Livetime)
pdata$Age       = as.numeric(pdata$Age)
pdata$MSI_Score = as.numeric(pdata$MSI_Score)
pdata$Purity    = as.numeric(pdata$Purity)
pdata$Ploidy    = as.numeric(pdata$Ploidy)
pdata$`Tumor_Size（cm）` = as.numeric(pdata$`Tumor_Size（cm）`)
str(pdata)
matMut <- matMut[, pdata$Tumor_Sample_Barcode]
###以上操作
alter_fun <- list(
  background = function(x, y, w, h) {                                 
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),             #grid.rect用来绘制矩形【x,y确定矩形位置；w,h确定矩形高度宽度
              gp = gpar(fill = "white", col = NA))                    #gp为图形参数，col=NA不设置矩形边框
  },
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Missense_Mutation"], col = NA))
  },
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),  
              gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
  },
  Multi_Hit = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Multi_Hit"], col = NA))
  },
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["In_Frame_Del"], col = NA))
  }
  #  Splice_Site = function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(0.5, "mm"),h-unit(0.5, "mm"),
  #            gp = gpar(fill = col["Splice_Site"], col = NA))
  #}
)
heatmap_legend_param <- list(title = "Alternations", 
                             at = c("Missense_Mutation", "Nonsense_Mutation" , "Multi_Hit","In_Frame_Del"), 
                             labels = c("Missense_Mutation", "Nonsense_Mutation" , "Multi_Hit","In_Frame_Del"))   
###配色环节
library(paletteer)     ##连续型配色
library(RColorBrewer)
library(circlize)
library(ggsci)
display.brewer.all()   ##查看所有配色
brewer.pal(9,"YlOrRd")

col <- c(Missense_Mutation = "#D73925", 
         Nonsense_Mutation = "#FF8E10", 
         Multi_Hit = "#66ABAF",
         Frame_Shift_Del="#FEC739",
         In_Frame_Del="#B47A93",
         Frame_Shift_Ins="#FFFFCC",
         In_Frame_Ins='#3D71A3')        # 指定颜色, 调整颜色代码即可 
#分类变量
color_group        = c(P = "#006093",R = "#CF1C29")
color_Pathology    = c(ULMS   ='#80C9CA',
                       `HG-ESS` ='#9575CD',
                       `LG-ESS` ='#81C784',
                       AS     ='#FFE0B2')
color_FIGO         = c(`IA`   ='#FDC593',
                       `IB`   ='#F57B2D',
                       `II`   ='#C2E7BB',
                       `IIB`  ='#5EB96B',
                       `IIIA` ='#84BADB',
                       `IIIB` ='#3484BE',
                       `IIIC` ='#105BA2',
                       `IVB`  ='#BE93C6')
color_Radiotherapy = c(Yes = "#B8D4E6",
                       No  = '#BFB7AB')
color_Chemotherapy = c(Yes = "#216FB0",
                       No  = '#BFB7AB')
color_Immunity     = c(Yes = "#08519C",
                       No  = '#BFB7AB')
color_status       = c(death = '#D3D3D3',
                       live =  '#A5D6A7')
col_OS       = colorRamp2(c(0, 100), c("white", "red3"))
col_age      = colorRamp2(c(20,80),c("white","#CC4763"))
col_size     = colorRamp2(c(0,40),c("white","skyblue"))
color_MSI    = colorRamp2(c(0,20),c("white","#A16179"))
color_Purity = colorRamp2(c(0,1),c("white","purple3"))
color_ploidy = colorRamp2(c(0,5.5),c("white","#B80EC4"))
color_TMB    = colorRamp2(c(0,152),c("white","#9D4097"))
# 定义注释信息 自定义颜色 连续性变量设置颜色（外）
ha_bottom <- HeatmapAnnotation(Pathology = pdata$Pathol_Subtype,
                               FIGO=pdata$FIGO,
                               Age = pdata$Age,
                               Tumor_size = pdata$`Tumor_Size（cm）`,
                               Radiotherapy=pdata$Radiotherapy,
                               Chemotherapy=pdata$Chemotherapy,
                               Immutherapy=pdata$Immunity,
                               Status=pdata$Prognosis,
                               OS = pdata$Livetime,
                               col = list(OS = col_OS,
                                          Age=col_age,
                                          Tumor_size=col_size,
                                          Patholoy=color_Pathology,
                                          Radiotherapy=color_Radiotherapy,
                                          Chemotherapy=color_Chemotherapy,
                                          Immutherapy=color_Immunity,
                                          Status=color_status,
                                          FIGO=color_FIGO), 
                               show_annotation_name = TRUE,
                               annotation_name_gp = gpar(fontsize = 10)
)
ha_top<- HeatmapAnnotation(Group = pdata$sequencing_site,
                           TMB = pdata$TMB,
                           MSI_Score = pdata$MSI_Score,
                           Purity = pdata$Purity,
                           Ploidy= pdata$Ploidy,
                           col = list(MSI_Score=color_MSI,Purity=color_Purity,Ploidy=color_ploidy,Group = color_group,TMB=color_TMB), 
                           show_annotation_name = TRUE,
                           annotation_name_gp = gpar(fontsize = 10)
)
column_title <- "TJUS_ONCOPRINT "
P2 <- oncoPrint(matMut,
                bottom_annotation = ha_bottom,
                top_annotation = ha_top, #注释信息在底部
                #   top_annotation=top_annotation,
                #right_annotation=NULL,
                alter_fun = alter_fun, 
                col = col,  
                column_title = column_title, 
                heatmap_legend_param = heatmap_legend_param,
                row_names_side = "left",
                pct_side = "right",
                # column_order=sample_order,
                #       column_split=3
                alter_fun_is_vectorized = FALSE
)
P2
##save result
pdf("../03.Output/filtered_oncoprint_R.pdf", width = 10.33, height = 8.16)
P2
dev.off()
png("../03.Output/filtered_oncoprint_R.png", width = 992, height = 783)
P2
dev.off()
