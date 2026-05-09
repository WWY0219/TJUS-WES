setwd("E:/#WWY_Project/TJUS_WES/TJUS_Analysis")
###Environment R packages
library(sigminer)
library(pheatmap)
library(maftools)
library(tidyverse)
library(NMF)
library(RColorBrewer) #color
library(paletteer)    #color
library(TCGAmutations)
library(ggpubr)
library(patchwork)

####Data input and clean####
TJUS <- read.maf(maf = './#TJUS_Data/TJUS.maf',clinicalData = './#TJUS_Data/clinicalall.tsv')
TJUS
##cat data
head(TJUS@data)
head(TJUS@maf.silent)
slotNames(TJUS)
pData <- getClinicalData(TJUS)
cancerGeneList <- read_delim("./#TJUS_data/cancerGeneList.tsv", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)
gene <- cancerGeneList$`Hugo Symbol`
clinmaf <- subsetMaf(TJUS,genes= gene)
clinmaf




####COSMIC Signature####
###TAlly Componenet
##SBS-signature Matrix
mt_tally_SBS <- sig_tally(
  TJUS,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
  mode = "SBS",
  add_trans_bias = F,
#SBS分析中，96 种突变类型会被扩展为 192 种(每种突变类型再细分为 “转录链” 和 “非转录链” 两种情况，以体现链偏倚差异)
  useSyn = TRUE                                #include all variant records in MAF object to generate sample matrix
)
mt_tally_SBS$nmf_matrix[1:5, 1:5]
str(mt_tally_SBS$all_matrices, max.level = 1)

output_SBS96 <- as.data.frame(t(mt_tally_SBS$nmf_matrix))
output_SBS96 <- output_SBS96 %>%
  tibble::rownames_to_column(var = "MutationType")
write.table(output_SBS96, "output_sigminer.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(output_SBS96,"output_SBS96.csv",row.names=F)

##DBS-signature Matrix
mt_tally_DBS <- sig_tally(
  TJUS,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
  mode = "DBS",
  add_trans_bias = F,                       
  useSyn = TRUE
)
output_DBS78 <- as.data.frame(t(mt_tally_DBS$nmf_matrix))
output_DBS78 <- output_DBS78 %>%
  tibble::rownames_to_column(var = "MutationType")
write.table(output_DBS78, "output_DBS78.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(output_DBS78,"output_DBS78.csv",row.names=F)

##ID-signature Matrix
mt_tally_ID <- sig_tally(
  TJUS,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
  useSyn = TRUE,
  mode = "ID",
  add_trans_bias = F
)
str(mt_tally_ID, max.level = 1)


####Extract Signatures####
##Estimate Signature Number
mt_est_SBS <- sig_estimate(mt_tally_SBS$all_matrices$SBS_96,
                       range = 2:6,                #测试的特征数量范围
                       nrun = 100,                  # increase this value if you wana a more stable estimation
                       use_random = FALSE,         # if TRUE, add results from randomized input
                       cores = 4,
                       verbose = TRUE              #是否输出运行过程信息
)
show_sig_number_survey2(mt_est_SBS$survey)
show_sig_number_survey(mt_est_SBS$survey, right_y = NULL)

##Extract Signatures
#Auto(bayesian NMF)
mt_sig2 <- sig_auto_extract(mt_tally_SBS$all_matrices$SBS_96,
                            K0 = 10, 
                            nrun = 100,
                            strategy = "stable")

mt_sig_SBS <- sig_extract(mt_tally_SBS$all_matrices$SBS_96,
                      n_sig = 4,
                      nrun = 100,
                      cores = 4
)
##储存sig文件
write.csv(mt_sig_SBS$Signature,"./TJUS_Sigminer/03.Output/mt_sig_SBS$Signature.csv",row.names =T)
write.csv(mt_sig_SBS$Signature.norm,"./TJUS_Sigminer/03.Output/mt_sig_SBS$Signature-norm.csv",row.names =T)
write.csv(mt_sig_SBS$Exposure,"./TJUS_Sigminer/03.Output/mt_sig_SBS$Exposure.csv",row.names =T)
write.csv(mt_sig_SBS$Exposure.norm,"./TJUS_Sigminer/03.Output/mt_sig_SBS$Exposure-norm.csv",row.names =T)




##Match Signatures
sim_SBS <- get_sig_similarity(mt_sig_SBS,
                          sig_db="SBS_hg19")
str(sim_SBS)
write.csv(sim_SBS$similarity,"./TJUS_Sigminer/03.Output/SBS_hg19_similarity.csv",row.names =T)
write.csv(sim_SBS$best_match,"./TJUS_Sigminer/03.Output/SBS_hg19_bestmatch.csv",row.names =T)
sim_SBS_cosmic <- get_sig_similarity(mt_sig_SBS,
                              sig_db="legacy")
write.csv(sim_SBS_cosmic$similarity,"./TJUS_Sigminer/03.Output/legacy_similarity.csv",row.names =T)
write.csv(sim_SBS_cosmic$best_match,"./TJUS_Sigminer/03.Output/legacy_bestmatch.csv",row.names =T)


pheatmap::pheatmap(sim_SBS$similarity)
sim_similarity <- sim_SBS$similarity
custom_order <- c("Sig1", "Sig2", "Sig3", "Sig4")
sim_similarity <- sim_similarity[custom_order, ]

red_palette <- rev(c("#780522", "#A81428", "#C6403D", "#F5AC8B", "#FAD4BF", "#FBE3D6", "#F8EAE1"))
blue_palette <- c("#134B87", "#256CAE", "#3685BB", "#74B0D2", "#84BDDA", "#C1DDE9", "#EDF2F6")
custom_colors <- c(blue_palette, "white", red_palette)
custom_colors_smooth <- colorRampPalette(c(blue_palette, "white", red_palette))(100)
P1 <- pheatmap(
  sim_similarity,
  cluster_rows = FALSE,  # 关闭行聚类（保持自定义顺序）
  cluster_cols = TRUE,   # 可保留列聚类
  main = "SBS_96 signature",
  color = custom_colors_smooth)
P1
pheatmap::pheatmap(sim_SBS_cosmic$similarity)
sim_similarity_cosmic <- sim_SBS_cosmic$similarity
custom_order <- c("Sig1", "Sig2", "Sig3", "Sig4")
sim_similarity_cosmic <- sim_similarity_cosmic[custom_order, ]

red_palette <- rev(c("#780522", "#A81428", "#C6403D", "#F5AC8B", "#FAD4BF", "#FBE3D6", "#F8EAE1"))
blue_palette <- c("#134B87", "#256CAE", "#3685BB", "#74B0D2", "#84BDDA", "#C1DDE9", "#EDF2F6")
custom_colors <- c(blue_palette, "white", red_palette)
custom_colors_smooth <- colorRampPalette(c(blue_palette, "white", red_palette))(100)
P2 <- pheatmap(
  sim_similarity_cosmic,
  cluster_rows = FALSE,  # 关闭行聚类（保持自定义顺序）
  cluster_cols = TRUE,   # 可保留列聚类
  main = "COSMIC signature",
  color = custom_colors_smooth)
P2



##Operate Signature
sig_SBS <- get_sig_exposure(mt_sig_SBS)




###Exposure Profile
custom_colors <- c("#fc7f71", "#7fb2d3", "#ffb55f", "#8dd3c9")
show_sig_exposure(mt_sig_SBS,
                  palette = custom_colors,
                  hide_samps = FALSE)
##按照exposure分组
grp_exposure <- get_groups(mt_sig_SBS,
                           method = c("exposure"))
grp_exposure_label <- grp_exposure$group
names(grp_exposure_label) <- grp_exposure$sample
grp_exposure
write.csv(grp_exposure,"./TJUS_Sigminer/03.Output/grp_exposure.csv")
P_exposure <- show_sig_exposure(mt_sig_SBS, 
                        style = "cosmic",
                        palette = custom_colors,
                        groups = grp_exposure$enrich_sig,
                        hide_samps = FALSE)
P_exposure




####Signature Fit: Sample Signature Exposure Quantification and Analysis####
###Fit Signatures from reference databases
TJUS_mat <- t(mt_tally_SBS$all_matrices$SBS_96)
##看样本与COSMIC的吻合性
sig_fit_SBS <- sig_fit(TJUS_mat, 
                       sig_index = 1:30, 
                       return_class = "data.table", 
                       rel_threshold = 0.05)
write.csv(sig_fit_SBS,"./TJUS_Sigminer/03.Output/sig_fit_SBS.csv",row.names =F)

##Fit Custom Signatures
sig_Custom_SBS <- sig_fit(TJUS_mat, sig = mt_sig_SBS)
write.csv(sig_Custom_SBS,"./TJUS_Sigminer/03.Output/sig_Custom_SBS.csv",row.names =T)



###Performance Comparison####
#NMF 方法的 RSS 
sum((apply(mt_sig_SBS$Signature, 2, function(x) x / sum(x)) %*% mt_sig_SBS$Exposure - t(mt_tally_SBS$all_matrices$SBS_96))^2)
#sig_fit的RSS
H_estimate <- apply(mt_sig_SBS$Signature, 2, function(x) x / sum(x)) %*% sig_fit(t(mt_tally_SBS$all_matrices$SBS_96), sig = mt_sig_SBS)
H_estimate <- apply(H_estimate, 2, function(x) ifelse(is.nan(x), 0, x))
H_real <- t(mt_tally_SBS$all_matrices$SBS_96)
sum((H_estimate - H_real)^2)

###Estimate Exposure Stability by Bootstrap
bt_result_SBS <- sig_fit_bootstrap_batch(TJUS_mat, 
                                         sig = mt_sig_SBS, 
                                         n = 10)    #run 10times sig_fit
bt_result_SBS
custom_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
#Visulize
show_sig_bootstrap_exposure(bt_result_SBS,
                            palette= custom_colors,
                            plot_fun = c("boxplot"))
ggplot(bt_result_SBS[["expo"]], aes(x = sig, y = exposure, color = sig)) +
  # 箱线图：展示分布的中位数、四分位距等
  geom_boxplot(width = 0.5, alpha = 0.5) +
  # 散点图：添加原始数据点（jitter 使点水平抖动）
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
  # 标题
  ggtitle("All samples") +
  # 轴标签
  xlab("Signature") +
  ylab("Signature exposure") +
  # 主题调整（可选，如去除背景网格）
  theme_bw() +
  theme(panel.grid = element_blank())

show_sig_bootstrap_error(bt_result_SBS)
show_sig_bootstrap_stability(bt_result_SBS)


####Signature Object####
###SBS Signature profile
show_sig_profile(mt_sig_SBS, 
                 mode = "SBS", 
                 paint_axis_text = FALSE, 
                 x_label_angle = 90)

P1<- show_sig_profile(mt_sig_SBS, 
                 mode = "SBS", 
                 style = "cosmic", 
                 x_label_angle = 90)

P1
add_labels(P1, x = 0.74, y = 0.25, y_end = 0.9, labels = sim_SBS, n_label = 4)
add_labels(P1,x = 0.74, y = 0.25, y_end = 0.9, labels = sim_SBS_cosmic, n_label = 4)


###COSMIC Signature Profile
show_cosmic_sig_profile(sig_db = "legacy")
show_cosmic_sig_profile(sig_db = "SBS")
show_cosmic_sig_profile(sig_index = c(1:30), style = "cosmic")



##按照consensus分组
grp_consensus <- get_groups(mt_sig_SBS,
                  method = c("consensus"))
grp_consensus_label <- grp_consensus$group
names(grp_consensus_label) <- grp_consensus$sample
grp_consensus
write.csv(grp_consensus,"./03.Output/grp_consensus.csv")
P3 <- show_sig_exposure(mt_sig_SBS, 
                  style = "cosmic",
                  palette = custom_colors,
                  groups = grp_label,
                  hide_samps = FALSE)
P3

show_sig_consensusmap(mt_sig_SBS,
                      annColors = list(cluster = custom_colors))

show_catalogue(t(mt_tally_SBS$all_matrices$SBS_96), style = "cosmic", x_label_angle = 90)

####Group Analysis####
groups_kmeans <- get_groups(mt_sig_SBS, method = "k-means")
##处理excel文件


groups.cmp <- get_group_comparison(grp_exposure[, -1],
                                   col_group = "group",
                                   cols_to_compare = c("Ploidy", "TMB"),
                                   type = c("co", "co"),                  #
                                   verbose = TRUE
)
ggcomp <- show_group_comparison(groups.cmp)
ggcomp$co_comb

###按照consensus分类进行命名
G1 <- grp_exposure[grp_exposure$group == 1,"sample"]
G2 <- grp_exposure[grp_exposure$group == 2,"sample"]
G3 <- grp_exposure[grp_exposure$group == 3,"sample"]
G4 <- grp_exposure[grp_exposure$group == 4,"sample"]
G1 <- as.character(G1$sample)
G2 <- as.character(G2$sample)
G3 <- as.character(G3$sample)
G4 <- as.character(G4$sample)
G1maf <- subsetMaf(clinmaf,tsb=G1)
G2maf <- subsetMaf(clinmaf,tsb=G2)
G3maf <- subsetMaf(clinmaf,tsb=G3)
G4maf <- subsetMaf(clinmaf,tsb=G4)
vc_cols = c(Missense_Mutation = "#D73925", 
                   Nonsense_Mutation = "#FF8E10", 
                   Multi_Hit = "#66ABAF",
                   Frame_Shift_Del="#FEC739",
                   In_Frame_Del="#B47A93",
                   Frame_Shift_Ins="#FFFFCC",
                   In_Frame_Ins='#3D71A3') 
oncoplot(G1maf,draw_titv = TRUE,colors=vc_cols)
oncoplot(G2maf)
oncoplot(G3maf)
oncoplot(G4maf)
