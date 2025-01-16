library(tidyverse)
library(DESeq2)
#import data
setwd("D:/")
mycounts<-read.table("mMERGE_without_SP_more_than_equal_to_5.txt",header = TRUE,row.names = 1,sep ="\t" )
# mycounts<-mycounts[,c(1,3,4,5,7,8)]
condition<-factor(c(rep("DN",3),rep("DP",3)),levels = c("DN","DP"))
colData<-data.frame(row.names = colnames(mycounts),condition)

dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
dds <- DESeq(dds)

res= results(dds)
res = res[order(res$pvalue),]
head(res)
summary(res)
write.csv(res,file="All_results.csv")

data <- read.csv("d:/Download/jointpa_matched_features.csv",row.names = 1,header = TRUE)
data <- cbind(rownames(data),data)
result <- data[grepl("cpd", data$matched_features) & grepl("mmu", data$matched_features), ]
df_long <- result %>%
  separate_rows(matched_features, sep = ";") %>%
  mutate(Type = case_when(
    grepl("cpd", matched_features) ~ "cpd",
    grepl("mmu", matched_features) ~ "mmu",
    TRUE ~ "other"
  ))
cpd_elements <- df_long %>% filter(Type == "cpd")
mmu_elements <- df_long %>% filter(Type == "mmu")
colnames(cpd_elements)[1] <- "name"
colnames(mmu_elements)[1] <- "name"

result2 <- cpd_elements %>%
  inner_join(mmu_elements, by = "name") %>%
  select("name","matched_features.x", "matched_features.y")
result2$new <- paste0(result2$matched_features.x,"-",result2$matched_features.y)

metabolites <- read.csv("d:/Download/name_map.csv",header = TRUE)
genes <- read.csv("d:/Download/gene_name_map.csv",header = TRUE)
colnames(metabolites)[1] <- "Query2"
metabolites <- na.omit(metabolites)
genes <- na.omit(genes)
result3 <- result2 %>%
  mutate(Second_Element_metabolites = sapply(strsplit(as.character(matched_features.x), ":"), function(x) x[2])) %>%
  mutate(Second_Element_genes = sapply(strsplit(as.character(matched_features.y), ":"), function(x) x[2])) 
result3$Second_Element_genes <- as.integer(result3$Second_Element_genes)
metabolites <- metabolites[!is.na(metabolites$KEGG),]
result4 <- result3 %>%
  inner_join(metabolites, by = c("Second_Element_metabolites" = "KEGG")) %>%
  mutate(Second_Element_metabolites = Query2) %>%
  inner_join(genes, by = c("Second_Element_genes" = "Entrez")) %>%
  mutate(Second_Element_genes = Query)

library(openxlsx)
RNAseq_counts <- read.table("d:/mMERGE_without_SP_more_than_equal_to_5.txt",sep = "\t",header = TRUE,row.names = 1)
RNAseq_counts <-RNAseq_counts[unique(result4$Second_Element_genes),]
library(DESeq2)
coldata <- data.frame(
  condition = c("DN","DN","DN","DP","DP","DP"),
  row.names = c("DN1","DN2","DN3","DP1","DP2","DP3")
)
dds <- DESeqDataSetFromMatrix(countData = RNAseq_counts, colData = coldata, design = ~1)
dds <- DESeq(dds)
RNAseq_normalized <- counts(dds, normalized = TRUE)
metabolomics_data <- read.xlsx("d:/Differentially Expressed Metabolites.xlsx",sheet = "DP - DN")
rownames(metabolomics_data) <- metabolomics_data$Name
metabolomics_data <- metabolomics_data[unique(result4$Second_Element_metabolites),c(11:16)]
metabolomics_data_log <- log2(metabolomics_data + 1) 
RNAseq_normalized <- as.data.frame(RNAseq_normalized)[,c("DP1","DP2","DP3","DN1","DN2","DN3")]

correlation_matrixAA <- cor(t(RNAseq_normalized), t(metabolomics_data_log), method = "spearman")
cor_and_pvalue <- function(x, y) {
  test <- cor.test(x, y, method = "spearman")
  return(c(correlation = test$estimate, p_value = test$p.value))
}

# 初始化矩阵以存储相关系数和p值
correlation_matrix <- matrix(NA, nrow = ncol(t(RNAseq_normalized)), ncol = ncol(t(metabolomics_data_log)))
pvalue_matrix <- matrix(NA, nrow = ncol(t(RNAseq_normalized)), ncol = ncol(t(metabolomics_data_log)))

# 计算相关性和p值
for (i in 1:ncol(t(RNAseq_normalized))) {
  for (j in 1:ncol(t(metabolomics_data_log))) {
    result <- cor_and_pvalue(as.data.frame(t(RNAseq_normalized))[, i], t(metabolomics_data_log)[, j])
    correlation_matrix[i, j] <- result[1]  # 相关系数
    pvalue_matrix[i, j] <- result[2]      # p值
  }
}

# 转换矩阵为数据框，方便查看
correlation_matrix <- as.data.frame(correlation_matrix)
pvalue_matrix <- as.data.frame(pvalue_matrix)

# 设置行列名
colnames(correlation_matrix) <- colnames(t(metabolomics_data_log))
rownames(correlation_matrix) <- colnames(t(RNAseq_normalized))
colnames(pvalue_matrix) <- colnames(t(metabolomics_data_log))
rownames(pvalue_matrix) <- colnames(t(RNAseq_normalized))

# 初始化一个空的结果数据框来存储变量对、相关系数和p值
correlation_results <- data.frame(
  var1 = character(0),
  var2 = character(0),
  correlation = numeric(0),
  p_value = numeric(0)
)

# 计算所有可能的变量对之间的相关性和p值
for (i in 1:ncol(t(RNAseq_normalized))) {
  for (j in 1:ncol(t(metabolomics_data_log))) {
    result <- cor_and_pvalue(t(RNAseq_normalized)[, i], t(metabolomics_data_log)[, j])
    correlation_results <- rbind(correlation_results, data.frame(
      var1 = colnames(t(RNAseq_normalized))[i],
      var2 = colnames(t(metabolomics_data_log))[j],
      correlation = result[1],
      p_value = result[2]
    ))
  }
}
correlation_results$new_new <- paste0(correlation_results$var1,"-",correlation_results$var2)
result4$new_new <- paste0(result4$Second_Element_genes,"-",result4$Second_Element_metabolites)
result5 <- inner_join(result4,correlation_results,by="new_new")
result6 <- result5 %>% filter(p_value < 0.05)

write.table(result6[,c(18:21)],file = "Network_result.tsv",sep = "\t")
result7 <- data.frame(c(unique(result6$var1),unique(result6$var2)),c(rep("gene",length(unique(result6$var1))),rep("metabolites",length(unique(result6$var2)))))
write.table(result7,file = "TYPE_result.tsv",sep = "\t")

data7 <- read.csv("d:/scorescore.csv")
data7 <- cbind(rownames(data7),data7)
colnames(data7) <- c(colnames(data7)[2:14],"new")
data7 <- data7[,-14]
colnames(result7) <- c("node_name","type")
data7 <- inner_join(data7,result7,by="node_name")
result8 <- data7 %>%
  group_by(type) %>%
  mutate(
    mean_Degree = median(Degree),            # 计算每组 B 列的均值
    sd_Degree = sd(Degree),                # 计算每组 B 列的标准差
    Degree_standardized = (Degree - mean_Degree) / sd_Degree # 标准化
  ) %>%
  ungroup()  # 取消分组
result8 <- result8 %>% arrange(desc(Degree_standardized))
pdf("rank.pdf",width = 6,height = 6)
# result8$rank <- rank(-result8$Degree_standardized)
result8$rank <- rownames(result8)
result8$rank <- as.numeric(result8$rank)
result8 <- result8 %>% arrange(rank)
result8 <- result8[c(1:50),]
library(showtext)

# 启用 showtext
showtext_auto()
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
# result8$color <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
#                    "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
#                    "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", rep("gray",30))  # 前20名颜色，其他为灰色
# result8$color <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
#                                      "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
#                                        "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072")
# 使用 ggplot2 绘制散点图，并自定义字体、去除网格线和设置白色背景
ggplot(result8, aes(x = rank, y = Degree_standardized, color = type)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  
  labs(x = "Rank of Degree", y = "Connecting score", title = "Rank of Connecting score") +
  theme(
    plot.title = element_text(size = 12, family = "Arial", face = "bold"),  # 设置标题字体大小、字体、加粗
    axis.title.x = element_text(size = 10, family = "Arial"),              # 设置 X 轴标题字体大小和字体
    axis.title.y = element_text(size = 10, family = "Arial"),              # 设置 Y 轴标题字体大小和字体
    axis.text = element_text(size = 8, family = "Arial"),                  # 设置坐标轴刻度字体大小和字体
    axis.line = element_line(size = 1),                                     # 设置坐标轴线条粗细
    panel.background = element_rect(fill = "white"),                       # 设置背景为白色
    panel.grid = element_blank(),                                          # 去除网格线
    plot.background = element_rect(fill = "white")                         # 设置图形背景为白色
  )
dev.off()
write.table(result8,file = "CJJ_result8.txt",sep = "\t")
