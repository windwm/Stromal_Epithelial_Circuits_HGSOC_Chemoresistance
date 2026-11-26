#######重新画点图
##
library(Seurat)
library(circlize)
library(viridis)
library(dplyr)
library(tidyr)

setwd("F:\\20261102zxx修图\\renjy\\renjy\\")

# ==========================================
# 1. 读取 CellPhoneDB 数据
# ==========================================

mypvals <- read.delim(
  "pvalues.txt",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

mymeans <- read.delim(
  "means.txt",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ==========================================
# 2. 提取 CAF 与 EC-Mix2 的六种细胞组合
# ==========================================

caf <- c("CAF-1", "CAF-2", "CAF-3")

cellcols <- c(
  paste0(caf, "|EC-Mix2"),
  paste0("EC-Mix2|", caf)
)

stopifnot(all(cellcols %in% names(mypvals)))
stopifnot(all(cellcols %in% names(mymeans)))

# ==========================================
# 3. 转换为长格式
# ==========================================

means_long <- mymeans %>%
  select(
    id_cp_interaction,
    interacting_pair,
    gene_a,
    gene_b,
    partner_a,
    partner_b,
    receptor_a,
    receptor_b,
    all_of(cellcols)
  ) %>%
  pivot_longer(
    cols = all_of(cellcols),
    names_to = "cell_pair",
    values_to = "means"
  )

pvals_long <- mypvals %>%
  select(
    id_cp_interaction,
    all_of(cellcols)
  ) %>%
  pivot_longer(
    cols = all_of(cellcols),
    names_to = "cell_pair",
    values_to = "pvals"
  )

# ==========================================
# 4. 合并 mean 和 pvalue
# 注意：此处不提前筛选 pvals < 0.01
# ==========================================

pldf <- means_long %>%
  inner_join(
    pvals_long,
    by = c("id_cp_interaction", "cell_pair")
  ) %>%
  separate(
    cell_pair,
    into = c("cell_a", "cell_b"),
    sep = "\\|",
    remove = FALSE
  ) %>%
  mutate(
    means = as.numeric(means),
    pvals = as.numeric(pvals)
  )

# ==========================================
# 5. 拆分 interacting_pair
# ==========================================

pldf <- pldf %>%
  mutate(
    receptor_a = as.logical(receptor_a),
    receptor_b = as.logical(receptor_b),
    
    n_underscore = str_count(
      interacting_pair, "_"
    ),
    
    # 1个下划线：在第1个下划线处分割
    # 2个下划线：在第2个下划线处分割
    member_a = case_when(
      n_underscore == 1 ~
        sub("_.*$", "", interacting_pair),
      
      n_underscore == 2 ~
        sub("_[^_]+$", "", interacting_pair),
      
      TRUE ~ NA_character_
    ),
    
    member_b = case_when(
      n_underscore == 1 ~
        sub("^[^_]*_", "", interacting_pair),
      
      n_underscore == 2 ~
        sub("^.*_", "", interacting_pair),
      
      TRUE ~ NA_character_
    )
  )

# 检查无法自动拆分的相互作用
unresolved_pairs <- pldf %>%
  filter(
    is.na(member_a) |
      is.na(member_b)
  ) %>%
  distinct(
    interacting_pair,
    n_underscore
  )

print(unresolved_pairs)

# 如果存在无法拆分的相互作用，停止运行
if (nrow(unresolved_pairs) > 0) {
  stop("存在无法自动拆分的 interacting_pair，请先检查")
}


# ==========================================
# 6. 判断 Sender、Receiver、Ligand、Receptor
# ==========================================

pldf <- pldf %>%
  mutate(
    # 一真一假才自动判断配受体方向
    valid_LR = !is.na(receptor_a) &
      !is.na(receptor_b) &
      xor(receptor_a, receptor_b),
    
    # Sender
    Sender = if_else(
      valid_LR & receptor_a,
      cell_b,
      cell_a
    ),
    
    # Receiver
    Receiver = if_else(
      valid_LR & receptor_a,
      cell_a,
      cell_b
    ),
    
    # Ligand
    Ligand = if_else(
      valid_LR & receptor_a,
      member_b,
      member_a
    ),
    
    # Receptor
    Receptor = if_else(
      valid_LR & receptor_a,
      member_a,
      member_b
    ),
    
    # 统一 Ligand_Receptor 名称
    LR_pair = paste0(
      Ligand, "_", Receptor
    ),
    
    direction_status = if_else(
      valid_LR,
      "Based_on_receptor_annotation",
      "Original_order_unverified"
    )
  )

# ==========================================
# 7. 人工校正 FGFR1_NCAM1
# ==========================================

pldf <- pldf %>%
  mutate(
    Sender = if_else(
      interacting_pair == "FGFR1_NCAM1",
      cell_b, Sender
    ),
    
    Receiver = if_else(
      interacting_pair == "FGFR1_NCAM1",
      cell_a, Receiver
    ),
    
    Ligand = if_else(
      interacting_pair == "FGFR1_NCAM1",
      "NCAM1", Ligand
    ),
    
    Receptor = if_else(
      interacting_pair == "FGFR1_NCAM1",
      "FGFR1", Receptor
    ),
    
    LR_pair = if_else(
      interacting_pair == "FGFR1_NCAM1",
      "NCAM1_FGFR1", LR_pair
    ),
    
    direction_status = if_else(
      interacting_pair == "FGFR1_NCAM1",
      "Manually_curated",
      direction_status
    )
  )

# ==========================================
# 8. 按真实 Sender/Receiver 方向分组
# ==========================================

pldf <- pldf %>%
  mutate(
    direction = case_when(
      Sender %in% caf &
        Receiver == "EC-Mix2" ~ "CAF_to_EC",
      
      Sender == "EC-Mix2" &
        Receiver %in% caf ~ "EC_to_CAF",
      
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(direction))


# ==========================================
# 9. 两个方向分别进行显著性筛选
# ==========================================

# 每个配受体对、每个方向：
# 至少一个细胞组合同时满足
# pvals < 0.01 且 means > 0

pldf <- pldf %>%
  group_by(
    id_cp_interaction,
    direction
  ) %>%
  mutate(
    logi = sum(
      !is.na(pvals) &
        pvals < 0.01 &
        !is.na(means) &
        means > 0
    )
  ) %>%
  filter(logi >= 1) %>%
  ungroup() %>%
  filter(
    !is.na(means),
    means > 0
  )

# ==========================================
# 10. 分别提取两个方向
# ==========================================

CAF_to_EC <- pldf %>%
  filter(direction == "CAF_to_EC") %>%
  mutate(
    CellPair = Sender
  )

EC_to_CAF <- pldf %>%
  filter(direction == "EC_to_CAF") %>%
  mutate(
    CellPair = Receiver
  )

# 固定横轴顺序
CAF_to_EC$CellPair <- factor(
  CAF_to_EC$CellPair,
  levels = caf
)

EC_to_CAF$CellPair <- factor(
  EC_to_CAF$CellPair,
  levels = caf
)


# ==========================================
# 12. 绘制 CAF -> EC-Mix2
# ==========================================

p1 <- ggplot(
  CAF_to_EC,
  aes(x = CellPair, y = LR_pair)
) +
  geom_point(
    aes(
      color = means,
      size = -log10(pvals + 0.0001)
    )
  ) +
  scale_size_continuous(
    range = c(1, 3)
  ) +
  scale_color_gradient2(
    high = "red",
    mid = "yellow",
    low = "darkblue",
    midpoint = 1
  ) +
  labs(
    title = "CAF to EC-Mix2",
    x = NULL,
    y = "Ligand_Receptor",
    color = "Means",
    size = "-log10(pvals+1e-04)"
  ) +
  theme_bw()+ theme(
    axis.text.x = element_text(
      colour = "black"
    ),
    axis.text.y = element_text(
      colour = "black"
    )
  )

print(p1)

# ==========================================
# 13. 绘制 EC-Mix2 -> CAF
# ==========================================

p2 <- ggplot(
  EC_to_CAF,
  aes(x = CellPair, y = LR_pair)
) +
  geom_point(
    aes(
      color = means,
      size = -log10(pvals + 0.0001)
    )
  ) +
  scale_size_continuous(
    range = c(1, 3)
  ) +
  scale_color_gradient2(
    high = "red",
    mid = "yellow",
    low = "darkblue",
    midpoint = 0.6
  ) +
  labs(
    title = "EC-Mix2 to CAF",
    x = NULL,
    y = "Ligand_Receptor",
    color = "Means",
    size = "-log10(pvals+1e-04)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      colour = "black"
    ),
    axis.text.y = element_text(
      colour = "black"
    )
  )

print(p2)

# ==========================================
# 14. 导出数据及图片
# ==========================================

write.csv(
  pldf,
  "CAF_ECMix2_filtered_all.csv",
  row.names = FALSE
)

write.csv(
  CAF_to_EC,
  "CAF_to_ECMix2.csv",
  row.names = FALSE
)

write.csv(
  EC_to_CAF,
  "ECMix2_to_CAF.csv",
  row.names = FALSE
)

ggsave(
  "CAF_to_ECMix2.pdf",
  p1,
  width = 7,
  height = 8
)

ggsave(
  "ECMix2_to_CAF.pdf",
  p2,
  width = 7,
  height = 8
)

####弦图
library(dplyr)
library(tidyr)
library(stringr)
library(circlize)

setwd("F:\\20261102zxx修图\\renjy\\renjy\\")

# ==========================================
# 1. 读取 CellPhoneDB 全部结果
# ==========================================

mymeans <- read.delim(
  "means.txt",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

mypvals <- read.delim(
  "pvalues.txt",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ==========================================
# 2. 提取所有细胞组合
# ==========================================

# CellPhoneDB 细胞组合列通常包含 |
cell_pairs <- intersect(
  grep("\\|", colnames(mymeans), value = TRUE),
  grep("\\|", colnames(mypvals), value = TRUE)
)

if (length(cell_pairs) == 0) {
  stop("没有找到细胞组合列，请检查文件")
}

# 提取全部细胞类型
all_cells <- unique(unlist(
  strsplit(cell_pairs, "|", fixed = TRUE)
))

print(all_cells)

cat("细胞类型数量：", length(all_cells), "\n")
cat("细胞组合数量：", length(cell_pairs), "\n")

# ==========================================
# 3. 转换成长格式
# ==========================================

annotation_cols <- c(
  "id_cp_interaction",
  "interacting_pair",
  "partner_a",
  "partner_b",
  "gene_a",
  "gene_b",
  "receptor_a",
  "receptor_b"
)

annotation_cols <- intersect(
  annotation_cols,
  colnames(mymeans)
)

means_long <- mymeans %>%
  select(
    all_of(annotation_cols),
    all_of(cell_pairs)
  ) %>%
  pivot_longer(
    cols = all_of(cell_pairs),
    names_to = "cell_pair",
    values_to = "means"
  )

pvals_long <- mypvals %>%
  select(
    id_cp_interaction,
    all_of(cell_pairs)
  ) %>%
  pivot_longer(
    cols = all_of(cell_pairs),
    names_to = "cell_pair",
    values_to = "pvals"
  )

# ==========================================
# 4. 合并 mean 和 pvalue
# ==========================================

pldf_all <- means_long %>%
  inner_join(
    pvals_long,
    by = c("id_cp_interaction", "cell_pair")
  ) %>%
  separate(
    cell_pair,
    into = c("cell_a", "cell_b"),
    sep = "\\|",
    remove = FALSE
  ) %>%
  mutate(
    means = as.numeric(means),
    pvals = as.numeric(pvals),
    receptor_a = as.logical(receptor_a),
    receptor_b = as.logical(receptor_b)
  )

# ==========================================
# 1. 筛选显著互作
# ==========================================

pldf_sig <- pldf_all %>%
  filter(
    !is.na(pvals),
    pvals < 0.01,
    !is.na(means),
    means > 0
  )

# ==========================================
# 2. 统一细胞组合顺序
# ==========================================

# 例如：
# CAF-1 | EC-Mix2
# EC-Mix2 | CAF-1
# 都统一为 CAF-1 | EC-Mix2

pldf_sig1 <- pldf_sig %>%
  mutate(
    Cell1 = pmin(cell_a, cell_b),
    Cell2 = pmax(cell_a, cell_b)
  )

# ==========================================
# 3. 统计每两个细胞类型之间的互作对数
# ==========================================

interaction_count <- pldf_sig1 %>%
  group_by(Cell1, Cell2) %>%
  summarise(
    n_interactions = n_distinct(id_cp_interaction),
    .groups = "drop"
  ) %>%
  arrange(Cell1, Cell2)

print(interaction_count)

write.csv(
  interaction_count,
  "All_cell_interaction_count_undirected.csv",
  row.names = FALSE
)


# ==========================================
# 2. 定义无方向互作统计函数
# ==========================================

count_interactions <- function(data, cells) {
  
  data %>%
    filter(
      cell_a %in% cells,
      cell_b %in% cells
    ) %>%
    mutate(
      Cell1 = pmin(cell_a, cell_b),
      Cell2 = pmax(cell_a, cell_b)
    ) %>%
    group_by(Cell1, Cell2) %>%
    summarise(
      n_interactions = n_distinct(id_cp_interaction),
      .groups = "drop"
    ) %>%
    arrange(Cell1, Cell2)
}


# ==========================================
# 3. 图1：指定7种细胞
# ==========================================

cells_fig1 <- c(
  "CAF-1",
  "CAF-2",
  "CAF-3",
  "EC-Mix1",
  "EC-Mix2",
  "EC-Mix3",
  "EC-Mix4"
)

# 统计互作对数
interaction_count1 <- count_interactions(
  pldf_sig,
  cells_fig1
)

print(interaction_count1)

write.csv(
  interaction_count1,
  "CAF_ECMix_interaction_count.csv",
  row.names = FALSE
)

# ==========================================
# 5. 定义需要合并的细胞类型
# ==========================================

ec_cells <- c(
  "EC-Mix1",
  "EC-Mix2",
  "EC-Mix3",
  "EC-Mix4"
)

t_cells <- c(
  "TC_CD4_FOXP3",
  "TC_CD4_CXCL13",
  "TC_CD4_S100A6",
  "TC_CD8_IL7R",
  "TC_CD8_IL4R"
)

b_cells <- c(
  "BC_CD38",
  "BC_MS4A1"
)

macrophage_cells <- c(
  "M_MRC1",
  "M_SPP1",
  "M_HLA-DRB5",
  "M_MKI67",
  "M_IGKC"
)

# ==========================================
# 6. 定义细胞合并函数
# ==========================================

merge_celltype <- function(x) {
  
  case_when(
    x %in% ec_cells ~ "EC-mixs",
    
    x %in% macrophage_cells ~ "Macrophage",
    
    x %in% t_cells ~ "T",
    
    x %in% b_cells ~ "B",
    
    TRUE ~ x
  )
}

# ==========================================
# 7. 对原始数据进行细胞类型合并
# ==========================================

pldf_merged <- pldf_sig %>%
  mutate(
    cell_a = merge_celltype(cell_a),
    cell_b = merge_celltype(cell_b)
  )

# ==========================================
# 8. 指定图2细胞类型
# ==========================================

cells_fig2 <- c(
  "CAF-1",
  "CAF-2",
  "CAF-3",
  "Endothelial cell",
  "EC-mixs",
  "Mesenchymal-fib",
  "T",
  "NK_GNLY",
  "B",
  "Macrophage",
  "DC",
  "Mesothelial cell",
  "Theca cell"
)

# ==========================================
# 9. 统计合并后的显著互作数量
# ==========================================

interaction_count2 <- count_interactions(
  pldf_merged,
  cells_fig2
)

print(interaction_count2)

write.csv(
  interaction_count2,
  "Merged_cell_interaction_count.csv",
  row.names = FALSE
)

##不要CAF-1，CAF-2,CAF-3之间的互作，只保留CAF和其他细胞群之间的互作

library(dplyr)
library(circlize)

# 定义CAF细胞类型
caf_cells <- c("CAF-1", "CAF-2", "CAF-3")

# 只保留CAF与非CAF之间的互作
filter_CAF_interactions <- function(data) {
  
  data %>%
    filter(
      (Cell1 %in% caf_cells & !Cell2 %in% caf_cells) |
        (Cell2 %in% caf_cells & !Cell1 %in% caf_cells)
    )
}

# 图1：CAF与EC-Mix1~4
interaction_count1_CAF <- filter_CAF_interactions(
  interaction_count1
)

# 图2：CAF与其他合并后的细胞类型
interaction_count2_CAF <- filter_CAF_interactions(
  interaction_count2
)

print(interaction_count1_CAF)
print(interaction_count2_CAF)

# 保存
write.csv(
  interaction_count1_CAF,
  "CAF_ECMix_interaction_count_CAFonly.csv",
  row.names = FALSE
)

write.csv(
  interaction_count2_CAF,
  "CAF_other_cells_interaction_count_CAFonly.csv",
  row.names = FALSE
)





# ==========================================
# 4. 图1：绘制弦图
# ==========================================

colors_fig1 <- c(
  "CAF-1"   = "#519259",
  "CAF-2"   = "#3E7C17",
  "CAF-3"   = "#064635",
  "EC-Mix1" = "#651fff",
  "EC-Mix2" = "#b39ddb",
  "EC-Mix3" = "#9575cd",
  "EC-Mix4" = "#7e57c2"
)


# ==========================================
# 10. 设置图2细胞颜色
# ==========================================

colors_fig2 <- c(
  "CAF-1"            = "#519259",
  "CAF-2"            = "#3E7C17",
  "CAF-3"            = "#064635",
  "Endothelial cell" = "#4DBBD5",
  "EC-mixs"          = "#00A087",
  "Mesenchymal-fib"  = "#7E6148",
  "T"                = "#3C5488",
  "NK_GNLY"               = "#8491B4",
  "B"                = "#91D1C2",
  "Macrophage"       = "#DC0000",
  "DC"               = "#F1C40F",
  "Mesothelial cell"      = "#7A68A6",
  "Theca cell"            = "#8FBC8F"
)


###图1 CAF与EC-mix
pdf("CAF_ECMix_chord_CAFonly.pdf", width = 11, height = 11)
circos.clear()

circos.par(
  start.degree = 90,
  gap.degree = 5,
  cell.padding = c(0, 0, 0, 0)
)

chordDiagram(
  interaction_count1_CAF %>%
    select(Cell1, Cell2, n_interactions),
  
  order = cells_fig1,
  grid.col = colors_fig1,
  col = colors_fig1[interaction_count1_CAF$Cell1],
  transparency = 0.45,
  directional = 0,
  
  # 增加数量刻度
  annotationTrack = c("name", "grid", "axis"),
  
  annotationTrackHeight = c(
    mm_h(5),
    mm_h(3),
    mm_h(5)
  )
)

circos.clear()
dev.off()




###图2 CAF与所有celltype
pdf("CAF_other_cells_chord_CAFonly.pdf", width = 11, height = 11)

circos.clear()

circos.par(
  start.degree = 90,
  gap.degree = 3,
  cell.padding = c(0, 0, 0, 0)
)

chordDiagram(
  interaction_count2_CAF %>%
    select(Cell1, Cell2, n_interactions),
  
  order = cells_fig2,
  grid.col = colors_fig2,
  col = colors_fig2[interaction_count2_CAF$Cell1],
  transparency = 0.55,
  directional = 0,
  
  # 增加数量刻度
  annotationTrack = c("name", "grid", "axis"),
  
  annotationTrackHeight = c(
    mm_h(5),
    mm_h(3),
    mm_h(5)
  )
)

circos.clear()
dev.off()


