# 資料處理與讀取
library(tidyverse)      # 含 dplyr, ggplot2, tidyr, stringr, purrr
library(haven)          # read_sav
library(labelled)       # zap_labels

# 描述統計與相關分析
library(psych)          # mardia
library(Hmisc)          # rcorr

# 視覺化
library(scales)         # percent, expansion
library(RColorBrewer)   # 配色
library(patchwork)      # 圖形拼接 (+, /)
library(viridis)        # 色盲友善配色
library(ggpubr)         # ggarrange
library(ggrepel)        # geom_text_repel

# 地圖
library(sf)             # st_read, st_centroid

# 統計檢定與 ANOVA
library(carData)
library(car)            # leveneTest, Anova
library(MASS)           # polr
library(rstatix)        # games_howell_test, dunn_test, box_m
library(emmeans)        # emmeans, pairs
library(effectsize)     # 效果量

# 迴歸模型評估
library(rms)            # lrm
library(pROC)           # roc, auc, ggroc

# 分群分析
library(cluster)        # daisy, silhouette
library(clustMixType)   # kproto
library(patchwork)

#MCA
library(FactoMineR)
library(factoextra)

# 遺失值處理
library(missMDA)        # imputeFAMD, estim_ncpFAMD
# ============================================
# 讀取資料
# ============================================

# 資料
my_data=read_sav("E:/碩一/資料漫步/data/data.sav")

# ====================================
# Part 0: FAMD 遺失值填補
# ====================================

cat("=== Part 0: FAMD 遺失值填補 ===\n\n")

# 移除前4欄（通常是 ID 或系統欄位）
my_data_raw = my_data[, -(1:4)]

# ====================================
# 【0.1】初始篩選（按題號）
# ====================================

cat("【0.1】初始篩選（按題號）\n")

# 取得所有欄位標籤
all_labels = sapply(my_data_raw, function(x) attr(x, "label"))

# 篩選目標題號
target_nums = c(1, 2, 3, 4, 16, 18, 19, 22, 23, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 38, 39)
pattern = paste0("^(", paste(target_nums, collapse = "|"), ")(\\.|\\()")
matched_indices = grep(pattern, all_labels)
new_data = my_data_raw[, matched_indices]

cat("  目標題號：", paste(target_nums, collapse = ", "), "\n")
cat("  匹配欄位數：", ncol(new_data), "\n\n")

# ====================================
# 【0.2】過濾文字填空題
# ====================================

cat("【0.2】過濾文字填空題\n")

# 檢查是否為字元型態
is_char_type = sapply(new_data, is.character)

# 檢查欄位名稱是否包含 txt 或 other
is_txt_name = grepl("txt|other", names(new_data), ignore.case = TRUE)

# 檢查標籤是否為開放式問題
is_open_ended = sapply(new_data, function(x) {
  lab = attr(x, "label")
  if (is.null(lab)) return(FALSE)
  return(grepl("文字|填寫|寫下|內容", lab, ignore.case = TRUE))
})

# 執行過濾
new_data_cleaned = new_data[, !(is_char_type | is_txt_name | is_open_ended)]

cat("✓ 文字題過濾完成\n")
cat("  剩餘欄位數：", ncol(new_data_cleaned), "\n\n")

# ====================================
# 【0.3】處理 Q28、Q29 複選題
# ====================================

cat("【0.3】處理 Q28、Q29 複選題\n")

# 複選題沒填代表「未勾選」(0)，而非隨機漏填
target_multi = grep("^q28_|^q29_", names(new_data_cleaned), value = TRUE)

for (col in target_multi) {
  new_data_cleaned[[col]][is.na(new_data_cleaned[[col]])] = 0
}

cat("✓ 複選題處理完成\n")
cat("  共處理", length(target_multi), "個欄位\n\n")

# ====================================
# 【0.4】資料轉型與品質檢查
# ====================================

cat("【0.4】資料轉型與品質檢查\n")

# 移除標籤轉為資料框
new_data_final = zap_labels(new_data_cleaned)
new_data_final = as.data.frame(new_data_final)

# 強制轉型（q2 數值，其餘因子）
for (col in names(new_data_final)) {
  if (col == "q2") {
    new_data_final[[col]] = as.numeric(as.character(new_data_final[[col]]))
  } else {
    new_data_final[[col]] = as.factor(as.character(new_data_final[[col]]))
  }
}

cat("✓ 資料型態轉換完成\n")

# 檢查變異數：移除選項過少的欄位（少於3人）
valid_cols = sapply(new_data_final, function(x) {
  obs_data = na.omit(x)
  counts = table(obs_data)
  return(length(counts) > 1 && min(counts) > 3)
})

new_data_final = new_data_final[, valid_cols]

cat("✓ 低變異欄位移除完成\n")
cat("  最終欄位數：", ncol(new_data_final), "\n")
cat("  樣本數：", nrow(new_data_final), "\n")
cat("  缺失值總數：", sum(is.na(new_data_final)), "\n\n")

# ====================================
# 【0.5】執行 FAMD 填補
# ====================================

cat("【0.5】執行 FAMD 填補\n")

set.seed(123)

# 估計最佳主成分數（含錯誤容錯機制）
cat("  估計最佳主成分數中...\n")

res_ncp = tryCatch({
  estim_ncpFAMD(new_data_final, ncp.max = 3, maxiter = 1000)
}, error = function(e) {
  cat("  ⚠ 自動估計失敗，使用預設 ncp = 2\n")
  list(ncp = 2)
})

cat("  最佳主成分數：", res_ncp$ncp, "\n")

# 執行填補
cat("  執行填補中...\n")
imputed_res = imputeFAMD(new_data_final, ncp = res_ncp$ncp, maxiter = 1000)

cat("✓ FAMD 填補完成\n\n")

# ====================================
# 【0.6】結果數值化處理
# ====================================

cat("【0.6】結果數值化處理\n")

# 清理填補結果
complete_data = as.data.frame(imputed_res$completeObs)

for (col in names(complete_data)) {
  clean_val = gsub(paste0("^", col, "_"), "", as.character(complete_data[[col]]))
  complete_data[[col]] = as.numeric(clean_val)
}

cat("✓ 數值化處理完成\n")
cat("  填補後 NA 數量：", sum(is.na(complete_data)), "\n\n")

# ====================================
# 【0.7】資料合併回原始資料
# ====================================

cat("【0.7】資料合併\n")

# 取得填補後的欄位名稱
famd_cols = names(complete_data)

# 用填補後的值替換原始資料中的對應欄位
for (col in famd_cols) {
  if (col %in% names(my_data)) {
    my_data[[col]] = complete_data[[col]]
  }
}

cat("✓ 資料合併完成\n")
cat("  已填補欄位數：", length(famd_cols), "\n")
cat("  合併後 NA 數量：", sum(is.na(my_data[, famd_cols])), "\n\n")

# ====================================
# 【0.8】最終檢查
# ====================================

cat("【0.8】最終檢查\n")

# 檢查 Q28、Q29 是否正確處理
q28_q29_cols = grep("^q28|^q29", names(complete_data), value = TRUE)
cat("  Q28、Q29 欄位：\n")
print(q28_q29_cols)

cat("\n  資料前 5 列預覽：\n")
print(head(complete_data, 5))

cat("\n✓ FAMD 遺失值填補流程完成！\n\n")
# ====================================
# 資料清理與變項建構
# ====================================
data1 = my_data %>%
  dplyr::select(-c(
    "q3_other", "q4_88_text", "q8_88_text", "q9_88_text", "q10_88_text", 
    "q11_88_text", "q12_8_text", "q13_02_8_text", "q13_03_5_text", 
    "q13_04_8_text", "q13_05_5_text", "q14_01_5_text", "q14_02_5_text", 
    "q14_03_4_text", "q15_02_13_text", "q21a_6_text", "q28_4_text", 
    "q29_5_text"
  ))
data = data1 %>%
  mutate(
    #人口變項 
    性別 = factor(q1, levels = 1:2, labels = c("男", "女")),
    出生年 = q2,
    年齡 = 110 - q2,
    年齡組 = cut(年齡, 
              breaks = c(0, 25, 35, 45, 55, 100),
              labels = c("18-25歲", "26-35歲", "36-45歲", "46-55歲", "56+"),
              right = FALSE),
    # 出生縣市
    出生縣市_raw = q3,
    
    出生縣市 = case_when(
      出生縣市_raw == 1 ~ "臺北市",
      出生縣市_raw == 2 ~ "新北市",
      出生縣市_raw == 3 ~ "基隆市",
      出生縣市_raw == 4 ~ "桃園市",
      出生縣市_raw == 5 ~ "新竹市",
      出生縣市_raw == 6 ~ "新竹縣",
      出生縣市_raw == 7 ~ "苗栗縣",
      出生縣市_raw == 8 ~ "臺中市",
      出生縣市_raw == 9 ~ "彰化縣",
      出生縣市_raw == 10 ~ "南投縣",
      出生縣市_raw == 11 ~ "雲林縣",
      出生縣市_raw == 12 ~ "嘉義市",
      出生縣市_raw == 13 ~ "嘉義縣",
      出生縣市_raw == 14 ~ "臺南市",
      出生縣市_raw == 15 ~ "高雄市",
      出生縣市_raw == 16 ~ "屏東縣",
      出生縣市_raw == 17 ~ "宜蘭縣",
      出生縣市_raw == 18 ~ "花蓮縣",
      出生縣市_raw == 19 ~ "臺東縣",
      出生縣市_raw == 20 ~ "澎湖縣",
      出生縣市_raw == 21 ~ "金門縣",
      出生縣市_raw == 22 ~ "連江縣",
      出生縣市_raw == 24 ~ "其他國家",
      TRUE ~ NA_character_
    ),
    出生縣市 = factor(出生縣市),
    
    # 出生地區分類（六大區域）
    出生地區 = case_when(
      出生縣市 %in% c("臺北市", "新北市", "基隆市", "桃園市", 
                  "新竹市", "新竹縣", "宜蘭縣") ~ "北部",
      出生縣市 %in% c("苗栗縣", "臺中市", "彰化縣", "南投縣", "雲林縣") ~ "中部",
      出生縣市 %in% c("嘉義市", "嘉義縣", "臺南市") ~ "南部",
      出生縣市 %in% c("高雄市", "屏東縣") ~ "高屏",
      出生縣市 %in% c("花蓮縣", "臺東縣") ~ "東部",  
      出生縣市 %in% c("澎湖縣", "金門縣", "連江縣") ~ "離島",
      出生縣市 == "其他國家" ~ "其他國家",
      TRUE ~ NA_character_
    ),
    出生地區 = factor(出生地區, 
                  levels = c("北部", "中部", "南部", "高屏", "東部", 
                             "離島",  "其他國家")),
    # 簡化分類：台灣本島 vs. 境外
    出生地_大分類 = case_when(
      出生地區 %in% c("北部", "中部", "南部", "高屏", "東部") ~ "台灣本島",
      出生地區 == "離島" ~ "台灣離島",
      出生地區 == "其他國家" ~ "其他國家",
      TRUE ~ NA_character_
    ),
    出生地_大分類 = factor(出生地_大分類,
                     levels = c("台灣本島", "台灣離島", "其他國家")),
    #教育程度
    教育程度 = case_when(
      q4 <= 3 ~ "國小以下",
      q4 >= 4 & q4 <= 6 ~ "國高中",
      q4 >= 7 & q4 <= 15 ~ "專科",
      q4 >= 16 & q4 <= 19 ~ "大學",
      q4 >= 20 ~ "研究所",
      TRUE ~ NA_character_
    ),
    教育程度 = factor(教育程度, 
                  levels = c("國小以下", "國高中", "專科", "大學", "研究所"),
                  ordered = TRUE),
    
    # ===== Q22: 被動攻擊曝露（觀察到的不文明言論）=====
    # 檢查並處理缺失值
    q22_01_clean = if_else(
      is.na(q22_01_1) | q22_01_1 < 1 | q22_01_1 > 4, 
      1, 
      as.numeric(q22_01_1)
    ),
    q22_02_clean = if_else(
      is.na(q22_02_1) | q22_02_1 < 1 | q22_02_1 > 4, 
      1, 
      as.numeric(q22_02_1)
    ),
    q22_03_clean = if_else(
      is.na(q22_03_1) | q22_03_1 < 1 | q22_03_1 > 4, 
      1, 
      as.numeric(q22_03_1)
    ),
    q22_04_clean = if_else(
      is.na(q22_04_1) | q22_04_1 < 1 | q22_04_1 > 4, 
      1, 
      as.numeric(q22_04_1)
    ),
    q22_05_clean = if_else(
      is.na(q22_05_1) | q22_05_1 < 1 | q22_05_1 > 4, 
      1, 
      as.numeric(q22_05_1)
    ),
    
    # 五個題項
    Q22_髒話 = q22_01_clean,
    Q22_兇人 = q22_02_clean,
    Q22_罵人 = q22_03_clean,
    Q22_不雅玩笑 = q22_04_clean,
    Q22_諷刺 = q22_05_clean,
    
    # 平均分數（1-4）
    被動攻擊分數 = (Q22_髒話 + Q22_兇人 + Q22_罵人 + Q22_不雅玩笑 + Q22_諷刺) / 5,
    
    # 總分（5-20）
    被動攻擊分數_總分 = Q22_髒話 + Q22_兇人 + Q22_罵人 + Q22_不雅玩笑 + Q22_諷刺,
    
    # 分類
    被動攻擊等級 = case_when(
      被動攻擊分數 < 1.5 ~ "極低",
      被動攻擊分數 < 2.5 ~ "低",
      被動攻擊分數 < 3.5 ~ "中",
      被動攻擊分數 >= 3.5 ~ "高",
      TRUE ~ NA_character_
    ),
    被動攻擊等級 = factor(被動攻擊等級, 
                    levels = c("極低", "低", "中", "高"),
                    ordered = TRUE),
    
    高被動攻擊 = if_else(被動攻擊分數 >= 3, "高", "低"),
    
    # ===== Q23: 主動攻擊（自己使用的不文明言論）=====
    q23_01_clean = if_else(
      is.na(q23_01_1) | q23_01_1 < 1 | q23_01_1 > 4, 
      1, 
      as.numeric(q23_01_1)
    ),
    q23_02_clean = if_else(
      is.na(q23_02_1) | q23_02_1 < 1 | q23_02_1 > 4, 
      1, 
      as.numeric(q23_02_1)
    ),
    q23_03_clean = if_else(
      is.na(q23_03_1) | q23_03_1 < 1 | q23_03_1 > 4, 
      1, 
      as.numeric(q23_03_1)
    ),
    q23_04_clean = if_else(
      is.na(q23_04_1) | q23_04_1 < 1 | q23_04_1 > 4, 
      1, 
      as.numeric(q23_04_1)
    ),
    q23_05_clean = if_else(
      is.na(q23_05_1) | q23_05_1 < 1 | q23_05_1 > 4, 
      1, 
      as.numeric(q23_05_1)
    ),
    
    Q23_髒話 = q23_01_clean,
    Q23_兇人 = q23_02_clean,
    Q23_罵人 = q23_03_clean,
    Q23_不雅玩笑 = q23_04_clean,
    Q23_諷刺 = q23_05_clean,
    
    # 平均分數
    主動攻擊分數 = (Q23_髒話 + Q23_兇人 + Q23_罵人 + Q23_不雅玩笑 + Q23_諷刺) / 5,
    
    # 總分
    主動攻擊分數_總分 = Q23_髒話 + Q23_兇人 + Q23_罵人 + Q23_不雅玩笑 + Q23_諷刺,
    
    # 分類
    主動攻擊等級 = case_when(
      主動攻擊分數 < 1.5 ~ "無",
      主動攻擊分數 < 2.5 ~ "低",
      主動攻擊分數 < 3.5 ~ "中",
      主動攻擊分數 >= 3.5 ~ "高",
      TRUE ~ NA_character_
    ),
    主動攻擊等級 = factor(主動攻擊等級, 
                    levels = c("無", "低", "中", "高"),
                    ordered = TRUE),
    
    高主動攻擊 = if_else(主動攻擊分數 >= 3, "高", "低"),
    
    # ===== Q28: 抵制行為 =====
    被動抵制 = if_else(q28_1 == 1 | q28_2 == 1, 1, 0, missing = 0),
    主動表達 = if_else(q28_3 == 1, 1, 0, missing = 0),
    
    # ===== Q16/Q19: 惡搞行為 =====
    無害惡搞 = if_else(q16 == 1, "有", "無", missing = "無"),
    有害惡搞 = if_else(q18 == 1, "有", "無", missing = "無"),
    
    # ===== 滿意度變項 =====
    生活滿意度 = q38_01_1,
    快樂程度 = q39_1,
    同理心程度 = q31_1,
    
    # ===== 四象限行為分類 =====
    行為類型 = case_when(
      # 未參與：沒有任何行為
      被動抵制 == 0 & 主動表達 == 0 & 主動攻擊分數 < 2.5 & 被動攻擊分數 < 2.5 ~ "未參與",
      
      # 被動觀察：只有被動行為（取消關注/拒看），攻擊言論低
      被動抵制 == 0 & 主動表達 == 0 & 主動攻擊分數 < 2.5 & 被動攻擊分數 >= 2.5 ~ "被動觀察",
      
      # 被動抵制：有被動行為，無主動表達，低攻擊
      被動抵制 == 1 & 主動表達 == 0 & 主動攻擊分數 < 2.5 ~ "被動抵制",
      
      # 主動問責：有主動表達，但攻擊言論低
      主動表達 == 1 & 主動攻擊分數 < 2.5 ~ "主動問責",
      
      # 主動攻擊：有主動表達且攻擊言論高，或攻擊言論非常高
      (被動抵制 == 0 & 主動表達 == 1 & 主動攻擊分數 >= 2.5) | 
        主動攻擊分數 >= 3 ~ "主動攻擊",
      
      TRUE ~ "混合/其他"
    ),
    行為類型 = factor(
      行為類型,
      levels = c("未參與", "被動觀察", "被動抵制", "主動問責", "主動攻擊", "混合/其他")
    )
  )

# ====================================
# 敘述統計
# ====================================
#性別分析
print(table(data$性別))
# 性別分布長條圖
gender_dist <- data %>%
  dplyr::count(性別) %>%
  dplyr::mutate(
    pct = n / sum(n) * 100,
    label = paste0(n, "\n(", round(pct, 1), "%)")
  )

p_gender <- ggplot(gender_dist, aes(x = 性別, y = n, fill = 性別)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = label), vjust = 0.1, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("男" = "skyblue", "女" = "pink")) +
  labs(
    title = "受訪者性別分布",
    subtitle = paste0("N = ", nrow(data)),
    x = "性別",
    y = "人數"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(size = 11, face = "bold")
  )

print(p_gender)
# 年齡分析
print(table(data$年齡組))
# 年齡組分布長條圖
age_dist <- data %>%
  dplyr::count(年齡組) %>%
  dplyr::mutate(
    pct = n / sum(n) * 100,
    label = paste0(n, "\n(", round(pct, 1), "%)")
  )

p_age <- ggplot(age_dist, aes(x = 年齡組, y = n, fill = 年齡組)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = label), vjust = 0.1, size = 4, fontface = "bold") +
  scale_fill_brewer(palette = "YlOrRd") +
  labs(
    title = "受訪者年齡分布",
    subtitle = paste0("N = ", nrow(data)),
    x = "年齡組",
    y = "人數"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(size = 11, face = "bold")
  )

print(p_age)


# 檢查出生地區分布
cat("出生地區分布：\n")
print(table(data$出生地區, useNA = "ifany"))

cat("\n出生縣市分布（前10名）：\n")
print(head(sort(table(data$出生縣市), decreasing = TRUE), 10))

# 出生地分析

# 計算各縣市樣本數（只計算台灣地區）
sample_dist_taiwan = data %>%
  filter(出生地_大分類 %in% c("台灣本島", "台灣離島")) %>%  # 只保留台灣
  count(出生縣市, name = "樣本數") %>%
  mutate(
    百分比 = 樣本數 / sum(樣本數) * 100,
    標籤 = paste0(樣本數, "(", round(百分比, 1), "%)")
  )

cat("=== 台灣地區樣本分布 ===\n")
print(sample_dist_taiwan %>% arrange(desc(樣本數)))

# 境外樣本統計
foreign_count <- data %>%
  filter(出生地_大分類 %in% "其他國家") %>%
  count(出生地_大分類)

cat("\n=== 境外樣本 ===\n")
print(foreign_count)

# 下載與合併地圖
taiwan_map = st_read("https://raw.githubusercontent.com/g0v/twgeojson/master/json/twCounty2010.geo.json")
head(taiwan_map)
# 修復地圖
taiwan_map = st_make_valid(taiwan_map)
taiwan_map = st_transform(taiwan_map, crs=3826)

# 建立縣市名稱
taiwan_map = taiwan_map %>%
  mutate(
    縣市名稱 = if("COUNTYNAME" %in% names(.)) {
      str_replace(COUNTYNAME, "台", "臺")
    } else if("name" %in% names(.)) {
      str_replace(name, "台", "臺")
    } else {
      # 如果都沒有，使用第一個文字欄位
      names(.)[1]
    }
  )
# 計算樣本數
sample_dist = data %>%
  filter(出生地_大分類 %in% c("台灣本島", "台灣離島")) %>%
  count(出生縣市, name = "樣本數") %>%
  mutate(百分比 = 樣本數 / sum(樣本數) * 100)

# 合併資料
map_data = taiwan_map %>%
  left_join(sample_dist, by = c("縣市名稱" = "出生縣市")) %>%
  mutate(樣本數 = replace_na(樣本數, 0))

map_data = map_data %>%
  mutate(
    樣本等級 = cut(
      樣本數,
      breaks = c(-1, 0, 50, 100, 200, 300, Inf),
      labels = c("0人", "1-50人", "51-100人", "101-200人", "201-300人", "300人以上")
    )
  )
# 計算縣市中心點（用於標註）
centroids = st_centroid(map_data)
# 繪製地圖
p_map = ggplot(map_data) +
  geom_sf(aes(fill = 百分比), color = "white", size = 0.5) +
  geom_text_repel(
    data = centroids %>% 
      mutate(
        lon = st_coordinates(.)[,1],
        lat = st_coordinates(.)[,2],
        標籤 = if_else(
          樣本數 > 0,
          paste0(縣市名稱, "\n", sprintf("%.1f", 百分比), "%"),
          縣市名稱
        )
      ),
    aes(x = lon, y = lat, label = 標籤),
    size = 3.5,
    fontface = "bold",
    color = "black",
    bg.color = "white",
    bg.r = 0.15,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size = 0.3,
    max.overlaps = Inf,
    min.segment.length = 0,
    lineheight = 0.85
  ) +
  
  scale_fill_gradient(
    low = "#EFF3FF", 
    high = "#08519C",
    name = "百分比 (%)",
    na.value = "grey95"
  ) +
  
  labs(
    title = paste0("樣本台灣分布 (N=", sum(sample_dist$樣本數), ")"),
    subtitle = "依出生縣市統計"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )
print(p_map)

# 完整出生地統計表
birth_summary = data %>%
  count(出生地_大分類, name = "樣本數") %>%
  mutate(百分比 = 樣本數 / sum(樣本數) * 100) %>%
  arrange(desc(樣本數))

cat("\n=== 完整出生地統計 ===\n")
print(birth_summary)

# 2. 台灣地區詳細統計
taiwan_detail = data %>%
  filter(出生地_大分類 %in% c("台灣本島", "台灣離島")) %>%
  count(出生地區, name = "樣本數") %>%
  mutate(百分比 = 樣本數 / sum(樣本數) * 100) %>%
  arrange(desc(樣本數))

cat("\n=== 台灣地區詳細統計 ===\n")
print(taiwan_detail)

# 3. 境外樣本詳細資訊
if(sum(data$出生地_大分類 == "其他國家", na.rm = TRUE) > 0) {
  other_countries_detail = data %>%
    filter(出生地_大分類 == "其他國家")
  
  cat("\n=== 其他國家詳細資訊 ===\n")
  print(other_countries_detail)
}


# ====================================
# EDA - 長條圖
# ====================================

cat("=== Part 4: 長條圖 ===\n")

# ------ 圖1：行為類型分布（核心圖）------
behavior_dist <- data %>%
  dplyr::count(行為類型) %>%
  dplyr::mutate(
    pct = n / sum(n) * 100,
    label = paste0(n, "\n(", round(pct, 1), "%)")
  )

p1 <- ggplot(behavior_dist, aes(x = 行為類型, y = n, fill = 行為類型)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = label), vjust = 0.1, size = 4, fontface = "bold") +
  scale_fill_manual(
    values = c(
      "未參與" = "gray95",
      "被動觀察" = "#FFA500",
      "被動抵制" = "gold",
      "主動問責" = "#00BA38",
      "主動攻擊" = "#F8766D",
      "混合/其他" = "lightblue"
    )
  ) +
  labs(
    title = "台灣網民網路正義行為分布",
    subtitle = paste0("N = ", nrow(data)),
    x = "行為類型",
    y = "人數"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 11, face = "bold")
  )

print(p1)

# ------ 圖2：Q22被動攻擊分數分布 ------
q22_long <- data %>%
  dplyr::select(Q22_髒話, Q22_兇人, Q22_罵人, Q22_不雅玩笑, Q22_諷刺) %>%
  tidyr::pivot_longer(everything(), names_to = "題項", values_to = "分數") %>%
  dplyr::mutate(題項 = stringr::str_remove(題項, "Q22_"))

# 計算百分比（用於標籤）
q22_summary <- q22_long %>%
  dplyr::count(題項, 分數) %>%
  dplyr::group_by(題項) %>%
  dplyr::mutate(
    pct = n / sum(n) * 100,
    label = paste0(round(pct, 1), "%")
  ) %>%
  dplyr::ungroup()

p2 <- ggplot(q22_summary, aes(x = factor(分數), y = pct, fill = factor(分數))) +
  geom_col(alpha = 0.85, width = 0.7) +
  geom_text(aes(label = label), vjust = -0.5, size = 3, fontface = "bold") +
  scale_fill_manual(
    values = c(
      "1" = "#52C7B8",
      "2" = "#FFD662",
      "3" = "#FFA07A",
      "4" = "#E76F51"
    ),
    labels = c("從來沒有", "很少", "有時", "經常"),
    name = "頻率"
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.1))
  ) +
  facet_wrap(~題項, ncol = 3) +
  labs(
    title = "被動攻擊 - 觀察到不文明言論的頻率",
    subtitle = paste0("N = ", nrow(data)),
    x = "頻率 (1=從來沒有, 4=經常)",
    y = "百分比"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "grey95", color = NA),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p2)

# ------ 圖3：Q23主動攻擊分布 ------
q23_long <- data %>%
  dplyr::select(Q23_髒話, Q23_兇人, Q23_罵人, Q23_不雅玩笑, Q23_諷刺) %>%
  tidyr::pivot_longer(everything(), names_to = "題項", values_to = "分數") %>%
  dplyr::mutate(題項 = stringr::str_remove(題項, "Q23_"))

# 計算百分比
q23_summary <- q23_long %>%
  dplyr::count(題項, 分數) %>%
  dplyr::group_by(題項) %>%
  dplyr::mutate(
    pct = n / sum(n) * 100,
    label = paste0(round(pct, 1), "%")
  ) %>%
  dplyr::ungroup()

p3 <- ggplot(q23_summary, aes(x = factor(分數), y = pct, fill = factor(分數))) +
  geom_col(alpha = 0.85, width = 0.7) +
  geom_text(aes(label = label), vjust = -0.5, size = 3, fontface = "bold") +
  scale_fill_manual(
    values = c(
      "1" = "#52C7B8",
      "2" = "#FFD662",
      "3" = "#FFA07A",
      "4" = "#E76F51"
    ),
    labels = c("從來沒有", "很少", "有時", "經常"),
    name = "頻率"
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.1))
  ) +
  facet_wrap(~題項, ncol = 3) +
  labs(
    title = "主動攻擊 - 自己使用不文明言論的頻率",
    subtitle = paste0("N = ", nrow(data)),
    x = "頻率 (1=從來沒有, 4=經常)",
    y = "百分比"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "grey95", color = NA),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p3)

# ====================================
# Part 6: 滿意度分析
# ====================================

cat("=== Part 6: 滿意度分析 ===\n")

# ANOVA
aov_life = aov(生活滿意度 ~ 行為類型, data = data)
cat("\n【生活滿意度 ANOVA】\n")
print(summary(aov_life))

aov_happy = aov(快樂程度 ~ 行為類型, data = data)
cat("\n【快樂程度 ANOVA】\n")
print(summary(aov_happy))

# 視覺化
satisfaction_long = data %>%
  dplyr::select(行為類型, 生活滿意度, 快樂程度) %>%
  tidyr::pivot_longer(
    cols = c(生活滿意度, 快樂程度),
    names_to = "變項",
    values_to = "分數"
  )

p7 = ggplot(satisfaction_long, 
            aes(x = 行為類型, y = 分數, fill = 行為類型)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", 
               shape = 23, size = 3, fill = "white") +
  facet_wrap(~變項) +
  scale_fill_manual(
    values = c("未參與" = "gray70",
               "被動觀察" = "#FFA500",
               "被動抵制" = "gold",
               "主動問責" = "#00BA38",
               "主動攻擊" = "#F8766D",
               "混合/其他" = "lightblue")
  ) +
  labs(
    title = "不同行為類型的滿意度比較",
    subtitle = "白色菱形 = 平均數",
    x = "行為類型",
    y = "分數 (1-5)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 12)
  )

print(p7)

# ====================================
# 五類行為類型二維定位示意圖
# ====================================

quadrant_5 = data.frame(
  x = c(1, 1, 2, 3, 3),
  y = c(1, 2.5, 2, 1.5, 3),
  類型 = c("未參與", "被動觀察", "被動抵制", "主動問責", "主動攻擊"),
  描述 = c(
    "不關注\n無行動",
    "常觀察\n但不行動",
    "取消關注\n拒看內容",
    "公開指責\n理性批評",
    "言語攻擊\n情緒宣洩"
  )
)

p_quadrant_5 = ggplot(quadrant_5, aes(x = x, y = y)) +
  # 背景區塊
  annotate("rect", xmin = 0.4, xmax = 1.7, ymin = 0.4, ymax = 2.3, 
           fill = "#E8F8F5", alpha = 0.4) +
  annotate("rect", xmin = 1.7, xmax = 2.4, ymin = 0.4, ymax = 2.3, 
           fill = "#FEF9E7", alpha = 0.4) +
  annotate("rect", xmin = 2.4, xmax = 3.6, ymin = 0.4, ymax = 2.3, 
           fill = "#EBF5FB", alpha = 0.4) +
  annotate("rect", xmin = 2.4, xmax = 3.6, ymin = 2.3, ymax = 3.6, 
           fill = "#FDEDEC", alpha = 0.4) +
  # 區塊標籤
  annotate("text", x = 1.05, y = 3.4, label = "無行動區", 
           size = 4, fontface = "italic", color = "gray50") +
  annotate("text", x = 2.05, y = 3.4, label = "被動區", 
           size = 4, fontface = "italic", color = "gray50") +
  annotate("text", x = 3, y = 3.4, label = "主動區", 
           size = 4, fontface = "italic", color = "gray50") +
  # 點
  geom_point(aes(color = 類型), size = 22, alpha = 0.7) +
  geom_text(aes(label = 類型), size = 5, fontface = "bold") +
  geom_text(aes(label = 描述), vjust = 3.2, size = 3.5, lineheight = 0.85) +
  # 分隔線
  geom_vline(xintercept = c(1.5, 2.4), linetype = "dotted", color = "gray50", linewidth = 0.5) +
  # 配色
  scale_color_manual(
    values = c(
      "未參與" = "#95A5A6",
      "被動觀察" = "#F39C12",
      "被動抵制" = "gold",
      "主動問責" = "#27AE60",
      "主動攻擊" = "#E74C3C"
    )
  ) +
  scale_x_continuous(
    limits = c(0.4, 3.6),
    breaks = c(1, 2, 3),
    labels = c("無", "被動", "主動")
  ) +
  scale_y_continuous(
    limits = c(0.4, 3.6),
    breaks = c(1, 2, 3),
    labels = c("低", "中", "高")
  ) +
  labs(
    title = "網路正義行為類型分布圖",
    x = "參與程度",
    y = "攻擊程度"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 13, face = "bold"),
    axis.title = element_text(size = 14, face = "bold")
  )

print(p_quadrant_5)


# 交叉分析：行為類型 × 有害惡搞
cat("\n\n【行為類型 × 有害惡搞 交叉表】\n")
cross_tab = table(data$行為類型, data$有害惡搞)
print(cross_tab)

# 視覺化
#圖一分組長條圖
# 計算百分比
q19_summary = data %>%
  group_by(行為類型, 有害惡搞) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(行為類型) %>%
  mutate(
    total = sum(n),
    pct = n / total * 100,
    label = paste0(n, "\n(", round(pct, 1), "%)")
  )
# 分組長條圖
p_q19_v1 = ggplot(q19_summary, aes(x = 行為類型, y = pct, fill = 有害惡搞)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = label), 
            position = position_dodge(width = 0.8),
            vjust = 0, size = 3.5) +
  scale_fill_manual(values = c("無" = "skyblue", "有" = "pink")) +
  labs(
    title = "不同行為類型參與有害惡搞的比例",
    x = "行為類型",
    y = "百分比 (%)",
    fill = "有害惡搞"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold"),
    legend.position = "top"
  )

print(p_q19_v1)

# ====================================
# 核心研究變項：4大熱力圖
# ====================================

# ===== 熱力圖1：攻擊敏感度機制 =====
cat("【1/4】攻擊敏感度機制熱力圖...\n")

# 計算攻擊敏感指數
data = data %>%
  mutate(
    攻擊敏感度指數 = (q20_01_1 + q20_02_1 + q25_01_1 + q25_02_1 + q25_03_1 + q25_04_1) / 6,
    攻擊敏感程度 = case_when(
      攻擊敏感度指數 < 2 ~ "低敏感度",
      攻擊敏感度指數 < 3 ~ "中敏感度",
      TRUE ~ "高敏感度"
    ),
    攻擊敏感程度 = factor(攻擊敏感程度, 
                     levels = c("低敏感度", "中敏感度", "高敏感度"))
  )

moral_data = data %>%
  count(攻擊敏感程度, 行為類型) %>%
  group_by(攻擊敏感程度) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

p1_moral = ggplot(moral_data, 
                   aes(x = 攻擊敏感程度, y = 行為類型, fill = pct)) +
  geom_tile(color = "white", linewidth = 2) +
  geom_text(aes(label = paste0(round(pct, 1), "%\n(n=", n, ")")), 
            color = "black", size = 5, fontface = "bold", lineheight = 0.9) +
  scale_fill_gradient2(
    low = "#2166AC", 
    mid = "white", 
    high = "#B2182B",
    midpoint = 20, 
    name = "比例 (%)"
  ) +
  labs(
    title = "攻擊敏感程度 × 行為類型",
    subtitle = "敏感程度越高，越傾向極端攻擊行為",
    x = "攻擊敏感程度",
    y = "行為類型"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "#8B0000", size = 11),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    legend.position = "right",
    panel.grid = element_blank()
  )

print(p1_moral)

# ===== 熱力圖2：情緒驅動機制 =====
cat("【2/4】情緒驅動機制熱力圖...\n")

emotion_data = data %>%
  mutate(
    生活滿意度_等級 = case_when(
      生活滿意度 <= 2 ~ "不滿意",
      生活滿意度 == 3 ~ "普通",
      TRUE ~ "滿意"
    ),
    生活滿意度_等級 = factor(生活滿意度_等級, 
                      levels = c("不滿意", "普通", "滿意")),
    
    快樂程度_等級 = case_when(
      快樂程度 <= 2 ~ "不快樂",
      快樂程度 == 3 ~ "普通",
      TRUE ~ "快樂"
    ),
    快樂程度_等級 = factor(快樂程度_等級, 
                     levels = c("不快樂", "普通", "快樂")),
    
    主動攻擊_等級 = case_when(
      主動攻擊分數 < 1.5 ~ "低攻擊",
      主動攻擊分數 < 2.5 ~ "中攻擊",
      TRUE ~ "高攻擊"
    ),
    主動攻擊_等級 = factor(主動攻擊_等級, 
                     levels = c("低攻擊", "中攻擊", "高攻擊"))
  ) %>%
  count(生活滿意度_等級, 快樂程度_等級, 主動攻擊_等級) %>%
  group_by(主動攻擊_等級) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

p2_emotion = ggplot(emotion_data, 
                     aes(x = 快樂程度_等級, y = 生活滿意度_等級, fill = pct)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(label = paste0(round(pct, 1), "%","\n(n=", n, ")")), 
            color = "black", size = 4.5, fontface = "bold") +
  facet_wrap(~主動攻擊_等級, ncol = 3) +
  scale_fill_gradient(
    low = "#FFF5F0", 
    high = "#8B0000", 
    name = "比例 (%)"
  ) +
  labs(
    title = "情緒狀態 × 主動攻擊行為",
    subtitle = "生活不滿意 + 不快樂 → 高攻擊傾向",
    x = "快樂程度",
    y = "生活滿意度"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "#8B0000", size = 10),
    strip.text = element_text(size = 11, face = "bold", color = "white"),
    strip.background = element_rect(fill = "#8B0000", color = NA),
    axis.text = element_text(size = 10, face = "bold"),
    legend.position = "right",
    panel.border = element_rect(color = "gray80", fill = NA)
  )

print(p2_emotion)


# ===== 熱力圖3：同理心保護因子 =====
cat("【3/4】同理心保護因子熱力圖...\n")

# 同理心 × 惡搞行為 × 行為類型 熱力圖
empathy_data = data %>%
  mutate(
    同理心等級 = case_when(
      同理心程度 <= 2 ~ "低同理心",
      同理心程度 == 3 ~ "中同理心",
      同理心程度 == 4 ~ "高同理心",
      TRUE ~ NA_character_
    ),
    同理心等級 = factor(同理心等級, 
                   levels = c("低同理心", "中同理心", "高同理心")),
    
    惡搞類型 = case_when(
      有害惡搞 == "有" ~ "有害惡搞",
      無害惡搞 == "有" ~ "無害惡搞",
      TRUE ~ "無惡搞"
    ),
    惡搞類型 = factor(惡搞類型, levels = c("無惡搞", "無害惡搞", "有害惡搞"))
  ) %>%
  filter(!is.na(同理心等級)) %>%
  count(同理心等級, 惡搞類型, 行為類型) %>%
  group_by(同理心等級) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

p3_empathy = ggplot(empathy_data, 
                    aes(x = 惡搞類型, y = 行為類型, fill = pct)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(label = paste0(round(pct, 1), "%\n(n=", n, ")")), 
            color = "black", size = 3, fontface = "bold", lineheight = 0.9) +
  facet_wrap(~同理心等級, ncol = 3) +
  scale_fill_gradient(
    low = "skyblue", 
    high = "pink", 
    name = "比例 (%)"
  ) +
  labs(
    title = "同理心 × 惡搞行為 × 行為類型",
    subtitle = "高同理心作為保護因子：降低有害惡搞參與",
    x = "惡搞類型",
    y = "行為類型"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "blue", size = 10),
    strip.text = element_text(size = 11, face = "bold", color = "white"),
    strip.background = element_rect(fill = "blue", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    legend.position = "right",
    panel.border = element_rect(color = "gray80", fill = NA)
  )

print(p3_empathy)

# ===== 熱力圖4：核心變項相關矩陣 =====
cat("【4/4】核心變項相關矩陣熱力圖...\n")

# Spearman 相關係數熱力圖

# 準備資料
cor_data <- data %>%
  dplyr::select(
    被動攻擊分數, 
    主動攻擊分數,
    攻擊敏感度指數,
    生活滿意度, 
    快樂程度,
    同理心程度,
    年齡
  ) %>%
  na.omit()

cat("樣本數：", nrow(cor_data), "\n")

# 計算 Spearman 相關係數與 p 值
cor_spearman <- rcorr(as.matrix(cor_data), type = "spearman")
r_spearman <- cor_spearman$r
p_spearman <- cor_spearman$P

# 轉換為長格式
spearman_long <- expand.grid(
  變項1 = rownames(r_spearman),
  變項2 = colnames(r_spearman),
  stringsAsFactors = FALSE
) %>%
  mutate(
    r = mapply(function(x, y) r_spearman[x, y], 變項1, 變項2),
    p = mapply(function(x, y) p_spearman[x, y], 變項1, 變項2),
    sig = case_when(
      is.na(p) ~ "",
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ ""
    ),
    label = ifelse(
      變項1 == 變項2, 
      "1.00", 
      paste0(format(round(r, 2), nsmall = 2), sig)
    )
  )

# 繪製熱力圖
p_spearman_heatmap <- ggplot(spearman_long, 
                             aes(x = 變項2, y = 變項1, fill = r)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(label = label), color = "black", size = 4.5, fontface = "bold") +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman ρ"
  ) +
  scale_y_discrete(limits = rev) +
  labs(
    title = "Spearman 相關係數矩陣",
    subtitle = "*** p < .001, ** p < .01, * p < .05",
    caption = paste0("N = ", nrow(cor_data)),
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    plot.caption = element_text(hjust = 1, size = 10, face = "italic"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    legend.position = "right",
    panel.grid = element_blank()
  )

print(p_spearman_heatmap)

# 輸出相關係數表格
cat("\n【Spearman 相關係數矩陣】\n")
print(round(r_spearman, 3))

cat("\n【顯著性矩陣】\n")
print(round(p_spearman, 4))

# ===== 關鍵發現摘要 =====

# 1. 攻擊敏感度效果
moral_summary = data %>%
  group_by(攻擊敏感程度) %>%
  summarise(
    主動攻擊比例 = round(mean(行為類型 == "主動攻擊", na.rm = TRUE),4)
  )

cat("【發現1】攻擊敏感度效果\n")
print(moral_summary)

# 2. 同理心保護效果
empathy_summary = data %>%
  group_by(同理心程度) %>%
  summarise(
    極端行為比例 = mean(q18 == 1 & (q19_01 == 1 | q19_02 == 1), na.rm = TRUE)
  )

cat("【發現2】同理心保護效果\n")
print(empathy_summary)

# ====================================
# 網路正義行為研究 - 統計分析
# ====================================
# 第1部分：邏輯迴歸（Logistic Regression）
cat("【1】邏輯迴歸：預測高主動攻擊行為")

# 準備邏輯迴歸數據
logistic_data = data %>%
  mutate(
    高主動攻擊 = if_else(主動攻擊分數 >= 2.5, 1, 0),
    性別_numeric = if_else(性別 == "男", 1, 0),
    年齡_標準化 = scale(年齡)[,1],
    教育 = factor(教育程度),
    出生地 = factor(出生縣市),
    被動攻擊_標準化 = scale(被動攻擊分數)[,1],
    攻擊敏感度_標準化 = scale(攻擊敏感度指數)[,1],
    同理心_標準化 = scale(同理心程度)[,1],
    滿意度_標準化 = scale(生活滿意度)[,1],
    快樂度_標準化 = scale(快樂程度)[,1]
  ) %>%
  dplyr::select(
    高主動攻擊, 年齡_標準化, 性別_numeric, 被動攻擊_標準化,教育,出生地, 
    攻擊敏感度_標準化, 同理心_標準化, 滿意度_標準化,快樂度_標準化
  ) %>%
  na.omit()

cat("樣本數：", nrow(logistic_data), "\n")
cat("高主動攻擊:", sum(logistic_data$高主動攻擊), "人 (", 
    round(sum(logistic_data$高主動攻擊)/nrow(logistic_data)*100, 1), "%)\n\n")

# 2.1 單變項邏輯迴歸
cat("【2.1】單變項邏輯迴歸\n")
univariate_vars = c("年齡_標準化", "性別_numeric", "被動攻擊_標準化","教育", "出生地",
                     "攻擊敏感度_標準化", "同理心_標準化", "滿意度_標準化","快樂度_標準化")

univariate_results = list()
for(var in univariate_vars) {
  formula = as.formula(paste("高主動攻擊 ~", var))
  model = glm(formula, data = logistic_data, family = "binomial")
  univariate_results[[var]] = model
}

print(univariate_results)

# 2.2 多變項邏輯迴歸
cat("\n【2.2】多變項邏輯迴歸\n")
model_logistic_full = glm(
  高主動攻擊 ~ 年齡_標準化+性別_numeric+被動攻擊_標準化+
  攻擊敏感度_標準化+ 快樂度_標準化,
  data = logistic_data,
  family = "binomial"
)
print(summary(model_logistic_full))

# 建立 logistic_results（在 model_logistic_full 後面加上）
coef_full = summary(model_logistic_full)$coefficients[-1, ]
logistic_results <- tibble(
  變項 = rownames(coef_full),
  係數 = coef_full[, 1],
  標準誤 = coef_full[, 2],
  z值 = coef_full[, 3],
  p值 = coef_full[, 4],
  OR = exp(係數),
  CI_下 = exp(係數 - 1.96 * 標準誤),
  CI_上 = exp(係數 + 1.96 * 標準誤),
  顯著 = case_when(
    p值 < 0.001 ~ "***",
    p值 < 0.01 ~ "**",
    p值 < 0.05 ~ "*",
    TRUE ~ ""
  )
)

# 2.3 模型擬合度
cat("\n【2.3】模型擬合度評估\n")
AIC_full = AIC(model_logistic_full)
BIC_full = BIC(model_logistic_full)
cat("AIC:", round(AIC_full, 2), "\n")
cat("BIC:", round(BIC_full, 2), "\n")

# Nagelkerke R²
cat("Nagelkerke R²:", round(lrm(高主動攻擊 ~ 年齡_標準化+性別_numeric+被動攻擊_標準化+
                                  攻擊敏感度_標準化+快樂度_標準化, data = logistic_data)$stats['R2'], 3), "\n")

# 2.4 ROC 曲線與 AUC
cat("\n【2.4】ROC 曲線分析\n")
pred_prob = predict(model_logistic_full, type = "response")
roc_obj = roc(logistic_data$高主動攻擊, pred_prob)
auc_value = auc(roc_obj)

cat("AUC:", round(auc_value, 3), "\n")
cat("解釋：模型正確分類的概率為", round(auc_value*100, 1), "%\n\n")

# ROC 圖
p_roc = ggroc(roc_obj, color = "#E31A1C", size = 1.2) +
  geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "gray50") +
  annotate("text", x = 0.6, y = 0.2, 
           label = paste0("AUC = ", round(auc_value, 3)), 
           size = 5, fontface = "bold") +
  labs(
    title = "ROC 曲線：預測主動攻擊行為",
    x = "1 - 特異度",
    y = "靈敏度"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 11)
  )

print(p_roc)

# 2.5 OR 值圖表
p_or = ggplot(logistic_results, aes(x = OR, y = reorder(變項, OR))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 1) +
  geom_point(color = "black", size = 3) +
  geom_errorbar(
    aes(xmin = CI_下, xmax = CI_上),
    width = 0.2, size = 1
  ) +
  geom_text(aes(label = paste0(round(OR, 3), " ", 顯著)), 
            hjust = -0.3, size = 5, fontface = "bold",nudge_y = -0.2) +
  scale_x_log10() +
  labs(
    title = "勝算比圖（OR）：主動攻擊的風險因素",
    subtitle = "95% 信心區間 | 紅線 = OR=1（無效應）",
    x = "勝算比（Odds Ratio）",
    y = "變項"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.text = element_text(size = 14, face = "bold")
  )

print(p_or)

# 第3部分：ANOVA 與事後檢定

# 【3.0】準備資料

cat("【3.0】準備 ANOVA 資料\n")

anova_data = data %>%
  dplyr::select(
    主動攻擊分數, 被動攻擊分數, 攻擊敏感度指數,
    同理心程度, 生活滿意度, 快樂程度,
    年齡組, 性別, 出生地區, 教育程度
  ) %>%
  na.omit()

cat("樣本數：", nrow(anova_data), "\n\n")

# 【3.1】ANOVA 前提假設檢定

cat("【3.1】ANOVA 前提假設檢定\n\n")

# 3.1.1 常態性檢定（Shapiro-Wilk）
cat("--- 3.1.1 常態性檢定（Shapiro-Wilk）---\n")

# 依變項常態性
dv_vars = c("主動攻擊分數", "被動攻擊分數", "攻擊敏感度指數", 
            "同理心程度", "生活滿意度", "快樂程度")

normality_results = data.frame()

for (var in dv_vars) {
    test_result = shapiro.test(anova_data[[var]])
  
  normality_results = rbind(normality_results, data.frame(
    變項 = var,
    W統計量 = round(test_result$statistic, 4),
    p值 = round(test_result$p.value, 4),
    結論 = ifelse(test_result$p.value < 0.05, "違反常態", "符合常態")
  ))
}

print(normality_results)

# 各組內常態性（以性別 × 主動攻擊分數為例）
cat("\n各組內常態性檢定（性別 × 主動攻擊分數）：\n")
group_normality = anova_data %>%
  group_by(性別) %>%
  dplyr::summarise(
    n = dplyr::n(),
    W = shapiro.test(主動攻擊分數)$statistic,
    p = shapiro.test(主動攻擊分數)$p.value,
    .groups = "drop"
  )
print(group_normality)

# 3.1.2 變異數同質性檢定（Levene's Test）
cat("\n 3.1.2 變異數同質性檢定（Levene's Test）\n")

levene_results = data.frame()

# 性別
levene_sex = leveneTest(主動攻擊分數 ~ 性別, data = anova_data)
levene_results = rbind(levene_results, data.frame(
  分組變項 = "性別",
  依變項 = "主動攻擊分數",
  F值 = round(levene_sex$`F value`[1], 3),
  p值 = round(levene_sex$`Pr(>F)`[1], 4),
  結論 = ifelse(levene_sex$`Pr(>F)`[1] < 0.05, "違反同質性", "符合同質性")
))

# 年齡組
levene_age = leveneTest(主動攻擊分數 ~ 年齡組, data = anova_data)
levene_results = rbind(levene_results, data.frame(
  分組變項 = "年齡組",
  依變項 = "主動攻擊分數",
  F值 = round(levene_age$`F value`[1], 3),
  p值 = round(levene_age$`Pr(>F)`[1], 4),
  結論 = ifelse(levene_age$`Pr(>F)`[1] < 0.05, "違反同質性", "符合同質性")
))

# 出生地區
levene_region = leveneTest(主動攻擊分數 ~ 出生地區, data = anova_data)
levene_results = rbind(levene_results, data.frame(
  分組變項 = "出生地區",
  依變項 = "主動攻擊分數",
  F值 = round(levene_region$`F value`[1], 3),
  p值 = round(levene_region$`Pr(>F)`[1], 4),
  結論 = ifelse(levene_region$`Pr(>F)`[1] < 0.05, "違反同質性", "符合同質性")
))

# 教育程度
levene_edu = leveneTest(主動攻擊分數 ~ 教育程度, data = anova_data)
levene_results = rbind(levene_results, data.frame(
  分組變項 = "教育程度",
  依變項 = "主動攻擊分數",
  F值 = round(levene_edu$`F value`[1], 3),
  p值 = round(levene_edu$`Pr(>F)`[1], 4),
  結論 = ifelse(levene_edu$`Pr(>F)`[1] < 0.05, "違反同質性", "符合同質性")
))

cat("\nLevene's Test 結果：\n")
print(levene_results)

cat("\n✓ 前提假設檢定完成\n")
cat("建議：若違反常態性，可使用無母數檢定（Kruskal-Wallis）\n")
cat("建議：若違反變異數同質性，可使用 Welch's ANOVA\n\n")

# 【3.2】單因子 ANOVA

cat("【3.2】單因子 ANOVA）\n\n")

# 3.2.1 性別差異
# 標準 ANOVA
anova_sex = aov(主動攻擊分數 ~ 性別, data = anova_data)
cat("\n標準 ANOVA：\n")
print(summary(anova_sex))

# Welch's ANOVA（不假設變異數同質）
welch_sex = oneway.test(主動攻擊分數 ~ 性別, data = anova_data, var.equal = FALSE)
cat("\nWelch's ANOVA（穩健版）：\n")
cat("F =", round(welch_sex$statistic, 3), 
    ", df1 =", round(welch_sex$parameter[1], 2),
    ", df2 =", round(welch_sex$parameter[2], 2),
    ", p =", format.pval(welch_sex$p.value), "\n")

# 事後檢定
cat("\n事後檢定（Tukey HSD）：\n")
print(TukeyHSD(anova_sex))

# Games-Howell（不假設變異數同質）
cat("\n事後檢定（Games-Howell，穩健版）：\n")
games_howell_sex = anova_data %>%
  games_howell_test(主動攻擊分數 ~ 性別)
print(games_howell_sex)

#  3.2.2 年齡組差異
cat("\n--- 3.2.2 年齡組差異檢驗 ---\n")

anova_age = aov(主動攻擊分數 ~ 年齡組, data = anova_data)
cat("\n標準 ANOVA：\n")
print(summary(anova_age))

# Welch's ANOVA
welch_age = oneway.test(主動攻擊分數 ~ 年齡組, data = anova_data, var.equal = FALSE)
cat("\nWelch's ANOVA：\n")
cat("F =", round(welch_age$statistic, 3), ", p =", format.pval(welch_age$p.value), "\n")

# 事後檢定
cat("\n事後檢定（Tukey HSD）：\n")
age_tukey = TukeyHSD(anova_age)
print(age_tukey)

# Games-Howell（不假設變異數同質）
cat("\n事後檢定（Games-Howell，穩健版）：\n")
games_howell_age = anova_data %>%
  games_howell_test(主動攻擊分數 ~ 年齡組)
print(games_howell_age)

# 3.2.3 出生地區差異 
cat("\n--- 3.2.3 出生地區差異檢驗 ---\n")

anova_region = aov(主動攻擊分數 ~ 出生地區, data = anova_data)
cat("\n標準 ANOVA：\n")
print(summary(anova_region))

cat("\n事後檢定（Tukey HSD）：\n")
print(TukeyHSD(anova_region))

# 3.2.4 教育程度差異
cat("\n--- 3.2.4 教育程度差異檢驗 ---\n")

anova_edu = aov(主動攻擊分數 ~ 教育程度, data = anova_data)
cat("\n標準 ANOVA：\n")
print(summary(anova_edu))

cat("\n事後檢定（Tukey HSD）：\n")
print(TukeyHSD(anova_edu))

# ====================================
# 【3.3】無母數替代：Kruskal-Wallis 檢定
# ====================================

cat("\n【3.3】無母數替代：Kruskal-Wallis 檢定\n")
cat("當違反常態性假設時使用\n\n")

# 性別
kw_sex = kruskal.test(主動攻擊分數 ~ 性別, data = anova_data)
cat("--- 性別 ---\n")
cat("χ² =", round(kw_sex$statistic, 3), 
    ", df =", kw_sex$parameter,
    ", p =", format.pval(kw_sex$p.value), "\n")

# 年齡組
kw_age = kruskal.test(主動攻擊分數 ~ 年齡組, data = anova_data)
cat("--- 年齡組 ---\n")
cat("χ² =", round(kw_age$statistic, 3), 
    ", df =", kw_age$parameter,
    ", p =", format.pval(kw_age$p.value), "\n")

epsilon_age = anova_data %>%
  kruskal_effsize(主動攻擊分數 ~ 年齡組)

print(epsilon_age)

# 事後檢定（Dunn's test）
cat("\n事後檢定（Dunn's test with Bonferroni）：\n")
dunn_age = anova_data %>%
  dunn_test(主動攻擊分數 ~ 年齡組, p.adjust.method = "bonferroni")
print(dunn_age)

# 出生地區
kw_region = kruskal.test(主動攻擊分數 ~ 出生地區, data = anova_data)
cat("\n--- 出生地區 ---\n")
cat("χ² =", round(kw_region$statistic, 3), 
    ", df =", kw_region$parameter,
    ", p =", format.pval(kw_region$p.value), "\n")
# 事後檢定（Dunn's test）
cat("\n事後檢定（Dunn's test with Bonferroni）：\n")
dunn_region = anova_data %>%
  dunn_test(主動攻擊分數 ~ 出生地區, p.adjust.method = "bonferroni")
print(dunn_region)

# 教育程度
kw_edu = kruskal.test(主動攻擊分數 ~ 教育程度, data = anova_data)
cat("\n--- 教育程度 ---\n")
cat("χ² =", round(kw_edu$statistic, 3), 
    ", df =", kw_edu$parameter,
    ", p =", format.pval(kw_edu$p.value), "\n")
# 事後檢定（Dunn's test）
cat("\n事後檢定（Dunn's test with Bonferroni）：\n")
dunn_edu = anova_data %>%
  dunn_test(主動攻擊分數 ~ 教育程度, p.adjust.method = "bonferroni")
print(dunn_edu)
# ====================================
# 【3.4】雙因子 ANOVA（Two-way ANOVA）
# ====================================

cat("\n\n【3.4】雙因子 ANOVA\n\n")

# 3.4.1 性別 × 年齡組
cat("--- 3.4.1 性別 × 年齡組 對主動攻擊分數的影響 ---\n")

anova_2way_1 = aov(主動攻擊分數 ~ 性別 * 年齡組, data = anova_data)
cat("\n雙因子 ANOVA 結果：\n")
print(summary(anova_2way_1))

# Type III SS（更適合不平衡設計）
cat("\nType III ANOVA（不平衡設計）：\n")
anova_type3_1 = Anova(anova_2way_1, type = 3)
print(anova_type3_1)

# 交互作用圖
cat("\n繪製交互作用圖...\n")

interaction_data_1 = anova_data %>%
  group_by(性別, 年齡組) %>%
  dplyr::summarise(
    平均分數 = mean(主動攻擊分數),
    標準誤 = sd(主動攻擊分數) / sqrt(dplyr::n()),
    n = dplyr::n(),
    .groups = "drop"
  )

p_interaction_1 = ggplot(interaction_data_1, 
                         aes(x = 年齡組, y = 平均分數, 
                             color = 性別, group = 性別)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c("男" = "#1F77B4", "女" = "#FF7F0E")) +
  labs(
    title = "性別 × 年齡組 交互作用圖",
    subtitle = "主動攻擊分數",
    x = "年齡組",
    y = "平均主動攻擊分數",
    color = "性別"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_interaction_1)

# 3.4.2 性別 × 教育程度
cat("\n--- 3.4.2 性別 × 教育程度 對主動攻擊分數的影響 ---\n")

anova_2way_2 = aov(主動攻擊分數 ~ 性別 * 教育程度, data = anova_data)
cat("\n雙因子 ANOVA 結果：\n")
print(summary(anova_2way_2))

# 交互作用圖
interaction_data_2 = anova_data %>%
  group_by(性別, 教育程度) %>%
  dplyr::summarise(
    平均分數 = mean(主動攻擊分數),
    標準誤 = sd(主動攻擊分數) / sqrt(dplyr::n()),
    n = dplyr::n(),
    .groups = "drop"
  )

p_interaction_2 = ggplot(interaction_data_2, 
                         aes(x = 教育程度, y = 平均分數, 
                             color = 性別, group = 性別)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c("男" = "#1F77B4", "女" = "#FF7F0E")) +
  labs(
    title = "性別 × 教育程度 交互作用圖",
    subtitle = "主動攻擊分數",
    x = "教育程度",
    y = "平均主動攻擊分數",
    color = "性別"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_interaction_2)

# 3.4.3 簡單主效果分析
cat("\n--- 3.4.3 簡單主效果分析 ---\n")
cat("當交互作用顯著時，進行簡單主效果分析\n\n")

# 使用 emmeans 進行簡單主效果
emm_1 = emmeans(anova_2way_1, ~ 性別 | 年齡組)
cat("各年齡組內的性別差異：\n")
print(pairs(emm_1))

emm_2 = emmeans(anova_2way_1, ~ 年齡組 | 性別)
cat("\n各性別內的年齡組差異：\n")
print(pairs(emm_2))

# ====================================
# 【3.6】MANOVA（多變量變異數分析）
# ====================================

cat("\n\n【3.6】MANOVA（多變量變異數分析）\n")
cat("同時檢驗多個依變項的群組差異\n\n")

# 準備 MANOVA 資料
manova_data = data %>%
  dplyr::select(
    主動攻擊分數, 被動攻擊分數, 攻擊敏感度指數,
    同理心程度, 生活滿意度,
    性別, 年齡組, 教育程度
  ) %>%
  na.omit()

cat("MANOVA 樣本數：", nrow(manova_data), "\n\n")

# ----- 3.6.5 MANOVA 前提假設檢定 -----
cat("\n--- 3.6.5 MANOVA 前提假設檢定 ---\n")

# Box's M 檢定（共變異數矩陣同質性）
cat("\nBox's M 檢定（共變異數矩陣同質性）：\n")
cat("注意：Box's M 對常態性非常敏感，大樣本時常顯著\n")

# 使用 rstatix 的 box_m
box_result = box_m(manova_data[, c("主動攻擊分數", "被動攻擊分數", 
                                   "攻擊敏感度指數", "同理心程度", 
                                   "生活滿意度")],
                   manova_data$性別)
print(box_result)

# 多變量常態性（Mardia's test）
cat("\n多變量常態性）：\n")
library(psych)

mvn_vars = manova_data[, c("主動攻擊分數", "被動攻擊分數", 
                           "攻擊敏感度指數", "同理心程度", 
                           "生活滿意度")]

# Mardia 檢定
mardia_result = mardia(mvn_vars, plot = T)

# ----- 3.6.1 性別對多變項的影響 -----
cat("--- 3.6.1 性別對多變項的影響 ---\n")

# 建立依變項矩陣
dv_matrix_sex = cbind(
  manova_data$主動攻擊分數,
  manova_data$被動攻擊分數,
  manova_data$攻擊敏感度指數,
  manova_data$同理心程度,
  manova_data$生活滿意度
)
colnames(dv_matrix_sex) = c("主動攻擊", "被動攻擊", "攻擊敏感度", "同理心", "生活滿意度")

# MANOVA
manova_sex = manova(dv_matrix_sex ~ 性別, data = manova_data)

cat("\nMANOVA 結果（Pillai's Trace）：\n")
print(summary(manova_sex, test = "Pillai"))

cat("\nMANOVA 結果（Wilks' Lambda）：\n")
print(summary(manova_sex, test = "Wilks"))

cat("\nMANOVA 結果（Hotelling-Lawley）：\n")
print(summary(manova_sex, test = "Hotelling-Lawley"))

cat("\nMANOVA 結果（Roy's Greatest Root）：\n")
print(summary(manova_sex, test = "Roy"))

# 單變項 ANOVA（後續檢定）
cat("\n單變項 ANOVA（分解 MANOVA）：\n")
print(summary.aov(manova_sex))

# ----- 3.6.2 年齡組對多變項的影響 -----
cat("\n--- 3.6.2 年齡組對多變項的影響 ---\n")

manova_age = manova(dv_matrix_sex ~ 年齡組, data = manova_data)

cat("\nMANOVA 結果（Pillai's Trace）：\n")
print(summary(manova_age, test = "Pillai"))

cat("\nMANOVA 結果（Wilks' Lambda）：\n")
print(summary(manova_age, test = "Wilks"))

cat("\n單變項 ANOVA：\n")
print(summary.aov(manova_age))

# ----- 3.6.3 教育程度對多變項的影響 -----
cat("\n--- 3.6.3 教育程度對多變項的影響 ---\n")

manova_edu = manova(dv_matrix_sex ~ 教育程度, data = manova_data)

cat("\nMANOVA 結果（Pillai's Trace）：\n")
print(summary(manova_edu, test = "Pillai"))

cat("\n單變項 ANOVA：\n")
print(summary.aov(manova_edu))

# ----- 3.6.4 雙因子 MANOVA -----
cat("\n--- 3.6.4 雙因子 MANOVA（性別 × 年齡組）---\n")

manova_2way = manova(dv_matrix_sex ~ 性別 * 年齡組, data = manova_data)

cat("\nMANOVA 結果（Pillai's Trace）：\n")
print(summary(manova_2way, test = "Pillai"))

cat("\n單變項 ANOVA：\n")
print(summary.aov(manova_2way))

# ====================================
# 第5部分：聚類分析（Cluster Analysis）使用 kproto 方法
# ====================================
cat("\n\n【5】聚類分析：參與者原型\n")

# 準備聚類數據
cluster_data = data %>%
  dplyr::select(
    被動攻擊分數,
    主動攻擊分數,
    攻擊敏感度指數,
    同理心程度,
    生活滿意度
  ) %>%
  na.omit()
str(cluster_data)
# 將同理心程度,生活滿意度轉為因子（類別變數），這樣 kproto 才能運作
cluster_data_kproto <- cluster_data %>%
  mutate(
    # 數值欄位標準化
    被動攻擊分數 = scale(被動攻擊分數)[,1],
    主動攻擊分數 = scale(主動攻擊分數)[,1],
    攻擊敏感度指數 = scale(攻擊敏感度指數)[,1],
    # 類別欄位轉為 factor
    同理心程度 = factor(同理心程度),
    生活滿意度 = factor(生活滿意度)
  )
str(cluster_data_kproto)
# 5.1 確定最佳聚類數
cat("【5.1】確定最佳聚類數\n")

# 使用 kproto 計算不同 k 值的 withinss
wss = c()
for(k in 1:10) {
  cat("計算 k =", k, "...\n")
  kp_temp = kproto(cluster_data_kproto, k, nstart = 10, verbose = F)
  wss[k] = sum(kp_temp$withinss)
}

p_elbow = ggplot(data.frame(k = 1:10, wss = wss), aes(x = k, y = wss)) +
  geom_point(size = 3, color = "blue") +
  geom_line(color = "blue", size = 1) +
  geom_vline(xintercept = 3, linetype = "dashed", color = "red", size = 1.2) +
  labs(
    title = "Elbow Method：確定最佳聚類數",
    x = "聚類數 (k)", y = "群內距離平方和 (WSS)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 11)
  )

print(p_elbow)



# 5.2 K-Prototype 聚類 (k=3)
cat("\n【5.2】K-Prototype 聚類 (k=3)\n")
set.seed(123)
kp = kproto(cluster_data_kproto, 3, nstart = 25, verbose = 0)

cat("聚類大小：", table(kp$cluster), "\n")
cat("群內距離平方和（WSS）：", round(sum(kp$withinss), 2), "\n\n")

# ====================================
# 5.3 聚類特徵描述 - 修正版
# ====================================

cat("\n【5.3】聚類特徵描述\n")

# 將聚類結果加回原始資料
cluster_result <- cluster_data_kproto %>%
  mutate(Cluster = kp$cluster,
         # 確保所有欄位都是數值
         同理心程度 = as.numeric(as.character(同理心程度)),
         生活滿意度 = as.numeric(as.character(生活滿意度))
  )

# 計算各聚類的數值變項平均值
cluster_profiles <- cluster_result %>%
  group_by(Cluster) %>%
  summarise(
    n = n(),
    佔比 = paste0(round(n() / nrow(cluster_result) * 100, 1), "%"),
    被動攻擊_M = round(mean(被動攻擊分數), 2),
    被動攻擊_SD = round(sd(被動攻擊分數), 2),
    主動攻擊_M = round(mean(主動攻擊分數), 2),
    主動攻擊_SD = round(sd(主動攻擊分數), 2),
    攻擊敏感度_M = round(mean(攻擊敏感度指數), 2),
    攻擊敏感度_SD = round(sd(攻擊敏感度指數), 2),
    同理心_M = round(mean(同理心程度), 2),
    同理心_SD = round(sd(同理心程度), 2),
    滿意度_M = round(mean(生活滿意度), 2),
    滿意度_SD = round(sd(生活滿意度), 2)
  )

cat("\n【各聚類描述統計】\n")
print(as.data.frame(cluster_profiles))

# 為聚類命名（根據特徵）
cluster_names <- cluster_profiles %>%
  mutate(
    聚類名稱 = case_when(
      主動攻擊_M == max(主動攻擊_M) ~ "高攻擊型",
      主動攻擊_M == min(主動攻擊_M) ~ "低攻擊型",
      TRUE ~ "中間型"
    )
  )

cat("\n【聚類命名建議】\n")
for(i in 1:nrow(cluster_names)) {
  cat(sprintf("聚類 %d (%s, n=%d):\n", 
              cluster_names$Cluster[i], 
              cluster_names$聚類名稱[i],
              cluster_names$n[i]))
  cat(sprintf("  - 被動攻擊: %.2f ± %.2f\n", 
              cluster_names$被動攻擊_M[i], cluster_names$被動攻擊_SD[i]))
  cat(sprintf("  - 主動攻擊: %.2f ± %.2f\n", 
              cluster_names$主動攻擊_M[i], cluster_names$主動攻擊_SD[i]))
  cat(sprintf("  - 攻擊敏感度: %.2f ± %.2f\n", 
              cluster_names$攻擊敏感度_M[i], cluster_names$攻擊敏感度_SD[i]))
  cat(sprintf("  - 同理心: %.2f ± %.2f\n", 
              cluster_names$同理心_M[i], cluster_names$同理心_SD[i]))
  cat(sprintf("  - 生活滿意度: %.2f ± %.2f\n\n", 
              cluster_names$滿意度_M[i], cluster_names$滿意度_SD[i]))
}

# ====================================
# 5.4 聚類特徵視覺化 - 雷達圖風格長條圖
# ====================================

# 準備視覺化資料（標準化後比較）
cluster_viz <- cluster_result %>%
  group_by(Cluster) %>%
  summarise(
    被動攻擊 = mean(被動攻擊分數),
    主動攻擊 = mean(主動攻擊分數),
    攻擊敏感度 = mean(攻擊敏感度指數),
    同理心 = mean(同理心程度),
    生活滿意度 = mean(生活滿意度)
  ) %>%
  tidyr::pivot_longer(-Cluster, names_to = "變項", values_to = "平均值") %>%
  group_by(變項) %>%
  mutate(標準化值 = scale(平均值)[,1]) %>%  # 跨聚類標準化方便比較
  ungroup()

# 繪製分面長條圖
p_cluster_profile <- ggplot(cluster_viz, 
                            aes(x = 變項, y = 標準化值, fill = factor(Cluster))) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(
    values = c("1" = "#2ECC71", "2" = "#F39C12", "3" = "#E74C3C"),
    name = "聚類"
  ) +
  labs(
    title = "各聚類特徵比較（標準化分數）",
    subtitle = "0 = 全體平均 | 正值 = 高於平均 | 負值 = 低於平均",
    x = NULL,
    y = "標準化分數 (Z-score)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold"),
    legend.position = "top"
  )

print(p_cluster_profile)

# ===========================
# Kproto 多 k 值評估指標函數
# ===========================

evaluate_kproto <- function(data, k_range = 1:10, nstart = 10, verbose = TRUE) {
  # 計算 kproto 在不同 k 值下的各種評估指標
  
  library(clustMixType)
  library(cluster)
  library(dplyr)
  
  results <- data.frame(
    k = integer(),
    wss = numeric(),
    silhouette = numeric()
  )
  
  # 計算距離矩陣（用 Gower 距離處理混合型資料）
  dist_matrix <- daisy(data, metric = "gower")
  n <- nrow(data)
  
  for (k in k_range) {
    if (verbose) cat("計算 k =", k, "...\n")
    
    tryCatch({
      # 執行 kproto
      kp <- kproto(data, k, nstart = nstart, verbose = 0)
      
      # 計算指標
      wss <- sum(kp$withinss)
      
      # Silhouette 係數
      if (k == 1) {
        silhouette_score <- NA
      } else {
        sil <- silhouette(kp$cluster, dist_matrix)
        silhouette_score <- mean(sil[, 3])
      }
      
      # 儲存結果
      results <- rbind(results, data.frame(
        k = k,
        wss = round(wss, 2),
        silhouette = round(silhouette_score, 3)
      ))
    }, error = function(e) {
      if (verbose) cat("  警告: k =", k, "計算失敗\n")
    })
  }
  
  if (verbose) cat("計算完成！\n\n")
  
  return(results)
}

cat("開始計算 kproto 評估指標...\n\n")

# 執行函數（計算 k = 1:10）
evaluation_results <- evaluate_kproto(
  data = cluster_data_kproto,
  k_range = 1:10,
  nstart = 10,
  verbose = TRUE
)

# 顯示結果
cat("【評估結果摘要】\n")
print(evaluation_results)

# 找出最佳 k 值
cat("\n【最佳 k 值建議】\n")

# 基於 Silhouette 的建議
if (sum(!is.na(evaluation_results$silhouette)) > 0) {
  best_silhouette <- evaluation_results$k[which.max(evaluation_results$silhouette)]
  cat("根據 Silhouette 係數：k =", best_silhouette, 
      "（Silhouette =", 
      evaluation_results$silhouette[evaluation_results$k == best_silhouette], "）\n")
}

# WSS 的 Elbow 點
cat("\nWSS 值（用於 Elbow 法則）：\n")
print(evaluation_results[, c("k", "wss")])

# ====================================
# 視覺化評估指標
# ====================================

# 準備視覺化資料
eval_long <- evaluation_results %>%
  pivot_longer(
    cols = c("wss", "silhouette"),
    names_to = "metric",
    values_to = "value"
  )


# Silhouette 圖
p_silhouette_detailed <- ggplot(
  evaluation_results[!is.na(evaluation_results$silhouette), ], 
  aes(x = k, y = silhouette)
) +
  geom_point(size = 3, color = "#3498DB") +
  geom_line(color = "#3498DB", linewidth = 1) +
  geom_vline(xintercept = 3, linetype = "dashed", color = "red", size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray", linewidth = 0.8) +
  labs(
    title = "Silhouette 係數：確定最佳聚類數",
    subtitle = "越接近 1 表示聚類越好，<0.25 表示聚類結構弱",
    x = "聚類數 (k)",
    y = "平均 Silhouette 係數"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.text = element_text(size = 11)
  )

print(p_silhouette_detailed)

cat("✓ 評估完成！\n")

# 評估指標
cat("\n【聚類評估指標】\n")
cat("每個聚類的 WSS：", round(kp$withinss, 2), "\n")
cat("總 WSS：", round(sum(kp$withinss), 2), "\n")

# 5.3 聚類特徵 - 準備資料
kp_out <- as.data.frame(kp$cluster)
colnames(kp_out) <- "Cluster"
final <- cbind(cluster_data, kp_out)

# 為了計算平均值，需要將因子轉回數值
cluster_data_numeric = cluster_data %>%
  mutate(同理心程度 = as.numeric(as.character(同理心程度)),
         生活滿意度 = as.numeric(as.character(生活滿意度))
        )

final_numeric = cbind(cluster_data_numeric, kp_out)

cluster_profiles = final_numeric %>%
  group_by(Cluster) %>%
  summarise(
    n = n(),
    被動攻擊_平均 = mean(被動攻擊分數),
    主動攻擊_平均 = mean(主動攻擊分數),
    攻擊敏感度_平均 = mean(攻擊敏感度指數),
    同理心_平均 = mean(同理心程度),
    滿意度_平均 = mean(生活滿意度)
  ) %>%
  mutate(
    聚類名稱 = case_when(
      Cluster == 1 ~ "聚類1",
      Cluster == 2 ~ "聚類2",
      Cluster == 3 ~ "聚類3",
      TRUE ~ paste0("聚類", Cluster)
    )
  )

cat("【5.3】聚類特徵描述：\n")
print(cluster_profiles)

# ====================================
# 第7部分：效應大小與信心區間
# ====================================
cat("\n\n【7】效應大小與信心區間\n")

# 計算性別間的效應大小
cat("【7.1】Cohen's d：性別間的主動攻擊差異\n")
male_data = anova_data %>% filter(性別 == "男") %>% pull(主動攻擊分數)
female_data = anova_data %>% filter(性別 == "女") %>% pull(主動攻擊分數)

cohens_d = (mean(male_data) - mean(female_data)) / 
  sqrt(((length(male_data)-1)*sd(male_data)^2 + 
          (length(female_data)-1)*sd(female_data)^2) / 
         (length(male_data) + length(female_data) - 2))

cat("Cohen's d:", round(cohens_d, 4), "\n")
cat("效應大小解釋：", 
    if(abs(cohens_d) < 0.2) "微小" 
    else if(abs(cohens_d) < 0.5) "小" 
    else if(abs(cohens_d) < 0.8) "中等" 
    else "大", "\n\n")

# 計算信心區間
ci_d = c(cohens_d - 1.96 * sqrt(2/length(male_data) + 2/length(female_data)),
          cohens_d + 1.96 * sqrt(2/length(male_data) + 2/length(female_data)))
cat("95% CI: [", round(ci_d[1], 4), ", ", round(ci_d[2], 4), "]\n\n")

# ====================================
# 第8部分：多重比較修正（Bonferroni）
# ====================================
cat("【7.2】多重比較修正\n")

# 年齡組間成對比較
age_groups = levels(anova_data$年齡組)
n_comparisons = choose(length(age_groups), 2)
bonferroni_alpha = 0.05 / n_comparisons

cat("年齡組配對數：", n_comparisons, "\n")
cat("Bonferroni 調整後 α:", round(bonferroni_alpha, 4), "\n\n")

# kproto分組組間成對比較
kproto_groups = length(unique(final$Cluster))
n1_comparisons = choose(kproto_groups, 2)
bonferroni_alpha1 = 0.05 / n1_comparisons

cat("年齡組配對數：", n1_comparisons, "\n")
cat("Bonferroni 調整後 α:", round(bonferroni_alpha1, 4), "\n\n")

# ====================================
# 第9部分：分群 MCA 分析 (變項與維度相關矩陣版)
# ====================================
cat("\n\n【9】分群 MCA 分析\n")

# --- 9.1 資料準備 ---
cat("【9.1】準備 MCA 資料\n")
pattern_to_remove <- "^(q22|q23|q28|q20|q25|q31_1|q39)"
cols_to_del <- grep(pattern_to_remove, names(complete_data), value = TRUE)

data_for_mca <- complete_data %>%
  dplyr::select(-all_of(cols_to_del)) %>%
  mutate(across(where(is.numeric), as.factor)) %>%
  mutate(
    Cluster = factor(kp$cluster, 
                     levels = c(1, 2, 3), 
                     labels = c("Cluster_1", "Cluster_2", "Cluster_3"))
  )

# --- 9.2 Benzécri 修正函數 ---
calculate_benzecri <- function(mca_obj) {
  eig_table <- as.data.frame(get_eigenvalue(mca_obj))
  eig <- eig_table[, 1]
  J <- ncol(mca_obj$call$X)
  threshold <- 1/J
  idx <- which(eig > threshold)
  if (length(idx) == 0) return(NULL)
  adjusted_lambda <- ((J / (J - 1)) * (eig[idx] - threshold))^2
  prop <- adjusted_lambda / sum(adjusted_lambda) * 100
  cum_prop <- cumsum(prop)
  data.frame(
    維度 = paste0("Dim.", idx),
    原始特徵值 = round(eig[idx], 4),
    修正慣量 = round(adjusted_lambda, 4),
    修正解釋比例 = paste0(round(prop, 2), "%"),
    累積解釋比例 = paste0(round(cum_prop, 2), "%")
  )
}

# --- 9.3 執行分群 MCA 並輸出數據 ---
cat("【9.3】執行 MCA 運算與相關性矩陣輸出\n")

mca_results_list <- list()
target_clusters <- c("Cluster_1", "Cluster_2", "Cluster_3")

for (cluster_name in target_clusters) {
  cat("\n", paste(rep("=", 40), collapse = ""), "\n")
  cat(">>> 正在處理：", cluster_name, " <<<\n")
  
  # 1. 資料過濾
  mca_input <- data_for_mca %>%
    filter(Cluster == cluster_name) %>%
    dplyr::select(-Cluster) %>%
    droplevels() %>%
    na.omit()
  
  keep_cols <- sapply(mca_input, function(x) nlevels(x) >= 2)
  mca_input <- mca_input[, keep_cols]
  
  if (ncol(mca_input) < 2) {
    cat("⚠ ", cluster_name, " 有效變項不足，跳過\n")
    next
  }
  
  # 2. 執行 MCA (不繪圖)
  n_dims <- min(23, ncol(mca_input) - 1)
  res_mca <- MCA(mca_input, ncp = n_dims, graph = FALSE)
  mca_results_list[[cluster_name]] <- res_mca
  
  # 3. 輸出 Benzécri 修正解釋量
  cat("\n[1] Benzécri 修正解釋量：\n")
  print(calculate_benzecri(res_mca))
  
  # 4. 輸出「變項與維度相關矩陣 (Eta2)」
  # 這反映了原始變項對各維度的解釋強度
  cat("\n[2] 原始變項與維度相關矩陣 (η²)：\n")
  eta2_matrix <- as.data.frame(res_mca$var$eta2)
  # 僅顯示前 5 個維度 (或根據需要調整)
  print(round(eta2_matrix[, 1:min(5, ncol(eta2_matrix))], 4))
  
  # 5. 核心變項篩選總結 (η² ≥ 0.3)
  cat("\n[3] 強相關變項篩選清單 (η² ≥ 0.3)：\n")
  cor_mat_df <- eta2_matrix
  cor_mat_df$Variable <- rownames(cor_mat_df)
  
  top_features <- cor_mat_df %>%
    dplyr::select(Variable, matches("Dim.1$|Dim.2$|Dim.3$")) %>%
    pivot_longer(cols = starts_with("Dim"), 
                 names_to = "維度", 
                 values_to = "Eta2") %>%
    filter(Eta2 >= 0.3) %>%
    arrange(維度, desc(Eta2))
  
  if (nrow(top_features) > 0) {
    print(as.data.frame(top_features))
  } else {
    cat("該群體無強相關變項 (η² < 0.3)\n")
  }
}

cat("\n", paste(rep("=", 40), collapse = ""), "\n")
cat("✓ MCA 所有群體數據分析完成！\n")
