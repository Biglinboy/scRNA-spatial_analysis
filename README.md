# 单细胞与空间转录组分析代码

本仓库收录论文 *<论文标题>* 中人类肿瘤单细胞（scRNA-seq）与空间转录组（Visium HD）
的分析代码。代码为整理后的代表性版本，按模块编号组织，可用于复现论文主要图件。

## 1. 分析内容

空间转录组（RPCA 整合）
- 空间对象聚类与空间域识别
- RCTD 反卷积（major celltype 与 subtype 两套参考）
- KEGG 通路模块打分与空间可视化
- 重点基因空间表达与共定位
- SCENIC 转录因子活性（上皮 CNAs 与髓系两套）

单细胞（Harmony 整合）
- 8 对样本合并、聚类与细胞类型注释
- 肿瘤与髓系亚群细分
- 配对样本差异表达与多面板火山图
- CellChat 细胞通讯

## 2. 目录结构

```
.
├── analysis.R                     综合分析入口
├── README.md
├── LICENSE
├── .gitignore
├── config/
│   └── config.R                   唯一的路径与参数入口
├── functions/
│   └── utils.R                    通用函数
├── scripts/
│   ├── 00_setup.R                 依赖检查与目录准备
│   ├── 01_spatial_merge_cluster.R 空间聚类与空间域
│   ├── 02_rctd_deconvolution.R    RCTD 反卷积
│   ├── 03_rctd_spatial_plot.R     RCTD 权重空间图
│   ├── 04_pathway_score.R         通路模块打分
│   ├── 05_genelist_spatial.R      基因空间表达与共定位
│   ├── 06_scenic_tf_activity.R    SCENIC 转录因子
│   ├── 07_sc_merge_annotate.R     单细胞聚类与注释
│   ├── 08_sc_subcluster.R         亚群细分
│   ├── 09_sc_diff_volcano.R       配对差异与火山图
│   └── 10_cellchat_interaction.R  细胞通讯
├── metadata/
│   └── sample_sheet.csv           样本与分组对照表
├── data/
│   └── README.md                  数据获取与放置说明
├── results/                       运行产物（不入库）
└── figures/                       成图（不入库）
```

## 3. 环境

- R >= 4.1.2
- Seurat >= 4.3 / SeuratObject >= 4.1
- spacexr >= 2.2.0（RCTD），STdeconvolve >= 1.3.1
- SCENIC >= 1.3.1，SCopeLoomR，AUCell >= 1.20
- CellChat >= 1.6.1，NMF，ggalluvial
- scRNAtoolVis，pheatmap，cowplot，patchwork
- viridis，RColorBrewer，paletteer，reshape2，argparser

版本明细见运行后生成的 `sessionInfo.txt`。

## 4. 快速开始

```bash
git clone https://github.com/<user>/<repo>.git
cd <repo>

# 默认取当前目录为项目根，数据目录为上一级；
# 也可以用环境变量指定：
export BJXW_PROJECT_ROOT=$(pwd)
export BJXW_DATA_ROOT=/path/to/data

Rscript analysis.R            # 全流程
Rscript analysis.R 04 05      # 只跑通路打分与基因空间图
```

## 5. 代码与图对应

| 图 | 脚本 |
|---|---|
| 空间 UMAP、空间域 | 01 |
| RCTD major / subtype 比例图 | 02, 03 |
| 通路得分空间图 | 04 |
| 重点基因空间表达、共定位 | 05 |
| TF rank 图与活性热图 | 06 |
| 单细胞 UMAP、注释、点图 | 07 |
| 肿瘤 / 髓系亚群 | 08 |
| 配对差异多火山图 | 09 |
| 细胞通讯网络与信号流 | 10 |

## 6. 说明

- 所有脚本不含绝对路径，换机器只需改 `config/config.R` 或设环境变量。
- 随机种子固定为 42，关键参数集中在配置文件。
- 参数按论文实际使用的设置整理；聚类编号到细胞类型的映射需按各自数据核对。

## 7. 许可与联系

代码采用 MIT 许可证，见 `LICENSE`。
通讯作者：<name>, <email>
代码维护：<name>, <email>
