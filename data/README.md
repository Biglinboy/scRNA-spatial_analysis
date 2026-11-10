# 数据说明

本目录不存放原始测序数据。请按下述方式获取并在 `config/config.R`
（或环境变量 `BJXW_DATA_ROOT`）中指向数据根目录。

## 1. 数据来源

- scRNA-seq 表达矩阵：GEO GSEXXXXXX
- Visium HD 空间转录组：GEO GSEXXXXXX
- 处理后 Seurat 对象与细胞元数据：Zenodo 10.5281/zenodo.XXXXXXX
- 公共参考：CellChatDB.human、KEGG 基因集、cisTarget 数据库

## 2. 目录约定

```
<data_root>/
├── Seurat/
│   ├── new16Harmony/
│   │   └── merge_rds/merge_data.rds        单细胞合并对象
│   └── RPCA/
│       ├── merge_spatial_cluster.rds       空间合并对象
│       ├── 5.RCTD/CNAs_result.rds          上皮亚群对象
│       ├── 5.RCTD/Myeloid_result.rds       髓系亚群对象
│       └── 14.Scenic/CNAs/result/CNAs_Finished.loom
└── script/
    └── getGoTerm/Human_kegg_gene_sets_symbols.xls
```

## 3. 关键对象说明

- 单细胞对象需含 `harmony` 降维，以及 `seurat_clusters`、`annot`、`group` 三列。
- 空间对象需含 `orig.ident`、`seurat_clusters` 与 `@images`，供 RCTD 取坐标。
- SCENIC loom 需含 `Regulons` 与 `RegulonsAUC` 两个属性。

## 4. 受控数据

若涉及受控访问，请写明申请渠道并提供下采样示例对象，保证代码可被独立检查。
