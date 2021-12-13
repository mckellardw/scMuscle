# Script to subset whole Seurat object, which will display by default on the Shiny app
library(Seurat)

load("/local/workdir/dwm269/muscle_data/scMuscle/web_app/data/scMuscle_mm10_slim_v1-1.RData")

sample.rate=0.1
scMuscle.sub.seurat <- subset(
  scMuscle.slim.seurat,
  cells =  sample(
    Cells(scMuscle.slim.seurat),
    size = sample.rate*ncol(scMuscle.slim.seurat)
    )
)

saveRDS(
  scMuscle.sub.seurat,
  file="data/scMuscle_subsampled_mm10_v1-1.rds",
  compress=F
)
