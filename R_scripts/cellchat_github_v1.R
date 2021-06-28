#
# Load libraries ----
library(Matrix)
library(dplyr)
library(Seurat)

library(NMF)
library(CellChat)

# Load data ----

# load in complete scMuscle.seurat object

# CellChat processing ----
# Vignette used to build analysis:
#     https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html

CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# create cellchat object
cells.include <- Cells(scMuscle.seurat)[!is.na(scMuscle.seurat$harmony_PHATE_IDs)]

scMuscle.cellchat <- createCellChat(
  object = GetAssayData(
    scMuscle.seurat,
    assay = "RNA",
    slot = "data"
  )[,cells.include],
  meta = data.frame(
    harmony_PHATE_IDs = as.factor(scMuscle.seurat$harmony_PHATE_IDs)[cells.include],
    row.names = cells.include
  ),
  group.by = "harmony_PHATE_IDs"
)

scMuscle.cellchat@idents <- factor(
  scMuscle.cellchat@idents,
  levels = c(
    "B Cells","NK Cells","T Cells", "Neutrophils",
    "Monocyte (Patrolling)","Monocyte (Inflammatory)", "Monocyte (Cxcl10+)", 
    "M2 Macro. (Cx3cr1_lo)","M2 Macro. (Cx3cr1_hi)","Dendritic",   
    "Endothelial (Capillary)","Endothelial (Artery)", "Endothelial (Vein)","Smooth Muscle & Pericytes",
    "FAPs (Stem)", "FAPs (Pro-remodeling)", "FAPs (Adipogenic)","Tenocytes",   
    "Neural",  
    "Quiescent_MuSCs", "Activated_MuSCs", "Committed_Myoblasts","Fusing_Myocytes",
    "Myonuclei (Type IIx)","Myonuclei (Type IIb)"
  )
)

# Set the database to use (whole CellChat database)
scMuscle.cellchat@DB <- CellChatDB.mouse

# subset the expression data of signaling genes for saving computation cost
scMuscle.cellchat <- subsetData(scMuscle.cellchat) # This step is necessary even if using the whole database

# Check parallelization settings to speed up runtime

# CellChat processing (memory intensive for large datasets)
scMuscle.cellchat <- identifyOverExpressedGenes(scMuscle.cellchat)
scMuscle.cellchat <- identifyOverExpressedInteractions(scMuscle.cellchat)
gc()

# project gene expression data onto PPI network (optional)
scMuscle.cellchat <- projectData(scMuscle.cellchat, PPI.mouse)

# scMuscle.cellchat@idents = droplevels(scMuscle.cellchat@idents, exclude = setdiff(levels(meta$labels),unique(meta$labels)))


scMuscle.cellchat <- computeCommunProb(scMuscle.cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
scMuscle.cellchat <- filterCommunication(scMuscle.cellchat, min.cells = 10)

scMuscle.cellchat <- computeCommunProbPathway(scMuscle.cellchat)

scMuscle.cellchat <- aggregateNet(scMuscle.cellchat)

# Compute the network centrality scores
scMuscle.cellchat <- netAnalysis_computeCentrality(
  scMuscle.cellchat, 
  slot.name = "netP"
) # the slot 'netP' means the inferred intercellular communication network of signaling pathways

plan("sequential") # clears memory usage from future parallelization
gc()

# Visualization ----
#   heatmap ----

sources.use = c(
  "Monocyte (Patrolling)",
  "Monocyte (Inflammatory)",
  "Monocyte (Cxcl10+)",
  "M2 Macro. (Cx3cr1_lo)",
  "M2 Macro. (Cx3cr1_hi)",
  "Endothelial (Capillary)",
  "Endothelial (Artery)",
  "Endothelial (Vein)",
  "Smooth Muscle & Pericytes",
  "FAPs (Stem)",
  "FAPs (Adipogenic)",
  "FAPs (Pro-remodeling)"
)
targets.use = c(
  "Quiescent_MuSCs",
  "Activated_MuSCs",
  "Committed_Myoblasts",
  "Fusing_Myocytes"
)

sub.cellchat <-subsetCellChat(
  scMuscle.cellchat,
  idents.use=c(sources.use,targets.use)
) 
sub.cellchat@idents <- factor(
  sub.cellchat@idents,
  levels= c(
    "Quiescent_MuSCs","Activated_MuSCs","Committed_Myoblasts","Fusing_Myocytes",
    "Monocyte (Patrolling)",
    "Monocyte (Inflammatory)",
    "Monocyte (Cxcl10+)",
    "M2 Macro. (Cx3cr1_lo)",
    "M2 Macro. (Cx3cr1_hi)",
    "Endothelial (Capillary)",
    "Endothelial (Artery)",
    "Endothelial (Vein)",
    "Smooth Muscle & Pericytes",
    "FAPs (Stem)",
    "FAPs (Adipogenic)",
    "FAPs (Pro-remodeling)"
  )
)
groupSize <- as.numeric(table(sub.cellchat@idents))

par(mfrow = c(1,1), xpd=TRUE)
netVisual_heatmap(
  sub.cellchat,
  slot.name="net",
  measure="weight",
  color.use = celltype.colors[sort(c(sources.use,targets.use))],
  sources.use = sources.use,
  targets.use = targets.use,
  # cluster.rows=T, cluster.cols=T,
  color.heatmap = "YlOrRd",
  font.size.title = big.font,
  font.size = small.font,
  width=width,
  height=height
)
dev.off()