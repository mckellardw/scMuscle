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
gc() #be sure to clear up memory usage...

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

# Focus on FAPs -> Myogenic cells ----

sources.use = c(
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

scMuscle.cellchat@meta$sample <- scMuscle.seurat$sample[colnames(scMuscle.cellchat@data)]
fap.cellchat <-subsetCellChat(
  scMuscle.cellchat,
  idents.use=c(sources.use,targets.use)
)
fap.cellchat@idents <- factor(
  fap.cellchat@idents,
  levels= c(
    "Quiescent_MuSCs","Activated_MuSCs","Committed_Myoblasts","Fusing_Myocytes",
    "FAPs (Stem)",
    "FAPs (Adipogenic)",
    "FAPs (Pro-remodeling)"
  )
)

# Numbers of different types of interactions between FAPs and myogenic cells
CellChatDB.mouse$interaction$annotation[CellChatDB.mouse$interaction$pathway_name %in% fap.cellchat@netP$pathways] %>%table()

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


#     FAPs->Myo chord plot ####
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
  "Tenocytes",
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

#plot
par(mfrow = c(1,1), xpd=TRUE)
netVisual_chord_gene(
  sub.cellchat, 
  sources.use = sources.use,
  targets.use = targets.use,
  slot.name = "netP",
  thresh=10^-100,
  small.gap=0.25,
  big.gap=5,
  color.use = celltype.colors,
  show.legend = F, #legend.pos.x = 10,
  lab.cex=0.6
  # width=width,
  # height=height
)
dev.off()
# Focus on FAPs ---------------------------------------------------------------
sources.use = c(
  "Tenocytes",
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

# scMuscle.cellchat@meta$sample <- scMuscle.seurat$sample[colnames(scMuscle.cellchat@data)]
fap.cellchat <-subsetCellChat(
  scMuscle.cellchat,
  # cells.use=Cells(scMuscle.seurat)[scMuscle.seurat$harmony_PHATE_IDs %in% c(sources.use,targets.use)]
  idents.use=c(sources.use,targets.use)
) 
fap.cellchat@idents <- factor(
  fap.cellchat@idents,
  levels= c(
    "Quiescent_MuSCs","Activated_MuSCs","Committed_Myoblasts","Fusing_Myocytes",
    "Tenocytes",
    "FAPs (Stem)",
    "FAPs (Adipogenic)",
    "FAPs (Pro-remodeling)"
  )
)

#
#       chord plots for different types of ligand-receptor interactions ----
tmp.pairLR.use <- list()
tmp.pairLR.use[["Cell-Cell Contact"]] <- data.frame( 
  interaction_name=CellChatDB$interaction$interaction_name[CellChatDB$interaction$annotation=="Cell-Cell Contact"]
)
tmp.pairLR.use[["Secreted Signaling"]] <- data.frame( 
  interaction_name=CellChatDB$interaction$interaction_name[CellChatDB$interaction$annotation=="Secreted Signaling"]
)
tmp.pairLR.use[["ECM-Receptor"]] <- data.frame( 
  interaction_name=CellChatDB$interaction$interaction_name[CellChatDB$interaction$annotation=="ECM-Receptor"]
)

par(mfrow = c(1,3), xpd=TRUE)
netVisual_chord_gene(
  fap.cellchat, 
  sources.use = sources.use,
  targets.use = targets.use,
  slot.name = "netP",
  pairLR.use = tmp.pairLR.use$`Cell-Cell Contact`,
  title.name="Cell-Cell Contact",
  thresh=0.01,
  small.gap=0.25,
  big.gap=5,
  color.use = celltype.colors,
  show.legend = F, 
  lab.cex=0.6
)
netVisual_chord_gene(
  fap.cellchat, 
  sources.use = sources.use,
  targets.use = targets.use,
  slot.name = "netP",
  pairLR.use = tmp.pairLR.use$`Secreted Signaling`,
  title.name="Secreted Signaling",
  thresh=0.01,
  small.gap=0.25,
  big.gap=5,
  color.use = celltype.colors,
  show.legend = F, 
  lab.cex=0.6
)
netVisual_chord_gene(
  fap.cellchat, 
  sources.use = sources.use,
  targets.use = targets.use,
  slot.name = "netP",
  pairLR.use = tmp.pairLR.use$`ECM-Receptor`,
  title.name="ECM-Receptor",
  thresh=0.01,
  small.gap=0.25,
  big.gap=5,
  color.use = celltype.colors,
  show.legend = F, 
  lab.cex=0.6
) 

#     Midkine (MK) signaling plots ####
#       LR pair contributions bar plot ----
pathways.show="MK"
netAnalysis_contribution(
  fap.cellchat, 
  signaling = pathways.show,
  title="Midkine Signaling" #pathways.show
)

#       Myo feature plot ----
FeaturePlot(
  myo.seurat,
  cells = sample(colnames(myo.seurat)),
  reduction='phate_harmony',
  features=c("Ncl","Lrp1","Sdc2"),
  combine=FALSE
) %>% lapply(
  FUN = function(x) x + #NoLegend() +
    umap.theme +
    lims(
      x=range(myo.seurat@reductions$phate_harmony@cell.embeddings[,1]),
      y=range(myo.seurat@reductions$phate_harmony@cell.embeddings[,2])
    ) +
    scale_color_viridis_c()+
    theme(
      plot.title=element_text(size=big.font,face="bold",hjust=0.5),
      legend.title = element_text(size=small.font,face="bold",hjust=0.5),
      legend.text = element_text(size=small.font),
      legend.position = "right",
      # legend.direction = "horizontal",
      plot.margin = unit(rep(0,4),"inches"),
      axis.line=element_blank(),axis.title=element_blank()
    )+
    labs(x="PHATE_Harmony_1",y="PHATE_Harmony_2",col="Log-Normalized\nExpression")
) %>%
  wrap_plots(ncol=1)

#       Myogenic cell MDK receptor vln plots ----
tmp.genes <- c("Ncl","Sdc4","Lrp1")

Idents(myo.seurat)<-"bin_IDs"
signal.vlns <- myo.seurat %>% VlnPlot(
  features = tmp.genes,
  idents = c(
    "Quiescent_MuSCs",
    "Activated_MuSCs",
    "Committed_Myoblasts",
    "Fusing_Myocytes"
  ),
  slot = "data",
  pt.size = 0,
  group.by="bin_IDs", 
  # cols = bin.colors,
  cols = celltype.colors[c(
    "Quiescent_MuSCs",
    "Activated_MuSCs",
    "Committed_Myoblasts",
    "Fusing_Myocytes")],
  combine = F
) %>% lapply(
  function(x) x + 
    vln.theme + 
    theme(
      plot.title = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x=element_blank(),
      axis.title.y=element_text(size=big.font, color="black", face="bold"),
      axis.text.y=element_text(size=small.font, color="black"),
      plot.margin = unit(rep(0,4),units="inches")
    )+
    NoLegend()+
    scale_y_continuous(expand=c(0,0))
)

for(i in 1:length(tmp.genes)){
  signal.vlns[[i]] <-  signal.vlns[[i]] + 
    labs(y=tmp.genes[i]) +
    theme(axis.title.y=element_text(face="bold.italic",size=small.font))
}

signal.vlns[[length(signal.vlns)]] <- signal.vlns[[length(signal.vlns)]] + 
  scale_y_continuous(expand = c(0.05,0))+
  theme(
    axis.line.x = element_line(color="black"),
    axis.ticks.x = element_line(color="black"),
    axis.text.x = element_text(color="black")
  )

wrap_plots(signal.vlns, ncol=1)


#       FAP vln plot for Mdk ----

tmp.genes <- c("Mdk")

lig.vlns <- fap.seurat %>% VlnPlot(
  features = tmp.genes,
  slot = "data",
  pt.size = 0,
  group.by="harmony_factorIDs", 
  cols = celltype.colors[celltypes %in% unique(fap.seurat$harmony_factorIDs)],
  combine = F
) %>% lapply(
  function(x) x + 
    vln.theme + 
    theme(
      plot.title = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x=element_blank(),
      axis.title.y=element_text(size=big.font, color="black", face="bold"),
      axis.text.y=element_text(size=small.font, color="black"),
      plot.margin = unit(rep(0,4),units="inches")
    )+
    NoLegend()+
    scale_y_continuous(expand=c(0,0))
)

for(i in 1:length(tmp.genes)){
  lig.vlns[[i]] <-  lig.vlns[[i]] + 
    labs(y=tmp.genes[i]) +
    theme(axis.title.y=element_text(face="bold.italic",size=small.font))
}

lig.vlns[[length(lig.vlns)]] <- lig.vlns[[length(lig.vlns)]] + 
  scale_y_continuous(expand = c(0.05,0))+
  theme(
    axis.line.x = element_line(color="black"),
    axis.ticks.x = element_line(color="black"),
    axis.text.x = element_text(color="black")
  )

wrap_plots(
  wrap_plots(signal.vlns, ncol=1)&theme(axis.text.x=element_blank()),
  plot_spacer(),
  wrap_plots(lig.vlns,ncol=1)&theme(axis.text.x=element_blank()), 
  ncol=1, heights = c(3,0.1,1)
)

#
#       Visium feature plots ----
# See 'helper_functions_v1.R' for this function
visListPlot(
  seu.list = seu.list,
  features = c("Mdk","Ncl","Sdc4","Lrp1"),
  pt.size = 0.2*2,
  font.size = small.font*2
)
