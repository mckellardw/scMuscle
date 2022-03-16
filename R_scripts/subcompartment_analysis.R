# Libs, setwd #####
library(reticulate) # Package that translates python modules to R
scanorama <- import(module = 'scanorama') # import the python module for scanorama
# library(CytoTRACE)
library(harmony)

library(Matrix)
library(dplyr)
library(Seurat)
# library(loomR)
library(future)
# library(sctransform)
# library(metacell)
library(DoubletFinder)
library(cluster)
library(data.table)

#Pseudotime
# library(monocle)

library(readr)
# library(Rmagic)
library(phateR)

#Plotting
library(ggplot2)
library(patchwork)
library(scales)
library(pals)
library(viridis)
# library(schex)


source("/home/dwm269/DWM_utils/sc_utils/seurat_helpers/seutils.R")
#
# Additional functions #####

FACSPlot <- 
  function(seu, features=NULL, assay='RNA', slot='data', 
           cells = NULL,
           col.by='orig.ident', size.by=0.05, 
           alpha=0.5, cols=NULL, size=0.5){
    if(! features[1] %in% rownames(seu)){
      cat(features[1],' not found!\n')
      return(NULL)
    }
    if(!features[2] %in% rownames(seu)){
      cat(features[2], ' not found!\n')
      return(NULL)
    }
    
    if(is.null(cells)){ # view all cells by default
      cells <- colnames(seu)
    }
    
    out =  data.frame(
      feature1=GetAssayData(seu, assay=assay, slot=slot)[features[1],cells],
      feature2=GetAssayData(seu, assay=assay, slot=slot)[features[2],cells],
      col.by=seu[[col.by]][cells,]
      # size.by = seu[[size.by]]
    ) %>% ggplot() + 
      geom_point(
        alpha=alpha,
        size=size,
        # col='black',
        aes(x=feature1,
            y=feature2,
            col=col.by
        )
      ) + 
      theme_minimal() + 
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
      guides(colour = guide_legend(override.aes = list(size=5)))+
      xlab(features[1]) + ylab(features[2])
    
    if(!is.null(cols)){
      out <- out + scale_color_manual(values=cols)
    }
    
    return(out)
  }
# Load data and metadata  ####

# load('../scMuscle_v20.RData')

meta <- fread("../scMuscle_metadata_v7.csv")
# 
# # Remove low quality/undesired datasets from the analysis
# meta <- meta[!meta$sample%in%c('Primary_MB','Yng_D0_C'),]

myo.seurat$injury <- factor(
  myo.seurat$injury,
  levels = c(
    "D0","D0.5","D1","D2","D2.5", "D3","D3.5","D4" ,"D5","D7","D10","D21","18hrs Culture"   
  )
)


myo.seurat$material <- myo.seurat$data.dir
for(samp in unique((myo.seurat$data.dir))){
  myo.seurat$material[myo.seurat$data.dir == samp] <- meta$material[meta$data.dir==samp]
}

# Clean up object for Shiny App ####
myo.slim.seurat <- myo.seurat

to.remove <- c( # clean up metadata
    "data.dir.mm10","harmony_snn_res.1.2","RNA_snn_res.0.8",
    "RNA_snn_res.1.2","harmony_res.1.0",
    "harmony_snn50_res.1.2","harmony_50_res.1.2",
    "bbknn_res.0.8","bbknn_res.1.0","bbknn_res.1.2",
    "scanorama_snn50_res.0.8","scanorama_50_res.0.8",
    "scanorama_res.1.2","scanorama_res.0.8",
    "millay.types" ,"millay.age","S.Score",
    "S.Score_G0v2","G2M.Score_G0v2",
    "G0.pos.Score_G0v2","G0.neg.Score_G0v2",
    "Phase_G0v2","sample_cont",
    colnames(myo.seurat@meta.data)[grep(colnames(myo.seurat@meta.data),pattern="mm39")]
  )
for(tmp in to.remove){
  if(tmp %in% colnames(myo.slim.seurat@meta.data)){
    myo.slim.seurat[[tmp]] <- NULL
  }
}

# clean up reductions
to.remove <- c( # clean up reductions
  "pca_RNA", "harmony", "scanorama","umap_RNA", "umap_harmony", "umap_bbknn", "umap_scanorama"
)
for(tmp in to.remove){
  if(!is.null(myo.slim.seurat[[tmp]])){
    myo.slim.seurat@reductions[[tmp]] <- NULL
  }
}

save(myo.slim.seurat,file="robjs/myo_slim_seurat_v3.RData")

#
# Subset to each celltype compartment and re-run UMAP on narrowed dataset #####
DefaultAssay(scMuscle.clean.seurat) <- 'RNA'
Idents(scMuscle.clean.seurat) <- 'harmony_factorIDs'

#     muscle ####
# myo.seurat <- 
#   subset(
#     x= scMuscle.seurat,
#     subset = harmony_superIDs %in% c("Myogenic", "Myonuclei") &
#       bbknn_superIDs %in% c("Myogenic", "Myonuclei") &
#       scanorama_superIDs %in% c("Myogenic", "Myonuclei")
#   )

myo.seurat <- 
  subset(
    x= scMuscle.seurat,
    subset = harmony_factorIDs %in% c(
      "MuSCs", #1 29
      "Myoblasts/Progenitors", #21
      "Myonuclei (Type IIx)", #11
      "Myonuclei (Type IIb)"
    ) &
      bbknn_factorIDs %in% c(
        "MuSCs", #11
        "Myoblasts_1", #7    
        "Myoblasts_2", #10
        "Committed Progenitors", #22
        "Myonuclei (Cell Prep)", #8 
        "Myonuclei (Type IIx)", #16
        "Myonuclei (Type IIb)" #13
      ) &
      scanorama_factorIDs %in% c(
        "MuSCs", 
        "Myoblasts", 
        "Myoblasts (Cultured)",
        "Committed Progenitors", 
        
        "Myonuclei (Nuclei Prep)",
        "Myonuclei (Cell Prep)", 
        "Myonuclei (Oprescu 2020)",
        "Myonuclei (Type IIx)",
        "Myonuclei (Type IIb)",
        
        "Mixed (Nuclei Prep)"
      )
  )

# n.pcs <- npcs(myo.seurat, reduction="harmony")
# myo.seurat <- RunUMAP(
#   object = myo.seurat, 
#   reduction = 'harmony',
#   dims = 1:n.pcs,
#   reduction.name='umap_harmony',
#   reduction.key = 'umapharmony_'
# )
# 
# n.pcs <- npcs(myo.seurat, reduction="scanorama")
# myo.seurat <- RunUMAP(
#   object = myo.seurat, 
#   reduction = 'scanorama',
#   dims = 1:n.pcs,
#   reduction.name='umap_scanorama',
#   reduction.key = 'umapscanorama_'
# )

#     FAPs ####

fap.seurat <- 
  subset(
    x= scMuscle.seurat,
    subset = harmony_factorIDs %in% c(
      "FAPs (Stem)", 
      "FAPs (Pro-remodeling)", 
      "FAPs (Adipogenic)"
    ) &
      bbknn_factorIDs %in% c(
        "FAPs (Stem)", 
        "FAPs (Pro-remodeling)", 
        "FAPs (Adipogenic)"
      ) &
      scanorama_factorIDs %in% c(
        "FAPs (Stem)", 
        "FAPs (Pro-remodeling)", 
        "FAPs (Adipogenic)",
        "FAPs (De Micheli 2020b)"
      )
  )

n.pcs <- npcs(fap.seurat, reduction="harmony")
fap.seurat <- RunUMAP(
  object = fap.seurat,
  reduction = 'harmony',
  dims = 1:n.pcs,
  reduction.name='umap_harmony',
  reduction.key = 'umapharmony_'
)
unique(fap.seurat$harmony_factorIDs)

#     endo ####
endo.seurat <- 
  subset(
    x= scMuscle.seurat,
    subset = harmony_factorIDs %in% c(
      "Endothelial (Capillary)", 
      "Endothelial (Artery)", 
      "Endothelial (Vein)"
    ) &
      bbknn_factorIDs %in% c(
        "Endothelial (Capillary)", 
        "Endothelial (Artery)", 
        "Endothelial (Vein)"
      ) &
      scanorama_factorIDs %in% c(
        "Endothelial (Capillary)", 
        "Endothelial (Mixed)", 
        "Endothelial (Vein)",
        "Endothelial (De Micheli 2020b)"
      )
  )

n.pcs <- npcs(endo.seurat, reduction="harmony")
endo.seurat <- RunUMAP(
  object = endo.seurat,
  reduction = 'harmony',
  dims = 1:n.pcs,
  reduction.name='umap_harmony',
  reduction.key = 'umapharmony_'
)
unique(endo.seurat$harmony_factorIDs)

#
#     myeloid ####
myeloid.seurat <- 
  subset(
    x= scMuscle.seurat,
    subset = harmony_factorIDs %in% c(
      "Monocyte (Patrolling)","Monocyte (Inflammatory)",
      "Monocyte (Cxcl10+)","M2 Macro. (Cx3cr1_lo)","M2 Macro. (Cx3cr1_hi)",
      "Dendritic"
    ) &
      bbknn_factorIDs %in% c(
        "Monocyte (Patrolling)","Monocyte (Inflammatory)",
        "Monocyte (Cxcl10+)","M2 Macro. (Cx3cr1_lo)","M2 Macro. (Cx3cr1_hi)",
        "Dendritic"
      ) &
      scanorama_factorIDs %in% c(
        "Monocyte (Patrolling)","Monocyte (Inflammatory)","Monocyte (Infiltrating)",
        "Monocyte (De Micheli 2020b)","Monocyte (Mixed)",
        "M2 Macro. (Cx3cr1_lo)","M2 Macro. (Cx3cr1_hi)",
        "Dendritic"
      )
  )

n.pcs <- npcs(myeloid.seurat, reduction="harmony")
myeloid.seurat <- RunUMAP(
  object = myeloid.seurat,
  reduction = 'harmony',
  dims = 1:n.pcs,
  reduction.name='umap_harmony',
  reduction.key = 'umapharmony_'
)
unique(myeloid.seurat$harmony_factorIDs)

#
# Comparison of integration methods ####

harmony.sil <- CalcSilhouette(
  myo.seurat,
  "umap_harmony",
  "source"#,return.df = F
)
rna.sil <- CalcSilhouette(
  myo.seurat,
  "umap_RNA",
  "source"
)
ggplot(as.data.frame(rna.sil)) + 
  geom_violin(
    aes(x=as.character(cluster),y=as.numeric(sil_width))
  ) +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    title = element_text(size=8, face='bold'),
    axis.title.x = element_blank(),
    # axis.title.y = element_blank(),#element_text(size = 8, face='bold', hjust = 0.5, vjust = 0.5),
    axis.text.x = element_text(size = 8, face='bold', hjust = 1, vjust = 0.5, angle = 90),
    axis.text.y = element_text(size = 8, face='bold', hjust = 1, vjust = 0.5),
    axis.ticks = element_blank()
  )

# PHATE ####
# https://github.com/KrishnaswamyLab/phateR/
#     harmony values ####
#       myo ####

# reorder harmony based on std dev
myo.seurat@reductions$harmony@cell.embeddings <- 
  myo.seurat@reductions$harmony@cell.embeddings[,myo.seurat@reductions$harmony@stdev %>% order(decreasing = T)]

# subset out cells that don't line up
remove.cells <- c(
  "oprescu_D2_ACTTTCACATAGGAGC", "oprescu_D2_CCTGCATAGCCTGGAA", "oprescu_D2_CTTCTCTCACATGAAA", 
  "oprescu_D2_GCACGTGGTACTTCCC", "oprescu_D2_TCTCACGGTGGACCTC", "oprescu_D2_TGGTAGTAGTGAACAT", 
  "oprescu_D2_TTACGTTAGCGTGAGT"
)
myo.seurat <- subset(
  myo.seurat,
  cells = Cells(myo.seurat)[!Cells(myo.seurat)%in%remove.cells]
)

# select high stddev dimensions

# n.pcs <- npcs(myo.seurat, reduction = "harmony")
# n.pcs = 32

harmony.dims <- myo.seurat@reductions$harmony@stdev %>% order(decreasing = T)
harmony.dims <- harmony.dims[1:35]

# run PHATE
tmp.phate <- phate(
  myo.seurat@reductions$harmony@cell.embeddings[,1:39], #cells x reduc dims
  gamma=1,
  n.jobs = 50
) 
colnames(tmp.phate$embedding) <- paste0("phate_harmony_", 1:2)
tmp.phate$params$data <- NULL

tmp.phate$embedding[,1] <- tmp.phate$embedding[,1]*-1 #flip x-axis

myo.seurat[["phate_harmony"]] <- CreateDimReducObject(
  embeddings = tmp.phate$embedding,
  key = "phate_harmony_", 
  assay = 'RNA',
  misc = tmp.phate$params
)

myo.seurat@reductions$phate_harmony@stdev <- 
  apply(myo.seurat@reductions$phate_harmony@cell.embeddings, 2, sd) #find std dev for phate vals


DimPlot(
  myo.seurat,
  cells=sample(colnames(myo.seurat)),
  reduction='phate_harmony',
  group.by = 'phate1.bins',
  pt.size = 0.5,
  # cols=unique(celltype.colors),
  label=F
) #+ ylim(c(-0.03,max(myo.seurat@reductions$phate_harmony@cell.embeddings[1,])))

#
#       FAPs ####

# reorder harmony based on std dev
fap.seurat@reductions$harmony@cell.embeddings <- 
  fap.seurat@reductions$harmony@cell.embeddings[,fap.seurat@reductions$harmony@stdev %>% order(decreasing = T)]


# select high stddev dimensions
n.pcs <- npcs(fap.seurat, reduction = "harmony")
# n.pcs = 32

harmony.dims <- fap.seurat@reductions$harmony@stdev %>% order(decreasing = T)
# harmony.dims <- harmony.dims[1:n.pcs]

# tmp<-Embeddings(fap.seurat[["harmony"]])[,harmony.dims]
# save(tmp, file="robjs/tmp_fap_harmony.RData")

# run PHATE
tmp.phate <- phate(
  Embeddings(fap.seurat[["harmony"]])[,harmony.dims], #cells x reduc dims
  gamma=1,
  n.jobs = 30
) 

# save(tmp.phate, file="robjs/tmp_fap_harmony_phate.RData")

colnames(tmp.phate$embedding) <- paste0("phate_harmony_", 1:2)
tmp.phate$params$data <- NULL

# tmp.phate$embedding[,1] <- tmp.phate$embedding[,1]*-1 #flip x-axis

fap.seurat[["phate_harmony"]] <- CreateDimReducObject(
  embeddings = tmp.phate$embedding,
  key = "phate_harmony_", 
  assay = 'RNA',
  misc = tmp.phate$params
)

fap.seurat@reductions$phate_harmony@stdev <- 
  apply(fap.seurat@reductions$phate_harmony@cell.embeddings, 2, sd) #find std dev for phate vals


DimPlot(
  fap.seurat,
  cells=sample(colnames(fap.seurat)),
  reduction='phate_harmony',
  group.by = 'harmony_factorIDs',
  pt.size = 0.5,
  # cols=unique(celltype.colors),
  label=F
) #+ ylim(c(-0.03,max(myo.seurat@reductions$phate_harmony@cell.embeddings[1,])))

save(fap.seurat,file="robjs/scMuscle_mm10_FAP_v1.RData")

#
#       endo ####

# reorder harmony based on std dev
endo.seurat@reductions$harmony@cell.embeddings <- 
  endo.seurat@reductions$harmony@cell.embeddings[,endo.seurat@reductions$harmony@stdev %>% order(decreasing = T)]


# select high stddev dimensions
n.pcs <- npcs(endo.seurat, reduction = "harmony")
# n.pcs = 32

harmony.dims <- endo.seurat@reductions$harmony@stdev %>% order(decreasing = T)
# harmony.dims <- harmony.dims[1:n.pcs]

tmp<-Embeddings(endo.seurat[["harmony"]])[,harmony.dims]
save(tmp, file="robjs/tmp_endo_harmony.RData")

# run PHATE
tmp.phate <- phate(
  Embeddings(endo.seurat[["harmony"]])[,harmony.dims], #cells x reduc dims
  gamma=1,
  n.jobs = 30
) 

# save(tmp.phate, file="robjs/tmp_fap_harmony_phate.RData")

colnames(tmp.phate$embedding) <- paste0("phate_harmony_", 1:2)
tmp.phate$params$data <- NULL

# tmp.phate$embedding[,1] <- tmp.phate$embedding[,1]*-1 #flip x-axis

endo.seurat[["phate_harmony"]] <- CreateDimReducObject(
  embeddings = tmp.phate$embedding,
  key = "phate_harmony_", 
  assay = 'RNA',
  misc = tmp.phate$params
)

endo.seurat@reductions$phate_harmony@stdev <- 
  apply(endo.seurat@reductions$phate_harmony@cell.embeddings, 2, sd) #find std dev for phate vals


DimPlot(
  endo.seurat,
  cells=sample(colnames(endo.seurat)),
  reduction='phate_harmony',
  group.by = 'harmony_factorIDs',
  pt.size = 0.5,
  # cols=unique(celltype.colors),
  label=F
) #+ ylim(c(-0.03,max(myo.seurat@reductions$phate_harmony@cell.embeddings[1,])))

save(endo.seurat,file="robjs/scMuscle_mm10_endo_v1.RData")

#
#       myeloid ####

# reorder harmony based on std dev
myeloid.seurat@reductions$harmony@cell.embeddings <- 
  myeloid.seurat@reductions$harmony@cell.embeddings[,myeloid.seurat@reductions$harmony@stdev %>% order(decreasing = T)]


# select high stddev dimensions
n.pcs <- npcs(myeloid.seurat, reduction = "harmony")
# n.pcs = 32

harmony.dims <- myeloid.seurat@reductions$harmony@stdev %>% order(decreasing = T)
# harmony.dims <- harmony.dims[1:n.pcs]

tmp<-Embeddings(myeloid.seurat[["harmony"]])[,harmony.dims]
save(tmp, file="robjs/tmp_myeloid_harmony.RData")

# run PHATE
tmp.phate <- phate(
  Embeddings(myeloid.seurat[["harmony"]])[,harmony.dims], #cells x reduc dims
  gamma=1,
  n.jobs = 30
) 

# save(tmp.phate, file="robjs/tmp_fap_harmony_phate.RData")

colnames(tmp.phate$embedding) <- paste0("phate_harmony_", 1:2)
tmp.phate$params$data <- NULL

# tmp.phate$embedding[,1] <- tmp.phate$embedding[,1]*-1 #flip x-axis

myeloid.seurat[["phate_harmony"]] <- CreateDimReducObject(
  embeddings = tmp.phate$embedding,
  key = "phate_harmony_", 
  assay = 'RNA',
  misc = tmp.phate$params
)

myeloid.seurat@reductions$phate_harmony@stdev <- 
  apply(myeloid.seurat@reductions$phate_harmony@cell.embeddings, 2, sd) #find std dev for phate vals


DimPlot(
  myeloid.seurat,
  cells=sample(colnames(myeloid.seurat)),
  reduction='phate_harmony',
  group.by = 'harmony_factorIDs',
  pt.size = 0.5,
  raster=F,
  # cols=unique(celltype.colors),
  label=F
) #+ ylim(c(-0.03,max(myo.seurat@reductions$phate_harmony@cell.embeddings[1,])))

save(myeloid.seurat,file="robjs/scMuscle_mm10_myeloid_v1.RData")

#
# Edit Cluster IDs ####
#   DotPlot ####
all.genes <- list(
  Nonmyogenic=c( #Non-myogenic
    'Ptprc',
    'C1qa', # Anti-inflammatory Macrophage
    'Ccl6', # Pro-inflammatory Macrophage
    'Cd74',
    'Pecam1', #Endothelial
    'Pdgfra', # FAPs
    'Gsn' # FAPs
  ),
 CellCycle= c( #Cell Cycle
    'Gas1',
    'Hes1',
    
    'Btg1',
    'Btg2',
    
    'Cdk1',
    'Cdkn1a',
    'Ccnd1',
    'Ccne1',
    'E2f1',
    'Brca1'
  ),
 Myogenic= c( #Myogenic
    'Cd34',
    'Itga7',
    'Vcam1',
    'Pax3',
    'Pax7', #MuSC
    'Myf5', #MB
    'Myod1',
    'Myog',
    'Tnnt2',
    'Acta1',
    'Myh1'
  )
)
tmp<- list()
for(i in 1:length(all.genes)){
  tmp[[i]] <- DotPlot(
  myo.seurat,
  group.by='harmony_res.0.5',
  assay='RNA',
  # col.min = 0.5,
  # dot.min = 0.25,
  #scale.by = 'size',
  cols = c('light gray', '#B31B1B'),
  # cols=c('light blue','dark blue'),
  # col.max = 4,
  dot.scale = 8,
  features=all.genes[[i]][length(all.genes[[i]]):1]
  )  +
    # coord_flip() + 
    labs(title=names(all.genes)[i]) +
    theme(
      title=element_text(hjust=0.5),
      axis.line = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.title.x = element_blank(), 
      axis.title.y = element_blank(),
      axis.text.x=element_text(angle = 90,hjust = 1,vjust= 0.5),
      legend.title = element_blank(),
      panel.grid.major = element_line(colour = "gray", size = 0.5)
    ) #+ DarkTheme()
}
wrap_plots(tmp)
#
#       harmony res.0.3 ####
myo.seurat <- AddCellTypeIdents(
  object = myo.seurat, #TODO rewrite this to speed it up...
  old.idents = 'harmony_res.0.3',
  newName = 'harmony_res.0.3_IDs',
  new.idents = c(
    "", #1
    "", #2
    "", #3
    "Myonuclei", #4
    "", #5
    "", #6
    "Myog. Progenitors", #7
    "Myoblasts", #8
    "", #9
    "Non-myogenic" #10
  )
)
# Order clusters
myo.seurat$harmony_res.0.3_IDs <- factor(
  x = myo.seurat$harmony_res.0.3_IDs, 
  levels = c(
    
  )
)
#       harmony res.0.5 ####
myo.seurat <- AddCellTypeIdents(
  object = myo.seurat, #TODO rewrite this to speed it up...
  old.idents = 'harmony_res.0.5',
  newName = 'harmony_res.0.5_IDs',
  new.idents = c(
      "", #1
      "", #2
      "Myonuclei", #3
      "", #4
      "", #5
      "", #6
      "", #7
      "Myoblasts", #8
      "Myog. Progenitors", #9
      "", #10
      "", #11
      "", #12
      "", #13
      "" #14
  )
)
# Order clusters
myo.seurat$harmony_res.0.5_IDs <- factor(
  x = myo.seurat$harmony_res.0.5_IDs, 
  levels = c(
    
  )
)
# (myo) DGE across clusters ####
plan("multiprocess", workers = 10)
plan()

options(future.globals.maxSize = 1000 * 1024^2)

Idents(myo.seurat) <- "harmony_factorIDs"
myo.markers <- FindAllMarkers(
  myo.seurat,
  verbose = T
)

myo.markers.filt = myo.markers[myo.markers$p_val_adj<10^-5 &
                                 myo.markers$avg_logFC>0.5,]

ggVolcano( myo.markers.filt[myo.markers.filt$cluster=="Myoblasts/Progenitors",])


tmp <-  myo.mon.markers[myo.mon.markers$cluster==20,]
tmp <- tmp[tmp$p_val_adj>0,]
# tmp <- tmp[tmp$pct.2<0.2,]
tmp <- tmp[tmp$avg_logFC<(-1),]
tmp <- tmp[1:25,"gene"]

myo.seurat %>% VlnPlot(
  idents = c("Myoblasts/Progenitors"),
  features = tspns,
  slot = "data",
  pt.size = 0.001,
  group.by="age.months",
  # cols = celltype.colors[celltypes %in% unique(myo.seurat$harmony_factorIDs)],
  combine = F
) %>% lapply(
  function(x) x + 
    vln.theme + 
    theme(
      # plot.title = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x=element_blank(),
      axis.title.y=element_text(size=big.font, color="black", face="bold"),
      axis.text.y=element_text(size=small.font, color="black"),
      plot.margin = unit(rep(0,4),units="inches")
    )+
    # NoLegend()+
    scale_y_continuous(expand=c(0,0))
) %>% wrap_plots( guides="collect")

####
FACSPlot(
  myo.seurat,
  features = c("Cd82","Pax7"),
  col.by = "harmony_factorIDs",
  alpha=1,
  cols = celltype.colors[celltypes %in% unique(myo.seurat$harmony_factorIDs)]
)


#       transcriptomic heterogeneity ####
ggplot(
  myo.seurat@meta.data[sample(colnames(myo.seurat)),],
  aes(
    x=pseudotime,
    y=nCount_RNA/nFeature_RNA
  )
)+
  geom_point(
    aes(col=harmony_factorIDs), #DF.pANN.individual
    alpha=0.5
  )+
  geom_smooth(col="black")+
  scale_color_manual(values=celltype.colors[celltypes %in% unique(myo.seurat$harmony_factorIDs)]) +
  # xlim(c(0,8.1))+
  labs(
    x="Pseudotime",
    y="nCount / nFeature"
  )

ggplot(
  myo.seurat@meta.data[sample(colnames(myo.seurat)),],
  aes(
    x=phate1.bins,
    y=nFeature_RNA/Sequencing.Saturation
  )
)+
  geom_boxplot(outlier.alpha = 0.5)+
  labs(
    x="harmony_PHATE_1 (bins)",
    y="nFeature / Sequencing Saturation"
  )+
  vln.theme+
  theme(
    axis.title = element_text(color="black", face="bold"),
    axis.ticks.x = element_line(color="black")
  )

#         binning ####

#     on phate_harmony values
nbins = 25
breaks = seq(
  min(myo.seurat@reductions$phate_harmony@cell.embeddings[,1]),
  max(myo.seurat@reductions$phate_harmony@cell.embeddings[,1]),
  diff(range(myo.seurat@reductions$phate_harmony@cell.embeddings[,1]))/nbins
)
binnames = as.factor(1:nbins)

myo.seurat$phate1.bins <- cut(
  myo.seurat@reductions$phate_harmony@cell.embeddings[,1]+10^-100,
  breaks = breaks,
  labels=binnames
)
myo.seurat$phate1.bins[is.na(myo.seurat$phate1.bins)] <- 1

ggplot(
  myo.seurat@meta.data,
  aes(x=phate1.bins)
)+
  geom_bar()+
  scale_y_log10(expand=c(0,0))+
  theme_minimal()

#         DGE along phate1 ####

plan("multiprocess", workers = 30)
options(future.globals.maxSize = 1000 * 1024^2)

Idents(myo.seurat) <- "phate1.bins"
myo.phate.bin.markers <- FindAllMarkers(
  myo.seurat, 
  verbose=T
)

myo.phate.bin.markers$is.surface.protein <- myo.phate.bin.markers$gene %in% cspa.surface.proteins$`ENTREZ gene symbol`
myo.phate.bin.markers$is.trascr.factor <- myo.phate.bin.markers$gene %in% tf.genes

write.csv(myo.phate.bin.markers, file="myo_phate_markers_v2.csv")

#
#         DGE between diss. and native differentiation ####

myo.seurat$bin_IDs <- levels(myo.seurat$harmony_factorIDs)[as.numeric(myo.seurat$harmony_factorIDs)]
# myo.seurat$bin_IDs[colnames(myo.slim.seurat)] <- paste0("myo_PHATE_", myo.slim.seurat$phate1.bins)


myo.seurat$bin_IDs[ colnames(myo.seurat)[myo.seurat$phate1.bins %in% 1:3] ] <- "Dissociated"
myo.seurat$bin_IDs[ colnames(myo.seurat)[myo.seurat$phate1.bins %in% 4:5] ] <- "Quiescent_MuSCs"
myo.seurat$bin_IDs[ colnames(myo.seurat)[myo.seurat$phate1.bins %in% 6:7] ] <- "Activated_MuSCs"
myo.seurat$bin_IDs[ colnames(myo.seurat)[myo.seurat$phate1.bins %in% 8:10] ] <- "Committed_Myoblasts"
myo.seurat$bin_IDs[ colnames(myo.seurat)[myo.seurat$phate1.bins %in% 11:18] ] <- "Fusing_Myocytes"
myo.seurat$bin_IDs[myo.seurat$bin_IDs %in% c(
  "MuSCs", "Myoblasts/Progenitors"
) ] <- "Spurious"

myo.seurat$bin_IDs <- factor(
  myo.seurat$bin_IDs,
  levels=c(
    "Quiescent_MuSCs",
    "Activated_MuSCs",
    "Committed_Myoblasts",
    "Fusing_Myocytes",
    "Myonuclei (Type IIx)",
    "Myonuclei (Type IIb)",
    "Dissociated",
    "Spurious"
  )
)

library(future)
plan("multiprocess", workers = 15)
# options(future.globals.maxSize = 1000 * 1024^2)

Idents(myo.seurat) <- bin_IDs
diss.markers <- FindMarkers(
  myo.seurat,
  group.by = "bin_IDs",
  ident.1 = "Dissociated",
  ident.2 = "Activated_MuSCs"
)
diss.markers <- diss.markers[
  diss.markers$p_val_adj<10^-20 &
    # diss.markers$p_val_adj > 0 &
    abs(diss.markers$avg_logFC)>0.75 
  # diss.markers$pct.1 > 0.5
  ,]

ggVolcano_v2(
  diss.markers,
  neg.log.pval.Thresh = 10,
  logFC_filter = 1
  )

VlnPlot(
  myo.seurat,
  group.by = "bin_IDs",
  pt.size = 0,
  # idents = c("Dissociated", "Activated MuSCs"),
  features = c(
    "percent.hsp",
    "Hspa1a",
    "Hsph1"
  )
)

#         binned vln plots ####
#             gene lists ####
# load TFs
# tf.table <- read.csv("/workdir/dwm269/muscle_data/tf_resources/GO_term_summary_20191104_112842_TFlist.csv")
# tf.genes <- unique(tf.table$Symbol)
# length(tf.genes) # 1029 genes

# http://www.informatics.jax.org/go/term/GO:0003700
# GO Term: "DNA-binding transcription factor activity"
tf.table <- fread("/workdir/dwm269/muscle_data/tf_resources/GO_term_summary_20201117_094707.txt")
tf.genes <- unique(tf.table$`MGI Gene/Marker ID`)
length(tf.genes) # 1362 genes

# Load cell surface markers list 
# https://wlab.ethz.ch/cspa/#abstract
cspa.surface.proteins <- fread("gene_list_resources/mouse_cspa_surfaceome_proteins.csv")
surface.genes <- unique(cspa.surface.proteins$`ENTREZ gene symbol`[cspa.surface.proteins$`CSPA category`=="1 - high confidence" ])
#
#             plotting ####
myo.phate.bin.markers.filt <- myo.phate.bin.markers[
  myo.phate.bin.markers$p_val_adj<10^-10 &
    # myo.phate.bin.markers$p_val_adj > 0 &
    # abs(myo.phate.bin.markers$avg_logFC)>0.5 
    myo.phate.bin.markers$avg_logFC>0.5   # only "up-regulated" genes
    # myo.phate.bin.markers$cluster %in% c(8:18)
    # !myo.phate.bin.markers$cluster %in% c(1:2) 
    # !myo.phate.bin.markers$cluster %in% c(19:25)
    # myo.phate.bin.markers$pct.2 < 0.3
  ,]

myo.phate.bin.markers.filt <- myo.phate.bin.markers.filt[order(myo.phate.bin.markers.filt$p_val_adj),]

# DWM
tmp.genes <- c("Purb", "Plaur", "Cd44", "Angptl4","Slc3a2", "Adm","Cd200","Anxa2", "Mest","Nctc1")

#intermediate cells
c("Ccnd1","S100a1", "Sepw1", "Mest", "Dclk1", "Lgals3", "Sparc", "Hes6", "Chrng", "Spg21", "Lgals1")

#canonical
tmp.genes <- c("Pax7", "Myf5","Myod1","Cdkn1c", "Myog","Gm7325", "Tmem8c", "Acta1", "Myh4")

# top genes by p value
tmp.genes <- unique(myo.phate.bin.markers.filt$gene)[1:20]

# genes that occur most often
tmp.genes <- names(sort(table(myo.phate.bin.markers.filt$gene),T))[1:10]

# TFs
tmp.genes <- as.vector(tf.genes[tf.genes %in% myo.phate.bin.markers.filt$gene])#[51:82]
tmp.genes <- c("Cebpb", "Purb", "Hes6", "Klf6", "Mef2a", "Mef2c") #committed
tmp.genes <- c("Atf5","Fosl1") #diss.

#surface markers
tmp.genes <- surface.genes[surface.genes %in% myo.phate.bin.markers.filt$gene]#[1:12]

# lincRNA
lincRNA.genes = read.csv("gene_list_resources/lincRNA_genes.csv")
tmp.genes <- lincRNA.genes$x[lincRNA.genes$x %in% myo.phate.bin.markers.filt$gene]#[1:10]

# all ncRNA
ncRNA.info = read.csv("gene_list_resources/ncRNA_info.csv")
ncRNA.genes = ncRNA.info$mgi_symbol %>% as.vector()
tmp.genes <- unique(ncRNA.genes[ncRNA.genes %in% myo.phate.bin.markers.filt$gene])#[11:20]

tmp.genes <- ncRNA.genes[ncRNA.genes %in% myo.phate.bin.markers.filt$gene[myo.phate.bin.markers.filt$cluster %in% c(7:12)]]#[1:10]

# diss. v. native markers
tmp.genes <- rownames(diss.markers)[1:10]

# quiescence markers 
tmp.genes <- unique(myo.phate.bin.markers.filt$gene[myo.phate.bin.markers.filt$cluster %in% c(4,5)])[1:7]

# tmp.genes<- c("Gm7325", "Tmem8c","Meg3","Gm10073", "Nctc1")

# plot
tmp.vlns <- myo.seurat %>% VlnPlot(
  features = tmp.genes,
  slot = "data",
  pt.size = 0,
  group.by="phate1.bins", #
  # cols = rainbow(25)%>%saturation(values=0.75) %>% brightness(values=0.9) %>% as.vector(),
  combine = F
) %>% lapply(
  function(x) x +
    vln.theme +
    theme(
      plot.title = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x=element_blank(),
      axis.title.y=element_text(size=big.font, color="black", face="bold", angle=0,vjust=0.5),
      axis.text.y=element_text(size=small.font, color="black"),
      plot.margin = unit(rep(0,4),units="inches")
    )+
    NoLegend()+
    scale_y_continuous(expand=c(0,0))
)
for(i in 1:length(tmp.vlns)){
  tmp.vlns[[i]] <-  tmp.vlns[[i]] + labs(y=tmp.genes[i])
}
tmp.vlns[[length(tmp.vlns)]] <- tmp.vlns[[length(tmp.vlns)]] +
  theme(
    axis.text.x = element_text(color="black"),
    axis.title.x = element_text(color="black", size=big.font, face="bold"),
    axis.line.x = element_line(color="black"),
    axis.ticks.x = element_line(color="black")
  )+
  scale_y_continuous(expand=c(0.1,0))+
  xlab("PHATE_Harmony_1")

wrap_plots(tmp.vlns, ncol=1, heights = c(rep(1,length(tmp.vlns)-1),1.1))

#
# (myo) sample diversity in each bin ####
tmp <- sc_entropy(
  myo.seurat, 
  group.by="phate1.bins", entropy.on="sample", #?
  out.name="H.bin.sample", weighted=T, norm2one=F, verbose=T
  )

simpson <- lapply(
  split(
    myo.seurat$sample,#[myo.seurat$is.whole.muscle], 
    myo.seurat$phate1.bins #[myo.seurat$is.whole.muscle]
  ), # split by sample
  FUN = function(X) vegan::diversity(
    as.numeric(as.factor(X)), # convert char to numeric
    index="simpson"
    )
) %>% unlist()

shannon <- lapply(
  split(
    myo.seurat$sample,#[myo.seurat$is.whole.muscle],
    myo.seurat$phate1.bins#[myo.seurat$is.whole.muscle]
    ),
  FUN = function(X) vegan::diversity(
    as.numeric(as.factor(X)), # convert char to numeric
    index="shannon"
  )
) %>% unlist()

# THIS!
left.col ="#3f49d9"
right.col="#b52e24"

# Simpson's plot
ggplot()+
  geom_bar(
    data=myo.seurat@meta.data,
    aes(x=phate1.bins),
    fill="light gray"
  )+
  geom_line(
    data=data.frame(
      x=1:25,
      y=10000^(10*(simpson-0.9))
    ),
    aes(x=x, y=y),
    col="black"
  )+
  scale_y_log10(
    name="# Cells",
    expand=c(0,0),
    breaks=c(1,10,100,1000,10000),
    labels = scientific,
    sec.axis = sec_axis( 
      trans = ~(log(., base=10000)/10)+0.9,
      breaks = c(0.9,0.925,0.95,0.975,1),
      name="Simpson Index (SampleID)"
      )
  )+
  vln.theme + 
  theme(
    axis.ticks.x = element_line(color="black"),
    axis.title.x = element_text(color="black", face="bold", size=14),
    axis.text.x =  element_text(color="black", size=12),
    axis.title.y.left = element_text(color="black", face="bold", size=14),
    axis.text.y.left = element_text(color="black", size=12),
    axis.title.y.right = element_text(color="black", face="bold", size=14),
    axis.text.y.right = element_text(color="black", size=12)
  ) +
  xlab("PHATE_Harmony Bins")



# Shannon plot
ggplot()+
  geom_bar(
    data=myo.seurat@meta.data,
    aes(x=phate1.bins),
    fill="gray"
  )+
  geom_line(
    data=data.frame(
      x=1:25,
      y=10000^((shannon-3)/10)
    ),
    aes(x=x, y=y),
    col="black"
  )+
  scale_y_log10(
    name="Number Cells",
    expand=c(0,0),
    breaks=c(1,10,100,1000,10000),
    labels = scientific,
    sec.axis = sec_axis( 
      trans = ~(log(., base=10000)/10)+3,
      # breaks = c(0.9,0.95,1),
      name="Shannon Index (SampleID)"
    )
  )+
  vln.theme + 
  theme(
    axis.ticks.x = element_line(color="black"),
    axis.title.x = element_text(color="black", face="bold", size=14),
    axis.text.x =  element_text(color="black", size=12),
    axis.title.y.left = element_text(color="black", face="bold", size=14),
    axis.text.y.left = element_text(color="black", size=12),
    axis.title.y.right = element_text(color="black", face="bold", size=14),
    axis.text.y.right = element_text(color="black", size=12)
  ) +
  xlab("PHATE_Harmony Bins")




ggplot(
  myo.seurat@meta.data,
  aes(
    x=as.numeric(phate1.bins),
    color=type
  )
)+
  geom_line(stat="count")+
  scale_y_log10() +
  theme_minimal() +
  theme(
    axis.line = element_line(color="black"),
    axis.text = element_text(color="black",size=12),
    axis.title = element_text(color="black",size=14,face="bold"),
    # legend.position="none",
    legend.title =  element_text(color="black",size=14,face="bold"),
    legend.text =  element_text(color="black",size=12)
  )+
  scale_color_manual(values=type.colors) +
  labs(
    x="PHATE_Harmony Bins",
    y="Count",
    color="SampleID"
  )


myo.seurat$source.label <- as.character(myo.seurat$source)
for(i in 1:length(meta$data.dir)){
  # print(meta$sample[i])
  myo.seurat$sample[myo.seurat$data.dir==meta$data.dir[i]]<- meta$sample[i]
}




#
# (myo) DGE along trajectories ####
myo.seurat$mon_clusters <- clusters(myo.cds, reduction_method = "UMAP")

Idents(myo.seurat) <-  "mon_clusters"
myo.mon.markers <- FindAllMarkers(
  myo.seurat
)
#
# (not myo) DGE across clusters ####

# load("/local/workdir/dwm269/muscle_data/robjs/scMuscle_mm10_myeloid_v1.RData")
# load("/local/workdir/dwm269/muscle_data/robjs/scMuscle_mm10_endo_v1.RData")
# load("/local/workdir/dwm269/muscle_data/robjs/scMuscle_mm10_FAP_v1.RData")

plan("multiprocess", workers = 16)
options(future.globals.maxSize = 5000 * 1024^2)

Idents(fap.seurat) <- "harmony_factorIDs"
fap.markers <- FindAllMarkers(fap.seurat)
# write.csv(fap.markers,file = "scMuscle/supplemental_data/FAP_markers_SupFile.csv")
fap.markers<-fread("scMuscle/supplemental_data/FAP_markers_SupFile.csv",header = T)
fap.markers[fap.markers$avg_log2FC>2 & fap.markers$pct.1>0.5 & fap.markers$pct.2<0.4,]

Idents(endo.seurat) <- "harmony_factorIDs"
endo.markers <- FindAllMarkers(endo.seurat)
# write.csv(endo.markers,file = "scMuscle/supplemental_data/endothelial_markers_SupFile.csv")
endo.markers<-fread("scMuscle/supplemental_data/endothelial_markers_SupFile.csv",header = T)

Idents(myeloid.seurat) <- "harmony_factorIDs"
myeloid.markers <- FindAllMarkers(myeloid.seurat)
# write.csv(myeloid.markers,file = "scMuscle/supplemental_data/myeloid_markers_SupFile.csv")
myeloid.markers<-fread("scMuscle/supplemental_data/myeloid_markers_SupFile.csv",header = T)

# Visualization (myogenic only) ####
#     meta data stuff ####
#       sample distribution for myogenic cells

ggplot(data = myo.seurat@meta.data) +
  geom_bar(aes(x=sample, fill=source))+
  theme_minimal()+ 
  scale_y_continuous(expand = c(0, 1))+
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.title.x = element_blank(),#element_text(size = 12, hjust = 0, vjust = 1), 
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 6,angle = 90,hjust = 1,vjust= 0.5 ),
    axis.text.y=element_text(size = 8),
    axis.ticks = element_blank(),
    legend.title=element_blank()
  )
#
#     VlnPlots ####
VlnPlot(
  # tmp,
  myo.seurat,
  # datasets.seurat.list[[1]],
  features = c(
    # "Pax7",
    # "Myod1",
    # "Myog",
    "percent.mt"
  ),
  group.by ='phate1.bins',split.by = "Phase",
  # split.by = 'DF.classifications',
  # group.by='harmony_res.0.3_IDs',
  # group.by='DF.classifications',
  # cols = as.vector(polychrome())[c(15,16,17,26,32)],
  cols=as.vector(polychrome())[6:36],
  pt.size=0,
  combine=T, ncol = 1
)  #NoLegend() + 
  theme(axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size = 14),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size = 10),
        axis.text.x=element_text(size=10,
                                 angle = 90 ,
                                 hjust = 1,
                                 vjust= 0.5
        ))


#         cell types ####
vln.genes <- c(
  # 'Itga7',
  # 'Vcam1',
  # 'Pax7', #MuSC
  # 'Myf5', 
  # 'Myod1',
  # 'Myog',
  # 'Acta1',
  # 'Myh1'
  "Hdac4"
)
myo.seurat %>% VlnPlot(
  assay = 'RNA',
  # features = vln.genes,
  features = c("Foxo1", "Foxo3", "Foxo4"),
  slot = "data",
  # split.by = 'DF.classifications',
  group.by='phate1.bins', #cols=as.vector(polychrome())[3:36],
  # idents=1:5,
  # group.by='DF.classifications',
  # cols = as.vector(polychrome())[c(15,16,17,26,32)],
  # split.by = "age.months",
  pt.size=0.0001,
  combine=F
)  %>%
  lapply(
    FUN = function(x) x + 
      # NoLegend() + 
      theme(
        # axis.line = element_blank(),
        # panel.border = element_rect(colour = "black", fill=NA, size=1),
        title = element_text(size=8, face='bold'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),#element_text(size = 8, face='bold', hjust = 0.5, vjust = 0.5),
        axis.text.x = element_blank(), #element_text(size = 8, face='bold', hjust = 1, vjust = 0.5, angle = 90),
        axis.text.y = element_text(size = 12, face='bold', hjust = 1, vjust = 0.5, angle = 90)
        # axis.ticks = element_blank()
        # legend.text = element_text(size=small.font)
      ) + 
      scale_y_continuous(expand = c(0, 0)) + 
      scale_color_viridis_d(option="D", aesthetics="fill")
    # ylab("Log-Normalized Expression") 
  ) %>% wrap_plots(ncol=1,guides = 'collect') + 
  theme(axis.text.x=element_text(size = 12, face='bold', hjust = 1, vjust = 0.5, angle = 90))

#     DimPlot ####
DimPlot(
  myo.seurat,
  cells=sample(colnames(myo.seurat)),
  reduction='phate_harmony',
  # group.by = 'harmony_factorIDs', cols=as.vector(polychrome())[3:30],
  # group.by = 'Phase', cols=as.vector(polychrome())[c(13,17,7,6)], #c(13,17,7,6)
  # group.by="cm.comparison", cols=c( "black","light gray"), #"steelblue4","tomato4",
  
  group.by = 'phate1.bins',
  # cols=rainbow(nbins)%>%saturation(values=0.75) %>% brightness(values=0.9) %>% as.vector(),
  # order="Pericytes",
  pt.size = 0.5,
  label=F,label.size = 4,repel = T
) +
  # scale_color_viridis_d(alpha = 0.6,option = "B") #+labs(col="Age") +
  # scale_color_manual(values= rainbow(42)%>%saturation(values = 0.75) %>% brightness(values = 0.9) %>% as.vector()) +
# NoLegend() +
  # umap_theme+
  theme(axis.title = element_text(size=14),legend.text =element_text(size=12), legend.title = element_text(face="bold",size=14) )+
  labs(x="PHATE_harmony_1",y="PHATE_harmony_2") #, color= "Source"

# save(myo.seurat, file="myo_v1.RData")
#
#     FeaturePlot ####
rownames(myo.seurat)[grep(x=rownames(scMuscle.seurat), pattern = 'Hb')] #search object for gene...

DefaultAssay(myo.seurat) <- 'RNA'

FeaturePlot(
  myo.seurat,
  cells = sample(colnames(myo.seurat)),
  reduction='phate_harmony',
  order=F,
  features="Spry1",#  c(tspns,"Myog"),
  # c(
  #   # "Cdkn1a",
  #   "Junb","Fosb","Jun", "Fos",
  #    "Egr1", "Mapk8", "Bcl2", "Trp53"
  #   # "Hspb1", "Cebpb", "Atf3"
  # ),
  combine=FALSE
  # cols = c('lightgrey', '#B31B1B')
) %>% lapply(
  FUN = function(x) x + #NoLegend() +
    # umap_theme +
    xlim(range(myo.seurat@reductions$phate_harmony@cell.embeddings[,1])) +
    ylim(range(myo.seurat@reductions$phate_harmony@cell.embeddings[,2])) +
    scale_color_viridis_c(option = "C",direction = 1)
) %>% wrap_plots()


seu_dotime(
  myo.seurat,
  cells = sample(colnames(myo.seurat)),
  genes=c(emt),
  col.by = "harmony_factorIDs"
)%>% lapply(
  FUN = function(x) x + #NoLegend() +
    # umap_theme +
    scale_color_manual(values=celltype.colors[celltypes %in% unique(myo.seurat$harmony_factorIDs)])
) %>% wrap_plots(nrow=2,guides="collect")

#
#     Dot Plot ####
dot.genes <- c(
  'Ptprc',
  'Ms4a1', #B Cell
  'Cd8a', # T cell 
  'Ccl5',# T cell 
  # 'Trac', #Pan T cell
  'Ctla4', # CD4 T cell
  'Cd4', # CD4 T cell
  # 'Lgals2', # Dendritic
  'C1qa', # Macrophage
  'Cd83', # Macrophage
  'Cd163',
  'Ccl6', #Activated Macrophage
  'Ccl9', # Activated Macrophage
  'Cd74', # Resident Macrophage or Dendritic
  # 'Cxcr3', # Dendritic
  # 'Cd209a', # Dendritic
  # 'Xcr1', # Dendritic
  # 'Tpsab1', # Mast Cell
  'Nkg7', # NK Cell
  'Igkc', # Plasma Cell
  'Hba-a1', # RBCs
  'Pdgfra', # Fibro/Adipogenic progenitor
  'Thy1', 
  'Vwf', # Endothelial
  'Pecam1', #Endothelial
  'Cdh5', # Endothelial
  'Dcn', #fibroblast
  'Tnmd', #Tendon
  'Scx',
  'Col1a1', #Tendon & Fibs
  'Myh11',
  #'Dmd',
  'Mpz' ,#Schwann
  'Ptn', # Neurall/Glial
  'Acta1',
  'Tnnt2',
  'Pax7', #MuSC
  'Myf5', #MB
  'Myod1',
  'Myog',
  # 'Runx1',
  # 'Tsc22d3',
  # 'Tsc22d4',
  # 'Tgif1',
  # 'Arid5b',
  # 'Fosl1',
  # 'Nr3c1',
  # 'Stat6', #End MCA
  # 'Cd34',
  'Itga7',
  # 'Hes1',
  'Vcam1'
  # 'Des'
)

#         scMuscle ####
#             final genes list ####
final.dot.genes <- c(
  'Ptprc',
  'Tnmd',
  'Myh11',
  'Hba-a1',
  'Cd8a',
  'Ccl5',
  'Cd4',
  'Nkg7',
  'Mpz',
  'Acta1',
  'Pax7',
  'Myf5',
  'Myod1',
  'Myog',
  'C1qa',
  'Cd163',
  'Pdgfra',
  'Col1a1',
  'Pecam1',
  'Ms4a1'
)
#             plotting ####
DotPlot(
  myo.seurat,
  group.by='harmony_res.0.3',
  # group.by='scanorama_tf_res.0.3',
  # group.by='harmonyres03IDs',
  # group.by='seurat_clusters',
  # group.by='sample',
  # group.by = 'scanorama_tf_res.0.3',
  # assay = 'scanorama',
  assay='RNA',
  # assay='integrated',
  # col.min = 0,
  # dot.min = 0.25,
  # split.by = 'DF.scanorama.tf',
  # split.by='DF.RNA.merged',
  #scale.by = 'size',
  cols = c('lightgrey', '#B31B1B'),
  # cols=c('red','blue'),
  # col.max = 4,
  # dot.scale = 8,
  features=dot.genes[length(dot.genes):1]
  # features=final.dot.genes[length(final.dot.genes):1]
) +
  theme(axis.line = element_line(colour = "black"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.x=element_text(angle = 90 ,
                                 hjust = 1, 
                                 vjust= 0.5
        )#,
        # axis.text.y=element_text(angle = 45, 
        #                          hjust = 1#, 
        #                          #vjust=-0.5
        # )
  )  + theme(panel.grid.major = element_line(colour = "gray", size = 0.5)) #+ DarkTheme()

#+ DarkTheme()

#
#     Doublet percent in each cluster ####
cell.doublets.df <- melt(data=scMuscle.seurat@meta.data, 
                         id.vars=c('sample', 'chemistry',
                                   'DF.classifications', 'DF.pANN','harmony_res.0.8', 'harmony_res.0.5', 'harmony_res.0.3'),
                         measure.vars=c('DF.classifications'),
                         na.rm=TRUE
)

doublet.perc <- ggplot(cell.doublets.df, aes(x=harmony_res.0.3, y=1, fill=DF.classifications)) + 
  geom_bar(width = 0.9, stat = "identity",position='fill') +
  theme_minimal() +
  # DarkTheme() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.x=element_text(angle = 45 ,
                                 hjust = 1,
                                 vjust= 1, 
                                 size=14),
        axis.text.y=element_text(size = 14)) +
  ylab('% sample') 
doublet.perc

#         RNA values ####
myo.rna.phate <- phate(t(GetAssayData(myo.seurat, assay= 'RNA', slot='data'))) #flip matrix b/c python...
colnames(myo.rna.phate$embedding) <- paste0("phate_rna", 1:2)

#add to seurat object
myo.seurat[["phate_rna"]] <- CreateDimReducObject(embeddings = myo.rna.phate$embedding,
                                                        key = "phate_rna_", 
                                                        assay = 'scanorama',
                                                        misc = myo.phate$params)

#     Volcano Plots ####
ggVolcano(DGEdata = , 
          FC_filter=0.05,
          pvalThresh=10^-50, 
          plotTitle='Volcano Plot', 
          xlim=NULL, ylim=NULL, 
          geneTextSize=0.5)


# end ####




