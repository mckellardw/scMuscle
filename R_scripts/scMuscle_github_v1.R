




# Load libs #####
library(reticulate) # Package that translates python modules to R
scanorama <- import(module = 'scanorama') # import the python module for scanorama
library(DoubletFinder)
library(KernSmooth)
library(Matrix)
library(dplyr)
library(Seurat)
library(harmony)
library(future)
library(SoupX)
library(phateR)

library(DropletUtils)
library(cluster)
library(parallel)
library(data.table)

# for plotting
library(ggplot2)
library(patchwork)
library(pals)
library(viridis)
library(shades)

source("helper_functions_v1.R") 

# Load metadata
meta <- fread("supplemental_data/sample_metadata.csv")

# Set number of cores to use
NCORES = 1 # 

# Preprocess objects individually ####
#     SoupX ####
# https://github.com/constantAmateur/SoupX
soup.list <- list()
soup.mat.list <- list()

cl <- makeCluster(NCORES)

#load souplist objects from cellranger outputs
data.dir <- list("/list/of/paths/to/cellranger/outputs")
soup.list <- mclapply( 
  data.dir,
  FUN = SoupX::load10X,
  keepDroplets=TRUE,
  mc.cores=NCORES
)

soup.list.est <- mclapply(
  soup.list,
  FUN = function(sc){
    return(tryCatch(autoEstCont(sc), error=function(e) NULL))
  },
  mc.cores=NCORES
)

adj.mat.list <- mclapply(
  soup.list.est,
  FUN = function(sc){
    return(tryCatch(adjustCounts(sc), error=function(e) NULL))
  },
  mc.cores=NCORES
)
stopCluster(cl)

#OPTIONAL: Save adjusted matrices to disk
for(i in 1:length(adj.mat.list)){
  if(!is.null(adj.mat.list[[i]]) & !file.exists(paste0(data.dir[i],"/soupx/matrix.mtx"))){
    DropletUtils:::write10xCounts(
      path="/path/to/cellranger/outputs/", #path to each sample's cellranger count output
      adj.mat.list[[i]]
    )
  }
}

#OPTIONAL: save Rho values
rhos <- list()
for(i in 1:length(soup.list.est)){
  rhos[[i]] <- mean(soup.list.est[[i]]$metaData$rho)
}
rhos <- do.call(rbind,rhos)

#     Read in the count matrices ####
mat.list <- list()
soupx.used <- list()
for(i in 1:length(meta$data.dir)){ 
  # only us SoupX values for samples from whole muscle (i.e., not FACS-sorted)
  if(file.exists(paste0(meta$data.dir[i], '/outs/soupx')) & meta$is.whole.muscle[i]){ 
    cat("Reading #",i, ": ", meta$data.dir[i], ' \n')
    mat.list[[i]] <- Read10X(data.dir = paste0(meta$data.dir[i],"/outs/soupx"))
    soupx.used[[i]] <- T
  }else if(file.exists(paste0(meta$data.dir[i], '/outs/filtered_feature_bc_matrix'))){ # use raw counts for FACS-sorted samples
    cat("Reading #",i, ": ", meta$data.dir[i], ' \n')
    mat.list[[i]] <- Read10X(data.dir = paste0(meta$data.dir[i], '/outs/filtered_feature_bc_matrix'))
    soupx.used[[i]] <- F
  }else{
    cat("Data not found for # ", i, " (", meta$data.dir[i], ")", "\n")
    soupx.used[[i]] <- NULL
  }
}

cat(sum(unlist(lapply(mat.list, ncol))),"cells (total) loaded...\n")

#     Seurat ####
seu.list <- list()

# Initialize seurat objects
cl <- makeCluster(NCORES)
seu.list <- mclapply( 
  mat.list,
  FUN = function(mat){
    return(CreateSeuratObject(
      counts = mat, 
      min.features = 500,
      min.cells = 3,
      project = 'muscle_data'
    ))
  },
  mc.cores = NCORES
)  
stopCluster(cl)

for(i in 1:length(seu.list)){
  cat(' #####################################\n',
      '### Processing dataset number ', i, '###\n',
      '#####################################\n')
  # Add meta data
  for(md in colnames(meta)){
    seu.list[[i]][[md]] <- meta[[md]][i]
  }
  
  # add %MT
  seu.list[[i]][["percent.mt"]]  <- PercentageFeatureSet(seu.list[[i]], pattern = "mt-") 
  
  # Filter out low quality cells according to the metrics defined above
  seu.list[[i]] <- subset(seu.list[[i]],
                           subset = percent.mt < 30 &
                             nCount_RNA > 1000
  )
  # Only mito and floor filtering; trying to find doublets
}
cat((sum(unlist(lapply(mat.list, ncol)))-sum(unlist(lapply(seu.list, ncol)))),"cells (total) removed...\n")

# Preprocess seurat objects
seuPreProcess <- function(seu, assay='RNA', n.pcs=50, res=0.8){
  # NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
  pca.name = paste0('pca_', assay)
  pca.key = paste0(pca.name,'_')
  umap.name = paste0('umap_', assay)
  
  seu = NormalizeData(
    seu
  ) %>% FindVariableFeatures(
    assay = assay,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = F
  ) %>% ScaleData(
    assay = assay
  ) %>% RunPCA(
    assay = assay,
    reduction.name = pca.name,
    reduction.key = pca.key,
    verbose = F,
    npcs = n.pcs
  )
  
  #find pcs to use (custom function, see helper functions script for details)
  n.pcs.use = npcs(seu, reduction = pca.name, var.toal = 0.95)
  
  # FindNeighbors %>% RunUMAP, FindClusters
  seu <- FindNeighbors(
    seu,
    reduction = pca.name,
    dims = 1:n.pcs.use,
    force.recalc = TRUE,
    verbose = FALSE
  ) %>% RunUMAP(
    reduction = pca.name,
    dims = 1:n.pcs.use,
    reduction.name=umap.name
  )
  
  seu@reductions[[umap.name]]@misc$n.pcs.used <- n.pcs.use
  
  seu <- FindClusters(object = seu,resolution = res)
  seu[[paste0('RNA_res.',res)]] <- as.numeric(seu@active.ident)
  
  return(seu)
}

# preprocess each dataset individually
seu.list <- lapply(seu.list, seuPreProcess)

#
#     DoubletFinder on RNA, individual datasets ####

# Mostly taken from the DoubletFinder vignette:
# https://github.com/chris-mcginnis-ucsf/DoubletFinder

# Made some minor adjustments to some of the DoubletFinder functions, see helper functions
#     for details.

# sweep.res.list <- list()
# sweep.stats <- list()
bcmvn <- list()
pK <- list()
homotypic.prop <- list()
nExp_poi <- list()
nExp_poi.adj <- list()

# Estimated Doublet Rate for each dataset
edr <- estimateDoubletRate.DWM(seu.list = seu.list)/100 #use your own known EDR here

for(i in 1:length(seu.list)){
  cat(' ############################################\n',
      '### DoubletFinder for dataset number ', i, '###\n',
      '############################################\n')
  
  ## pK Identification (no ground-truth)
  bcmvn[[i]]<- paramSweep_v3_DWM(
    seu=seu.list[[i]], 
    PCs = 1:seu.list[[i]]@reductions$umap_RNA@misc$n.pcs.used, 
    NCORES = NCORES
  ) %>% summarizeSweep(
    GT = FALSE
  ) %>% find.pK() 
  
  # Pull out max of bcmvn
  pK[[i]] <- as.numeric(as.character(bcmvn[[i]]$pK[bcmvn[[i]]$BCmetric==max(bcmvn[[i]]$BCmetric)])) # ugly, but functional...
  
  ## Homotypic Doublet Proportion Estimate
  homotypic.prop[[i]] <- modelHomotypic(seu.list[[i]]$seurat_clusters) 
  
  nExp_poi[[i]] <- round(edr[[i]]*length(colnames(seu.list[[i]])))  
  nExp_poi.adj[[i]] <- round(nExp_poi[[i]]*(1-homotypic.prop[[i]]))
}

# Run DoubletFinder
for(i in 1:length(seu.list)){
  seu.list[[i]] <- 
    doubletFinder_V3.DWM_v2( # just changed it so the output metadata column name is customizable
      seu=seu.list[[i]], 
      PCs = 1:seu.list[[i]]@reductions$umap_RNA@misc$n.pcs.used, 
      pN = 0.25, #default value
      pK= pK[[i]], 
      nExp = nExp_poi.adj[[i]],  
      reuse.pANN = F, 
      classification.name='DF.individual', 
      pANN.name='DF.pANN.individual'
    )
}
#
# Merging and Integration ####
#     Merge datasets ####
# set default assay to RNA
tmp.list <- list()
for(i in 1:length(seu.list)){
  DefaultAssay(seu.list[[i]]) <- "RNA"
  tmp.list[[i]] <- DietSeurat(seu.list[[i]], assays = "RNA")
}

# merging in batches seems to reduce runtime/memory usage
tmp.merged <- list()

tmp.merged[[1]] <- merge(
  tmp.list[[1]],
  y = tmp.list[2:20],
  add.cell.ids = meta$sample[1:20]
)
tmp.merged[[2]] <- merge(
  tmp.list[[21]],
  y = tmp.list[22:40],
  add.cell.ids = meta$sample[21:40]
)
tmp.merged[[3]] <- merge(
  tmp.list[[41]],
  y = tmp.list[42:60],
  add.cell.ids = meta$sample[41:60]
)
tmp.merged[[4]] <- merge(
  tmp.list[[61]],
  y = tmp.list[62:80],
  add.cell.ids = meta$sample[61:80]
)
tmp.merged[[5]] <- merge(
  tmp.list[[81]],
  y = tmp.list[82:length(seu.list)],
  add.cell.ids = meta$sample[81:length(tmp.list)]
)

# merge tmp count matrices
scMuscle.pref.seurat <- merge(
  tmp.merged[[1]],
  y = tmp.merged[2:length(tmp.merged)]
)

# clean up memory usage
if(!is.null(scMuscle.pref.seurat)){
  rm(tmp.merged)
  rm(tmp.list)
  gc()
}

VlnPlot(
  scMuscle.pref.seurat,
  features = c(
    'nCount_RNA',
    'nFeature_RNA',
    'percent.mt'
  ),
  group.by = 'source',
  pt.size = 0
)
#
#     Seurat preprocessing on merged data ####
DefaultAssay(scMuscle.pref.seurat) <- 'RNA'

scMuscle.pref.seurat <- 
  NormalizeData(
    scMuscle.pref.seurat,
    assay = 'RNA'
  ) %>% FindVariableFeatures(
    selection.method = 'vst',
    nfeatures = 2000,
    verbose = TRUE
  ) %>% ScaleData(
    assay = 'RNA',
    verbose = TRUE
  ) %>% RunPCA(
    assay = 'RNA',
    reduction.name = 'pca_RNA',
    reduction.key = 'pca_RNA_',
    verbose = TRUE,
    npcs = 50
  ) 

ElbowPlot(
  scMuscle.pref.seurat,
  reduction = 'pca_RNA',
  ndims = 50
)

n.pcs = npcs(scMuscle.pref.seurat, reduction ="pca_RNA", var.toal = 0.95)

scMuscle.pref.seurat <-
  RunUMAP(
    scMuscle.pref.seurat,
    reduction = 'pca_RNA', 
    dims = 1:n.pcs,
    reduction.name='umap_RNA'
  ) %>% FindNeighbors(
    reduction = 'pca_RNA',
    dims = 1:n.pcs,
    force.recalc = TRUE,
    verbose = F
  ) 

scMuscle.pref.seurat <- FindClusters(object = scMuscle.pref.seurat, resolution = 0.8)
scMuscle.pref.seurat[['RNA_res.0.8']] <- as.numeric(scMuscle.pref.seurat@active.ident)
#
# Harmony integration (RNA, w/ doublets) ####
# https://github.com/immunogenomics/harmony
scMuscle.pref.seurat <- 
  scMuscle.pref.seurat %>% RunHarmony(
    group.by.vars=c('sample'),
    reduction='pca_RNA',
    assay='RNA', 
    plot_convergence = TRUE,
    verbose=TRUE
  ) 

n.pcs = npcs(scMuscle.pref.seurat,reduction="harmony")

scMuscle.pref.seurat <- 
  scMuscle.pref.seurat %>% RunUMAP(
    reduction = 'harmony',
    dims = 1:n.pcs,
    reduction.name='umap_harmony'
  )  
scMuscle.pref.seurat@reductions$umap_harmony@misc$n.pcs.used <- n.pcs

scMuscle.pref.seurat <- 
  scMuscle.pref.seurat %>% FindNeighbors(
    reduction = 'harmony',
    dims = 1:n.pcs, 
    graph.name = 'harmony_snn',
    force.recalc = TRUE,
    verbose = FALSE
  )

scMuscle.pref.seurat <- FindClusters(
  object = scMuscle.pref.seurat, 
  resolution = 1.0,
  graph.name='harmony_snn'
)
scMuscle.pref.seurat[['harmony_res.1.0']] <- as.numeric(scMuscle.pref.seurat@active.ident) # save idents as numerics

scMuscle.pref.seurat <- FindClusters(
  object = scMuscle.pref.seurat,
  resolution = 1.2,
  graph.name='harmony_snn'
)
scMuscle.pref.seurat[['harmony_res.1.2']] <- as.numeric(scMuscle.pref.seurat@active.ident)
#
# BBKNN integration ####
# https://github.com/Teichlab/bbknn
use_condaenv("scRNA3")

anndata = import("anndata",convert=FALSE)
bbknn <- import(module = 'bbknn')
sc = import("scanpy.api",convert=FALSE)
np = import("numpy")
scipy = import("scipy")

tmp.bbknn <- bbknn$bbknn_pca_matrix( # ~35min for ~350k cells
  use_faiss = T,
  pca = scMuscle.pref.seurat@reductions$pca_RNA@cell.embeddings,
  batch_list = scMuscle.pref.seurat$sample
)

scMuscle.pref.seurat@graphs$bbknn <-Seurat::as.Graph(
  sparseMatrix(
    j = tmp.bbknn[[2]]$indices,
    p = tmp.bbknn[[2]]$indptr,
    x = as.vector(tmp.bbknn[[2]]$data),
    index1 = F, #b/c coming from python
    giveCsparse = T,
    dimnames = list(colnames(scMuscle.pref.seurat),colnames(scMuscle.pref.seurat))
  )
)

scMuscle.pref.seurat@graphs$bbknn@assay.used <- "RNA"

#       Clustering and dim reduction on bbknn output ####
scMuscle.pref.seurat <- 
  scMuscle.pref.seurat %>% RunUMAP(
    graph = 'bbknn',
    reduction.name='umap_bbknn',
    reduction.key = "umap_bbknn_"
  ) 

scMuscle.pref.seurat <- FindClusters(
  object = scMuscle.pref.seurat,
  resolution = 1.2,
  graph.name='bbknn'
)
scMuscle.pref.seurat[['bbknn_res.1.2']] <- as.numeric(scMuscle.pref.seurat@active.ident)

#
# scanorama  integration (w/ doublets) ####
# https://github.com/brianhie/scanorama

# Prepare a list of lists of gene names 
DefaultAssay(scMuscle.pref.seurat) <- 'RNA'
genes_list <- strLoop(rownames(scMuscle.pref.seurat),length(meta$sample))

# Split up the object and pipe into scanorama 
integrated.corrected.data <- #~36hrs runtime, ~300-400Gb RAM needed
  DietSeurat( #pull out just the stuff you need
    scMuscle.pref.seurat,
    counts=F, data=T, scale.data=F,
    assays='RNA'
  ) %>% SplitObject(
    split.by = "sample" 
  ) %>% scanoramaPrep( # custom function, see helper function script for details
  ) %>% scanorama$correct(
    genes_list,
    return_dimred=TRUE, 
    return_dense=TRUE, #maybe false?
    ds_names = meta$sample,
    verbose=TRUE
  ) 
gc()

# Combine all of the corrected dim_red outputs, transposing so cells are cols again
corrected_scanorama_dimred <- do.call(rbind, integrated.corrected.data[[1]])
rownames(corrected_scanorama_dimred) <- colnames(scMuscle.pref.seurat)
colnames(corrected_scanorama_dimred) <- paste0("scanorama.dr_",1:nrow(corrected_scanorama))

# Add scanorama data to the Seurat object
scMuscle.pref.seurat[["scanorama"]] <- CreateDimReducObject(
  embeddings = corrected_scanorama_dimred,
  key="scanoramadr_",
  assay=DefaultAssay(scMuscle.pref.seurat)
)

if(!is.null(scMuscle.pref.seurat[["scanorama"]])){
  rm(corrected_scanorama_dimred)
  rm(integrated.corrected.data)
}
gc()

# 
#     Seurat preprocessing of Scanorama integration ####
scMuscle.pref.seurat@reductions$scanorama@stdev <- 
  apply(scMuscle.pref.seurat@reductions$scanorama@cell.embeddings, 2, sd) #find std dev for scanorama vals

n.pcs = npcs(scMuscle.pref.seurat, reduction = "scanorama", var.toal = 0.95)

scMuscle.pref.seurat <- scMuscle.pref.seurat %>%
  RunUMAP(
    reduction = "scanorama",
    dims = 1:n.pcs,
    reduction.name = 'umap_scanorama',
    reduction.key = 'umap_scanorama_'
  )
scMuscle.pref.seurat <- scMuscle.pref.seurat %>%
  FindNeighbors(
    reduction = "pca_scanorama",
    dims = 1:n.pcs,
    force.recalc = TRUE,
    verbose = FALSE
  )

scMuscle.pref.seurat <- FindClusters(
  object = scMuscle.pref.seurat, 
  resolution = 0.8
)
scMuscle.pref.seurat[['scanorama_res.0.8']] <- as.numeric(scMuscle.pref.seurat@active.ident) 

#
# Find cell type IDs ####
#     Dot Plot ####
#       cell type marker genes ####
all.genes <- c(
  'Ptprc',
  'Ms4a1', #B Cell
  'Igkc', # Plasma Cell
  'Cd8a', # T cell
  # 'Cd27', # Cd8 memory T Cell
  # 'Ccl5',# T cell
  # 'Ctla4', # CD4 T cell
  # 'Cd4', # CD4 T cell
  # 'Tbx21',
  # 'Gata3',
  'Nkg7', # NK Cell
  'Klra7', # NK Cell
  'Klra4', # NK Cell
  'Gzma', # NK Cell
  
  'S100a9',# Neutrophil
  # 'S100a8', # Neutrophil
  # 'Cxcr4', # Neutrophil
  # "Csf3r",
  "Mmp9",
  "Cd14",
  
  #monocyte
  "Ctsa" ,
  "Ctsb",
  "Ctsl",
  "Ctsz",
  "Cx3cr1",
  "Fabp5",
  "Cxcl10",
  "Ly6c1",
  "Ly6c2",
  "Ccr2",
  "Treml4",
  "Nr4a1",
  
  ### Macrophages
  ## M1
  'Itgam', #Cd11b
  "Ccl2",
  'Ccl6', # Pro-infl. Macrophage
  'Ccl9', # Pro-infl. Macrophage
  "Tgfbi",
  
  ## M2 
  # "Mapk14", #p38, Loiben's work (Cx3cr1_hi macro)
  "Ctss", #antigen presentation
  "Fcgr3", #Cd16
  "Csf1r",
  'C1qa', # Anti-infl. Macrophage
  'Cd74', # Resident Macrophage or Dendritic
  'Cd83', # Macrophage
  
  ### Dendritic
  "Flt3",
  'Cd209a', # Dendritic
  'Xcr1', # Dendritic
  
  ### Endothelial
  # 'Pecam1', #pan-Endothelial
  # 'Cdh5', # pan-Endothelial
  # 'Eng', #Endoglin - large vessel?
  # 
  # "Fabp4",
  # 
  # 'Ly6a', #large-artery, Sca1 
  # "Fbln5", #Large artery 
  # "Hey1", #cap/arterial 
  # 
  # "Lrg1", # Large Vein 
  # "Icam1", 
  # "Hif1a", 
  # 'Vwf', # Endothelial, large vein 
  # 
  # "Fmo1",#cap venous
  # "Fmo2",# cap venous
  # 'Kdr',#Flk1, capillary 
  # "Cxcl12",#capillary 
  # "Lpl", # capillary
  # "Apln", #angiogenic - capillary
  # "Aplnr", #angiogenic - capillary
  
  ## pericyte
  # "Mcam", 
  # "Alpl",
  # "Pdgfrb", 
  # "Cspg4", 
  # "Angpt1", 
  # "Igf1",
  
  ### Smooth Muscle
  # "Acta2",
  # 'Myh11', 
  # "Myl9",
  
  ### Neural
  # 'Mpz' ,# Schwann
  # 'Ptn', # Neurall/Glial
  
  ### FAPs
  # 'Pdgfra', # Fibro/Adipogenic progenitor
  
  ## Cd34+ FAPs
  "Cd34", #stem markers
  # 'Thy1',
  # "Tnfaip6",
  # "Igfbp5",
  # "Dpp4",
  # "Il33",
  # "Il6",
  
  ### Pro-remodeling FAPs
  # "Egfr",
  # "Mmp2", #Pro-remodeling/fibrotic FAPs
  # "Mmp14", #Pro-remodeling/fibrotic FAPs
  # "Ctsk", #Pro-remodeling/fibrotic FAPs
  # "Adam12",
  # "Bmp1",
  # "Bmp5",
  # "Myoc",
  # 'Col1a1', #
  # "Col3a1",
  # 'Dcn', #fibroblast
  
  ### Adipo FAPs
  # 'Gsn', # Adipo FAPs
  # "Bgn", # Adipo FAPs
  # "Apod", # Adipo FAPs
  # "Hdlbp", # Adipo FAPs
  # "Dio2", #Adipocyte
  # "Cd164",
  # "Rspo3",
  
  ### Tenocytes
  # 'Tnmd', #Tendon
  # 'Scx',
  
  ### Myogenic
  # 'Cd9',
  # 'Itga7',
  # 'Vcam1',
  # "Sdc4",
  # 'Pax3',
  # 'Pax7', #MuSC
  # 'Myf5', #MB
  # 'Myod1',
  # "Cdkn1a",
  # 'Myog',
  # 'Tnnt2',
  # 'Dmd',
  # "Tnnt3",
  # 'Acta1',
  # "Vegfa", #Myh1+
  # "Mybpc1",
  # "Mybpc2",
  # 'Myh1', # Type IIx
  # "Myh2", #Type IIa
  "Myh4" # Type IIb
)
#
#       plotting ####
scMuscle.pref.seurat%>%DotPlot(
  group.by='harmony_res.1.0',
  assay='RNA',
  dot.scale = 5,
  features=rev(all.genes)
) + 
  NoLegend() +
  scale_color_viridis()+
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.title = element_blank(), 
    axis.text.x=element_text(angle = 90,hjust = 1,vjust= 0.5),
    panel.grid.major = element_line(colour = "gray", size = 0.5
    )
  ) 

#
#     Harmony res.1.0 -- DF  ####
scMuscle.pref.seurat <- AddCellTypeIdents(
  object = scMuscle.pref.seurat, 
  old.name = 'harmony_res.1.0',
  new.name = 'prefilter_IDs',
  new.idents = c(
    "FAPs (Pro-remodeling)", #1
    "MuSCs", #2
    "Monocyte (Patrolling)", #3
    "Monocyte (Inflammatory)", #4
    "Endothelial (Capillary)", #5
    "Myonuclei (Type IIb)", #6
    "M2 Macro.", #7
    NA, #8
    "FAPs (Stem)",#9
    "Myonuclei (Type IIx)",  #10
    "T Cells", #11
    "Dendritic",#12
    "Smooth Muscle", #13
    "Endothelial (Artery)", #14
    "Tenocytes",  #15
    "Neutrophils", #16
    "B Cells", #17 
    "NK Cells", #18
    "Endothelial (Vein)", #19
    "Neural", #20
    "Myoblasts/Progenitors", #21
    "Dendritic",  #22
    "M1 Macro.", #23
    NA, #24
    "Endothelial (Stem)", #25
    "MuSCs" #26 
  )
)

# order the cell types
scMuscle.pref.seurat$prefilter_factorIDs <- factor(
  x = scMuscle.pref.seurat$prefilter_IDs, 
  levels = c(
    "B Cells", #17 
    "T Cells", #11
    "NK Cells", #18
    
    "Neutrophils", #16
    
    "Monocyte (Patrolling)", #3
    "Monocyte (Inflammatory)", #4
    "M1 Macro.", #23
    "M2 Macro.", #7
    "Dendritic",#12 22
    
    "Endothelial (Stem)", #25
    "Endothelial (Capillary)", #5
    "Endothelial (Artery)", #14
    "Endothelial (Vein)", #19
    "Smooth Muscle", #13
    "Neural", #20
    
    "FAPs (Stem)",#9
    "FAPs (Pro-remodeling)", #1
    "Tenocytes",  #15
    
    "MuSCs", #2 26 
    "Myoblasts/Progenitors", #21
    "Myonuclei (Type IIx)",  #10
    "Myonuclei (Type IIb)" #6
  )
)
#
#     Harmony res.1.0 -- DF -- superclusters  ####
scMuscle.pref.seurat <- AddCellTypeIdents(
  object = scMuscle.pref.seurat, 
  old.idents = 'harmony_res.1.0',
  newName = 'group.1.0',
  new.idents = c(
    "FAPs", #1
    "Myogenic", #2
    "Monocyte", #3
    "Macrophage", #4
    "Endothelial", #5
    "Myonuclei", #6
    "Macrophage", #7
    NA, #8
    "FAPs",#9
    "Myonuclei",  #10
    "T Cells", #11
    "Dendritic",#12
    "Smooth Muscle", #13
    "Endothelial", #14
    "Tenocytes",  #15
    "Neutrophils", #16
    "B Cells", #17 
    "NK Cells", #18
    "Endothelial", #19
    "Neural", #20
    "Myogenic", #21
    "Dendritic",  #22
    "Macrophage", #23
    NA, #24
    "Endothelial", #25
    "MuSCs" #26 
  )
)

# order the cell types
scMuscle.pref.seurat$group.1.0 <- factor(
  x = scMuscle.pref.seurat$group.1.0, 
  levels = c(
    "B Cells", #17 
    "T Cells", #11
    "NK Cells", #18
    "Neutrophils", #16
    
    "Monocyte", #3
    "Macrophage", 
    "Dendritic",#12
    
    "Endothelial", 
    
    "Smooth Muscle", #13
    "FAPs", 
    "Tenocytes",  #15
    "Neural", #20
    
    "Myogenic", 
    "Myonuclei", 
    
    "NOISY"
  )
)

#
#     BBKNN res.1.2 ####
scMuscle.pref.seurat <- AddCellTypeIdents(
  seu = scMuscle.pref.seurat, 
  old.name = 'bbknn_res.1.2',
  new.name = 'bbknn_res.1.2_IDs',
  new.idents = c(
    "Endothelial (Capillary)", #1
    "FAPs (Adipogenic)", #2
    "Monocyte (Patrolling)", #3
    "Monocyte (Infiltrating)", #4  
    "FAPs (Pro-remodeling)", #5
    "M2 Macro.", #6   
    "Quiescent MuSCs", #7    
    "Monocyte (Inflmmatory)", #8 
    "Myonuclei (Acta1_hi)", #9
    "FAPs (Stem)", #10
    "Myoblasts", #11
    "Dendritic", #12
    "Myonuclei (Type IIb)", #13
    "Activated MuSCs", #14
    "T Cells", #15
    "Myonuclei (Type IIx)", #16
    "Smooth Muscle", #17
    "M1 Macro.", #18
    "Neutrophils", #19
    "Endothelial (Artery)", #20
    "Monocytes (Cxcl10+)", #21
    "Tenocytes", #22
    "Committed Progenitors", #23
    "Committed Progenitors (Noisy)", #24
    "B Cells", #25
    "Monocytes (Noisy)", #26
    "FAPs (Noisy)",  #27
    "Neural", #28
    "NK Cells", #29
    "Endothelial (Vein)" #30
  )
)
Monocytes (Patrolling)
scMuscle.seurat$bbknn_factorIDs <- factor(
  scMuscle.seurat$bbknn_res.1.2_IDs,
  levels = c(
    
  )
)
#
#     scanorama_res.0.8 ####
scMuscle.pref.seurat <- AddCellTypeIdents(
  seu = scMuscle.pref.seurat, 
  old.name = 'scanorama_res.0.8',
  new.name = 'scanorama_res.0.8_IDs',
  new.idents = c(
    "M2 Macro. (1)", #1
    "Endothelial (Capillary)", #2
    "FAPs (3)", #3
    "FAPs (4)", #4
    "Monocytes (Infiltrating)", #5
    "Quiescent MuSCs", #6
    "Monocytes (7)", #7
    "Myoblasts (8)", #8
    "Monocytes (Inflammatory)", #9
    "FAPs (Stem)", #10
    "Myonuclei (Type IIx)", #11
    "Activated MuSCs", #12
    "Dendritic", #13
    "Myonuclei (Type IIb)", #14
    "Endothelial (Capillary)", #15
    "T Cells", #16
    "Macrophage (17)", #17
    "Myonuclei (Acta1_hi)", #18
    "Tenocytes", #19
    "Committed Progenitors", #20
    "Monocytes (21)", #21
    "Myogenic (Noisy)", #22
    "Neutrophils", #23
    "B Cells", #24
    "Smooth Muscle", #25
    "Myonuclei (Type IIb)", #26
    "FAPs (Noisy_27)", #27
    "Endothelial (28)", #28
    "NK Cells", #29
    "Monocyte (Inflammatory)_30", #30
    "Neural", #31
    "Endothelial (32)", #32
    "Myonuclei (33)", #33
    "Endothelial (Capillary)", #34  
    "FAPs (Pro-remodeling)_35", #35
    "Monocytes (Noisy)", #36   
    "Neutrophils", #37    
    "FAPs (Noisy_38)", #38 
    "Endothelial (Capillary)", #39
    "FAPs (Noisy_40)", #40
    "T Cells", #41
    "Endothelial (Capillary)", #42 
    "MuSCs_43", #43
    "Myonuclei (Acta1_hi)", #44
    "MuSCs (Noisy)", #45
    "FAPs (46)", #46
    "Macrophages (47)" #47
  )
)

scMuscle.pref.seurat$scanorama_factorIDs <- factor(
  scMuscle.pref.seurat$scanorama_res.0.8_IDs,
  levels = c(
    
  )
)

#
# Analysis of noisy cells ####
#     noisy cells- violin plot ####
ggplot(
  scMuscle.pref.seurat@meta.data,
  aes(
    x=is.noisy,
    y=nFeature_RNA
  )
)+
  geom_violin(fill="gray") +
  theme_minimal() +
  theme(
    axis.line = element_line(size=0.5, color="black"),
    panel.grid.major = element_line(size=1, color="gray"),
    panel.grid.minor = element_line(size=1, color="light gray"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(color="black", size=small.font),
    axis.title =  element_text(color="black", size=big.font),
    axis.ticks.x = element_blank()
  ) + 
  scale_y_log10()

#
# Remove doublets  ####
Idents(scMuscle.pref.seurat)<-'sample' 
scMuscle.seurat <- subset(
  scMuscle.pref.seurat,
  subset = DF.individual == 'Singlet' 
)  

#     Seurat preprocessing on merged data (w/out doublets)####
DefaultAssay(scMuscle.seurat) <- 'RNA'

scMuscle.seurat <- 
  NormalizeData(
    scMuscle.seurat,
    assay = 'RNA'
  ) %>% FindVariableFeatures(
    selection.method = 'vst',
    nfeatures = 2000,
    verbose = F
  ) %>% ScaleData(
    assay = 'RNA',
    verbose = F
  ) %>% RunPCA(
    assay = 'RNA',
    reduction.name = 'pca_RNA',
    reduction.key = 'pca_RNA_',
    verbose = F,
    npcs = 50
  ) 

ElbowPlot(
  scMuscle.seurat,
  reduction = 'pca_RNA',
  ndims = 50
)

n.pcs = npcs(scMuscle.clean.seurat, reduction="pca_RNA")

scMuscle.seurat <-
  RunUMAP(
    scMuscle.seurat,
    reduction = 'pca_RNA', 
    dims = 1:n.pcs,
    reduction.name='umap_RNA'
  ) %>% FindNeighbors(
    reduction = 'pca_RNA',
    dims = 1:n.pcs,
    force.recalc = TRUE,
    verbose = F
  ) 
scMuscle.seurat@reductions[["umap_RNA"]]@misc$n.pcs.used <- n.pcs

scMuscle.seurat <- FindClusters(
  object = scMuscle.seurat,
  resolution = 0.8
  )
scMuscle.seurat[['RNA_res.0.8']] <- as.numeric(scMuscle.seurat@active.ident)
#
# harmony integration (after QC/doublet removal) ####
scMuscle.seurat <- 
  scMuscle.seurat %>% RunHarmony(
    group.by.vars=c('sample'),
    reduction='pca_RNA',
    assay='RNA', 
    plot_convergence = TRUE,
    verbose=TRUE
  ) 

#     Find n.pcs
n.pcs = npcs(scMuscle.clean.seurat, reduction = "harmony")

scMuscle.clean.seurat <- 
  scMuscle.clean.seurat %>% RunUMAP(
    reduction = 'harmony',
    dims = 1:n.pcs, 
    reduction.name='umap_harmony'
  )  

scMuscle.clean.seurat@reductions[["umap_harmony"]]@misc$n.pcs.used <- n.pcs

scMuscle.clean.seurat <- 
  scMuscle.clean.seurat %>% FindNeighbors(
    reduction = 'harmony',
    dims = 1:n.pcs, 
    graph.name = 'harmony_snn',
    force.recalc = TRUE,
    verbose = T
  )

# scMuscle.seurat <- FindClusters(object = scMuscle.seurat, resolution = 1.0,graph.name='harmony_snn')
# scMuscle.seurat[['harmony_res.1.0']] <- as.numeric(scMuscle.seurat@active.ident)

scMuscle.clean.seurat <- FindClusters(object = scMuscle.clean.seurat, resolution = 1.2,graph.name='harmony_snn')
scMuscle.clean.seurat[['harmony_res.1.2']] <- as.numeric(scMuscle.clean.seurat@active.ident)
#
# BBKNN integration ####
# https://github.com/Teichlab/bbknn
use_condaenv("scRNA3")

anndata = import("anndata",convert=FALSE)
bbknn <- import(module = 'bbknn')
sc = import("scanpy.api",convert=FALSE)
np = import("numpy")
scipy = import("scipy")

tmp.bbknn <- bbknn$bbknn_pca_matrix( # ~35min for 350k cells
  use_faiss = T,
  pca = scMuscle.noisy.seurat@reductions$pca_RNA@cell.embeddings,
  batch_list = scMuscle.noisy.seurat$sample
)

# Optional: save bbknn output
# scipy$sparse$save_npz(file = "tmp_bbknn",matrix = tmp.dimred[[2]])
# tmp.bbknn <- np$load("tmp_bbknn.npz") 
# tmp.bbknn$files
# [1] "indices" "indptr"  "format"  "shape"   "data"  

tmp <-Seurat::as.Graph(
  sparseMatrix(
    j = tmp.bbknn[[2]]$indices,
    p = tmp.bbknn[[2]]$indptr,
    x = as.vector(tmp.bbknn[[2]]$data),
    index1 = F, #b/c coming from python
    giveCsparse = T,
    dimnames = list(colnames(scMuscle.noisy.seurat),colnames(scMuscle.noisy.seurat))
    # dims = unlist(tmp.dimred[[2]]$shape)
  )
)
scMuscle.noisy.seurat@graphs$bbknn <- tmp
scMuscle.noisy.seurat[["bbknn"]] <- tmp
scMuscle.noisy.seurat@graphs$bbknn@assay.used <- "RNA"

#       Clustering and dim reduction on bbknn output ####
scMuscle.noisy.seurat <- 
  scMuscle.noisy.seurat %>% RunUMAP(
    graph = 'bbknn',
    reduction.name='umap_bbknn',
    reduction.key = "umap_bbknn_"
  ) 

scMuscle.noisy.seurat <- FindClusters(
  object = scMuscle.noisy.seurat,
  resolution = 1.2,
  graph.name='bbknn'
)
scMuscle.noisy.seurat[['bbknn_res.1.2']] <- as.numeric(scMuscle.noisy.seurat@active.ident)

#
# scanorama  integration (after QC/doublet removal) ####
# Prepare a list of lists of gene names 
DefaultAssay(scMuscle.noisy.seurat) <- 'RNA'
genes_list <- strLoop(rownames(scMuscle.noisy.seurat),length(meta$sample))

# Split up the object and pipe into scanorama 
timestamp()
integrated.corrected.data <- #~36hrs runtime, ~300-400Gb RAM needed, max
  DietSeurat( #pull out just the stff you need
    scMuscle.noisy.seurat,
    counts=F, data=T, scale.data=F,
    assays='RNA'
  ) %>% SplitObject(
    split.by = "sample" 
  ) %>% scanoramaPrep(
  ) %>% scanorama$correct(
    genes_list,
    # batch_size=as.integer(10000), #`int`, optional (default: `5000`)
    return_dimred=TRUE, 
    return_dense=TRUE, #maybe false?
    ds_names = meta$sample,
    verbose=TRUE,
    dimred = as.integer(50) #n.dim of dim_red
    # hvg: `int`, optional (default: None)
  ) 
gc()
timestamp()

# Combine all of the corrected dim_red outputs, transposing so cells are cols again
corrected_scanorama_dimred <- do.call(rbind, integrated.corrected.data[[1]])

rownames(corrected_scanorama_dimred) <- colnames(scMuscle.noisy.seurat)
colnames(corrected_scanorama_dimred) <- paste0("scanorama_",1:ncol(corrected_scanorama_dimred))

# Add scanorama data to the Seurat object
scMuscle.noisy.scanorama <- CreateDimReducObject(
  embeddings = corrected_scanorama_dimred,
  key="scanorama_",
  assay=DefaultAssay(scMuscle.noisy.seurat)
)
save(scMuscle.noisy.scanorama,file="scMuscle.noisy.scanorama.RData")
# 
#     Seurat preprocessing of Scanorama integration ####
scMuscle.noisy.seurat@reductions$scanorama@stdev <- 
  apply(scMuscle.noisy.seurat@reductions$scanorama@cell.embeddings, 2, sd) #find std dev for scanorama vals

n.pcs = npcs(scMuscle.noisy.seurat, reduction = "scanorama")

scMuscle.noisy.seurat <- 
  scMuscle.noisy.seurat %>% RunUMAP(
    reduction = "scanorama",
    dims = 1:n.pcs,
    reduction.name = 'umap_scanorama',
    reduction.key = 'umap_scanorama_'
  ) 

scMuscle.noisy.seurat <- 
  scMuscle.noisy.seurat %>% FindNeighbors(
    reduction = "scanorama",
    graph.name="scanorama_snn",
    dims = 1:n.pcs,
    force.recalc = TRUE,
    verbose = FALSE
  )

scMuscle.noisy.seurat <- FindClusters(object = scMuscle.noisy.seurat, resolution = 0.8)
scMuscle.noisy.seurat[['scanorama_res.0.8']] <- as.numeric(scMuscle.noisy.seurat@active.ident)

# scMuscle.seurat <- FindClusters(object = scMuscle.seurat, resolution = 1.0)
# object.seurat[['scanorama_res.1.0']] <- as.numeric(scMuscle.seurat@active.ident) 

# scMuscle.seurat <- FindClusters(object = scMuscle.seurat, resolution = 1.2)
# scMuscle.seurat[['scanorama_res.1.2']] <- as.numeric(scMuscle.seurat@active.ident) 
#
# Edit cluster IDs - w/ out noisy cells ####
#TODO     harmony res.1.2 ####
scMuscle.noisy.seurat <- AddCellTypeIdents(
  seu = scMuscle.noisy.seurat, 
  old.name = 'harmony_res.1.2',
  new.name = 'harmony_res.1.2_IDs',
  new.idents = c(
    "FAPs (Pro-remodeling)", #1
    "MuSCs/Myoblasts", #2
    "Endothelial (Capillary)", #3
    "Monocyte (Infiltrating)", #4  
    "Monocyte (Patrolling)", #5
    "M2 Macro.", #6 
    "Myonuclei (Type IIx)", #7    
    "Myonuclei (Type IIb)", #8  
    "Monocyte (Inflammatory)", #9
    "FAPs (Stem)", #10
    "noisy_11", #11
    "T Cells", #12
    "Dendritic", #13
    "Endothelial (Artery)", #14
    "Smooth Muscle", #15
    "Neutrophils", #16
    "Tenocytes", #17
    "B Cells", #18
    "Committed Progenitors", #19
    "Monocyte (Cxcl10+)", #20
    "Endothelial (Vein)", #21 
    "NK Cells", #22
    "Neural", #23
    "FAPs (Adipogenic)", #24
    "M1 Macro.", #25
    "noisy_26", #26
    "Dendritic", #27
    "noisy_28", #28
    "noisy_29", #29
    "MuSCs/Myoblasts" #30
  )
)

scMuscle.noisy.seurat$harmony_factorIDs <- factor(
  scMuscle.noisy.seurat$harmony_res.1.2_IDs,
  levels = c(
    
  )
)

#
#     harmony_res.1.2 -- super-clusters ####
scMuscle.seurat <- AddCellTypeIdents(
  object = scMuscle.seurat, 
  old.idents = 'harmony_res.1.2',
  newName = 'harmony_superIDs',
  new.idents = c(
    
  )
)

scMuscle.seurat$harmony_superFactor <- factor(
  scMuscle.seurat$harmony_superIDs,
  levels = c(
    "Lymphocytes",
    "Myeloid",
    "Endothelial",
    "Smooth Muscle",
    "FAPs", 
    "Tenocytes",
    "Neural",
    "Myogenic", 
    "Myonuclei"
  )
)

#
#     BBKNN res.1.2 ####
scMuscle.noisy.seurat <- AddCellTypeIdents(
  seu = scMuscle.noisy.seurat, 
  old.name = 'bbknn_res.1.2',
  new.name = 'bbknn_res.1.2_IDs',
  new.idents = c(
    "Monocytes (Patrolling)", #1
    "Endothelial (Capillary)", #2
    "Monocytes (Infiltrating)", #3
    "FAPs (Adipogenic)", #4  
    "FAPs (Pro-remodeling)", #5
    "Monocytes (Inflmmatory)", #6   
    "Quiescent MuSCs", #7    
    "Myoblasts", #8 
    "FAPs (Stem)", #9
    "Activated MuSCs", #10
    "Dendritic", #11
    "Myonuclei (Type IIb)", #12
    "Myonuclei (Acta1_hi)", #13
    "Myonuclei (Type IIx)", #14
    "T Cells", #15
    "Endothelial (Artery)", #16
    "Neutrophils", #17
    "Smooth Muscle", #18
    "M1/M2 Macro.", #19
    "Monocytes (Cxcl10+)", #20
    "Tenocytes", #21
    "Committed Progenitors_1", #22
    "Committed Progenitors_2", #23
    "B Cells", #24
    "Monocytes (25)", #25
    "NK Cells", #26
    "Endothelial (Vein)",  #27
    "Neural", #28
    "FAPs (Adipogenic)", #29
    "Monocytes (Inflmmatory)" #30
  )
)

scMuscle.seurat$bbknn_factorIDs <- factor(
  scMuscle.seurat$bbknn_res.1.2_IDs,
  levels = c(
    
  )
)
#
#     BBKNN res.1.2 - super clusters####
scMuscle.seurat <- AddCellTypeIdents(
  object = scMuscle.seurat, 
  old.idents = 'bbknn_res.1.2',
  newName = 'bbknn_superIDs',
  new.idents = c(
    "Myeloid", #1
    "Endothelial", #2
    "FAPs", #3
    "Myeloid", #4  
    "FAPs", #5
    "Myeloid", #6   
    "Myogenic", #7    
    "Myonuclei", #8 
    "FAPs", #9
    "Myogenic", #10
    "Myogenic", #11
    "Myeloid", #12
    "Myonuclei", #13
    "Myeloid", #14
    "Myeloid", #15
    "Myonuclei", #16
    "Lymphocytes", #17
    "Endothelial", #18
    "Myeloid", #19
    "Smooth Muscle", #20
    "Tenocytes", #21
    "Myogenic", #22
    "Endothelial", #23
    "Lymphocytes", #24
    "Lymphocytes", #25
    "Endothelial", #26
    "Neural",  #27
    "Myeloid", #28
    "Myeloid", #29
    "Endothelial" #30
  )
)

scMuscle.seurat$bbknn_superFactor <- factor(
  scMuscle.seurat$bbknn_superIDs,
  levels = c(
    "Lymphocytes",
    "Myeloid",
    "Endothelial",
    "Smooth Muscle",
    "FAPs", 
    "Tenocytes",
    "Neural",
    "Myogenic", 
    "Myonuclei"
  )
)
#
#     scanorama_res.0.8 ####
scMuscle.noisy.seurat <- AddCellTypeIdents(
  object = scMuscle.noisy.seurat, 
  old.idents = 'scanorama_res.0.8',
  newName = 'scanorama_res.0.8_IDs',
  new.idents = c(
    "M2 Macro. (Cx3cr1_hi)_1", #1
    "M2 Macro. (Cx3cr1_hi)_2", #2
    "Endothelial (Capillaries)", #3
    "MuSCs", #4  
    "FAPs (Pro-remodeling)", #5
    "Monocytes", #6   
    "M1 Macro.", #7    
    "Myoblasts", #8 
    "Myonuclei (Type IIx)", #9
    "FAPs (Adipogenic)", #10
    "Endothelial (Vein)", #11
    "FAPs (Stem)", #12
    "FAPs (Stem)", #13
    "Activated MuSCs", #14
    "NK Cells", #15
    "Myonuclei (Type IIb)", #16
    "Endothelial (Capillaries)", #17
    "FAPs (Pro-remodeling)", #18
    "Dendritic", #19
    "Neutrophils", #20
    "Endothelial (Capillaries)", #21
    "Committed Progenitors", #22
    "M1 Macro.", #23
    "Myonuclei (P10)", #24 - Vegfa_hi
    "Myonuclei_Dmd", #25
    "B Cells", #26
    "Tenocytes", #27
    "Dendritic", #28
    "Pericytes", #29!!!
    "T Cells", #30
    "Tenocytes", #31
    "Monocytes", #32
    "Neural", #33
    "Myonuclei (Type IIb)", #34  
    "Committed Progenitors", #35
    "Smooth Muscle", #36   
    "Myonuclei (Type IIx)", #37    
    "Myonuclei_Dmd", #38 
    "M2 Macro. (Cx3cr1_lo)", #39
    "Myoblasts", #40
    "Pericytes", #41
    "Ddx17_hi_42", #42 nuclei
    "Monocytes", #43
    "FAPs (Pro-remodeling)", #44
    "M1 Macro.", #45
    "Myoblasts", #46
    "Endothelial (Vein)", #47
    "Myonuclei_Dmd", #48
    "B Cells" #49
  )
)

scMuscle.seurat$scanorama_factorIDs <- factor(
  scMuscle.seurat$scanorama_res.0.8_IDs,
  levels = c(
    
  )
)

#
#     scanorama_res.0.8 - super clusters ####
scMuscle.seurat <- AddCellTypeIdents(
  object = scMuscle.seurat, 
  old.idents = 'scanorama_res.0.8',
  newName = 'scanorama_superIDs',
  new.idents = c(
    "Myeloid", #1
    "Myeloid", #2
    "Endothelial", #3
    "Myogenic", #4  
    "FAPs", #5
    "Myeloid", #6   
    "Myeloid", #7    
    "Myogenic", #8 
    "Myonuclei", #9
    "FAPs", #10
    "Endothelial", #11
    "FAPs", #12
    "FAPs", #13
    "Myogenic", #14
    "Lymphocytes", #15
    "Myonuclei", #16
    "Endothelial", #17
    "FAPs", #18
    "Myeloid", #19
    "Myeloid", #20
    "Endothelial", #21
    "Myogenic", #22
    "Myeloid", #23
    "Myonuclei", #24 - Vegfa_hi
    "Myonuclei", #25
    "Lymphocytes", #26
    "Tenocytes", #27
    "Myeloid", #28
    "Endothelial", #29
    "Lymphocytes", #30
    "Tenocytes", #31
    "Myeloid", #32
    "Neural", #33
    "Myonuclei", #34  
    "Myogenic", #35
    "Smooth Muscle", #36   
    "Myonuclei", #37    
    "Myonuclei", #38 
    "Myeloid", #39
    "Myogenic", #40
    "Endothelial", #41
    "Myonuclei", #42 nuclei
    "Myeloid", #43
    "FAPs", #44
    "Myeloid", #45
    "Myogenic", #46
    "Endothelial", #47
    "Myonuclei", #48
    "Lymphocytes", #49
    "Tenocytes", #50
    "Neural", #51
    "Myogenic", #52
    "Myogenic", #53
    "Myonuclei", #54
    "Lymphocytes", #55
    "Myeloid", #56
    "Lymphocytes", #57
    "Myogenic", #58
    "Myeloid" #59
  )
)

scMuscle.seurat$scanorama_superFactor <- factor(
  scMuscle.seurat$scanorama_superIDs,
  levels = c(
    "Lymphocytes",
    "Myeloid",
    "Endothelial",
    "Smooth Muscle",
    "FAPs", 
    "Tenocytes",
    "Neural",
    "Myogenic", 
    "Myonuclei"
  )
)

#

