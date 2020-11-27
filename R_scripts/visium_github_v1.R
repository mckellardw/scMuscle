# Libs ####
library(TED) 

library(Matrix)
library(dplyr)
library(Seurat)
library(future)

library(cluster)
library(data.table)
library(reshape2)

# Plotting
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(pals)
library(viridis)
library(shades)
library(ggplotify)

source("helper_functions_v1.R") 

# Metadata ####
meta <- list()

meta$data.dir <- c(
  "/path/to/Vis5A_spaceranger_count",
  "/path/to/Vis7B_spaceranger_count",
  "/path/to/Vis9A_spaceranger_count"
)

# Sample ID
meta$sample <- c(
  "Vis5A", "Vis7B", "Vis9A"
)

# Days post injury
meta$injury <- c(
  "D2", "D5", "D7"
)

# inputs for flipping the coordinates to make DimPlots easier with Visium data 
meta$coord_flip<- c(
  T, T,T
)
meta$h_flip<- c(
  -1, -1, 1
)
meta$v_flip<- c(
  -1, -1, -1
)

meta <- as.data.frame(meta)

#
# Read in count mats & initialize Seurat objects ####

seu.list <- list()
for(i in 1:length(meta$data.dir)){ 
  if(file.exists(paste0(meta$data.dir[i], '/outs/filtered_feature_bc_matrix'))){ 
    cat("Reading #",i, ": ", meta$data.dir[i], ' \n')
    seu.list[[i]] <- Seurat::Load10X_Spatial(
      data.dir = paste0(meta$data.dir[i], '/outs'),
      filter.matrix = T
    )
  }else{
    cat("Data not found for # ", i, " (", meta$data.dir[i], ")", "\n")
  }
}

# Add in meta data features
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
  
  # hemaglobin scores
  seu.list[[i]][["percent.Hb"]]  <- PercentageFeatureSet(
    seu.list[[i]], 
    features = c(
      rownames(seu.list[[i]])[grep(rownames(seu.list[[i]]),pattern="Hba")],
      rownames(seu.list[[i]])[grep(rownames(seu.list[[i]]),pattern="Hbb")]
    )
  )
  
  # ribosomal protein scores
  seu.list[[i]][["percent.Rp"]]  <- PercentageFeatureSet(
    seu.list[[i]], 
    features = c(
      rownames(seu.list[[i]])[grep(rownames(seu.list[[i]]),pattern="Rps")],
      rownames(seu.list[[i]])[grep(rownames(seu.list[[i]]),pattern="Rpl")]
    )
  )
  
  # Filter out low quality spots according to the metrics defined above- 
  #   Note, no percent.mt filter because these samples were not dissociated
  seu.list[[i]] <- subset(seu.list[[i]],
                          subset = nCount_Spatial > 1000 &
                            nFeature_Spatial > 500
  )
}

#       Preprocess seurat objects ####
seuPreProcess <- function(seu, assay='Spatial', n.pcs=50, res=0.8){
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
  
  #find pcs to use
  n.pcs.use = npcs(seu,var.toal = 0.95,reduction = pca.name)
  
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

seu.list <- lapply(seu.list, seuPreProcess)

# Add spatial locations as a reduction (for easy plotting with DimPlot)
seu.list <- mapply(
  FUN = function(SEU, coord_flip, h_flip, v_flip){
    if(coord_flip){
      tmp <- as.matrix(cbind(
        SEU@images$slice1@coordinates$row * h_flip,
        SEU@images$slice1@coordinates$col * v_flip
      ))
    }else{
      tmp <- as.matrix(cbind(
        SEU@images$slice1@coordinates$col * h_flip,
        SEU@images$slice1@coordinates$row * v_flip
      ))
    }
    colnames(tmp) <- paste0("space_", 1:2)
    rownames(tmp) <- colnames(SEU)
    
    SEU[["space"]] <- CreateDimReducObject(
      embeddings=as.matrix(tmp),
      assay="Spatial",
      key = "space_"
    )
    
    return(SEU)
  },
  seu.list,
  meta$coord_flip,
  meta$h_flip,
  meta$v_flip
)

# BayesPrism/TED deconvolution ####
# Load reference transcriptomic data

# Coming soon!
load("/path/to/scMuscle_seurat.RData") #integrated reference from Fig. 1
load("/path/to/myo__seurat.RData") #subset myogenic cells, after running PHATE, and binning (see Fig. 2)

# remove mito and ribo genes; important for deconvolution - 
#     See vignette at https://github.com/Danko-Lab/TED
genes.remove <-  unique(c(
  rownames(scMuscle.seurat)[grep("Rps",x = rownames(scMuscle.seurat))], # ribo, short
  rownames(scMuscle.seurat)[grep("Rpl",x = rownames(scMuscle.seurat))], # ribo long 
  rownames(scMuscle.seurat)[grep("mt-",x = rownames(scMuscle.seurat))] # mito
))

#       Deconvolution with PHATE subtypes (Raw Counts, DEGs) ####

n.cores=1

scMuscle.seurat$decon_IDs <- levels(scMuscle.seurat$harmony_factorIDs)[as.numeric(scMuscle.seurat$harmony_factorIDs)]

# Use PHATE bins to label myogenic states
scMuscle.seurat$decon_IDs[ colnames(myo.slim.seurat)[myo.slim.seurat$phate1.bins %in% 4] ] <- "Quiescent_MuSCs"
scMuscle.seurat$decon_IDs[ colnames(myo.slim.seurat)[myo.slim.seurat$phate1.bins %in% 5:7] ] <- "Activated_MuSCs"
scMuscle.seurat$decon_IDs[ colnames(myo.slim.seurat)[myo.slim.seurat$phate1.bins %in% 8:10] ] <- "Committed_Myoblasts"
scMuscle.seurat$decon_IDs[ colnames(myo.slim.seurat)[myo.slim.seurat$phate1.bins %in% 11:17] ] <- "Fusing_Myocytes"
scMuscle.seurat$decon_IDs[scMuscle.seurat$harmony_IDs %in% c(
  "MuSCs", "Myoblasts/Progenitors"
) ] <- NA

# Check labels from PHATE binning in the original UMAP
DimPlot(
  scMuscle.seurat,
  group.by="decon_IDs",
  reduction="umap_harmony"
)

Idents(scMuscle.seurat) <- "decon_IDs"

#Parallelize - optional, but runtime is multiple days without parallelization
plan("multiprocess", workers = n.cores)
options(future.globals.maxSize = 4500 * 1024^2)

# Output is in "supplemental_data/scMuscle_harmonytypes_plus_phatebins_markers.csv"
all.markers <- FindAllMarkers(
  scMuscle.seurat, 
  verbose=T
)
gc()


#Build cell type gene expression profile - see helper_functions_v1.R for details
celltype.gep <- scMuscle.seurat %>% get_cell_type_model(
  slot="counts",
  cluster.names = "decon_IDs",
  genes.ignore = genes.remove # remove mito and ribo genes
)
genes.keep = unique(all.markers$gene[all.markers$pct.1>0.5 & abs(all.markers$avg_logFC)>1])
genes.keep = genes.keep[!genes.keep %in% genes.remove]

celltype.gep <- celltype.gep[rownames(celltype.gep) %in% genes.keep, ]
dim(celltype.gep)

gc()

celltype.bp <- list()
for(i in 1:length(seu.list)){ 
  cat("Sample: ", meta$sample[i],"\n")
  celltype.bp[[i]] <-run.Ted( 
    ref.dat = t(celltype.gep), #cell types x genes
    X=as.matrix(t(GetAssayData(seu.list[[i]],slot="counts"))), # spot to deconvolve
    input.type = "GEP", 
    n.cores=n.cores
  )
  gc()
}

#
#     Add outputs to Seurat objects ####
for(i in 1:length(seu.list)){ # add theta (composition) values as an Assay
  colnames(celltype.bp[[i]]$res$final.gibbs.theta) <- stringr::str_replace_all(rownames(celltype.bp[[i]]$para$input.phi),pattern = " ", replacement="_")
  seu.list[[i]][["celltype.bp"]] <- CreateAssayObject(data=t(celltype.bp[[i]]$res$final.gibbs.theta))
}


#      calculate co-occurrence ####
# Adapted from:
#   https://github.com/madhavmantri/chicken_heart/blob/master/scripts/anchor_integeration.R

# Order cell types...
celltypes.order = c(
  "Quiescent-MuSCs", "Activated-MuSCs", "Committed-Myoblasts", 
  "Fusing-Myocytes",
  "Myonuclei-(Type-IIb)","Myonuclei-(Type-IIx)",
  
  "Endothelial-(Stem)","Endothelial-(Capillary)","Endothelial-(Artery)","Endothelial-(Vein)",
  
  "Smooth-Muscle", 
  
  "FAPs-(Stem)", "FAPs-(Adipogenic)","FAPs-(Pro-remodeling)",
  "Tenocytes","Neural",
  
  "Monocyte-(Cxcl10+)", "Monocyte-(Patrolling)","Monocyte-(Inflammatory)",
  "M1-Macro.",
  "M2-Macro.-(Cx3cr1-lo)",
  "M2-Macro.-(Cx3cr1-hi)",
  
  "Dendritic", "Neutrophils", "B-Cells", "NK-Cells", "T-Cells" 
)

prediction.scores <- lapply(
  seu.list[2:4],
  FUN = function(SEU){
    prediction.scores <- t(GetAssayData(SEU, assay="sub_raw_deg_ted"))[,celltypes.order]
    return(prediction.scores)
  }
)

interaction.list <- lapply(
  prediction.scores,
  FUN=function(PRED){
    tmp.feat = colnames(PRED)
    
    interaction_matrix = matrix(
      0, 
      ncol = ncol(PRED), 
      nrow = ncol(PRED)
    )
    rownames(interaction_matrix) <- tmp.feat
    colnames(interaction_matrix) <- tmp.feat
    
    celltype.filter <- 10 # filter for number of celltypes to look for in each spot
    
    for(i in 1:nrow(PRED)){
      tmp <- names(sort(PRED[i,PRED[i,] > 0], decreasing = T))[1:celltype.filter] # filter based on hard number of cell types
      
      if(length(tmp) == 2){
        interaction_matrix[tmp[1], tmp[2] ] <- interaction_matrix[tmp[1], tmp[2] ] + 1
      } else if(length(tmp) == 3){
        interaction_matrix[tmp[1], tmp[2] ] <- interaction_matrix[tmp[1], tmp[2] ] + 1
        interaction_matrix[tmp[2], tmp[3] ] <- interaction_matrix[tmp[2], tmp[3] ] + 1
        interaction_matrix[tmp[1], tmp[3] ] <- interaction_matrix[tmp[1], tmp[3] ] + 1
      } else if(length(tmp) == 4){
        interaction_matrix[tmp[1], tmp[2] ] <- interaction_matrix[tmp[1], tmp[2] ] + 1
        interaction_matrix[tmp[2], tmp[3] ] <- interaction_matrix[tmp[2], tmp[3] ] + 1
        interaction_matrix[tmp[3], tmp[4] ] <- interaction_matrix[tmp[3], tmp[4] ] + 1
        interaction_matrix[tmp[1], tmp[3] ] <- interaction_matrix[tmp[1], tmp[3] ] + 1
        interaction_matrix[tmp[1], tmp[4] ] <- interaction_matrix[tmp[1], tmp[4] ] + 1
        interaction_matrix[tmp[2], tmp[4] ] <- interaction_matrix[tmp[2], tmp[4] ] + 1
      } else if(length(tmp) == 5){
        interaction_matrix[tmp[1], tmp[2] ] <- interaction_matrix[tmp[1], tmp[2] ] + 1
        interaction_matrix[tmp[1], tmp[3] ] <- interaction_matrix[tmp[1], tmp[3] ] + 1
        interaction_matrix[tmp[1], tmp[4] ] <- interaction_matrix[tmp[1], tmp[4] ] + 1
        interaction_matrix[tmp[1], tmp[5] ] <- interaction_matrix[tmp[1], tmp[5] ] + 1
        
        interaction_matrix[tmp[2], tmp[3] ] <- interaction_matrix[tmp[2], tmp[3] ] + 1
        interaction_matrix[tmp[2], tmp[4] ] <- interaction_matrix[tmp[2], tmp[4] ] + 1
        interaction_matrix[tmp[2], tmp[5] ] <- interaction_matrix[tmp[2], tmp[5] ] + 1
        
        interaction_matrix[tmp[3], tmp[4] ] <- interaction_matrix[tmp[3], tmp[4] ] + 1
        interaction_matrix[tmp[3], tmp[5] ] <- interaction_matrix[tmp[3], tmp[5] ] + 1
        
        interaction_matrix[tmp[4], tmp[5] ] <- interaction_matrix[tmp[4], tmp[5] ] + 1
      } else if(length(tmp) == 6){
        interaction_matrix[tmp[1], tmp[2] ] <- interaction_matrix[tmp[1], tmp[2] ] + 1
        interaction_matrix[tmp[1], tmp[3] ] <- interaction_matrix[tmp[1], tmp[3] ] + 1
        interaction_matrix[tmp[1], tmp[4] ] <- interaction_matrix[tmp[1], tmp[4] ] + 1
        interaction_matrix[tmp[1], tmp[5] ] <- interaction_matrix[tmp[1], tmp[5] ] + 1
        interaction_matrix[tmp[1], tmp[6] ] <- interaction_matrix[tmp[1], tmp[6] ] + 1
        
        interaction_matrix[tmp[2], tmp[3] ] <- interaction_matrix[tmp[2], tmp[3] ] + 1
        interaction_matrix[tmp[2], tmp[4] ] <- interaction_matrix[tmp[2], tmp[4] ] + 1
        interaction_matrix[tmp[2], tmp[5] ] <- interaction_matrix[tmp[2], tmp[5] ] + 1
        interaction_matrix[tmp[2], tmp[6] ] <- interaction_matrix[tmp[2], tmp[6] ] + 1
        
        interaction_matrix[tmp[3], tmp[4] ] <- interaction_matrix[tmp[3], tmp[4] ] + 1
        interaction_matrix[tmp[3], tmp[5] ] <- interaction_matrix[tmp[3], tmp[5] ] + 1
        interaction_matrix[tmp[3], tmp[6] ] <- interaction_matrix[tmp[3], tmp[6] ] + 1
        
        interaction_matrix[tmp[4], tmp[5] ] <- interaction_matrix[tmp[4], tmp[5] ] + 1
        interaction_matrix[tmp[4], tmp[6] ] <- interaction_matrix[tmp[4], tmp[6] ] + 1
        
        interaction_matrix[tmp[5], tmp[6] ] <- interaction_matrix[tmp[5], tmp[6] ] + 1
      } else if(length(tmp) == 7){
        interaction_matrix[tmp[1], tmp[2] ] <- interaction_matrix[tmp[1], tmp[2] ] + 1
        interaction_matrix[tmp[1], tmp[3] ] <- interaction_matrix[tmp[1], tmp[3] ] + 1
        interaction_matrix[tmp[1], tmp[4] ] <- interaction_matrix[tmp[1], tmp[4] ] + 1
        interaction_matrix[tmp[1], tmp[5] ] <- interaction_matrix[tmp[1], tmp[5] ] + 1
        interaction_matrix[tmp[1], tmp[6] ] <- interaction_matrix[tmp[1], tmp[6] ] + 1
        interaction_matrix[tmp[1], tmp[7] ] <- interaction_matrix[tmp[1], tmp[7] ] + 1
        
        interaction_matrix[tmp[2], tmp[3] ] <- interaction_matrix[tmp[2], tmp[3] ] + 1
        interaction_matrix[tmp[2], tmp[4] ] <- interaction_matrix[tmp[2], tmp[4] ] + 1
        interaction_matrix[tmp[2], tmp[5] ] <- interaction_matrix[tmp[2], tmp[5] ] + 1
        interaction_matrix[tmp[2], tmp[6] ] <- interaction_matrix[tmp[2], tmp[6] ] + 1
        interaction_matrix[tmp[2], tmp[7] ] <- interaction_matrix[tmp[2], tmp[7] ] + 1
        
        interaction_matrix[tmp[3], tmp[4] ] <- interaction_matrix[tmp[3], tmp[4] ] + 1
        interaction_matrix[tmp[3], tmp[5] ] <- interaction_matrix[tmp[3], tmp[5] ] + 1
        interaction_matrix[tmp[3], tmp[6] ] <- interaction_matrix[tmp[3], tmp[6] ] + 1
        interaction_matrix[tmp[3], tmp[7] ] <- interaction_matrix[tmp[3], tmp[7] ] + 1
        
        interaction_matrix[tmp[4], tmp[5] ] <- interaction_matrix[tmp[4], tmp[5] ] + 1
        interaction_matrix[tmp[4], tmp[6] ] <- interaction_matrix[tmp[4], tmp[6] ] + 1
        interaction_matrix[tmp[4], tmp[7] ] <- interaction_matrix[tmp[4], tmp[7] ] + 1
        
        interaction_matrix[tmp[5], tmp[6] ] <- interaction_matrix[tmp[5], tmp[6] ] + 1
        interaction_matrix[tmp[5], tmp[7] ] <- interaction_matrix[tmp[5], tmp[7] ] + 1
        
        interaction_matrix[tmp[6], tmp[7] ] <- interaction_matrix[tmp[6], tmp[7] ] + 1
      } else if(length(tmp) == 8){
        interaction_matrix[tmp[1], tmp[2] ] <- interaction_matrix[tmp[1], tmp[2] ] + 1
        interaction_matrix[tmp[1], tmp[3] ] <- interaction_matrix[tmp[1], tmp[3] ] + 1
        interaction_matrix[tmp[1], tmp[4] ] <- interaction_matrix[tmp[1], tmp[4] ] + 1
        interaction_matrix[tmp[1], tmp[5] ] <- interaction_matrix[tmp[1], tmp[5] ] + 1
        interaction_matrix[tmp[1], tmp[6] ] <- interaction_matrix[tmp[1], tmp[6] ] + 1
        interaction_matrix[tmp[1], tmp[7] ] <- interaction_matrix[tmp[1], tmp[7] ] + 1
        interaction_matrix[tmp[1], tmp[8] ] <- interaction_matrix[tmp[1], tmp[8] ] + 1
        
        interaction_matrix[tmp[2], tmp[3] ] <- interaction_matrix[tmp[2], tmp[3] ] + 1
        interaction_matrix[tmp[2], tmp[4] ] <- interaction_matrix[tmp[2], tmp[4] ] + 1
        interaction_matrix[tmp[2], tmp[5] ] <- interaction_matrix[tmp[2], tmp[5] ] + 1
        interaction_matrix[tmp[2], tmp[6] ] <- interaction_matrix[tmp[2], tmp[6] ] + 1
        interaction_matrix[tmp[2], tmp[7] ] <- interaction_matrix[tmp[2], tmp[7] ] + 1
        interaction_matrix[tmp[2], tmp[8] ] <- interaction_matrix[tmp[2], tmp[8] ] + 1
        
        interaction_matrix[tmp[3], tmp[4] ] <- interaction_matrix[tmp[3], tmp[4] ] + 1
        interaction_matrix[tmp[3], tmp[5] ] <- interaction_matrix[tmp[3], tmp[5] ] + 1
        interaction_matrix[tmp[3], tmp[6] ] <- interaction_matrix[tmp[3], tmp[6] ] + 1
        interaction_matrix[tmp[3], tmp[7] ] <- interaction_matrix[tmp[3], tmp[7] ] + 1
        interaction_matrix[tmp[3], tmp[8] ] <- interaction_matrix[tmp[3], tmp[8] ] + 1
        
        interaction_matrix[tmp[4], tmp[5] ] <- interaction_matrix[tmp[4], tmp[5] ] + 1
        interaction_matrix[tmp[4], tmp[6] ] <- interaction_matrix[tmp[4], tmp[6] ] + 1
        interaction_matrix[tmp[4], tmp[7] ] <- interaction_matrix[tmp[4], tmp[7] ] + 1
        interaction_matrix[tmp[4], tmp[8] ] <- interaction_matrix[tmp[4], tmp[8] ] + 1
        
        interaction_matrix[tmp[5], tmp[6] ] <- interaction_matrix[tmp[5], tmp[6] ] + 1
        interaction_matrix[tmp[5], tmp[7] ] <- interaction_matrix[tmp[5], tmp[7] ] + 1
        interaction_matrix[tmp[5], tmp[8] ] <- interaction_matrix[tmp[5], tmp[8] ] + 1
        
        interaction_matrix[tmp[6], tmp[7] ] <- interaction_matrix[tmp[6], tmp[7] ] + 1
        interaction_matrix[tmp[6], tmp[8] ] <- interaction_matrix[tmp[6], tmp[8] ] + 1
        
        interaction_matrix[tmp[7], tmp[8] ] <- interaction_matrix[tmp[7], tmp[8] ] + 1
      } else if(length(tmp) == 9){
        interaction_matrix[tmp[1], tmp[2] ] <- interaction_matrix[tmp[1], tmp[2] ] + 1
        interaction_matrix[tmp[1], tmp[3] ] <- interaction_matrix[tmp[1], tmp[3] ] + 1
        interaction_matrix[tmp[1], tmp[4] ] <- interaction_matrix[tmp[1], tmp[4] ] + 1
        interaction_matrix[tmp[1], tmp[5] ] <- interaction_matrix[tmp[1], tmp[5] ] + 1
        interaction_matrix[tmp[1], tmp[6] ] <- interaction_matrix[tmp[1], tmp[6] ] + 1
        interaction_matrix[tmp[1], tmp[7] ] <- interaction_matrix[tmp[1], tmp[7] ] + 1
        interaction_matrix[tmp[1], tmp[8] ] <- interaction_matrix[tmp[1], tmp[8] ] + 1
        interaction_matrix[tmp[1], tmp[9] ] <- interaction_matrix[tmp[1], tmp[9] ] + 1
        
        interaction_matrix[tmp[2], tmp[3] ] <- interaction_matrix[tmp[2], tmp[3] ] + 1
        interaction_matrix[tmp[2], tmp[4] ] <- interaction_matrix[tmp[2], tmp[4] ] + 1
        interaction_matrix[tmp[2], tmp[5] ] <- interaction_matrix[tmp[2], tmp[5] ] + 1
        interaction_matrix[tmp[2], tmp[6] ] <- interaction_matrix[tmp[2], tmp[6] ] + 1
        interaction_matrix[tmp[2], tmp[7] ] <- interaction_matrix[tmp[2], tmp[7] ] + 1
        interaction_matrix[tmp[2], tmp[8] ] <- interaction_matrix[tmp[2], tmp[8] ] + 1
        interaction_matrix[tmp[2], tmp[9] ] <- interaction_matrix[tmp[2], tmp[9] ] + 1
        
        interaction_matrix[tmp[3], tmp[4] ] <- interaction_matrix[tmp[3], tmp[4] ] + 1
        interaction_matrix[tmp[3], tmp[5] ] <- interaction_matrix[tmp[3], tmp[5] ] + 1
        interaction_matrix[tmp[3], tmp[6] ] <- interaction_matrix[tmp[3], tmp[6] ] + 1
        interaction_matrix[tmp[3], tmp[7] ] <- interaction_matrix[tmp[3], tmp[7] ] + 1
        interaction_matrix[tmp[3], tmp[8] ] <- interaction_matrix[tmp[3], tmp[8] ] + 1
        interaction_matrix[tmp[3], tmp[9] ] <- interaction_matrix[tmp[3], tmp[9] ] + 1
        
        interaction_matrix[tmp[4], tmp[5] ] <- interaction_matrix[tmp[4], tmp[5] ] + 1
        interaction_matrix[tmp[4], tmp[6] ] <- interaction_matrix[tmp[4], tmp[6] ] + 1
        interaction_matrix[tmp[4], tmp[7] ] <- interaction_matrix[tmp[4], tmp[7] ] + 1
        interaction_matrix[tmp[4], tmp[8] ] <- interaction_matrix[tmp[4], tmp[8] ] + 1
        interaction_matrix[tmp[4], tmp[9] ] <- interaction_matrix[tmp[4], tmp[9] ] + 1
        
        interaction_matrix[tmp[5], tmp[6] ] <- interaction_matrix[tmp[5], tmp[6] ] + 1
        interaction_matrix[tmp[5], tmp[7] ] <- interaction_matrix[tmp[5], tmp[7] ] + 1
        interaction_matrix[tmp[5], tmp[8] ] <- interaction_matrix[tmp[5], tmp[8] ] + 1
        interaction_matrix[tmp[5], tmp[9] ] <- interaction_matrix[tmp[5], tmp[9] ] + 1
        
        interaction_matrix[tmp[6], tmp[7] ] <- interaction_matrix[tmp[6], tmp[7] ] + 1
        interaction_matrix[tmp[6], tmp[8] ] <- interaction_matrix[tmp[6], tmp[8] ] + 1
        interaction_matrix[tmp[6], tmp[9] ] <- interaction_matrix[tmp[6], tmp[9] ] + 1
        
        interaction_matrix[tmp[7], tmp[8] ] <- interaction_matrix[tmp[7], tmp[8] ] + 1
        interaction_matrix[tmp[7], tmp[9] ] <- interaction_matrix[tmp[7], tmp[9] ] + 1
        
        interaction_matrix[tmp[8], tmp[9] ] <- interaction_matrix[tmp[8], tmp[9] ] + 1
      } else if(length(tmp) >= 10){
        interaction_matrix[tmp[1], tmp[2] ] <- interaction_matrix[tmp[1], tmp[2] ] + 1
        interaction_matrix[tmp[1], tmp[3] ] <- interaction_matrix[tmp[1], tmp[3] ] + 1
        interaction_matrix[tmp[1], tmp[4] ] <- interaction_matrix[tmp[1], tmp[4] ] + 1
        interaction_matrix[tmp[1], tmp[5] ] <- interaction_matrix[tmp[1], tmp[5] ] + 1
        interaction_matrix[tmp[1], tmp[6] ] <- interaction_matrix[tmp[1], tmp[6] ] + 1
        interaction_matrix[tmp[1], tmp[7] ] <- interaction_matrix[tmp[1], tmp[7] ] + 1
        interaction_matrix[tmp[1], tmp[8] ] <- interaction_matrix[tmp[1], tmp[8] ] + 1
        interaction_matrix[tmp[1], tmp[9] ] <- interaction_matrix[tmp[1], tmp[9] ] + 1
        interaction_matrix[tmp[1], tmp[10] ] <- interaction_matrix[tmp[1], tmp[10] ] + 1
        
        interaction_matrix[tmp[2], tmp[3] ] <- interaction_matrix[tmp[2], tmp[3] ] + 1
        interaction_matrix[tmp[2], tmp[4] ] <- interaction_matrix[tmp[2], tmp[4] ] + 1
        interaction_matrix[tmp[2], tmp[5] ] <- interaction_matrix[tmp[2], tmp[5] ] + 1
        interaction_matrix[tmp[2], tmp[6] ] <- interaction_matrix[tmp[2], tmp[6] ] + 1
        interaction_matrix[tmp[2], tmp[7] ] <- interaction_matrix[tmp[2], tmp[7] ] + 1
        interaction_matrix[tmp[2], tmp[8] ] <- interaction_matrix[tmp[2], tmp[8] ] + 1
        interaction_matrix[tmp[2], tmp[9] ] <- interaction_matrix[tmp[2], tmp[9] ] + 1
        interaction_matrix[tmp[2], tmp[10] ] <- interaction_matrix[tmp[2], tmp[10] ] + 1
        
        interaction_matrix[tmp[3], tmp[4] ] <- interaction_matrix[tmp[3], tmp[4] ] + 1
        interaction_matrix[tmp[3], tmp[5] ] <- interaction_matrix[tmp[3], tmp[5] ] + 1
        interaction_matrix[tmp[3], tmp[6] ] <- interaction_matrix[tmp[3], tmp[6] ] + 1
        interaction_matrix[tmp[3], tmp[7] ] <- interaction_matrix[tmp[3], tmp[7] ] + 1
        interaction_matrix[tmp[3], tmp[8] ] <- interaction_matrix[tmp[3], tmp[8] ] + 1
        interaction_matrix[tmp[3], tmp[9] ] <- interaction_matrix[tmp[3], tmp[9] ] + 1
        interaction_matrix[tmp[3], tmp[10] ] <- interaction_matrix[tmp[3], tmp[10] ] + 1
        
        interaction_matrix[tmp[4], tmp[5] ] <- interaction_matrix[tmp[4], tmp[5] ] + 1
        interaction_matrix[tmp[4], tmp[6] ] <- interaction_matrix[tmp[4], tmp[6] ] + 1
        interaction_matrix[tmp[4], tmp[7] ] <- interaction_matrix[tmp[4], tmp[7] ] + 1
        interaction_matrix[tmp[4], tmp[8] ] <- interaction_matrix[tmp[4], tmp[8] ] + 1
        interaction_matrix[tmp[4], tmp[9] ] <- interaction_matrix[tmp[4], tmp[9] ] + 1
        interaction_matrix[tmp[4], tmp[10] ] <- interaction_matrix[tmp[4], tmp[10] ] + 1
        
        interaction_matrix[tmp[5], tmp[6] ] <- interaction_matrix[tmp[5], tmp[6] ] + 1
        interaction_matrix[tmp[5], tmp[7] ] <- interaction_matrix[tmp[5], tmp[7] ] + 1
        interaction_matrix[tmp[5], tmp[8] ] <- interaction_matrix[tmp[5], tmp[8] ] + 1
        interaction_matrix[tmp[5], tmp[9] ] <- interaction_matrix[tmp[5], tmp[9] ] + 1
        interaction_matrix[tmp[5], tmp[10] ] <- interaction_matrix[tmp[5], tmp[10] ] + 1
        
        interaction_matrix[tmp[6], tmp[7] ] <- interaction_matrix[tmp[6], tmp[7] ] + 1
        interaction_matrix[tmp[6], tmp[8] ] <- interaction_matrix[tmp[6], tmp[8] ] + 1
        interaction_matrix[tmp[6], tmp[9] ] <- interaction_matrix[tmp[6], tmp[9] ] + 1
        interaction_matrix[tmp[6], tmp[10] ] <- interaction_matrix[tmp[6], tmp[10] ] + 1
        
        interaction_matrix[tmp[7], tmp[8] ] <- interaction_matrix[tmp[7], tmp[8] ] + 1
        interaction_matrix[tmp[7], tmp[9] ] <- interaction_matrix[tmp[7], tmp[9] ] + 1
        interaction_matrix[tmp[7], tmp[10] ] <- interaction_matrix[tmp[7], tmp[10] ] + 1
        
        interaction_matrix[tmp[8], tmp[9] ] <- interaction_matrix[tmp[8], tmp[9] ] + 1
        interaction_matrix[tmp[8], tmp[10] ] <- interaction_matrix[tmp[8], tmp[10] ] + 1
        
        interaction_matrix[tmp[9], tmp[10] ] <- interaction_matrix[tmp[9], tmp[10] ] + 1
      }
      
      # add self-interactions (just count how many times each cell type shows up...)
      for(tmp.type in tmp){
        interaction_matrix[tmp.type,tmp.type] <- interaction_matrix[tmp.type,tmp.type] + 1
      }
    }
    
    interaction_matrix <- interaction_matrix + t(interaction_matrix)
    colnames(interaction_matrix)
    
    interaction_matrix[lower.tri(interaction_matrix)] <- NA 
    
    return(interaction_matrix)
  }
)

# Normalize each cell type by how often it is detected
interaction.list <- lapply(
  interaction.list,
  FUN = function(X){
    for(i in 1:nrow(X)){
      X[i,] <- X[i,]/diag(X)[i] # normalize each row by the diagonal value
    }
    return(X)
  }
)

# Set self-interactions to NA - helps with plotting scales
interaction.list <- lapply(
  interaction.list,
  FUN = function(X){
    diag(X) <- NA
    return(X)
  }
)

#      interaction heatmaps, Fig. 4e ####
injury=c("D2", "D5", "D7")

inter.heat <- list()
for(i in 1:length(injury)){
  draw.legend = injury[i]=="D7"
  inter.heat[[i]] <- ggplotify::as.ggplot(
    pheatmap(
      interaction.list[[i]][rev(rownames(interaction.list[[i]]))[-1],c(-1)], 
      color=viridis(1000, option="magma", direction=-1),
      border_color = NA,
      main=injury[i],
      legend = draw.legend, 
      
      show_rownames = draw.legend, show_colnames = draw.legend,
      labels_row = stringr::str_replace_all(rev(rownames(interaction.list[[i]])),"-"," ")[-1],
      labels_col = stringr::str_replace_all(colnames(interaction.list[[i]]),"-"," ")[-1],
      fontsize_row = 5, fontsize_col= 5,
      
      fontsize=small.font, 
      cluster_rows = F, cluster_cols = F,
      na_col = "white"
    )
  )+
    theme(plot.margin = unit(rep(0,4),"inches"))
}

wrap_plots(inter.heat, nrow=1, widths=c(1,1,1.7))

wrap_plots(
  wrap_plots(inter.heat[c(1,2)], ncol=1),
  inter.heat[[3]],
  widths=c(1,2.5)
)
#
#     line plot, Fig. 4f ####

# prep df
df <- list()
tmp.injury <- c(2, 5, 7)
for(i in 1: length(interaction.list)){
  df[[i]] <- melt(interaction.list[[i]])
  df[[i]]$injury <- tmp.injury[i]
}
df <- do.call(rbind, df)
df <- df[!is.na(df$value),]
colnames(df) <- c("row.type", "col.type", "score", "injury")
df$cat.type <- paste0(df$row.type,"_", df$col.type)


# plots
tmp.interactions <- c(
  paste0(
    c("Quiescent-MuSCs", "Activated-MuSCs", "Committed-Myoblasts", "Fusing-Myocytes"),
    "_Monocyte-(Patrolling)"
  )%>%c(),
  paste0(
    c("Quiescent-MuSCs", "Activated-MuSCs", "Committed-Myoblasts", "Fusing-Myocytes"),
    "_M2-Macro.-(Cx3cr1-lo)"
  )%>%c(),
  paste0(
    c("Quiescent-MuSCs", "Activated-MuSCs", "Committed-Myoblasts", "Fusing-Myocytes"),
    "_M2-Macro.-(Cx3cr1-hi)"
  )%>%c()
)

interaction.macro.line <- ggplot(
  df[df$cat.type %in% tmp.interactions,],
  aes(
    x=injury, y=score, 
    col=row.type
  )
) +
  geom_line() + 
  geom_point(size=0.75 ) +
  theme_minimal()+
  fig1bcd.theme +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(size = 0.25),
    panel.grid.minor = element_blank()
  )+
  ylim(c(0,0.45))+
  scale_color_manual(
    values= bin.colors[c(3,5,10,15)]
  )+
  labs(
    x="Days Post-Injury",
    y="Normalized Co-occurrence",
    color="Myogenic Cell Type"
  )+
  facet_wrap(facets="col.type")
interaction.macro.line

tmp.interactions <- c(
  paste0(
    c("Quiescent-MuSCs", "Activated-MuSCs", "Committed-Myoblasts", "Fusing-Myocytes"),
    "_FAPs-(Stem)"
  ),
  paste0(
    c("Quiescent-MuSCs", "Activated-MuSCs", "Committed-Myoblasts", "Fusing-Myocytes"),
    "_FAPs-(Adipogenic)"
  ),
  paste0(
    c("Quiescent-MuSCs", "Activated-MuSCs", "Committed-Myoblasts", "Fusing-Myocytes"),
    "_FAPs-(Pro-remodeling)"
  )
)

interaction.faps.line <- ggplot(
  df[df$cat.type %in% tmp.interactions,],
  aes(
    x=injury, y=score, 
    col=row.type
  )
) +
  geom_line() + 
  geom_point(size=0.75 ) +
  theme_minimal()+
  fig1bcd.theme +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(size = 0.25),
    panel.grid.minor = element_blank()
  )+
  ylim(c(0,0.45))+
  scale_color_manual(
    values= bin.colors[c(3,5,10,15)]
  )+
  labs(
    x="Days Post-Injury",
    y="Normalized Co-occurrence",
    color="Myogenic Cell Type"
  )+
  facet_wrap(facets="col.type")
interaction.faps.line


wrap_plots(
  interaction.macro.line,
  interaction.faps.line,
  ncol=1,
  guides="collect"
)

#    ####







