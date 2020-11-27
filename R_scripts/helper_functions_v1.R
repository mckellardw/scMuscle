
# Helper functions used in each of the other scripts


# Calculate the number of PCs that contain some proportion (95%) of the variance
#     Output: number of dimensions that contain the desired percentage of variance
npcs <- function(
  seu, 
  var.toal=0.95, 
  reduction="pca"
){
  if(is.null(seu@reductions[[reduction]])){
    cat("Reduction", reduction, "not found!")
    return(NULL)
  }
  
  tmp.var <- (seu@reductions[[reduction]]@stdev)^2
  var.cut <- var.toal*sum(tmp.var)
  n.pcs=0
  var.sum = 0
  while(var.sum < var.cut){
    n.pcs = n.pcs + 1
    var.sum <- var.sum + tmp.var[n.pcs]
  }
  
  return(n.pcs)
}


# Function to generate mean expression profiles for individual cell types
#     Output: genes-by-celltype matrix of mean expression profiles
get_cell_type_model <- 
  function(
    seur, #seurat object
    assay='RNA',
    slot='data',
    cluster.names, #cell type model; pass in metadata column name for desired cell types
    ignore.clusters=NULL, # vector of chars; cluster IDs to ignore in model generation
    cells.use=NULL, # vector of cell barcodes to use; works in combination with ignore.clusters, doublet.classifications
    genes.ignore=NULL, # genes to remove from model
    doublet.classifications=NULL, # name of doublet classification metadata column; default is to ignore classifications
    # nCores=1, #TODO: parallelize with dofor
    verbose=T
  ){
    require(Seurat)
    
    #build singlet_cell type_gene expression matrix
    ref.mat <- list() #initialize reference expression matrix
    exp.mat <- GetAssayData(seur, assay=assay, slot=slot) # Pull out expression data
    
    # subset based on desired cells
    if(!is.null(cells.use)){ 
      exp.mat <- exp.mat[,cells.use]
    }
    
    # remove ignored genes
    if(!is.null(genes.ignore)){
      exp.mat = exp.mat[!rownames(seur)%in%genes.ignore,]
    }
    
    #Grab meta.data, filtered for singlets
    if(is.null(cells.use)){ #TODO- subset based on cells.use
      meta <- seur@meta.data  #ignore sinlget/doublet classifications
    }else{
      meta <- seur@meta.data[cells.use,] 
    }
    
    cell.types <- sort(unlist(unique(meta[[cluster.names]]))) #Get cell type names
    if(!is.null(ignore.clusters)){
      cell.types <- cell.types[(!cell.types%in%ignore.clusters)]  
    }
    print(cell.types)
    
    # Remove NAs from exp.mat and meta
    exp.mat <- exp.mat[,!is.na(meta[[cluster.names]])]
    meta <- meta[!is.na(meta[[cluster.names]]),]
    
    # Generate cell type references, one cell type at a time
    if(verbose){cat("Building reference profiles... \n")}
    if(verbose & !is.null(cells.use)){cat("   Using the provided cells... \n")}
    
    for(i in 1:length(cell.types)){
      tmp.cells = rownames(meta)[ as.vector(meta[[cluster.names]])==cell.types[i] ]
      
      if(verbose){cat(cell.types[i]," (",i,"/",length(cell.types),")","... ", sep="")}
      
      if(length(tmp.cells)==0){
        if(verbose){cat(" NO CELLS", sep="")}
      }
      
      if(length(tmp.cells)==1){
        if(verbose){cat(" Only 1 cell!", sep="")}
        ref.mat[[i]] <- exp.mat[,rownames(meta)%in%tmp.cells]  # grab expression data
      }else{
        if(verbose){cat(length(tmp.cells)," cells.", sep="")}
        ref.mat[[i]] <- 
          Matrix::rowMeans( #TODO: should these expression profiles be built on row means?
            exp.mat[,rownames(meta)%in%tmp.cells]  # grab expression data
            
          )
      }
      
      if(verbose){cat("\n")}
    }
    ref.mat <- do.call(cbind, ref.mat)
    colnames(ref.mat) <- cell.types
    
    if(verbose){cat("Done!\n")}
    
    return(ref.mat)
  }



