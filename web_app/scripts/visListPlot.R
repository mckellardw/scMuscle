
# Generate feature plots given
visListPlot <- function(
  seu.list,
  features=NULL,
  alt.titles=NULL, # alternative titles for genes/features being passed
  assay='Spatial',
  reduction="space",
  pt.size=1,
  font.size=8
){
  require(Seurat)
  require(ggplot2)
  require(viridis)
  
  if(is.null(alt.titles)){
    alt.titles=features
  }
  
  seu.list <- lapply(
    seu.list,
    FUN = function(SEU){
      SEU@active.assay <- assay
      return(SEU)
    }
  )
  
  # Get expression limits for each gene, across all datasets 
  gene.lims <- lapply(
    features,
    FUN = function(FEAT){
      out.max <- lapply(
        seu.list,
        FUN = function(SEU) max(GetAssayData(SEU,assay=assay)[FEAT,])
      ) %>% unlist() %>% max()
      return(c(0,out.max))
    }
  )

  plot.list <- list()
  for(i in 1:length(features)){
    tmp <- lapply(
      seu.list,
      FUN = function(SEU)
        FeaturePlot(
          SEU,
          slot ="data",
          features = features[i],
          pt.size = pt.size,
          reduction=reduction
        ) +
        scale_color_viridis(limits=unlist(gene.lims[i]), na.value = gray(0.42))+ #,trans="log10"  limits=c(10^-2,6),
        # scale_color_gradientn(colors=RColorBrewer::brewer.pal(11,"Spectral")[11:1],limits=c(10^-2,1),na.value = gray(0.42)) +
        theme(
          plot.margin = unit(rep(0,4), "inches"),
          axis.ticks = element_blank(),
          axis.text=element_blank(),
          axis.title = element_blank(),
          axis.line=element_blank(),
          plot.title = element_blank(),
          legend.position="bottom",
          legend.title = element_text(size=font.size,face="bold", hjust=0.5),
          legend.text = element_text(size=font.size,face="bold")
        )
      # labs(color="Log-Normalized\nExpression")
    )
    tmp[[1]] <- tmp[[1]] +
      theme(
        plot.title = element_text(size=font.size,face="bold.italic",  vjust=1)
      ) +
      labs(title=alt.titles[i])
    plot.list[[i]] <- tmp
  }

  injury=c("D2", "D5", "D7")
  for(i in 1:length(plot.list[[1]]) ){
    plot.list[[1]][[i]] <- plot.list[[1]][[i]] +
      labs(y=injury[i]) +
      theme(axis.title.y = element_text(size=font.size, face="bold", color="black"))
  }

  plot.list <- lapply(
    plot.list,
    FUN = function(X)
      wrap_plots(X, ncol=1, guides="collect")&theme(legend.position="bottom",legend.margin = margin(0,0,0,0,"inches"))
  )

  return(
    wrap_plots(plot.list,nrow=1)
  )
}
