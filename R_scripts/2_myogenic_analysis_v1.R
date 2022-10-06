# Libs, setwd #####
library(Matrix)
library(dplyr)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(phateR)

library(cluster)
library(data.table)

#Plotting
library(ggplot2)
library(patchwork)
library(pals)
library(viridis)
library(shades)

source("helper_functions_v1.R") 

# Load table of sample metadata  ####

meta <- fread("supplemental_data/sample_metadata.csv")

# Figure settings ####
small.font = 6
big.font = 8
line.width = 0.5
pt.size=0.1
label.size=2
mon.cell.size=0.1
mon.graph.size=0.5
mon.label.size=1

celltypes = c(
  "B Cells","NK Cells","T Cells",
  "Neutrophils","Monocyte (Patrolling)","Monocyte (Inflammatory)",
  "Monocyte (Cxcl10+)","M1 Macro.","M2 Macro. (Cx3cr1_lo)",
  "M2 Macro. (Cx3cr1_hi)","Dendritic","Endothelial (Stem)",
  "Endothelial (Capillary)","Endothelial (Artery)","Endothelial (Vein)",
  "Smooth Muscle","Neural","FAPs (Stem)",
  "FAPs (Pro-remodeling)","FAPs (Adipogenic)","Tenocytes",
  "MuSCs","Myoblasts/Progenitors","Myonuclei (Type IIx)",
  "Myonuclei (Type IIb)","Monocyte (Infiltrating)" ,"Myoblasts_1",
  "Myoblasts_2","Committed Progenitors","Myonuclei (Cell Prep)",
  "Pericytes","Activated MuSCs","Myoblasts",
  "Myonuclei (P10)"  
)
celltype.colors = as.vector(polychrome())[c(3:9,17,11:13,21,14:16,32,18:20,28,10,30,29,31,33,34,35,36,22:27)] 
names(celltype.colors) <- celltypes

source.colors = as.vector(kelly())[5:22]
doublet.colors = c("#8A9197FF","#C80813FF")
age.colors = viridis(n=length(unique(myo.seurat$age.months)),option = "B")
type.colors = c("#24325F99","#82491E99","#526E2D99","#FB646799","#E762D799","#E8924299","#917C5D")
chemistry.colors = c("#FD7446","#709AE1FF")  
sex.colors = c("#9c4872", "#67a0d6","white") 
injury.agent.colors = c("salmon4", "peru")
injury.colors = c(viridis(length(unique(myo.seurat$injury))-1, option="plasma"), "gray")
is.injured.colors =c("snow4","seagreen")
bin.colors = rainbow(nbins)%>%saturation(values=0.75) %>% brightness(values=0.9) %>% as.vector()


# Plot themes ####
umap.theme <- theme(
  axis.line = element_line(    color = "black", size = line.width),
  axis.title = element_text(face='bold',size = small.font, hjust = 0, vjust = 1), 
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  legend.text = element_text(size=small.font),
  legend.title = element_text(size=big.font, face="bold", hjust=0.5)
)

fig1bcd.theme <- theme_minimal()+
  theme(
    axis.line = element_blank(),
    panel.grid.major = element_line(color='gray'),
    panel.grid.minor = element_line(color='light gray'),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.title.x = element_text(face = 'bold',size = small.font, hjust = 0.5, vjust = 0.5), 
    axis.title.y = element_text(face = 'bold',size = small.font, hjust = 0.5, vjust = 0.5),
    axis.text.x=element_text(
      size = small.font, color = "black",hjust=0.5,
      margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm")
    ),
    axis.text.y=element_text(
      size = small.font, angle=90,color = "black",hjust=0.5,
      margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm")
    ),
    axis.ticks.length = unit(-2, "mm"),
    axis.ticks = element_line(size = .5, color="black"),
    
    legend.background = element_rect(color = "black", fill='white', size=0.5),
    # legend.key.size = unit(0.5,'mm'),
    legend.text = element_text(size = small.font, hjust = 0, vjust = 0.5),
    legend.title = element_text(face = 'bold',size = small.font, hjust = 0.5, vjust = 0.5),
    legend.position = c(0.2,0.6)
  ) 
fig1bcd.labs <- labs(
  col='Sample Prep',
  shape='Chemistry', 
  size = "Median UMI Counts per Cell"# "Sequencing\nSaturation (%)"
) 

bar.theme <- theme_minimal() +
  theme(
    title=element_text(hjust=0.5,face="bold",size = big.font),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color='black',size = 1),
    axis.title.x = element_text(face = 'bold',size = small.font, hjust = 0.5, vjust = 1), 
    axis.title.y = element_text(face = 'bold',size = small.font, hjust = 0.5, vjust = 0.5),
    axis.text.x=element_text(
      size = small.font, color = "black",angle = 90,hjust=0.5,
      margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm")
    ),
    axis.text.y=element_text(
      size = small.font, color = "black",angle = 90,hjust=0.5,
      margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm")
    ),
    axis.ticks.length = unit(-2, "mm"),
    axis.ticks = element_line(size = .5, color="black")
  )+ NoLegend()

pie.theme <- 
  theme_minimal() + 
  theme(
    panel.border = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    title=element_text(face="bold", hjust=0.5,vjust = 0.5, size=big.font)
  ) + NoLegend()

vln.theme <- theme(
  panel.background = element_blank(),
  axis.line=element_line(color="black"),
  # panel.grid.major = element_line(color='gray'),
  # panel.grid.minor = element_line(color='light gray'),
  legend.text = element_text(size=small.font, color="black"),
  panel.grid = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y=element_text(size=big.font, color="black"),
  axis.text = element_text(color="black",size=small.font),
  axis.ticks.x = element_blank()
)

dot.theme <- theme(
  axis.line = element_blank(),
  panel.border = element_rect(colour = "black", fill=NA, size=1),
  axis.title.x = element_blank(), 
  axis.title.y = element_blank(),
  axis.text.x=element_text(angle = 90,hjust = 1,vjust= 0.5, size=small.font),
  axis.text.y=element_text(angle=0,hjust = 1,vjust= 0.5, size=small.font),
  legend.text = element_text(size=small.font,color="black"),
  legend.title = element_text(size=big.font,color="black"),
  panel.grid.major = element_line(colour = "gray", size = 0.5)
)
#
# Subset to each celltype compartment and re-run UMAP on narrowed dataset #####
DefaultAssay(scMuscle.clean.seurat) <- 'RNA'
Idents(scMuscle.clean.seurat) <- 'harmony_factorIDs'

myo.seurat <- 
  subset(
    x= scMuscle.clean.seurat,
    subset = harmony_superIDs %in% c("Myogenic", "Myonuclei") &
      bbknn_superIDs %in% c("Myogenic", "Myonuclei") &
      scanorama_superIDs %in% c("Myogenic", "Myonuclei")
  )

# re-run UMAP for cleaner visualization
n.pcs <- npcs(myo.seurat, reduction="harmony")
myo.seurat <- RunUMAP(
  object = myo.seurat, 
  reduction = 'harmony',
  dims = 1:n.pcs,
  reduction.name='umap_harmony',
  reduction.key = 'umapharmony_'
)

# Run PHATE on harmony values ####

tmp.phate <- phate(
  t(GetAssayData(myo.seurat,slot="data")),
  gamma=1,
  n.jobs = 25
) 
colnames(tmp.phate$embedding) <- paste0("phate_", 1:2)
tmp.phate$params$data <- NULL

myo.seurat[["phate"]] <- CreateDimReducObject(
  embeddings = tmp.phate$embedding,
  key = "phate_", 
  assay = 'RNA',
  misc = tmp.phate$params
)

#         binning ####
#     on pseudotime values
nbins = 25
breaks = seq(
  min(myo.seurat$pseudotime),
  max(myo.seurat$pseudotime),
  max(myo.seurat$pseudotime)/nbins
)
binnames = as.factor(1:nbins)

myo.seurat$pseudo.bins <- cut(
  myo.seurat$pseudotime+10^-100,
  breaks = breaks,
  labels=binnames
)

#     binning on phate_harmony values ####

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

#         DGE along phate1 ####
plan("multiprocess", workers = 25)
plan()

options(future.globals.maxSize = 1000 * 1024^2)

# Output is saved in "supplemental_data/myo_phate_bins_markers.csv"
Idents(myo.seurat) <- "phate1.bins"
myo.phate.bin.markers <- FindAllMarkers(
  myo.seurat, 
  verbose=T
)

#           gene biotype lists ####

# http://www.informatics.jax.org/go/term/GO:0003700
# GO Term: "DNA-binding transcription factor activity"
tf.table <- fread("/workdir/dwm269/muscle_data/tf_resources/GO_term_summary_20201117_094707.txt")
tf.genes <- unique(tf.table$`MGI Gene/Marker ID`)
length(tf.genes) # 1362 genes

# Load cell surface markers list 
# https://wlab.ethz.ch/cspa/#abstract
cspa.surface.proteins <- fread("../mouse_cspa_surfaceome_proteins.csv")
surface.genes <- unique(cspa.surface.proteins$`ENTREZ gene symbol`[cspa.surface.proteins$`CSPA category`=="1 - high confidence" ])

#   add gene biotypes to differential expression output
myo.phate.bin.markers$is.surface.protein <- myo.phate.bin.markers$gene %in% surface.genes
myo.phate.bin.markers$is.trascr.factor <- myo.phate.bin.markers$gene %in% tf.genes

# Plotting ####
#     harmony + phate, cell types (Fig. 2a) ####

DimPlot(
  myo.seurat,
  reduction='phate_harmony',
  cells = sample(colnames(myo.seurat)),
  group.by="harmony_factorIDs",
  cols=celltype.colors[celltypes %in% unique(myo.seurat$harmony_factorIDs)],
  pt.size = pt.size,
  label.size = label.size,
  repel = T,label=T
) +
  NoLegend() +
  labs(x="PHATE_harmony_1",y="PHATE_harmony_2")+
  umap.theme

#     phate bins (Fig. 2b) ####

DimPlot(
  myo.seurat,
  reduction='phate_harmony',
  cells = sample(colnames(myo.seurat)),
  group.by="phate1.bins",
  cols=bin.colors,
  pt.size = pt.size,
  label.size = label.size,
  repel = T,label=T
) +
  NoLegend() +
  labs(x="PHATE_Harmony_1",y="PHATE_Harmony_2")+
  umap.theme

#     sample diversity in each bin ####
simpson <- lapply(
  split(
    myo.seurat$sample,
    myo.seurat$phate1.bins 
  ), # split by sample
  FUN = function(X) vegan::diversity(
    as.numeric(as.factor(X)), # convert char to numeric
    index="simpson"
  )
) %>% unlist()

# plot with Simpson diversity
simpson.bins <- ggplot()+
  geom_bar(
    data=myo.seurat@meta.data,
    aes(x=phate1.bins,fill=phate1.bins)
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
    labels = scales::scientific,
    sec.axis = sec_axis( 
      trans = ~(log(., base=10000)/10)+0.9,
      breaks = c(0.9,0.95,1),
      name="Simpson Index\n(SampleID)"
    )
  )+
  scale_fill_manual(
    values=rainbow(nbins,end = 0.9)%>%saturation(values=0.75) %>% brightness(values=0.9) %>% as.vector()
  )+
  vln.theme + 
  theme(
    axis.ticks.x = element_line(color="black"),
    axis.ticks.y = element_line(color="black"),
    axis.title.x = element_text(color="black", face="bold", size=big.font),
    axis.text.x =  element_text(color="black", size=small.font),
    axis.title.y.left = element_text(color="black", face="bold", size=big.font),
    axis.text.y.left = element_text(color="black", size=small.font),
    axis.title.y.right = element_text(color="black", face="bold", size=big.font),
    axis.text.y.right = element_text(color="black", size=small.font)
  ) +
  xlab("PHATE_Harmony_1 (Binned)")

#     PHATE_bins - canonical markers ####
tmp.genes <- c(
  "Pax7", "Myf5","Myod1","Cdkn1c", "Myog","Gm7325", "Tmem8c", "Acta1", "Ckm", "Myh1", "Myh4"
)

can.vlns <- myo.seurat %>% VlnPlot(
  features = tmp.genes,
  slot = "data",
  pt.size = 0,
  group.by="phate1.bins", #
  cols = rainbow(nbins)%>%saturation(values=0.75) %>% brightness(values=0.9) %>% as.vector(),
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
for(i in 1:length(can.vlns)){
  can.vlns[[i]] <-  can.vlns[[i]] + labs(y=tmp.genes[i])
}
can.vlns[[length(can.vlns)]] <- can.vlns[[length(can.vlns)]] +
  theme(
    axis.text.x = element_text(color="black"),
    axis.title.x = element_text(color="black", size=big.font, face="bold"),
    axis.line.x = element_line(color="black"),
    axis.ticks.x = element_line(color="black")
  )+
  scale_y_continuous(expand=c(0.1,0))+
  xlab("PHATE_Harmony_1")

can.vlns[[6]] <- can.vlns[[6]] + ylab("Mymx")
can.vlns[[7]] <- can.vlns[[7]] + ylab("Mymk")

wrap_plots(can.vlns, ncol=1, heights = c(rep(1,length(can.vlns)-1),1.1))

#     Gene diversity boxplots (Fig. 2c) ####
ggplot(
  myo.seurat@meta.data[sample(colnames(myo.seurat)),],
  aes(
    x=phate1.bins,
    y=nFeature_RNA/Sequencing.Saturation
  )
)+
  geom_boxplot(outlier.alpha = 0.5, outlier.size = 0.1)+
  labs(
    x="PHATE_Harmony_1 (Binned)",
    y="# Genes Expressed /\nSequencing Saturation"
  )+
  vln.theme+
  theme(
    axis.title.x = element_text(size=big.font,color="black", face="bold"),
    axis.title.y = element_text(size=big.font,color="black", face="bold"),
    axis.ticks.x = element_line(color="black"),
    axis.ticks.y = element_line(color="black")
  )

#     Injury timepoints of intermediate cell states ####
# subset out intermediate cell states, and plot cell/nucleus counts across each sample (Fig. 3a)
trans.seu <- subset(myo.seurat, subset=phate1.bins %in% c(8:17))

ggplot(
  trans.seu@meta.data
) +
  geom_bar(
    stat="count",
    size=0.2,
    col="black",
    aes(
      x=source.label,
      fill=phate1.bins
    )
  )+
  labs(x="Source", y="# Cells\n(PHATE bins #8-17)", col="Injury")+
  scale_color_viridis_d(option = "plasma")+
  scale_fill_manual(
    values=bin.colors[8:17]
  ) + 
  scale_y_continuous(expand=expansion(c(0.01,0))) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line=element_line(color="black", size=line.width),
    axis.ticks=element_line(color="black", size=line.width),
    axis.title = element_text(color="black",size=big.font, face="bold"),
    axis.text = element_text(color="black",size=small.font),
    legend.title = element_text(color="black",size=big.font),
    legend.text = element_text(color="black",size=small.font),
    axis.text.x=element_text(angle=30, hjust=1, vjust=1)
  )+ 
  NoLegend()

#     DotPlots of transcription factors/surface markers highly expressed in intermediate cell states ####
#           transient dotplot ####
tmp.genes <- list(
  c( #Surface Markers
    "Itgb1", "Vcam1","Cd34","Cxcr4","Itga7",  # canonical markers
      
    # intermediate - surface
    "Cd164","Ppap2a", "Jam3", "Cd97",
    "Fndc5", "Erbb3", "Fam171a2"
    
  ),
  c( #TFs
    "Pax7","Myod1", "Myog",
    
    #intermediate - TFs
    "Id1","Mycl", "Scx",
    "Zbtb18", "Purb",
    "Hes6", "Mef2a",  
    "E2f8",  "Tead4"
  )
)

trans.dot <- lapply(
  tmp.genes,
  FUN=function(GENES)
    myo.seurat %>% DotPlot(
      features = GENES,
      group.by="phate1.bins",
      dot.min = 10^-10, # remove zeroes
      dot.scale=5
    ) +
    coord_flip() +
    scale_color_viridis() +
    dot.theme +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.justification = "center",legend.box = "vertical",
      legend.title=element_text(color="black", face="bold", size=small.font),
      axis.text.x = element_text(angle=0, hjust=0.5),
      axis.text.y = element_text(face="italic"),
      plot.margin = unit(rep(0,4),units="inches"),
      axis.title.x = element_text(color="black", face="bold", size=big.font)
    )+
    labs(y="PHATE_Harmony_1 (Binned)")
)

trans.dot[[1]]<- trans.dot[[1]] + 
  theme(axis.text.x = element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank())

wrap_plots(
  trans.dot, 
  ncol=1, 
  heights=lapply(tmp.genes,length), 
  guides = "collect"
)
