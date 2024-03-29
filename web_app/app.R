# scMuscle Explorer
#   written by Leo Song & David McKellar


# Load libraries and helper functions ----
library(shiny)
library(Seurat)
library(ggplot2)
library(shinythemes)
library(pals)
library(dplyr)
library(patchwork)
library(shiny)
library(shinycssloaders)
library(shades) 


source("scripts/visListPlot.R")

# Settings & version info ----
CURR_DATA_VERSION = 1.1
cornell_red = "#B31B1B"
spatial_theta_colors <- c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4", "#E6F598", "#FFFFBF", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142")
vis.height = "600px"

# Gene lists used to draw the app (prevents need for loading data immediately)
genes_all <- read.csv("resources/genes_all.csv") %>% unlist() %>% unique() %>% sort()
genes_myo <- read.csv("resources/genes_myo.csv") %>% unlist() %>% unique() %>% sort()
genes_vis <- read.csv("resources/genes_vis.csv") %>% unlist() %>% unique() %>% sort()
types_vis <- read.csv("resources/types_vis.csv") %>% unlist() %>% unique() %>% sort()

# User interface ----
ui <- fluidPage(
  theme = shinytheme("paper"),
  navbarPage(
    title="scMuscle",
    id="NAVBAR",
    # Each tabPanel makes a tab on the upper navbar
    # 'About' tab ----
    tabPanel(
      title="About",
      # Title of the app
      fluidRow(
        style = 'border-bottom: 4px solid gray', 
        br(),
        h1(
          "scMuscle",
          align = "center",
          style = "color: #B31B1B; face: bold"
        ),
        h4(
          "The Single-Cell Muscle Project (scMuscle) aims to collect, analyze, and provide to the research community skeletal muscle transcriptomic data.",
          align = "center",
          style = "color: #6B6B6B"
        ),
        h4(
          "Click one of the tabs above to begin viewing the data...",
          align = "center",
          style="color: #B31B1B"
        ),
        align = "center",
        br()
      ),
      # Who we are
      fluidRow(
        style="border-top: 4px",
        column(
          width=6, 
          br(), br(), 
          div(
            img(src = "img3.png", width = 500, height = 100), 
            align = "center"
          ),
          br(),
          h4(
            p(
              paste("Current data version:", CURR_DATA_VERSION, "(details found "),
              a(href = "https://github.com/mckellardw/scMuscle/tree/main/data_versions", "here"),
              ")",
              align = "center"
            ),
            style="color: #6B6B6B"
          )
        ),
        column(
          width=6,
          h4(
            p("This resource was developed by the",
              a(href = "https://cosgrovelab.bme.cornell.edu/", "Cosgrove"),
              "and",
              a(href = "http://devlaminck.bme.cornell.edu/", "De Vlaminck"), 
              "Labs at Cornell University."),
            style="color: #6B6B6B"
          ),
          h4(
            p("Code for this web app is available ",
              a(href = "https://github.com/mckellardw/scMuscle/tree/main/web_app", "here"), 
              ". To report errors or request new features, please post an issue on the ",
              a(href="https://github.com/mckellardw/scMuscle", "Github"), 
              "page or email", span(" scmuscle@cornell.edu", style = "color:#B31B1B")),
            style="color: #6B6B6B"
          ),
          h4(
            p("Please cite this ",
              a(href = "https://www.biorxiv.org/content/10.1101/2020.12.01.407460v1", "preprint"),
              " when using this resource:"),
            style="color: #6B6B6B"
          ),
          h6(
            "McKellar, D. W. et al. Strength in numbers: Large-scale integration of single-cell transcriptomic data reveals rare, transient muscle progenitor cell states in muscle regeneration. bioRxiv 2020.12.01.407460 (2020). doi:10.1101/2020.12.01.407460",
            style="color: #000000"
          ),
          br(), br()
        )
      )
    ),
    
    # 'All Cell Types' tab----
    tabPanel(
      title="All Cell Types",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          # UMAP Panel Inputs----
          conditionalPanel(
            condition = "input.tabselected==1",
            # inputs for reduction to be shown
            br(),
            # helpText("Select a cell type clustering, gene, or metadata feature to visualize in UMAP space"),
            # inputs gene for feature plot
            selectizeInput(
              "gene1",
              label = "Plot 1: Gene expression UMAP... Select gene",
              choices = genes_all,
              selected = "Myod1"
            ),
            actionButton(
              "action3", label = "Flip Color Scale",
              style="color: #B31B1B; background-color: #F7F7F7; border-color: #B31B1B"
            ),
            br(),
            selectizeInput(
              "umap.reduction",
              label = "Plot 2: Cell type UMAP... Select reduction",
              choices = c("Harmony", "BBKNN", "Scanorama"),
              selected = "Harmony"
            ),
            selectizeInput(
              "variables.umap",
              label = "Plot 3: Metadata feature UMAP... Select feature",
              choices = c(
                "source","sample", "chemistry", "injury days", "injury agent",
                "age", "type", "tissue", "mouse strain", "sex", "mice per sample", 
                "sequencing instrument"
              ),
              selected = "source"
            )
          ),
          # Single Violin Panel Inputs----
          conditionalPanel(
            condition = "input.tabselected==2",
            # selects gene to be examined on violin plot
            br(),
            selectizeInput(
              "gene2",
              label = "Select gene",
              choices = genes_all,
              selected = "Myod1"
            ),
            selectizeInput(
              "variables.singleVln",
              label = "Select feature",
              choices = c(
                "source","sample", "chemistry", "injury days", "injury agent", 
                "Harmony cell types", "BBKNN cell types", "Scanorama cell types", 
                "age", "type", "tissue", "mouse strain", "sex", "mice per sample", 
                "sequencing instrument"
              ),
              selected = "Harmony cell types"
            )
          ),
          # Split Violin Panel Inputs----
          conditionalPanel(
            condition = "input.tabselected==3",
            # selects gene to be examined on violin plots
            br(),
            selectizeInput(
              "gene3",
              label = "Select a gene to plot",
              choices = genes_all,
              selected = "Myod1"
            ),
            # selects cell types to label on x axis
            br(),
            selectizeInput(
              "splitviolincelltype",
              label = "Select cell type IDs:",
              choices = c("Harmony", "BBKNN", "Scanorama"),
              selected = "Harmony"
            ),
            selectizeInput(
              "variables.splitVln",
              label = "Select feature",
              choices = c(
                "source","sample", "chemistry", "injury days", "injury agent", 
                "Harmony cell types", "BBKNN cell types", "Scanorama cell types", 
                "age", "type", "tissue", "mouse strain", "sex", "mice per sample", 
                "sequencing instrument"
              ),
              selected = "sex"
            )
          ),
          # Dot Panel Inputs----
          conditionalPanel(
            condition = "input.tabselected==4",
            # dot plot action button
            br(),
            helpText("Click to apply input changes below"),
            actionButton(
              "action1", label = "Generate Dot Plot",
              style="color: #B31B1B; background-color: #F7F7F7; border-color: #B31B1B"
            ),
            # inputs for dot plot
            br(), br(), br(),
            selectizeInput(
              "dot",
              label = "Select gene(s) to display:",
              choices = genes_all,
              selected = "Myod1",
              multiple = TRUE
            ),
            selectizeInput(
              "variables.dot",
              label = "Select feature",
              choices = c(
                "source","sample", "chemistry", "injury days", "injury agent", 
                "Harmony cell types", "BBKNN cell types", "Scanorama cell types", 
                "age", "type", "tissue", "mouse strain", "sex", "mice per sample", 
                "sequencing instrument"
              ),
              selected = "Harmony cell types"
            )
          ),
          # Inputs present in all tabs----
          br(),
         
            # downloadable plot type and dimensions for "All Cell Types" NavBar Tab
          br(),
          helpText("Download Specifications"),
          selectizeInput(
            "downloadable1",
            label = "file type:",
            choices = c("pdf", "png", "eps"),
            selected = "pdf"
          ),
          numericInput(
            "plotsizex1",
            label = "horizontal dimension (inches)",
            value = 5
          ),
          numericInput(
            "plotsizey1",
            label = "vertical dimension (inches)",
            value = 5
          )
        ), #sidebarPanel
        
        # Establishes spaces for plots in the main panel ----
        mainPanel(
          tabsetPanel(
            # umap panel----
            tabPanel(
              "UMAP", value = 1,
              # umap grouped by cell types
              br(),
              downloadButton("down3", label = "Download"),
              br(), br(),
              plotOutput("feature") %>% withSpinner(type = 1, color = cornell_red),
              br(),
              # umap grouped by metadata features
              downloadButton("down2", label = "Download"),
              br(), br(),
              plotOutput("umap") %>% withSpinner(type = 1, color = cornell_red),
              br(),
              # feature plot
              downloadButton("down1", label = "Download"),
              br(), br(),
              plotOutput("grouping") %>% withSpinner(type = 1, color = cornell_red), 
              br() 
            ),
            
            # single violinplot panel----
            tabPanel(
              "Single Violin", value = 2,
              # by different metadata variables
              br(),
              downloadButton("down4", label = "Download"),
              br(), br(),
              plotOutput("violin1") %>% withSpinner(type = 1, color = cornell_red),
              br()
            ),
            # split violinplot panel----
            tabPanel(
              title="Split Violin", value = 3,
              # makes different violin plots for each unique instance of a metadata variable
              # grouped by cell types IDs of different reductions
              br(),
              downloadButton("down5", label = "Download"),
              br(), br(),
              imageOutput("violin2") %>% withSpinner(type = 1, color = cornell_red),
              br()
            ),
            # Dot Plot panel----
            tabPanel(
              title="Dot Plot", value = 4,
              # DotPlot
              br(),
              downloadButton("down6", label = "Download"),
              br(), br(),
              plotOutput("dotplot") %>% withSpinner(type = 1, color = cornell_red)
            ),
            
            id = "tabselected"
          )
        ) #mainPanel
      )
    ), #tabPanel
    
    # 'Myogenic Cells' tab----
    tabPanel(
      title="Myogenic Cells",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          # PHATE Panel Inputs
          # inputs for reduction to be shown
          br(),
          helpText("Visualize gene expression in myogenic cells alone..."),
          selectizeInput(
            "reduction3",
            label = "Choose a dimensional reduction to display:",
            choices = c(
              "PHATE + Harmony"#, 
              # "PHATE + Scanorama"
            ),
            selected = "PHATE + Harmony"
          ),
          br(),
          helpText("Group cells by metadata features..."),
          selectizeInput(
            "bins",
            label = "Metadata Features:",
            choices = c("PHATE bins"),
            selected = "PHATE bins"
          ),
          # violin plot action button
          br(),
          helpText("Click to generate new violin plots"),
          actionButton(
            "action2", label = "Generate",
            style="color: #B31B1B; background-color: #F7F7F7; border-color: #B31B1B"
          ),
          # inputs for violin plots
          br(), br(), br(),
          selectizeInput(
            "gene4",
            label = "Select gene(s) to plot",
            choices = genes_myo,
            selected = "Myod1",
            multiple = TRUE
          ),
          # inputs present on all tabs
          # downloadable plot type and dimensions for "Myogenic Cells" Tab
          br(),
          helpText("Download Specifications"),
          selectizeInput(
            "downloadable2",
            label = "file type:",
            choices = c("pdf", "png", "eps"),
            selected = "pdf"
          ),
          numericInput(
            "plotsizex2",
            label = "horizontal dimension (inches)",
            value = 5
          ),
          numericInput(
            "plotsizey2",
            label = "vertical dimension (inches)",
            value = 5
          )
        ), 

        # Establishes spaces for plots in the main panel
        mainPanel(
          # PHATE grouped by variables
          br(),
          downloadButton("down7", label = "Download"),
          br(), br(),
          plotOutput("PHATE") %>% withSpinner(type = 1, color = cornell_red),
          br(),
          downloadButton("down8", label = "Download"),
          br(), br(),
          imageOutput("phateviolin") %>% withSpinner(type = 1, color = cornell_red),
          br()
        )
      )
    ),
    
    # Spatial Panel----
    tabPanel(
      title="Spatial",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          #Spatial Panel Inputs
          br(),
          helpText("Visualize spatial gene expression across injury response"),
          selectizeInput(
            inputId = "vis_genes_selector",
            label = "Select gene(s) to plot:",
            choices = genes_vis,
            selected = "Myod1"
          ),
          br(),
          helpText("Visualize cell type distributions across injury response"),
          selectizeInput(
            inputId = "vis_types_selector",
            label = "Select cell type(s) to plot:",
            choices = types_vis,
            selected = "Fusing-Myocytes"
          )
        ),
        # Establishes spaces for plots in the main panel
        mainPanel(
          # Gene Expression----
          column(
            width=6,
            downloadButton("down9", label = "Download"),
            br(),
            plotOutput("spatialgene", height=vis.height) %>% withSpinner(type = 1, color = cornell_red),
            br()
          ),
          # Theta Values----
          column(
            width=6,
            downloadButton("down10", label = "Download"),
            br(),
            plotOutput("spatialtheta", height=vis.height) %>% withSpinner(type = 1, color = cornell_red),
            br()
          )
        )
      )
    ),
    # 'Downloads' tab----
    tabPanel(
      title="Downloads",
      # h4("Click on the links below to download Seurat objects, metadata files, etc."),
      br(),
      h4(
        p(
          "Supplemetal resources from our preprint may be found ",
          a(href = "https://github.com/mckellardw/scMuscle/tree/main/supplemental_data", "here.")
        )
      ),
      h2("Seurat object downloads coming soon!")
      
      # TODO: add functionality here:
      # downloadLink(
      #   "down11", 
      #   label = paste0("All Cells (Seurat object; ", round(file.size("data/scMuscle_mm10_slim_v5.RData")/10^6), "Mb)")
      # ),
      # br(),
      # downloadLink(
      #   'down12', 
      #   label = paste0("Myogenic Cells (seurat obj; ", round(file.size("data/myo_slim_seurat_v3.RData")/10^6), "Mb)")
      # )
    )
    
  )
)

# Server logic ----
server <- function(input, output){
  # Loading data----
  load("./data/scMuscle_mm10_slim_v1-1.RData") # All Cells
  load("./data/myo_slim_seurat_v1-1.RData") # Myo Cells
  load("./data/vis_slim_v1.RData") # Visium
 
  # Plot themes and colors----
  # Figure settings
  small.font = 10
  big.font = 12
  line.width = 0.8
  pt.size = 0.01
  vis.pt.size=1.5
  pt.stroke = 0.3
  label.size = 4
  feature_color_direction = reactive(1)
  
  # Color Palette to cell delegation
  colors.celltypes <- reactive({
    as.vector(polychrome())[c(3:9,17,11:13,21,14:16,32,18:20,28,10,30,29,31,33,34,35,36,22:27)]
  })
  
  # Plot themes
  umap.theme <- theme(
    panel.border = element_rect(color="black", size = line.width),
    axis.line = element_blank(), #element_line(color = "black", size = line.width),
    axis.title = element_text(face="bold",size = small.font, hjust = 0.5, vjust = 1),
    axis.text = element_text(size=small.font,color="black"),
    axis.ticks = element_line(color = "black", size = line.width),
    legend.text = element_text(size=small.font,color="black"),
    legend.title = element_text(size=small.font,color="black", face="bold"),
    title=element_text(face="bold.italic",size = small.font, hjust = 0.5)
  )
  dot.theme <- theme(
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    legend.text = element_text(size=small.font, color="black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x=element_text(angle = 90,hjust = 1,vjust= 0.5),
    panel.grid.major = element_line(colour = "gray", size = 0.5)
  )
  vln.theme <- theme(
    panel.background = element_blank(),
    axis.line = element_line(color = "black", size = line.width),
    legend.text = element_text(size=small.font, color="black"),
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y=element_text(size=small.font, face="bold",color="black"),
    axis.text = element_text(color="black",size=small.font),
    axis.ticks.x = element_line(color="black")
  )
  
  # Passes UMAP Panel Inputs----
  
  # passes the reduction type to be plotted
  umap.reduction <- reactive({
    switch(
      input$umap.reduction,
      "Harmony" = "umap_harmony",
      "BBKNN" = "umap_bbknn",
      "Scanorama" = "umap_scanorama"
    )
  })
  umapxlabel <- reactive(paste(umap.reduction(), "1", sep = "_"))
  umapylabel <- reactive(paste(umap.reduction(), "2", sep = "_"))
  
  # creates variable to group by for UMAPs by having the user just decide the reduction
  umap.idents <- reactive(
    if(umap.reduction()== "umap_harmony"){"harmony_factorIDs"}
    else if(umap.reduction()=="umap_bbknn"){"bbknn_factorIDs"}
    else if(umap.reduction()=="umap_scanorama"){"scanorama_factorIDs"}
  )
  
  # passes gene to feature
  gene1 <- reactive({input$gene1})
  
  # Passes Single Violin Panel Inputs----
  
  # passes gene for Violin Plot
  gene2 <- reactive({input$gene2})
  
  # Passes Split Violin Panel Inputs----
  
  # passes gene for Multiple Violin Plots
  gene3 <- reactive({input$gene3})
  
  # allows violin plots to be grouped by chosen cell types
  splitviolincelltype <- reactive(
    if (input$splitviolincelltype == "Harmony"){"harmony_factorIDs"}
    else if(input$splitviolincelltype =="BBKNN"){"bbknn_factorIDs"}
    else if(input$splitviolincelltype == "Scanorama"){"scanorama_factorIDs"}
  )
  
  # Passes Dot Plot Panel Inputs----
  
  # reverses vector to have dotplot plot genes list left to right
  dot <- reactive(rev(input$dot))
  
  # Passes Inputs for Seurat plots----
  variables.umap <- reactive({
    switch(
      input$variables.umap,
      "source" = "source.label",
      "sample" = "sample",
      "chemistry" = "chemistry",
      "injury days" = "injury",
      "injury agent" = "injury.agent",
      "age" = "age",
      "sample prep." = "type",
      "tissue" = "tissue",
      "mouse strain" = "mouse.strain",
      "sex" = "sex",
      "mice per sample" = "mice.per.sample",
      "sequencing instrument" = "Sequencing.Instrument"
    )
  })
  
  variables.singleVln <- reactive({
    switch(
      input$variables.singleVln,
      "source" = "source.label",
      "Harmony cell types" = "harmony_factorIDs",
      "BBKNN cell types" = "bbknn_factorIDs",
      "Scanorama cell types" = "scanorama_factorIDs",
      "sample" = "sample",
      "chemistry" = "chemistry",
      "injury days" = "injury",
      "injury agent" = "injury.agent",
      "age" = "age",
      "sample prep." = "type",
      "tissue" = "tissue",
      "mouse strain" = "mouse.strain",
      "sex" = "sex",
      "mice per sample" = "mice.per.sample",
      "sequencing instrument" = "Sequencing.Instrument"
    )
  })
  
  variables.splitVln <- reactive({
    switch(
      input$variables.splitVln,
      "source" = "source.label",
      "Harmony cell types" = "harmony_factorIDs",
      "BBKNN cell types" = "bbknn_factorIDs",
      "Scanorama cell types" = "scanorama_factorIDs",
      "sample" = "sample",
      "chemistry" = "chemistry",
      "injury days" = "injury",
      "injury agent" = "injury.agent",
      "age" = "age",
      "sample prep." = "type",
      "tissue" = "tissue",
      "mouse strain" = "mouse.strain",
      "sex" = "sex",
      "mice per sample" = "mice.per.sample",
      "sequencing instrument" = "Sequencing.Instrument"
    )
  })
  
  variables.dot <- reactive({
    switch(
      input$variables.dot,
      "source" = "source.label",
      "Harmony cell types" = "harmony_factorIDs",
      "BBKNN cell types" = "bbknn_factorIDs",
      "Scanorama cell types" = "scanorama_factorIDs",
      "sample" = "sample",
      "chemistry" = "chemistry",
      "injury days" = "injury",
      "injury agent" = "injury.agent",
      "age" = "age",
      "sample prep." = "type",
      "tissue" = "tissue",
      "mouse strain" = "mouse.strain",
      "sex" = "sex",
      "mice per sample" = "mice.per.sample",
      "sequencing instrument" = "Sequencing.Instrument"
    )
  })
  
  # Passes PHATE Tab Inputs----
  # passes the reduction type to be plotted
  reduction3 <- reactive({
    switch(
      input$reduction3,
      "PHATE + Harmony" = "phate_harmony"
    )
  })
  # passes variable to group by
  bins <- reactive({switch(input$bins, "PHATE bins" = "phate1.bins")})
  
  # passes gene for Multiple Violin Plots
  gene4 <- reactive({input$gene4})
  
  # generates first plot (All Cells - FeaturePlot/UMAP) ----
  output$feature <- renderPlot({
    if (input$action3 %% 2 == 0) {
      FeaturePlot(
        scMuscle.slim.seurat,
        cells = sample(Cells(scMuscle.slim.seurat)),
        #plot cells in random order
        features = gene1(),
        reduction = umap.reduction()
      ) +
        scale_colour_viridis_c() +
        labs(color = "Log-Normalized\nExpression") +
        umap.theme + theme(
          legend.position = "top",
          plot.title=element_text(face="bold.italic")
        )
    } else {
      FeaturePlot(
        scMuscle.slim.seurat,
        cells = sample(Cells(scMuscle.slim.seurat)),
        #plot cells in random order
        features = gene1(),
        reduction = umap.reduction()
      ) +
        scale_colour_viridis_c(direction = -1) +
        labs(color = "Log-Normalized\nExpression") +
        umap.theme + 
        theme(
          legend.position = "top",
          plot.title=element_text(face="bold.italic")
        )
    }
  })
  
  
  # generates second plot (All Cells - Cell Type UMAP) ----
  output$umap <- renderPlot({
    DimPlot(
      scMuscle.slim.seurat,
      cells = sample(colnames(scMuscle.slim.seurat)), #plot cells in random order
      reduction=umap.reduction(),
      group.by=umap.idents(),
      cols=colors.celltypes(), #adds colors for just the cell types present in this clustering
      na.value = NA, # removes noisy cells from plot
      pt.size = pt.size, # see value above
      label.size = label.size, # see value above
      repel = T,label= TRUE
    ) +
      NoLegend() +
      aes(stroke=pt.stroke)+
      umap.theme+
      theme(plot.title=element_blank())
  })
  
  # generates third plot (All Cells - metadata UMAP) ----
  output$grouping <- renderPlot({
    DimPlot(
      scMuscle.slim.seurat,
      cells = sample(colnames(scMuscle.slim.seurat)), #plot cells in random order
      reduction=umap.reduction(),
      group.by=variables.umap(),
      cols=colors.celltypes(), #adds colors for just the cell types present in this clustering
      na.value = NA, # removes noisy cells from plot
      pt.size = pt.size, # see value above
      label.size = label.size, # see value above
      repel = T,label= TRUE
    ) +
      aes(stroke=pt.stroke)+
      umap.theme+
      theme(plot.title=element_blank())+
      NoLegend()
  })
  
  # generates fourth plot (All Cells - Violin Plot)----
  output$violin1 <- renderPlot({
    VlnPlot( #TODO - add multiple gene plotting
      scMuscle.slim.seurat,
      features = gene2(),
      group.by = variables.singleVln(),
      cols = colors.celltypes(),
      combine = F,
      pt.size = 0
    ) %>% lapply(
      FUN=function(X) X +
        NoLegend() +
        scale_y_continuous(expand=c(0,0)) +
        scale_colour_viridis_c() +
        vln.theme
    ) %>% wrap_plots(ncol=1)
  })
  
  # renders image of fifth plot (All Cells - Split Violin Plot)----
  
  # scales the height of image (300 px of height given per violin plot)
  scaler <- reactive({200*length(unique(scMuscle.slim.seurat@meta.data[[variables.splitVln()]]))})
  
  # generates image of split violin plot
  output$violin2 <- renderImage({
    
    # A temp file to save the output. It will be deleted after renderImage
    # sends it, because deleteFile=TRUE.
    outfile <- tempfile(fileext='.png')
    
    # Generate a png (1100 px width)
    png(outfile, width = 1100, height = scaler())
    print(
      VlnPlot(
        scMuscle.slim.seurat,
        features = gene3(),
        group.by = splitviolincelltype(),
        pt.size = 0
      ) +
        NoLegend() +
        scale_y_continuous(expand=c(0,0.5))+
        facet_grid(rows = vars(scMuscle.slim.seurat@meta.data[[variables.splitVln()]])) +
        scale_colour_viridis_c() +
        vln.theme
    )
    dev.off()
    
    # Return a list
    list(
      src = outfile,
      contentType = 'image/png',
      alt = "This is alternate text"
    )
  }, deleteFile = TRUE)
  
  
  # generates sixth plot (All Cells - DotPlot)----
  output$dotplot <- renderPlot({
    input$action1
    isolate(
      DotPlot(
        scMuscle.slim.seurat,
        features = dot(),
        group.by = variables.dot()
      )+
        scale_colour_viridis_c() +
        dot.theme +
        labs(title = "Metadata Features")
    )
  })
  
  # generates seventh plot (Myogenic Cells - PHATE) ----
  output$PHATE <- renderPlot({
    DimPlot(
      myo.slim.seurat,
      cells = sample(colnames(myo.slim.seurat)), #plot cells in random order
      reduction=reduction3(),
      group.by=bins(),
      cols=rainbow(25)%>%saturation(values=0.75) %>% brightness(values=0.9) %>% as.vector(), #adds colors for just the cell types present in this clustering
      na.value = NA, # removes noisy cells from plot
      pt.size = pt.size, # see value above
      label.size = label.size, # see value above
      repel = T,label= TRUE
    ) +
      NoLegend() +
      aes(stroke=pt.stroke)+
      xlab("PHATE_Harmony_1") +
      ylab("PHATE_Harmony_2") +
      umap.theme+
      theme(plot.title=element_blank())
  })
  # generates eigth plot (Myogenic Cells - PHATE violins)----
  # scales the height of image (300 px of height given per violin plot)
  scaler2 <- reactive({200*length(gene4())})
  
  # generates image of PHATE violin plot
  output$phateviolin <- renderImage({
    input$action2
    isolate({
      # A temp file to save the output. It will be deleted after renderImage
      # sends it, because deleteFile=TRUE.
      outfile <- tempfile(fileext='.png')
      
      # Generate a png (1100 px width)
      png(outfile, width = 1100, height = scaler2())
      print(
        VlnPlot(
          myo.slim.seurat,
          features = gene4(),
          group.by = bins(),
          cols=rainbow(25)%>%saturation(values=0.75) %>% brightness(values=0.9) %>% as.vector(),
          combine=F,
          pt.size = 0
        ) %>% lapply(
          FUN = function(X) X +
            NoLegend() +
            scale_y_continuous(expand=c(0,.5))+
            scale_colour_viridis_c() +
            vln.theme
        ) %>% wrap_plots(ncol=1)
      )
      dev.off()
      
      # Return a list
      list(src = outfile,
           contentType = 'image/png',
           alt = "This is alternate text")
    })}, deleteFile = TRUE)
  
  
  # generates Spatial Plot (gene expression) ####
  
  # passes genes for spatial gene plots
  vis_genes_selector <- reactive({input$vis_genes_selector})
  
  #TODO- add adjustable point sizes & color scale reversal here
  output$spatialgene <- renderPlot({
      visListPlot(
        seu.list=vis.list,
        features=vis_genes_selector(),
        assay='Spatial',
        reduction="space",
        legend.position="right",
        pt.size=vis.pt.size,
        font.size=big.font
      ) 
    }
  )
  
  # generates Spatial Plot (theta values) ####
  
  # passes cell types for spatial cell type plots
  vis_types_selector <- reactive({input$vis_types_selector})
  
  output$spatialtheta <- renderPlot({
    visListPlot(
        seu.list=vis.list,
        features=vis_types_selector(),
        assay='sub_raw_deg_ted',
        reduction="space",
        legend.position="right",
        pt.size=vis.pt.size,
        font.size=big.font
      )&scale_color_gradientn(colors=spatial_theta_colors, na.value=gray(0.42), limits=c(10^-2,1)) #%>% print()
  })
  
  # DownloadHandler----
  # Allows plots to be donwloaded in specificed file type
    # All Cells - Cell Type UMAP----
  output$down1 <- downloadHandler(
    # specify file name
    filename = function() {
      paste("scMuscle_plot", input$downloadable1, sep = ".")
    },
    # creates the plot
    content = function(file){
      print(
        DimPlot(
          scMuscle.slim.seurat,
          cells = sample(colnames(scMuscle.slim.seurat)), #plot cells in random order
          reduction=umap.reduction(),
          group.by=umap.idents(),
          cols=colors.celltypes(), #adds colors for just the cell types present in this clustering
          na.value = NA, # removes noisy cells from plot
          pt.size = pt.size, # see value above
          label.size = label.size, # see value above
          repel = T,label= TRUE
        ) +
          NoLegend() +
          aes(stroke=pt.stroke)+
          xlab(umapxlabel()) +
          ylab(umapylabel()) +
          umap.theme+
          theme(plot.title=element_blank())
      )
      
      ggsave(
        file = file,#paste0(file,".",input$downloadable),
        device = input$downloadable1,
        width = input$plotsizex1, 
        height = input$plotsizey1,
        units = "in", 
        dpi=300
      )
    }
  )
  
    # All Cells - metadata UMAP----
  output$down2 <- downloadHandler(
    # specify file name
    filename = function() {
      paste("name", input$downloadable1, sep = ".")
    },
    # creates the plot
    content = function(file){
      
      print(
        DimPlot(
          scMuscle.slim.seurat,
          cells = sample(colnames(scMuscle.slim.seurat)), #plot cells in random order
          reduction=umap.reduction(),
          group.by=variables.umap(),
          cols=colors.celltypes(), #adds colors for just the cell types present in this clustering
          na.value = NA, # removes noisy cells from plot
          pt.size = pt.size, # see value above
          label.size = label.size, # see value above
          repel = T,label= TRUE
        ) +
          NoLegend() +
          aes(stroke=pt.stroke)+
          xlab(umapxlabel()) +
          ylab(umapylabel()) +
          umap.theme+
          theme(plot.title=element_blank())
      )
      
      ggsave(
        file = file,#paste0(file,".",input$downloadable),
        device = input$downloadable1,
        width = input$plotsizex1, 
        height = input$plotsizey1,
        units = "in", 
        dpi=300
      )
    }
  )
  
    # All Cells - FeaturePlot/UMAP----
  output$down3 <- downloadHandler(
    # specify file name
    filename = function() {
      paste("name", input$downloadable1, sep = ".")
    },
    # creates the plot
    content = function(file){
      if (input$action3 %% 2 == 0) {
        print(FeaturePlot(
          scMuscle.slim.seurat,
          cells = sample(Cells(scMuscle.slim.seurat)),
          #plot cells in random order
          features = gene1(),
          reduction = umap.reduction()
        ) +
          scale_colour_viridis_c() +
          labs(color = "Log-Normalized\nExpression") +
          umap.theme) + theme(legend.position = "top")
      } else {
        print(FeaturePlot(
          scMuscle.slim.seurat,
          cells = sample(Cells(scMuscle.slim.seurat)),
          #plot cells in random order
          features = gene1(),
          reduction = umap.reduction()
        ) +
          scale_colour_viridis_c(direction = -1) +
          labs(color = "Log-Normalized\nExpression") +
          umap.theme) + theme(legend.position = "top")
      }
      
      ggsave(
        file = file,#paste0(file,".",input$downloadable),
        device = input$downloadable1,
        width = input$plotsizex1, 
        height = input$plotsizey1,
        units = "in", 
        dpi=300
      )
    }
  )
  
    # All Cells - Single Violin Plot----
  output$down4 <- downloadHandler(
    # specify file name
    filename = function() {
      paste("name", input$downloadable1, sep = ".")
    },
    # creates the plot
    content = function(file) {
      print(
        VlnPlot( #TODO - add multiple gene plotting
          scMuscle.slim.seurat,
          features = gene2(),
          group.by = variables.singleVln(),
          cols = colors.celltypes(),
          combine = F,
          pt.size = 0
        ) %>% lapply(
          FUN=function(X) X +
            NoLegend() +
            scale_y_continuous(expand=c(0,0)) +
            scale_colour_viridis_c() +
            vln.theme
        ) %>% wrap_plots(ncol=1)
      )
      ggsave(
        file = file,#paste0(file,".",input$downloadable),
        device = input$downloadable1,
        width = input$plotsizex1, 
        height = input$plotsizey1,
        units = "in", 
        dpi=300
      )
    }
  )
  
    # All Cells - Split Violin Plot----
  output$down5 <- downloadHandler(
    # specify file name
    filename = function() {
      paste("name", input$downloadable1, sep = ".")
    },
    # creates the plot
    content = function(file){
      print(
        VlnPlot(
          scMuscle.slim.seurat, 
          features = gene3(), 
          group.by = splitviolincelltype(), 
          pt.size = 0
        ) + 
          NoLegend() +
          scale_y_continuous(expand=c(0,.5))+
          facet_grid(rows = vars(scMuscle.slim.seurat@meta.data[[variables.splitVln()]]))+
          scale_colour_viridis_c()+
          vln.theme
      )
      
      ggsave(
        file = file,#paste0(file,".",input$downloadable),
        device = input$downloadable1,
        width = input$plotsizex1, 
        height = input$plotsizey1,
        units = "in", 
        dpi=300
      )
    }
  )
    # All Cells - DotPlot----
  output$down6 <- downloadHandler(
    # specify file name
    filename = function() {
      paste("name", input$downloadable1, sep = ".")
    },
    # creates the plot
    content = function(file){
      print(
        DotPlot(
          scMuscle.slim.seurat,
          features = dot(),
          group.by = variables.dot()
        )+
          scale_colour_viridis_c()+
          dot.theme
      )
      ggsave(
        file = file,#paste0(file,".",input$downloadable),
        device = input$downloadable1,
        width = input$plotsizex1,
        height = input$plotsizey1,
        units = "in",
        dpi=300
      )
    }
  )

    # Myogenic Cells - PHATE----
  output$down7 <- downloadHandler(
    # specify file name
    filename = function() {
      paste("name", input$downloadable2, sep = ".")
    },
    # creates the plot
    content = function(file){
      print(
        DimPlot(
          myo.slim.seurat,
          cells = sample(colnames(myo.slim.seurat)), #plot cells in random order
          reduction=reduction3(),
          group.by=bins(),
          cols=rainbow(25)%>%saturation(values=0.75) %>% brightness(values=0.9) %>% as.vector(), #adds colors for just the cell types present in this clustering
          na.value = NA, # removes noisy cells from plot
          pt.size = pt.size, # see value above
          label.size = label.size, # see value above
          repel = T,label= TRUE
        ) +
          NoLegend() +
          aes(stroke=pt.stroke)+
          xlab("PHATE_Harmony_1") +
          ylab("PHATE_Harmony_2") +
          umap.theme+
          theme(plot.title=element_blank())
      )
      
      ggsave(
        file = file,#paste0(file,".",input$downloadable),
        device = input$downloadable2,
        width = input$plotsizex2, 
        height = input$plotsizey2,
        units = "in", 
        dpi=300
      )
    }
  )
  
    # Myogenic Cells - PHATE violins----
  output$down8 <- downloadHandler(
    # specify file name
    filename = function() {
      paste("name", input$downloadable2, sep = ".")
    },
    # creates the plot
    content = function(file){
      print(VlnPlot(
        myo.slim.seurat,
        features = gene4(),
        group.by = bins(),
        cols=rainbow(25)%>%saturation(values=0.75) %>% brightness(values=0.9) %>% as.vector(),
        combine=F,
        pt.size = 0
      ) %>% lapply(
        FUN = function(X) X +
          NoLegend() +
          scale_y_continuous(expand=c(0,.5))+
          scale_colour_viridis_c() +
          vln.theme
      ) %>% wrap_plots(ncol=1))
      
      ggsave(
        file = file,#paste0(file,".",input$downloadable),
        device = input$downloadable2,
        width = input$plotsizex2, 
        height = input$plotsizey2,
        units = "in", 
        dpi=300
      )
    }
  )
    # All Cells - .RData----
  output$down11 <- downloadHandler(
    filename = function() {
      paste("name", "Rdata", sep=".")
    },
    content = function(file){
      saveRDS(scMuscle.slim.seurat, file=file)
    }
  )
}

# Run app ----
shinyApp(ui, server)
