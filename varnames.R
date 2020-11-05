##########CpDAA app varnames#############################
# menuItem
intro
protocol
accounting #Data Summary
#Protein-level
geneConstraint
enrichr
#AA-level
aalevel1 #Detected vs Undetected
aalevel2 #Missense Analysis
#######################################
#######################################
#aalevel1...
selectGene1
aasubmit1
aa1Table
  downloadaa1Table # delay showing
proaadf()

caddPlotly1 # plot output
  input$caddTypeInput
  lineTypeCADD

fathmmPlotly1
  lineTypeFathmm

dannPlotly1
  lineTypeDann

## DT tables
  dataObs$clinvar
  dataObs$clinvarMean
  dataObs$comboCKmean #output$comboCKmean

# Enrichment variables
selectENRsymbol #"All Detected", "Cys Detected", "Lys Detected"
  queryENR
  DEToutput #DT
  downloadDataENR # button



  ENRoutput #DT
  queryENRresults
  selectENRresults choices = list("GObp", "GOcc", "GOmf", "KEGG", "ChEA", "HPM"), selected = "GObp"),


# data file
enrichdf 
  UniProt.ID, HGNC.symbol, sumCpDC, sumCpDK

detectedDF # subset based on user input




#---plots---
ENRGObp
sortByGObp

ENRGOcc
sortByGOcc

ENRGOmf
sortByGOmf

ENRkegg
sortByKEGG

ENRchea16
sortByChea

ENRhpm
sortByHpm

InterPro_Domains_2019

Kinase_Perturbations_from_GEO_up

OMIM_Disease

PPI_Hub_Proteins

ProteomicsDB_2020

Rare_Diseases_AutoRIF_ARCHS4_Predictions

Rare_Diseases_AutoRIF_Gene_Lists

Reactome_2016

RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO

SILAC_Phosphoproteomics


  # "WikiPathways_2019_Human", # 3
  # "MSigDB_Oncogenic_Signatures", # 4
  # "NCI-60_Cancer_Cell_Lines", # 5
  # "DepMap_WG_CRISPR_Screens_Sanger_CellLines_2019", # 6
  # "Achilles_fitness_increase", # 7
  # "Achilles_fitness_decrease", # 8
  # "GTEx_Tissue_Sample_Gene_Expression_Profiles_up", #9
  # "GTEx_Tissue_Sample_Gene_Expression_Profiles_down" # 10
  


#####################################
-------aalevel 2 missense------------

selectGene2
aasubmit2
aa2Table
  downloadaa2Table
caddPlotly2
caddTypeInput2 # Raw PHRED
assemblyCADD2 # GRCh37, GRCh38

fathmmPlotly2

dannPlotly2

clinvarTable











---------server------------------
aadf
missdf

geneDF <- read_csv(paste(getwd(),"data/AAlevel2_allCK_labeled5_shiny_189627.csv",sep="/"))


  

df <- read_csv(paste(getwd(),"data/AAlevel2_allCK_labeled5_shiny_189627.csv",sep="/")) %>%
  mutate(AAgroup = factor(AAgroup, levels = c(
    "Detected C",
    "Undetected C",
    "Detected K",
    "Undetected K"
  ))) %>%
  mutate(aaref = factor(aaref, levels = c(
    "C",
    "K"
  ))) %>%
  mutate(group = factor(group, levels = c(
    "Detected",
    "Undetected"
  ))) %>%
  mutate(overlap.sum.CVIDs.COSIDs = factor(overlap.sum.CVIDs.COSIDs)) %>%
  mutate(detected.group = factor(detected.group, levels = c(
    "High","High Target","High panReactive",
    "Medium","Medium Target","Medium panReactive",
    "Low","Low Target","Low panReactive",
    "Target","panReactive",
    "Undetected"
  ))) %>%
  mutate(matched.aapos = as.integer(matched.aapos))



  # select gene name
  updateSelectizeInput(session, "selectGene1",
                       choices = as.vector(unique(df$HGNCname)), server = TRUE)

selectizeInput(inputId = "selectGene1", label = "Select", choices = NULL,
  options = list(placeholder = "Type gene symbol",
    maxOptions = 4000)
  )















#######################################
#############DEApp functional section##########################
#######################################

  tabItem(tabName= "introTab",
          fluidRow(
            box(title="Functional Analysis of DE Results", 
                solidHeader = T, collapsible = F, 
                width = 12,
                sidebarLayout(
                  sidebarPanel(width=7,
                               h3("A. Use Methods Comparison Results"),
                               p("This option uses the DE gene list from 'Methods Comparison' section. The venn-diagram section overlapping all compared packages represents the DE gene list used by this option. Human or Mouse Ensembl IDs along with Mouse Entrez IDs will be converted to Human Entrez IDs first for best GSEA results. Successfully mapped IDs along with their associated and sorted log2FC values will serve as input for GSEA and Enrichment analysis. A preview of the input data and ID mapping results will be displayed after clicking 'Submit'."),
                               fluidRow(
                                 column(2,
                                        radioButtons(inputId="gseageneIDA", 
                                                     label="Gene ID type",
                                                     choices=c(ENSEMBL ='ENSEMBL',
                                                               ENTREZID ='ENTREZID'
                                                     ),selected='ENSEMBL')
                                 ),
                                 column(2, 
                                        radioButtons(inputId="organismTypeA", 
                                                     label="Organism",
                                                     choices=c(Human ='org.Hs.eg.db',
                                                               Mouse ='org.Mm.eg.db'),
                                                     selected='org.Hs.eg.db')
                                 ),
                                 column(3,
                                        radioButtons(inputId="gseaLogFCA", 
                                                     label="Source of log2FC column",
                                                     choices=c(edgeR ='edgeR',
                                                               voom ='voom',
                                                               DESeq2 = 'DESeq2'
                                                     ),selected='DESeq2')
                                 ),
                                 column(2),
                            
                                 column(3,
                                        br(),
                                        br(),
                                        actionButton(inputId = "goA", 
                                                     label = "Submit"),
                                        tags$style("button#goA {margin-left:auto;margin-right:auto;display:block;background-color:#EAEAEA; font-family:Andika, Arial, sans-serif; font-size:0.8em; letter-spacing:0.05em; text-transform:uppercase ;color:#fff; text-shadow: 0px 1px 10px #000;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}")
                                 )
                               ),
                               br(),
                               hr(),
                               br(),
                               h3("B. Upload File"),
                               p("This option accepts uploaded DE results for functional analysis. Similar to option A, Human or Mouse Ensembl IDs along with Mouse Entrez IDs will be converted to Human Entrez IDs first for best GSEA results. Successfully mapped IDs along with their associated and sorted log2FC values will serve as input for GSEA and Enrichment analysis. A preview of the input data and ID mapping results will be displayed after clicking 'Submit'."),
                               fluidRow(
                                 fileInput(inputId="gseaFile", 
                                           label=NULL, 
                                           accept=c('text/tab-separated-values',
                                                    'text/csv',
                                                    'text/comma-separated-values',
                                                    'text/tab-separated-values',
                                                    '.txt',
                                                    '.csv',
                                                    '.tsv')
                                 )
                               ), 
                               fluidRow(
                                 column(2,
                                        radioButtons(inputId="gseaFileSep", 
                                                     label="Separator",
                                                     choices=c(Comma=',',
                                                               Tab='\t'
                                                     ),selected=',')
                                 ),
                                 column(2,
                                        radioButtons(inputId="gseageneIDB", 
                                                     label="Gene ID type",
                                                     choices=c(ENSEMBL ='ENSEMBL',ENTREZID ='ENTREZID'
                                                     ),selected='ENSEMBL')
                                 ),
                                 column(2, 
                                        radioButtons(inputId="organismTypeB", 
                                                     label="Organism",
                                                     choices=c(Human ='org.Hs.eg.db',Mouse ='org.Mm.eg.db'),
                                                     selected='org.Hs.eg.db')
                                 ),
                                 column(3, 
                                        numericInput("geneidcol", 
                                                     label = "Gene ID column number", 
                                                     value = 1),
                                        numericInput("gseaLogFCB", 
                                                     label=" Log2FC column number",
                                                     value=2)
                                 ),
                                 column(3,
                                        br(),
                                        br(),
                                        actionButton(inputId = "goB", 
                                                     label = "Submit"),
                                        tags$style("button#goB {margin-left:auto;margin-right:auto;display:block;background-color:#EAEAEA; font-family:Andika, Arial, sans-serif; font-size:0.8em; letter-spacing:0.05em; text-transform:uppercase ;color:#fff; text-shadow: 0px 1px 10px #000;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}")
                                        
                                 )
                               ) # end sidebar row
                  ), # end side panel bar
                  
                  mainPanel(width=5,
                            DT::dataTableOutput("GOupload"),
                            br(),
                            verbatimTextOutput("mapText"),
                            br(),
                            uiOutput("downloadMAP",align="center")
                            
                  )
                ) # end side bar layout
            ) # end box
          )#, # end upload row
          #uiOutput("step4toPlots", align="center")
  ), # end tabItem intro


  tabItem(tabName="gsea",
          fluidPage(title = "GSEA options",
                    tabsetPanel(     
                      tabPanel(title = h4("GO"),  
                               
# ----------------- GO parameter setting ------------------
                               fluidRow(
                                 box(collapsible = F,
                                     width = 12,
                                 sidebarLayout(
                                   sidebarPanel(
                                     h3("GO GSEA Parameters"),
                                     fluidRow(
                                       column(6,
                                              # GO ontology
                                              selectInput("selectGO1", 
                                                          label = h4("Gene Ontology"),
                                                          choices = list("Biological Process" = "BP", 
                                                                         "Molecular Function" = "MF", 
                                                                         "Cellular Component" = "CC", 
                                                                         "ALL"),
                                                          selected = "MF"), 
                                              # PADJ
                                              selectInput("selectPADJ1", 
                                                          label = h4("P adj. method"), 
                                                          choices = list("holm", "hochberg", "hommel", 
                                                                         "bonferroni", "BH", "BY", "fdr", "none"), 
                                                          selected = "BH"), # BH
                                              # PCUTOFF
                                              numericInput("selectPCUT1", 
                                                           label = h4("P value cutoff"), 
                                                           value = 0.05)
                                              ),
                                       column(6,
                                              # min
                                              numericInput("minG1", 
                                                           label = h4("Min. gene set size"), value = 10),
                                              
                                              # max
                                              numericInput("maxG1", 
                                                           label = h4("Max gene set size"), value = 500),
                                              # nPerm
                                              numericInput("nperm1", 
                                                           label = h4("Permutations"), value = 1000)
                                              ) # end col
                                     ), # end row sidebarPanel
                                     
                                     # submit center button
                                     fluidRow(
                                       column(12, align="center",
                                              actionButton(inputId = "go2", # event reactive
                                                           label = "Submit"),
                                              tags$style("button#go2 {margin-left:auto;margin-right:auto;display:block;background-color:#EAEAEA; font-family:Andika, Arial, sans-serif; font-size:0.8em; letter-spacing:0.05em; text-transform:uppercase ;color:#fff; text-shadow: 0px 1px 10px #000;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}")
                                              )
                                     ) # end submit row
                                     
                                   ), # end sidebarPanel

                                   mainPanel(
                                           DT::dataTableOutput("gseaTable"),
                                           fluidRow(
                                             column(width=5),
                                             column(width=7,
                                                    br(),
                                                    uiOutput("downloadDataGO"), #  delay showing download button
                                                    br()
                                                    )
                                           ) # end fluid row
                                   ) # main panel
                                   ) # end side bar layout
                                 ) # end box
                                 ), # end upload row 1
                               
# ---------------------- GO PLOTS --------------------

                      fluidRow(
                        # ridge plot
                        box(title = "Ridge",
                            solidHeader = T, collapsible = T, collapsed = T,
                            width = 12, background="teal",
                        
                        column(9, 
                               plotOutput("GOridge")
                        ),
                        column(3,
                               h2("Ridge plot"),
                               p("The ridgeplot shows frequency of fold change values per core genes grouped within each set category (y axis). This plot is helpful for interpreting up-/down-regulated pathways."),
                               p("Figure details: displays only core enriched gene expression values with color based on p.adjust values."),
                               sliderInput("GOridgeslider", 
                                           label=h5("Categories displayed"),
                                           min=5, max=20, value=10),
                               actionButton(inputId = "renderGOridge", 
                                            label = "Render Ridge")
                        )
                        ) # end box
                      ) # end fluidRow

                    ), # end of GO fluid tab


# ---------------------- KEGG fluid tab  ---------------------------
                    tabPanel(title = h4("KEGG"),
                             fluidRow(
                                box(collapsible = F,
                                    width = 12,
                                    sidebarLayout(
                                      sidebarPanel(
                                        h3("KEGG GSEA Parameters"),
                                        fluidRow(
                                          column(6,
                                                 # PADJ
                                                 selectInput("selectPADJK", 
                                                             label = h4("P adj. method"), 
                                                             choices = list("holm", "hochberg", "hommel", 
                                                                            "bonferroni", "BH", "BY", "fdr", "none"), 
                                                             selected = "BH"),
                                                 # PCUTOFF
                                                 numericInput("selectPCUTK", 
                                                              label = h4("P value cutoff"), 
                                                              value = 0.05)
                                              ),
                                            column(6,
                                                   # min
                                                   numericInput("minGK", 
                                                                label = h4("Min. gene set size"), value = 10),
                                                   
                                                   # max
                                                   numericInput("maxGK", 
                                                                label = h4("Max gene set size"), value = 500),
                                                   # nPerm
                                                   numericInput("npermK", 
                                                                label = h4("Permutations"), value = 1000)
                                            ) # end col
                                          ), # end row sidebarPanel
                                          
                                            # submit center button
                                        fluidRow(
                                          column(12, align="center",
                                                 actionButton(inputId = "kegg2", # event reactive
                                                              label = "Submit"),
                                                 tags$style("button#kegg2 {margin-left:auto;margin-right:auto;display:block;background-color:#EAEAEA; font-family:Andika, Arial, sans-serif; font-size:0.8em; letter-spacing:0.05em; text-transform:uppercase ;color:#fff; text-shadow: 0px 1px 10px #000;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}")
                                          )
                                        ) # end submit row
                                            
                                      ), # end sidebarPanel
                                        
                                    mainPanel(
                                          DT::dataTableOutput("keggTable"),
                                          fluidRow(
                                            column(width=5),
                                            column(width=7,
                                                   br(),
                                                   uiOutput("downloadDataKEGG"),
                                                   br()
                                                   )
                                          ) # end fluid row
                                    ) # main panel
                                    
                                  ) # end side bar layout
                              ) # end box
                            ), # end  row 
                                  
# ---------------------- KEGG PLOTS  ------------------
                                  fluidRow(
                                    # pathview plot
                                    box(title = "Pathview",
                                        solidHeader = T, collapsible = T, collapsed = T,
                                        width = 12, background="aqua",
                                        
                                        column(9,
                                               plotOutput()
                                        ),
                                        column(3, 
                                               h2("Pathview plot"), 
                                               p("This plot displays gene expression data mapped to the native KEGG view plot 
                                                 for the KEGG ID entered below."),
                                               tags$b("To save the Pathview plot, please right click the rendered image and select 'Save image as...'"),
                                               textInput("KEGGpathID", 
                                                         label=h3("KEGG ID"), 
                                                         value="hsa..."),
                                               actionButton(inputId = 
                                                              "renderKEGGpath",
                                                            label = "Render Pathview")
                                        )
                                    ) # end box
                                  ) # end row

                ) # end fluid tab
          ) # end tabset panel w/ OPTIONS 
  )#, # end fluid page 
   # uiOutput("step4GSEAtoENR", align="center") # GSEA page
  ), # end tabItem-gsea


  tabItem(tabName= "enrichr",
          fluidPage(title="enrichPage",
                    tabsetPanel(
                      tabPanel(title = h4("Enrichr"), 
                        fluidRow(box(collapsible = F, width = 12,
                          sidebarLayout(
                             sidebarPanel(width=4,
                                          h3("Select Genes for Enrichment Analysis"),
                                          fluidRow(
                                            column(12,
                                                   # Query all genes from DGE results
                                                   p("Gene IDs from DE results are converted into gene symbol IDs and then sent to the Enrichr server to calculate enrichment scores using 10 pre-defined gene sets. Either All, only Up-regulated, or only Down-regulated genes can be selected for enrichment analysis."),
                                                   selectInput("selectENRsymbol", 
                                                               label = h4("Query differentially expressed genes"), 
                                                               choices = list("All", "Up", "Down"), 
                                                               selected = "All"),
                                                   actionButton(inputId="queryENR",
                                                                label = "Submit"),
                                                   tags$style("button#queryENR {margin-left:auto;margin-right:auto;display:block;background-color:#EAEAEA; font-family:Andika, Arial, sans-serif; font-size:0.8em; letter-spacing:0.05em; text-transform:uppercase ;color:#fff; text-shadow: 0px 1px 10px #000;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}")
                                            ))
                             ), # end sidebarpanel
                             
                             mainPanel(width=8,
                                       #box(width=12, ,
                                       p("Preview Enrichment Results in Human Proteome Map"),
                                           DT::dataTableOutput("ENRoutput"),
                                           fluidRow(
                                             column(width=5),
                                             column(width=7,
                                                    br(),
                                                    uiOutput("downloadDataENR"),
                                                    br()
                                             )
                                           )
                                       #) # end box
                             ) # main panel
                             
                           ) # sidebar layout
                       ) # end box
                     ), # end row query
              
                     # ---------------------- ENRICHR PLOTS  ------------------
                     fluidRow(
                       # gseaplot
                       box(title = "ChEA_2016",
                           solidHeader = T, collapsible = T, collapsed =T,
                           width = 12, background="purple",
                           column(9,
                                  plotOutput("ENRchea16")
                           ),
                           column(3,
                                  h2("ChEA 2016"),
                                  p("Target genes of transcription factors from published ChIP-chip, ChIP-seq, and other transcription factor binding site profiling studies."),
                                  p("Stats: 21,220 genes; 199 transcription factors; 386, 625 gene-transcription factor associations"),
                                  tags$a(href="https://amp.pharm.mssm.edu/Harmonizome/dataset/CHEA+Transcription+Factor+Targets", "website"),
                                  radioButtons(inputId="sortByChea",
                                               label="Sort bars by",
                                               choices=c('P.value', 'Combined.Score'),
                                               selected='P.value')
                           ) # end col
                       )
                     ), # end row chea
                     fluidRow(
                       box(title = "Tissue_Protein_Expression_from_Human_Proteome_Map",
                           solidHeader = T, collapsible = T, collapsed = T,
                           width = 12, background="purple",
                           column(9,
                                  plotOutput("ENRhpm")
                           ),
                           column(3,
                                  h2("Human Proteome Map"),
                                  p("HPM contains direct evidence of translation of a number of protein products derived from over 17,000 human genes covering >84% of the annotated protein-coding genes in humans. The database is based on >290,000 non-redundant peptide identifications from proteomic studies of multiple organs/tissues and cell types from individuals with clinically defined healthy tissues."),
                                  p("Stats: includes 17 adult tissues, 6 primary hematopoietic cells and 7 fetal tissues"),
                                  tags$a(href="https://www.humanproteomemap.org/", "website"),
                                  radioButtons(inputId="sortByHpm",
                                               label="Sort bars by",
                                               choices=c('P.value', 'Combined.Score'),
                                               selected='P.value')
                           ) # end col
                       )
                     ), # end row hpm
                     fluidRow(
                       box(title = "WikiPathways_2019_Human",
                           solidHeader = T, collapsible = T, collapsed = T,
                           width = 12, background="purple",
                           column(9,
                                  plotOutput("ENRwiki")
                           ),
                           column(3,
                                  h2("Human Wikipathways 2019"),
                                  p("WikiPathways is a community resource for contributing and maintaining content dedicated to biological pathways. Any registered WikiPathways user can contribute and while contributions are monitored by a group of admins, the bulk of peer review, editorial curation, and maintenance is done by the user community."),
                                  tags$a(href="https://www.wikipathways.org/index.php/WikiPathways", "website"),
                                  radioButtons(inputId="sortByWiki",
                                               label="Sort bars by",
                                               choices=c('P.value', 'Combined.Score'),
                                               selected='P.value')
                           ) # end col
                       )
                     ), # end row wiki
                     fluidRow(
                       box(title = "MSigDB_Oncogenic_Signatures",
                           solidHeader = T, collapsible = T, collapsed = T,
                           width = 12, background="purple",
                           column(9,
                                  plotOutput("ENRmsigdb")
                           ),
                           column(3,
                                  h2("MSigDB Oncogenic Signatures"),
                                  p("The Molecular Signatures Database (MSigDB) Oncogenic Signatures ontology contains gene sets relevant to cancer."),
                                  tags$a(href="https://www.gsea-msigdb.org/gsea/msigdb/index.jsp", "website"),
                                  radioButtons(inputId="sortByMsigdb",
                                               label="Sort bars by",
                                               choices=c('P.value', 'Combined.Score'),
                                               selected='P.value')
                           ) # end col
                       )
                     ), # end row msigdb
                     
                     fluidRow(
                       box(title = "NCI-60_Cancer_Cell_Lines",
                           solidHeader = T, collapsible = T, collapsed = T,
                           width = 12, background="purple",
                           column(9,
                                  plotOutput("ENRnci")
                           ),
                           column(3,
                                  h2("NCI 60 Cancer Cell Lines"),
                                  p("National Cancer Institute (NCI) 60 human tumor cell line anticancer drug discovery project was intended to screen thousands of pharmacological compounds for anticancer activity, thus replacing earlier cumbersome en vivo models."),
                                  p("Stats: 60 created cell lines represent nine human cancers- breast, central nervous system, colon, kidney, leukemia, lung, melanoma, ovary, and prostate"),
                                  tags$a(href="http://biogps.org/downloads/", "website"),
                                  radioButtons(inputId="sortByNci",
                                               label="Sort bars by",
                                               choices=c('P.value', 'Combined.Score'),
                                               selected='P.value')
                           ) # end col
                       )
                     ), # end row nci60
                     fluidRow(
                       box(title = "DepMap_WG_CRISPR_Screens_Sanger_CellLines_2019",
                           solidHeader = T, collapsible = T, collapsed = T,
                           width = 12, background="purple",
                           column(9,
                                  plotOutput("ENRdepmap")
                           ),
                           column(3,
                                  h2("Sanger Dependency Map - CRISPR screens"),
                                  p("The consequences of alterations in the DNA of cancer cells and subsequent vulnerabilities are not fully understood. This project aims to assign a dependency to every cancer cell in a patient which could be exploited to develop new therapies."),
                                  p("Using large-scale drug and genetic perturbation platforms, the project aims to exhaustively define dependencies which are operative in cancer cells. The gene set comes from large-scale CRISPR-Cas9 synthetic-lethal genetic screens to map gene function and identify dependencies across a diverse collection of cancer cell models."),
                                  tags$a(href="https://depmap.sanger.ac.uk/", "website"),
                                  radioButtons(inputId="sortByDemap",
                                               label="Sort bars by",
                                               choices=c('P.value', 'Combined.Score'),
                                               selected='P.value')
                           ) # end col
                       )
                     ), # end row depmap
                     
                     # ach in ui  
                     fluidRow(
                       box(title = "Achilles_fitness_increase",      
                           solidHeader = T, collapsible = T, collapsed = T,
                           width = 12, background="purple",
                           column(9,
                                  plotOutput("ENRfitIncrease")
                           ),
                           column(3,
                                  h2("Achilles Fitness Increase"),     
                                  p("Achilles systematically identifies and catalogs gene essentiality across hundreds of genomically characterized cancer cell lines. They use genome-scale RNAi and, more recently, CRISPR-Cas9 genetic perturbation reagents to silence or knockout individual genes and identify those genes that affect cell survival. This gene set is specific to increased fitness results from gene-knockout screens."),
                                  tags$a(href="http://www.broadinstitute.org/achilles", "website"),
                                  radioButtons(inputId="sortByFitIn",          
                                               label="Sort bars by",
                                               choices=c('P.value', 'Combined.Score'),
                                               selected='P.value')
                           ) # end col
                       )
                     ), # end row achilles increase
                     
                     # ach do ui  
                     fluidRow(
                       box(title = "Achilles_fitness_decrease",     
                           solidHeader = T, collapsible = T, collapsed = T,
                           width = 12, background="purple",
                           column(9,
                                  plotOutput("ENRfitDecrease")      
                           ),
                           column(3,
                                  h2("Achilles Fitness Decrease"),   
                                  p("Using Project Achilles, the ‘oncogene addictions’ of many known oncogenes such as PIK3CA, KRAS, and BRAF are identified by showing that cell lines harboring mutations of these genes exhibit higher sensitivity to their suppression. This gene set is specific to decreased fitness results from gene-knockout screens."),
                                  tags$a(href="http://www.broadinstitute.org/achilles", "website"),
                                  radioButtons(inputId="sortByFitDe",  
                                               label="Sort bars by",
                                               choices=c('P.value', 'Combined.Score'),
                                               selected='P.value')
                           ) # end col
                       )
                     ), # end row ach decrease
                     
                     
                     
                     # gtex up ui
                     fluidRow(
                       box(title = "GTEx_Tissue_Sample_Gene_Expression_Profiles_up",     
                           solidHeader = T, collapsible = T, collapsed = T,
                           width = 12, background="purple",
                           column(9,
                                  plotOutput("ENRgtexUp")     
                           ),
                           column(3,
                                  h2("GTEx Tissue Gene Expression Profiles - Up"),   
                                  p("The Genotype-Tissue Expression (GTEx) project is a resource to study human gene expression/regulation and its relationship to genetic variation. Samples were collected from 54 non-diseased tissue sites across nearly 1000 individuals, primarily for molecular assays including WGS, WES, and RNA-Seq. This gene set is specific to highly expressed genes in tissue samples."),
                                  tags$a(href="http://www.gtexportal.org/", "website"),
                                  radioButtons(inputId="sortByGtexUp",         
                                               label="Sort bars by",
                                               choices=c('P.value', 'Combined.Score'),
                                               selected='P.value')
                           ) # end col
                       )
                     ), # end row gtex up
                     
                     # gtex down ui  
                     fluidRow(
                       box(title = "GTEx_Tissue_Sample_Gene_Expression_Profiles_down",    
                           solidHeader = T, collapsible = T, collapsed = T,
                           width = 12, background="purple",
                           column(9,
                                  plotOutput("ENRgtexDown")    
                           ),
                           column(3,
                                  h2("GTEx Tissue Gene Expression Profiles - Down"),   
                                  p("This GTEx gene set is specific to lowly expressed genes in tissue samples."),
                                  tags$a(href="http://www.gtexportal.org/", "website"),
                                  radioButtons(inputId="sortByGtexDo",          
                                               label="Sort bars by",
                                               choices=c('P.value', 'Combined.Score'),
                                               selected='P.value')
                           ) # end col
                       )
                     ), # end row gtex down
# --------- Download -------------------------
                     fluidRow(
                       box(title = "Download",
                           solidHeader = T, collapsible = T, collapsed = T,
                           width = 12, background="purple",
                           column(9),
                           column(3,
                                  h2("Saving plots"), 
                                  p("To successfully download, the figure must first be rendered. All plots are saved as PDF."),
                                  selectInput("pickplotENR", 
                                              label = h4("Select Enrichr plot(s)"), 
                                              choices = list("All","Chea16", "HPM", 
                                                             "Wiki", "MSigDB","NCI60", "DeMap", "FitIncrease", "FitDecrease", "GtexUp", "GtexDown"), selected = "All"),
                                  downloadButton("downloadENRpval", 
                                                 label = "Download",
                                                 class = NULL),
                                  tags$style("#downloadENRpval {left:right; }"),
                                  br()
                           ) # end col
                       ) # end box
                     ) # end row
                    
            ) # end Enrichr fluid tab
          ) # end tabset panel
          )#, # end fluid page
          #uiOutput("step4ENRtoGSEA", align="center") # enrichr page
) # end tabItem enrichr