## CpDAA mapping results shiny app
## mfpfox 9.09.20
#######################################
library(shiny)
library(shinydashboard)
library(dashboardthemes)
library(DT)
library(ggplot2)
library(gplots)
library(grDevices)
library(plotly)
# run in r console
# library(rsconnect)
# rsconnect::deployApp('path/to/your/app')

#######################################
## trim function used for comparision group names' processing
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
#######################################
header <- dashboardHeader(
  title = "Chemoproteomic App"
)
#######################################
sidebar <- dashboardSidebar( 
  sidebarMenu(id="menu1",
              menuItem("Introduction", icon = icon("user"),
                       menuSubItem("Introduction", tabName = "intro"),
                       menuSubItem("Mapping Protocol", tabName = "protocol"),
                       menuSubItem("Data Summary", tabName ="accounting")
                       ),
              menuItem("Protein-level", icon = icon("file"),
                       menuSubItem("Gene Constraint", tabName = "geneConstraint"),
                       menuSubItem("Enrichment Analysis", tabName ="enrichr")
                       ),
              menuItem("AA-level", icon = icon("table"),
                       menuSubItem("Detected vs Undetected", tabName="aalevel1"),
                       menuSubItem("Missense Analysis", tabName="aalevel2")
                       )
              ) # end sidebar menu

  ) # end dashboard

#menuItem("Feedback", icon = icon("comment", lib="glyphicon"), 

# intro
# dashboard body
#######################################
body <- dashboardBody(
  ### changing theme
  shinyDashboardThemes(
    theme = "blue_gradient" #blue_gradient grey_light
  ),
  tabItems(
    ## Introduction tab panel
    tabItem(tabName="intro",
            fluidRow(
              box(title = "Abstract",
                  solidHeader = T, collapsible = T, collapsed = F,
                         width = 12,
                  h2("From Chemoproteomic-Detected Amino Acids to Genomic Coordinates: Insights into Precise Multi-omic Data Integration"),
                         p("Maria F. Palafox1, Valerie A. Arboleda1,2,5,7,8*, Keriann M. Backus3,4,5,6,7,8*", style="padding-left: 1em"),
                         p("*Corresponding Authors: varboleda@mednet.ucla.edu and kbackus@mednet.ucla.edu. Our preprint can be found", a("here", href="https://www.biorxiv.org/content/10.1101/2020.07.03.186007v1"), style="padding-left: 1em"),
                         tags$ol(
                           tags$li("Department of Human Genetics, David Geffen School of Medicine, UCLA, Los Angeles, CA, 90095, USA."),
                           tags$li("Department of Pathology and Laboratory Medicine, David Geffen School of Medicine, UCLA, Los Angeles, CA, 90095, USA."), 
                           tags$li("Department of Biological Chemistry, David Geffen School of Medicine, UCLA, Los Angeles, CA, 90095, USA."),
                           tags$li("Department of Chemistry and Biochemistry, College of Arts and Sciences, UCLA, Los Angeles, CA, 90095, USA."),
                           tags$li("Molecular Biology Institute, UCLA, Los Angeles, CA, 90095, USA."),
                           tags$li("DOE Institute for Genomics and Proteomics, UCLA, Los Angeles, CA, 90095, USA."),
                           tags$li("Jonsson Comprehensive Cancer Center, UCLA, Los Angeles, CA, 90095, USA. "),
                           tags$li("Eli and Edythe Broad Center of Regenerative Medicine and Stem Cell Research, UCLA, Los Angeles, CA, 90095, USA."), style="padding-left: 3em"),
                         h3("Abstract"), 
                  p("The integration of proteomic, transcriptomic, and genetic-variant annotation data will improve our understanding of genotype-phenotype associations. Due, in part, to challenges associated with accurate inter-database mapping, such multi-omic studies have not extended to chemoproteomics, a method that measures the intrinsic reactivity and potential 'druggability' of nucleophilic amino acid side chains. Here, we evaluated two mapping approaches to match chemoproteomic-detected cysteine and lysine residues with their genetic coordinates. Our analysis reveals that databases update cycles and reliance on stable identifiers can lead to pervasive misidentification of labeled residues. Enabled by this examination of mapping strategies, we then integrated our chemoproteomic data with in silico generated predictions of genetic variant pathogenicity, which revealed that codons of highly reactive cysteines are enriched for genetic variants that are predicted to be more deleterious. Our study provides a roadmap for more precise inter-database comparisons and points to untapped opportunities to improve the predictive power of pathogenicity scores and to advance prioritization of putative druggable sites through integration of predictions of pathogenicity with chemoproteomic datasets.", style="padding-left: 2em; text-align: justify")
            )# end box
            ), # end row
            fluidRow(
              box(title = "Graphical abstract",
                  solidHeader = T, collapsible = T, collapsed = T,
                  width = 12,
                  column(12, 
                         div(img(img(src ="dark_abstract.png", height = "600px", width = "800px")), 
                         style="text-align: center"), 
                         br()
                         )
                  )
              ),
            fluidRow(
              box(title = "Mapping terms",
                  solidHeader = T, collapsible = T, collapsed = T,
                  width = 12,
                  column(12, 
                         div(img(img(src ="vocab.png", height = "800px", width = "1000px")), 
                             style="text-align: center"), 
                         br()
                  )
              )
            ),
            fluidRow(
              box(title = "Figure 1",
                  solidHeader = T, collapsible = T, collapsed = T,
                  width = 12,
                  column(8, 
                         div(img(src ="Figure1_V7.png", height = "600px", width = "1000px"), style="text-align: center;")),
                         br(),
                  column(1),
                  column(3, 
                         br(),
                         h3("1. Stable identifier mapping across database releases."),
                         br(),
                         p(tags$b("A) Overall objective of this study is to use stable identifiers to map chemoproteomic detected amino acids (CpDAAs), denoted here as red/blue ‘X,’ to variant pathogenicity scores. B) Shows the number of stable UniProtKB protein IDs from cysteine and lysine chemoproteomics studies in original legacy chemoproteomics dataset (4,119 Uniprot stable IDs in aggregate)(Hacker et al. 2017; Backus et al. 2016; Weerapana et al. 2010) that fail to map to IDs in more recent releases of Ensembl and UniProtKB. C) Average update cycle length across major gene annotation databases. D) Timeline of gene annotation database releases, including Ensembl releases tested for compatibility (Figure 3) to CpDAA coordinates based on canonical UniProtKB protein sequences."), style="text-align: justify")
              )
            )
            ), 
            fluidRow(
              box(title = "Figure 2",
                  solidHeader = T, collapsible = T, collapsed = T,
                  width = 12,
                  column(3, 
                         br(),
                         h3("2. Stable and versioned identifier mapping between UniProtKB releases."),
                         br(),
                         p(tags$b("A) Residue mismapping due to UniProtKB database updates to sequences associated with stable identifiers. Shown are two representative mapping outcomes. For PRMT1, addition of 10 amino acids at the protein N-terminus results in incorrect mapping of all 13 CpDAAs. For FKBP7, CpDAA Lys83 is correctly mapped in the 2018 release, as the sequence update occurs C-terminal to the modified residue. B) The number of isoforms associated with CpDAA UniProtKB stable IDs. C) Canonical isoform ID number associated with UniProtKB entries. Only CpDAA stable IDs with >1 isoform sequence were used for analysis."), style="text-align: justify")),
                 # column(1),
                  column(8, 
                         div(img(img(src ="Figure2_v2-01.png", height = "600px", width = "900px"), style="text-align: center;")), br()
                  )
              )
            ), 
            fluidRow(
              box(title = "Figure 3",
                  solidHeader = T, collapsible = T, collapsed = T,
                  width = 12,
                  column(7, 
                         div(img(img(src ="Figure3_V4-01.png", height = "600px", width = "800ps"), style="text-align: center;")), br()
                  ),
                  column(1),
                  column(4, 
                         br(),
                         h3("3. Mapping between UniProtKB stable IDs and Ensembl stable and versioned IDs across releases."),
                         br(),
                         p(tags$b("A) CpDAA UniProtKB stable IDs map to multiple Ensembl stable protein IDs. The protein glucose-6-phosphate dehydrogenase (G6PD, UniProtKB ID P11413) maps to multiple Ensembl stable protein IDs including identical and non-identical sequences across all five Ensembl releases investigated. B) Number of stable and versioned Ensembl gene, transcript and protein IDs for G6PD across all five Ensembl releases shown in 'A.' C) Cumulative sequence re-annotations for Ensembl gene, transcript, and protein IDs since the v85 release. D-E) Average number of Ensembl gene, transcript, and protein IDs for (D) single isoform (n=1,466) and (E) multi-isoform (n=2,487) CpDAA UniProt entries. F) Heatmap depicts protein alignment scores by normalized Hamming distance, with 0 indicating no difference, comparing CpDAA UniProtKB protein sequences (n=3,953) to cross-referenced IDs of Ensembl protein sequences (n= 29,450) across the five Ensembl releases. For D-E, bar plots represent mean values ± SD for the number of Ensembl IDs per stable UniProt ID. Statistical significance was calculated using an unpaired Student’s T-test, **** p-value <0.0001."), style="text-align: justify"))
              )
            ), 
            fluidRow(
              box(title = "Figure 4",
                  solidHeader = T, collapsible = T, collapsed = T,
                  width = 12,
                  column(3, 
                         h3("4. Mapping CpDAA genomic coordinates to predictions of pathogenicity."),
                         br(),
                         p(tags$b("A-B) Aggregate number of cysteines (A) and lysines (B) in CpDAA-containing proteins (n=3,840), including both detected and undetected residues. C-D) Spearman’s correlation (r) of scores for all possible non-synonymous SNVs at cysteine (C) and Lysine (D) CpDAA codons. E) Odds ratio (OR) comparing the predicted pathogenicity of amino acid substitutions at detected versus undetected residue positions. Cys>Trp (yellow) and Lys>Glu (purple) missense scores were compared using a two-tailed Fisher’s Exact test with a 95% two sided confidence interval at the indicated score thresholds (y axis). OR> 1.0 indicates enrichment and OR< 1.0 indicates depletion of deleterious nonsynonymous SNVs at CpDAA relative to undetected residues. P-value cut-off = 5.0e-05, *** for values < 5.0e-21."),style="text-align: justify")),
                  #column(1),
                  column(8, 
                         div(img(img(src ="Figure4_V2-01.png", height = "800px", width = "900px"), style="text-align: center;")), br()
                  )
              )
            ),
            fluidRow(
              box(title = "Figure 5",
                  solidHeader = T, collapsible = T, collapsed = T,
                  width = 12,
                  column(8, 
                         div(img(img(src ="Figure5_V6-01.png", height = "600px", width = "900px"), style="text-align: center;")), br()
                  ),
                  #column(1),
                  column(3, 
                         br(),
                         h3("5. Association between amino acid reactivity and CADD score."),
                         br(),
                         p(tags$b("A-B)  Distribution of the max CADD38 (model for GRCh38) PHRED score/codon for (A) cysteine (n=1,401) and (B) lysine (n=4,363) CpDAAs of low, medium, and high intrinsic reactivities, defined by isoTOP-ABPP ratios, low (R10:1>5), Medium (2<R10:1<5), High (R<2) (Weerapana et al. 2010; Hacker et al. 2017). Kruskal-Wallis test was used for multiple pairwise-comparisons and Wilcox test was used for pairwise comparisons. Significant FDR adjusted p-values marked, *, p < 0.1, **, p < 0.01, ***, p < 0.0001. C) Shows CADD38 max codon missense scores for residues 1-300 of G6PD (UniProt ID P11413). D) Crystal structure of G6PD (PDB ID: 2BH9) shows K205 and K171 located within the enzyme active site. NADP+ cofactor shown in yellow. Surface colored by CADD38 max codon missense scores. Image generated in PyMOL (R. H. B. Smith, Dar, and Schlessinger, n.d.; DeLano and Others 2002)."), style="text-align: justify"))
              )
            )
      
    ), # end tabItem
# protocol   
#######################################
  tabItem(tabName="protocol",
    ##Workflow section
    h3("protocol section"),
    p(strong("Step 1:"), "Upload.", style="padding-left: 2em; padding-top: 1em")
  ), # end tabitem
# accounting
# accounting
#######################################
  tabItem(tabName="accounting",
          fluidRow(
            box(title = "Overview of mapping route",
                solidHeader = T, collapsible = T, collapsed = F,
                width = 12,
                column(12, 
                       h2("Data accounting during residue-codon mapping to pathogenicity scores."),
                       p(
                         tags$b("A) Shows workflow. B) Number of IDs retained after each mapping step shown in ‘A.’ C) Number of CpD cysteine and lysine residues retained during each mapping step shown in 'A'."),
                         style="padding-left: 2em"),
                       div(img(img(src ="Figure S6_V2.png", height = "500px", width = "700px")), 
                           style="text-align: center"), 
                       br(),
                       h4("Sources of Chemoproteomic Data: "),
                       p("Hacker, S. M., K. M. Backus, M. R. Lazear, S. Forli, B. E. Correia, and B. F. Cravatt. 2017. 'Global Profiling of Lysine Reactivity and Ligandability in the Human Proteome.' Nature Chemistry 9 (12): 1181–90.", style="padding-left: 2em"),
                       p("Weerapana, E., C. Wang, G. M. Simon, F. Richter, S. Khare, M. B. Dillon, D. A. Bachovchin, K. Mowen, D. Baker, and B. F. Cravatt. 2010. 'Quantitative Reactivity Profiling Predicts Functional Cysteines in Proteomes.' Nature 468 (7325): 790–95.", style="padding-left: 2em")
                )
            )
          )
  ), # end tabItem




# PROTEIN LEVEL geneContraint and enrichr for detected subsets
# gene constraint
#######################################
  tabItem(tabName="geneConstraint",
    h3("genecontraint section")
  ), # end tabitem
# enrichr
#######################################
  tabItem(tabName="enrichr",
          fluidPage(title="enrichrPage",
                    tabsetPanel(
                      tabPanel(
                        title=h4("Select Geneset"),
                        fluidRow(
                          box(collapsible = F,
                              width = 12,
                              sidebarLayout(
                                sidebarPanel(
                                  width=4,
                                  h3("Select Geneset for Enrichment Analysis Input"),
                                  fluidRow(
                                    column(12,
                                           p("Gene symbol IDs are sent to the Enrichr server to calculate enrichment scores using 10 pre-defined gene sets. Either All Detected, Cys Detected, or Lys Detected genes can be selected for enrichment analysis."),
                                           selectInput("selectENRsymbol",
                                                       label = h4("Query"),
                                                       choices = list("All Detected", "Cys Detected", "Lys Detected"), selected = "All Detected"),
                                           actionButton(inputId="queryENR",
                                                        label = "Submit"),
                                           tags$style("button#queryENR {margin-left:auto;margin-right:auto;display:block;background-color:#EAEAEA; font-family:Andika, Arial, sans-serif; font-size:0.8em; letter-spacing:0.05em; text-transform:uppercase ;color:#fff; text-shadow: 0px 1px 10px #000;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}")
                                    )
                                  )
                                ), # end sidebarpanel
                                mainPanel(width=8,
                                          p("Selected Geneset Data"),
                                          DT::dataTableOutput("DEToutput")
                                ) # main panel
                              ) # end sidebarlayout
                          ) # end box
                        ) # row
                      ),
                      tabPanel(
                        title=h4("Query EnrichR"),
                        fluidRow(
                          box(collapsible = F,
                              width = 12,
                              sidebarLayout(
                                sidebarPanel(
                                  width=4,
                                  #h3("Select Genes for Enrichment Analysis"),
                                  fluidRow(
                                    column(12,
                                           p("Select which Enrichr results you would like to view"),
                                           selectInput("selectENRresults",
                                                       label = h4("Select EnrichR Results"),
                                                       choices = list("GObp"=1, "GOcc"=2, "GOmf"=3, "KEGG"=4, "ChEA"=5, "HPM"=6), selected = 3),
                                           # actionButton(inputId="queryENRresults",
                                           #              label = "View"),
                                           # tags$style("button#queryENRresults {margin-left:auto;margin-right:auto;display:block;background-color:#EAEAEA; font-family:Andika, Arial, sans-serif; font-size:0.8em; letter-spacing:0.05em; text-transform:uppercase ;color:#fff; text-shadow: 0px 1px 10px #000;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}"),
                                           h5("Download all EnrichR results as excel file"),
                                           uiOutput("downloadDataENR")
                                    )
                                  )
                                ), # end sidebarpanel
                                mainPanel(width=8,
                                          p("Selected EnrichR Results Table"),
                                          DT::dataTableOutput("ENRoutput")
                                ) # main panel
                              ) # end sidebarlayout
                          ) # end box
                        ),
                        fluidRow(
                          box(title = "GO_Biological_Process_2018", solidHeader = T,
                              collapsible = T, collapsed =T,
                              width = 12, background="purple",
                              column(9,
                                     plotOutput("ENRGObp")
                              ),
                              column(3,
                                     h2("GO Biological Process 2018"),
                                     p(""),
                                     p(""),
                                     tags$a(href="http://www.geneontology.org/", "website"),
                                     radioButtons(inputId="sortByGObp",
                                                  label="Sort bars by",
                                                  choices=c('P.value', 'Combined.Score'),
                                                  selected='P.value')
                              ) # end col
                          )
                        ),
                        fluidRow(
                          box(title = "GO_Cellular_Component_2018", solidHeader = T,
                              collapsible = T, collapsed =T,
                              width = 12, background="purple",
                              column(9,
                                     plotOutput("ENRGOcc")
                              ),
                              column(3,
                                     h2("GO Cellular Component 2018"),
                                     p(""),
                                     p(""),
                                     tags$a(href="http://www.geneontology.org/", "website"),
                                     radioButtons(inputId="sortByGOcc",
                                                  label="Sort bars by",
                                                  choices=c('P.value', 'Combined.Score'),
                                                  selected='P.value')
                              ) # end col
                          )
                        ),
                        fluidRow(
                          box(title = "GO_Molecular_Function_2018", solidHeader = T,
                              collapsible = T, collapsed =T,
                              width = 12, background="purple",
                              column(9,
                                     plotOutput("ENRGOmf")
                              ),
                              column(3,
                                     h2("GO Molecular Function 2018"),
                                     p(""),
                                     p(""),
                                     tags$a(href="http://www.geneontology.org/", "website"),
                                     radioButtons(inputId="sortByGOmf",
                                                  label="Sort bars by",
                                                  choices=c('P.value', 'Combined.Score'),
                                                  selected='P.value')
                              ) # end col
                          )
                        ),
                        fluidRow(
                          box(title = "KEGG_2019_Human", solidHeader = T,
                              collapsible = T, collapsed =T,
                              width = 12, background="purple",
                              column(9,
                                     plotOutput("ENRkegg")
                              ),
                              column(3,
                                     h2("KEGG Pathways 2019"),
                                     p(""),
                                     p(""),
                                     tags$a(href="https://www.kegg.jp/", "website"),
                                     radioButtons(inputId="sortByKEGG",
                                                  label="Sort bars by",
                                                  choices=c('P.value', 'Combined.Score'),
                                                  selected='P.value')
                              ) # end col
                          )
                          
                        ),
                        
                        fluidRow(
                          box(title = "ChEA_2016", solidHeader = T,
                              collapsible = T, collapsed =T, width = 12, background="purple",
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
                        ),
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
                        )
                      ) # end tabpanel
                    ) # end tabset
          ) # end fluid page
   ), # end tabItem

# AA LEVEL
#######################################
  tabItem(tabName="aalevel1",
    fluidPage(title = "Detected vs Undetected AA Analysis",
      tabsetPanel(     
        tabPanel(title = h4("All CK"),  
                 fluidRow(
                   box(title = "AA group score summaries",
                       solidHeader = T, collapsible = T, collapsed = F,
                       width = 12,
                   h3("Mean scores for all possible missense variants overlapping Cysteine and Lysine in 3,840 UniProtKB proteins"),
                   br(), 
                   column(12,
                          DT::dataTableOutput("comboCKmean")
                   )
                   ) # end box
                 ),
                 fluidRow(
                   box(title = "Line plot analysis",
                       solidHeader = T, collapsible = T, collapsed = F,
                       width = 12,
                   sidebarLayout(
                     sidebarPanel(
                       h3("Analysis of detected vs undetected cysteine and lysine residues in chemoproteomic-detected UniProt proteins"),
                 fluidRow(
                   column(6,
                      selectizeInput("selectGene1", "Select", choices=NULL, 
                                     options = list(placeholder = "Type gene symbol",
                                                    maxOptions = 4000))
                          ) # end col
                        ), # end row sidebarPanel
                 fluidRow(
                   column(12, align="center",
                          actionButton(inputId = "aasubmit1", # event reactive
                                       label = "Submit"),
                          tags$style("button#aasubmit1 {margin-left:auto;margin-right:auto;display:block;background-color:#EAEAEA; font-family:Andika, Arial, sans-serif; font-size:0.8em; letter-spacing:0.05em; text-transform:uppercase ;color:#fff; text-shadow: 0px 1px 10px #000;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}")
                          ) # end col
                 ) # end submit row
               ), # end sidebarPanel

               mainPanel(
                       DT::dataTableOutput("aa1Table"),
                       fluidRow(
                         column(width=5),
                         column(width=7,
                                br(),
                                uiOutput("downloadaa1Table"), #  delay showing download button
                                br()
                                )
                       ) # end fluid row
               ) # main panel
               ) # end side bar layout
             ) # end box
             ), # end row 1

# AAlevel1 line plots
#######################################
# CADD
                  fluidRow(
                    box(title = "CADD",
                        solidHeader = T, collapsible = T, collapsed = T,
                        width = 12, #background="teal",
                    column(9, 
                           plotlyOutput("caddPlotly1") # plotly
                    ),
                    column(3,
                           h2("CADD Scores"),
                           p("The plot shows."),
                           p("Figure details: "),
                           # select Raw or PHRED scores
                           radioButtons(
                             "caddTypeInput",
                             "Select CADD Score Type: ",
                             choices = c("Raw", "PHRED"),
                             selected = "PHRED",
                             inline = FALSE,
                             width = NULL,
                             choiceNames = NULL,
                             choiceValues = NULL
                           ),
                           radioButtons(
                             "lineTypeCADD",
                             "Select Combo or CK specific line plot: ",
                             choices = c("Combo", "CK-specific"),
                             selected = "Combo",
                             inline = FALSE,
                             width = NULL,
                             choiceNames = NULL,
                             choiceValues = NULL
                           )
                           
                           # select combined or CK specific 
                    )
                    ) # end box
                  ), # end fluidRow
# FATHMM
                  fluidRow(
                  box(title = "FATHMM",
                      solidHeader = T, collapsible = T, collapsed = T,
                      width = 12, #background="teal",
                      column(9, 
                             plotlyOutput("fathmmPlotly1") # plotly
                      ),
                      column(3,
                             h2("FATHMMM-mkl Coding Scores"),
                             p("The plot shows."),
                             p("Figure details: "),
                             radioButtons(
                               "lineTypeFathmm",
                               "Select Combo or CK specific line plot: ",
                               choices = c("Combo", "CK-specific"),
                               selected = "Combo",
                               inline = FALSE,
                               width = NULL,
                               choiceNames = NULL,
                               choiceValues = NULL
                             )
                             # select combined or CK specific 
                      )
                  ) # end box
                  ), # end fluidRow
# DANN
                  fluidRow(
                    box(title = "DANN",
                        solidHeader = T, collapsible = T, collapsed = T,
                        width = 12, #background="teal",
                        column(9, 
                               plotlyOutput("dannPlotly1") # plotly
                        ),
                        column(3,
                               h2("DANN Scores"),
                               p("The plot shows."),
                               p("Figure details: "),
                               radioButtons(
                                 "lineTypeDann",
                                 "Select Combo or CK specific line plot: ",
                                 choices = c("Combo", "CK-specific"),
                                 selected = "Combo",
                                 inline = FALSE,
                                 width = NULL,
                                 choiceNames = NULL,
                                 choiceValues = NULL
                               )
                               # select combined or CK specific 
                        )
                    ) # end box
                  ) # end fluidRow

        ),  # end fluid tab panel

        tabPanel(title = h4("ClinVar Patho"),  
                 fluidRow(
                   box(title = "Mean scores for AA groups overlapping ClinVar \n
                       Pathogenic & Likely Pathogenic variants",
                       solidHeader = T, collapsible = T, collapsed = F,
                       width = 12,
                       h3("Summary of mean scores for clinvar missense variants at \n
                          Detected vs Undetected Cysteine and Lysine residues"),
                       br(), 
                       column(12,
                              DT::dataTableOutput("clinvarMean")
                       )
                   ) # end box
                 ),
                 fluidRow(
                   box(title = "ClinVar Pathogenic Variants & CpDAA sites",
                       solidHeader = T, collapsible = T, collapsed = F,
                       width = 12,
                       column(12,
                              DT::dataTableOutput("clinvarTable")
                       )
                   )
                 )
        ) # end Clinvar tabPanel
#######################################
      ) # end tabset panel w/ OPTIONS 
    )  # end fluid page 
  ), 
# end tabItem

# AAlevel2
#######################################
  tabItem(tabName="aalevel2",
    fluidPage(title= "Detected Cysteine vs Detected Lysine Analysis",
              tabsetPanel(
                tabPanel(title=h4("All CpDAA in 3,840 proteins"),
                         fluidRow(
                           box(collapsible = T,
                               width=12,
                               sidebarLayout(
                                 sidebarPanel(
                                   h3("Analysis of detected cysteines (CpDCys) vs detected lysines (CpDLys)"),
                                   fluidRow(
                                     column(6,
                                             selectizeInput("selectGene2", "Select", choices=NULL, options = list(placeholder="Type gene symbol", maxOptions=4000))
                                             ) # end col
                                   ),
                                   fluidRow(
                                     column(12, align="center",
                                            actionButton(inputId = "aasubmit2",
                                                         label="Submit"),
                                            tags$style("button#aasubmit2 {margin-left:auto;margin-right:auto;display:block;background-color:#EAEAEA; font-family:Andika, Arial, sans-serif; font-size:0.8em; letter-spacing:0.05em; text-transform:uppercase ;color:#fff; text-shadow: 0px 1px 10px #000;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}")
                                            )  # end  col
                                     ) #  end submit  row
                                 ),           # end  row sidebarpanel
                                 mainPanel = (
                                   column(8, 
                                          DT::dataTableOutput("aa2Table")
                                          )
                                   

                                 ) # end main
                               ) # end sidebarLayout
                               ) # end  box
                         ),
                           # CADD
                           fluidRow(
                             box(title = "CADD",
                                 solidHeader = T, collapsible = T, collapsed = T,
                                 width = 12, #background="teal",
                                 column(9,
                                        plotlyOutput("caddPlotly2") # plotly
                                 ),
                                 column(3,
                                        h2("CADD Scores"),
                                        p("The plot shows."),
                                        p("Figure details: "),
                                        # select Raw or PHRED scores
                                        radioButtons(
                                          "caddTypeInput2",
                                          "Select CADD Score Type: ",
                                          choices = c("Raw", "PHRED"),
                                          selected = "PHRED",
                                          inline = FALSE,
                                          width = NULL,
                                          choiceNames = NULL,
                                          choiceValues = NULL
                                        ),
                                        radioButtons(
                                          "assemblyCADD2",
                                          "Select CADD model for y axis: ",
                                          choices = c("GRCh37", "GRCh38"),
                                          selected = "GRCh38",
                                          inline = FALSE,
                                          width = NULL,
                                          choiceNames = NULL,
                                          choiceValues = NULL
                                        )
                                 ) #  end col
                             ) # end box
                           ), # end fluidRow
                           # FATHMM
                           fluidRow(
                             box(title = "FATHMM",
                                 solidHeader = T, collapsible = T, collapsed = T,
                                 width = 12, #background="teal",
                                 column(9,
                                        plotlyOutput("fathmmPlotly2") # plotly
                                 ),
                                 column(3,
                                        h2("FATHMMM-mkl Coding Scores"),
                                        p("The plot shows."),
                                        p("Figure details: ") #,
                                        # radioButtons(
                                        #   "lineTypeFathmm2",
                                        #   "Select Combo or CK specific line plot: ",
                                        #   choices = c("Combo", "CK-specific"),
                                        #   selected = "Combo",
                                        #   inline = FALSE,
                                        #   width = NULL,
                                        #   choiceNames = NULL,
                                        #   choiceValues = NULL
                                        # )
                                        # select combined or CK specific
                                 )
                             ) # end box
                           ), # end fluidRow
                           # DANN
                           fluidRow(
                             box(title = "DANN",
                                 solidHeader = T, collapsible = T, collapsed = T,
                                 width = 12, #background="teal",
                                 column(9,
                                        plotlyOutput("dannPlotly2") # plotly
                                 ),
                                 column(3,
                                        h2("DANN Scores"),
                                        p("The plot shows."),
                                        p("Figure details: ") #,
                                        # radioButtons(
                                        #   "lineTypeDann2",
                                        #   "Select Combo or CK specific line plot: ",
                                        #   choices = c("Combo", "CK-specific"),
                                        #   selected = "Combo",
                                        #   inline = FALSE,
                                        #   width = NULL,
                                        #   choiceNames = NULL,
                                        #   choiceValues = NULL
                                        # )
                                        # select combined or CK specific
                                 )
                             ) # end box
                           ) # end fluidRow
                         ) # tabpanel
                ) # tabset panel
                ) # end fluid page
  ) # end tabItem

########## end #############################
  ) # end tab items
) # end dashboard
ui <- dashboardPage(header, sidebar, body, skin = "blue")
