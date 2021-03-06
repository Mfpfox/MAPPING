## CpDAA shiny app 
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
addResourcePath("vid", directoryPath = './www')
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

header <- dashboardHeader(
  title = "CpDAA", 
  titleWidth = 200
)

sidebar <- dashboardSidebar( 
  width=150,
  collapsed=TRUE,
  sidebarMenu(id="menu1",
              menuItem("Welcome", icon = icon("user"), tabName="Welcome"),
              menuItem("Analysis", icon = icon("table"), tabName="aalevel1")
              ) # end sidebar menu
  ) # end dashboard

body <- dashboardBody(
  ### changing theme
  shinyDashboardThemes(
    theme = "grey_dark" #blue_gradient grey_light
  ),

  tabItems(
    ## Introduction tab panel
    tabItem(tabName="Welcome",
            fluidRow(
              box(solidHeader = F, collapsible = F, collapsed = F,
                         width = 12,
                  column(12, align="center",
                  br(),
                  h1("Welcome to the CpDAA App!"),
                  h3("For the exploration of ChemoProteomic-Detected Amino Acids (CpDAA)\nwith genetic-based functional annotations"),
                  br()
                  ) # end col
                )# end box
            ), # end row
            
            fluidRow(
              box(solidHeader = F, collapsible = F, collapsed = F,
                  width = 6,
                  column(12, align="center",
                         h2("Manual Curation Workflow for Chemoproteomic-Detected Proteins"), 
                         br(),
                         div(img(img(src ="workflow_shiny.png", height = "700", width = "800")), 
                         style="text-align: center"), 
                         br()
                  ), # endcol centered
                  column(12, 
                         h3("Chemoproteomics Data Sources"),
                         tags$ol(
                           tags$li("Weerapana, E., C. Wang, G. M. Simon, F. Richter, S. Khare, M. B. Dillon, D. A. Bachovchin, K. Mowen, D. Baker, and B. F. Cravatt.",a("2010", href="https://www.nature.com/articles/nature09472"), "Quantitative Reactivity Profiling Predicts Functional Cysteines in Proteomes. Nature 468 (7325): 790–95."),
                           tags$li("Backus KM, Correia BE, Lum KM, Forli S, Horning BD, Gonzalez-Paez GE, Chatterjee S, Lanning BR, Teijaro JR, Olson AJ, et al.", a("2016", href="https://www.nature.com/articles/nature18002"), "Proteome-wide covalent ligand discovery in native biological systems. Nature 534: 570–574"), 
                           tags$li("Hacker, S. M., K. M. Backus, M. R. Lazear, S. Forli, B. E. Correia, and B. F. Cravatt.", a("2017", href="https://www.nature.com/articles/nchem.2826"), "Global Profiling of Lysine Reactivity and Ligandability in the Human Proteome. Nature Chemistry 9 (12): 1181–90."), 
                           tags$li("Cysteine 2019 isoTOP-ABPP data: 10.6019/PXD022151"), style="padding-left: 2em"),
                         br()
                         ) # end col
                  ), # end box
              
              box(solidHeader = F, collapsible = F, collapsed = F,
                  width = 6,
                  column(12,
                         br(),
                         br(),
                         h2("Annotated CpDAA Counts"),
                         DT::dataTableOutput("countTable"),
                         br(),
                         br(),
                         br(),
                         br(),
                         h2("Genetic-Annotation Sources"),
                         DT::dataTableOutput("AnnoTable"),
                         br(),
                         br(),
                         br(),
                         br()
                  ) # end col
                  ) # end box
              ), # end row
            
            # fluidRow(
            #   br()
            # ),
            
            fluidRow(
              box(title = "PDB Structures", solidHeader = T, collapsible = T, collapsed = F, width=12, 
                  column(6, align="center",
                         h2("Caspase-8 structure analysis of CpD Cysteines"),
                         h4("Color set as max missense CADD score per codon"),
                         br(),
                         tags$video(src="vid/CASP8_CADD38max_blues_hires.mp4",
                                    width="750px", height="500px",
                                    type="video/mp4", controls="controls"), 
                         br(),
                         br()
                         
                  ),
                  column(6, align="center",
                         h2("Glucose-6-phosphate dehydrogenase structure analysis of CpD Lysines"),
                         h4("Color set as max missense CADD score per codon"),
                         br(),
                         tags$video(src="vid/G6PD_CADD38max_blue_hiRes.mp4",
                                    width="750px", height="500px",
                                    type="video/mp4", controls="controls"), 
                         br(),
                         br()
                  ) # end col
              )# outer most box
            ), # end row
            
            # fluidRow(
            #   br()
            # ),
            
            fluidRow(
              box(title="Annotated Datasets", solidHeader = T, collapsible = T, collapsed = F, width=12, 
              tabBox(width = 12,
                id = "databaseTabs", height = "250px",
                tabPanel(h2("Proteins"), 
                         DT::dataTableOutput("proTable"),
                         fluidRow(
                           column(width=11),
                           column(width=1,align="center",
                                  br(),
                                  uiOutput("downloadProTable"), 
                                  br()
                           ))
                          ),
                tabPanel(h2("Amino Acids"), 
                         DT::dataTableOutput("aa1Table"),
                         fluidRow(
                           column(width=11),
                           column(width=1,align="center",
                                  br(),
                                  uiOutput("downloadAaTable"), 
                                  br()
                           )) # end fluid row
                         ),
                tabPanel(h2("Pathogenic Mutations"), 
                         DT::dataTableOutput("clinvarTable"),
                         fluidRow(
                           column(width=11),
                           column(width=1,align="center",
                                  br(),
                                  uiOutput("downloadCVTable"), 
                                  br()
                           ))
                         )
              ) # tab box
              ) # end box
            ), # end row
            
            fluidRow(
              box(solidHeader = F, collapsible = F, collapsed = F,
                  width=6, 
                  column(12, 
                         h4("CpDAA App is currently in development (Version 1.0)"),
                         h4("Contact email: ", a("mfpalafox@ucla.edu", href="mfpalafox@ucla.edu")),
                         h4("If you find CpDAA App useful please cite our work, ", a("From Chemoproteomic-Detected Amino Acids to Genomic Coordinates: Insights into Precise Multi-omic Data Integration", href="https://www.biorxiv.org/content/10.1101/2020.07.03.186007v2"))
                  )
              ) # end box
            ),
            br(),
            br(),
            uiOutput("welcomeLink", align="center"),
            tags$style("button#switchTab2 {margin-left:auto;margin-right:auto;display:block;background-color:#EAEAEA; font-family:Andika, Arial, sans-serif; font-size:1em; letter-spacing:0.05em; text-transform:uppercase ;color:#093044; text-shadow: 0px 1px 10px #00FFFF;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}"),
            br()
    ), # end tabItem

############## ANALYSIS #####################################
  tabItem(tabName="aalevel1",
    fluidPage(title = "Detected vs Undetected Amino Acids",
             fluidRow(
               box(
                   solidHeader = F, collapsible = T, collapsed = F,
                   width = 12,
                   column(4,
                          h2("Analysis of Cysteines and Lysines in 4,526 ChemoProteomic-Detected Proteins", style = "color: white"),
                          selectizeInput("selectGene1", "Select Gene or UniProtKB Protein ID", choices=NULL, 
                                         options = list(placeholder = "Search gene",
                                                        maxOptions = 4600)),
                          br(),
                          actionButton(inputId = "aasubmit1", # event reactive
                                       label = "Submit"),
                          tags$style("button#aasubmit1 {margin-left:auto;margin-right:auto;display:block;background-color:#EAEAEA; font-family:Andika, Arial, sans-serif; font-size:1em; letter-spacing:0.05em; text-transform:uppercase ;color:#000; text-shadow: 0px 1px 5px #00FFFF;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}"),
                          br(),
                          br()
                          ), # end col
               column(width=8,
                 tabBox(
                  width = 12,
                   id = "aaTabs", height = "250px",
                   tabPanel(h3("Protein-level"), 
                            DT::dataTableOutput("proteinInfo"),
                            fluidRow(
                              column(width=11),
                              column(width=1, align="center",
                                     br(),
                                     uiOutput("downloadProTable2"), #  delay showing download button
                                     br()
                              ))
                   ),
                   tabPanel(h3("AA-level"), 
                            DT::dataTableOutput("aa2Table"),
                            fluidRow(
                              column(width=11),
                              column(width=1, align="center",
                                     br(),
                                     uiOutput("downloadAaTable2"), #  delay showing download button
                                     br()
                              )) # end fluid row
                   )
                 ) # tabbox
               ) # col 8
             ) # end box
             ), # end row 1
             
      fluidRow(
        br()
      ),
#  line plots
                fluidRow(
                  box(title = "Amino acid summary of missense deleteriousness scores",
                      solidHeader = T, collapsible = T, collapsed = F,
                      width = 12,
                      column(12, 
                             plotlyOutput("linePlots") 
                      )
                      # column(3,
                      #        h2("Scores"),
                      #        p("info")
                      # )
                  ) # end box
                ), # end fluidRow
            br(),
            br(),
            uiOutput("lineplotsLink", align="center"),
            tags$style("button#switchTab1 {margin-left:auto;margin-right:auto;display:block;background-color:#EAEAEA; font-family:Andika, Arial, sans-serif; font-size:1em; letter-spacing:0.05em; text-transform:uppercase ;color:#093044; text-shadow: 0px 1px 10px #00FFFF;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}"),
br()
    )  # end fluid page 
  ) # end tabItem

########## end #############################
  ) # end tab items
) # end dashboard
ui <- dashboardPage(header, sidebar, body) # skin = "blue"
