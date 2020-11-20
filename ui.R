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

addResourcePath("vid", directoryPath = './www')

#######################################
## trim function used for comparision group names' processing
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
#######################################
header <- dashboardHeader(
  title = "CpDAA App", 
  titleWidth = 200
)
#######################################
sidebar <- dashboardSidebar( 
  width=200,
  sidebarMenu(id="menu1",
              menuItem("Welcome", icon = icon("user"), tabName="Welcome"),
              menuItem("Analysis", icon = icon("table"), tabName="aalevel1")
              # menuItem("Feedback", icon = icon("comment", lib="glyphicon")
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
                  h1("Welcome to the CpDAA database!"),
                  h3("For the exploration of ChemoProteomic-Detected Amino Acids\nwith Genetic-based Functional Annotations"),
                  br()
                  ) # end col
                )# end box
            ), # end row
            fluidRow(
              box(solidHeader = F, collapsible = F, collapsed = F,
                  width = 5,
                  column(12, align="center",
                         h3("Manual curation of CpDAA in 4,526 UniProtKB proteins"), 
                         br(),
                         div(img(img(src ="workflow_shiny.png", height = "600", width = "650")), 
                         style="text-align: center"), 
                         br(),
                         br(),
                         DT::dataTableOutput("countTable"),
                         br(),
                         br()
                         ) # end col
                  ), # end box
              box(solidHeader = F, collapsible = F, collapsed = F,
                  width = 7,
                  column(12,
                         br(),
                         br(), 
                         br(),
                         DT::dataTableOutput("AnnoTable"),
                         br(),
                         br(), 
                         br(),
                         h3("CpDAA Sources"),
                         tags$ol(
                           tags$li("Weerapana, E., C. Wang, G. M. Simon, F. Richter, S. Khare, M. B. Dillon, D. A. Bachovchin, K. Mowen, D. Baker, and B. F. Cravatt.",a("2010", href="https://www.nature.com/articles/nature09472"), "Quantitative Reactivity Profiling Predicts Functional Cysteines in Proteomes. Nature 468 (7325): 790–95."),
                           tags$li("Backus KM, Correia BE, Lum KM, Forli S, Horning BD, Gonzalez-Paez GE, Chatterjee S, Lanning BR, Teijaro JR, Olson AJ, et al.", a("2016", href="https://www.nature.com/articles/nature18002"), "Proteome-wide covalent ligand discovery in native biological systems. Nature 534: 570–574"), 
                           tags$li("Hacker, S. M., K. M. Backus, M. R. Lazear, S. Forli, B. E. Correia, and B. F. Cravatt.", a("2017", href="https://www.nature.com/articles/nchem.2826"), "Global Profiling of Lysine Reactivity and Ligandability in the Human Proteome. Nature Chemistry 9 (12): 1181–90."), 
                           tags$li("Cysteine 2019 isoTOP-ABPP data: 10.6019/PXD022151"), style="padding-left: 2em"),
                         br(),
                         br()
                  )
                  )
              ), # end row
            fluidRow(
              br(), 
              br()
            ),
            fluidRow(
              tabBox(
                title = "CpDAA database files", width = 12,
                id = "databaseTabs", height = "250px",
                tabPanel("Protein-level", 
                         DT::dataTableOutput("proTable"),
                         fluidRow(
                           column(width=5),
                           column(width=7,
                                  br(),
                                  uiOutput("downloadProTable"), #  delay showing download button
                                  br()
                           ))
                          ),
                tabPanel("AA-level", 
                         DT::dataTableOutput("aa1Table"),
                         fluidRow(
                           column(width=5),
                           column(width=7,
                                  br(),
                                  uiOutput("downloadAaTable"), #  delay showing download button
                                  br()
                           )) # end fluid row
                         ),
                tabPanel("ClinVar Missense", 
                         DT::dataTableOutput("clinvarTable"),
                         fluidRow(
                           column(width=5),
                           column(width=7,
                                  br(),
                                  uiOutput("downloadCVTable"), #  delay showing download button
                                  br()
                           ))
                         )
              ) # tabbox
            ), # end row
            fluidRow(
              br(), 
              br(), 
              br()
            ),
            fluidRow(
              box( 
                  solidHeader = F, collapsible = F, collapsed = F,
                  width=6, 
                  column(12, align="center",
                         h3("CpD Cysteine of Caspase-8 protein, structure colored by max missense CADD score"),
                         br(),
                         tags$video(src="vid/CASP8_CADD38max_blues_hires.mp4",
                                    width="750px", height="500px",
                                    type="video/mp4", controls="controls")
                         
                  )
              ),
              box(
                solidHeader = F, collapsible = F, collapsed = F,
                width=6, 
                column(12, align="center",
                       h3("CpD Lysine of G6PD protein, structure colored by max missense CADD scores"),
                       br(),
                       tags$video(src="vid/G6PD_CADD38max_blue_hiRes.mp4",
                                  width="750px", height="500px",
                                  type="video/mp4", controls="controls")
                       
                ) # end col
              )# end box
            ), # end row
            fluidRow(
              box(solidHeader = F, collapsible = F, collapsed = F,
                  width=6, 
                  column(12, 
                         h5("CpDAA is currently in development (Version 0.1)"),
                         h5("Contact email: ", a("mfpalafox@ucla.edu", href="mfpalafox@ucla.edu")),
                         h5("If you find CpDAA App helpful please cite our work, ", a("From Chemoproteomic-Detected Amino Acids to Genomic Coordinates: Insights into Precise Multi-omic Data Integration", href="https://www.biorxiv.org/content/10.1101/2020.07.03.186007v2"))
                  )
              ) # end box
            )
    ), # end tabItem

############## ANALYSIS #####################################
  tabItem(tabName="aalevel1",
    fluidPage(title = "Detected vs Undetected AA Analysis",
              #tags$head(tags$style('h4 {color:white;}')),
      #         tags$style(HTML("
      #      .navbar-nav li a:hover, .navbar-nav > .active > a {
      #       color: #fff !important;
      #       background-color:#2b8cc4 !important;
      #       background-image: none !important;
      #   }")),
      # tabsetPanel(     
      #   tabPanel(title = h4("Score mapping"),  type="pills",
             fluidRow(
               box(
                   solidHeader = F, collapsible = T, collapsed = F,
                   width = 12,
                   sidebarLayout(
                     sidebarPanel(
                       h3("Pathogenicity score analyis of Chemo-Proteomic Detected and Undetected  Cysteine and Lysine in 4,526 proteins", style = "color: black"),
                 fluidRow(
                   column(6,
                      selectizeInput("selectGene1", "Select", choices=NULL, 
                                     options = list(placeholder = "Search gene",
                                                    maxOptions = 4600))
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
                 tabBox(
                  width = 12,
                   id = "aaTabs", height = "250px",
                   tabPanel("Protein-level", 
                            DT::dataTableOutput("proteinInfo"),
                            fluidRow(
                              column(width=5),
                              column(width=7,
                                     br(),
                                     uiOutput("downloadProTable2"), #  delay showing download button
                                     br()
                              ))
                   ),
                   tabPanel("AA-level", 
                            DT::dataTableOutput("aa2Table"),
                            fluidRow(
                              column(width=5),
                              column(width=7,
                                     br(),
                                     uiOutput("downloadAaTable2"), #  delay showing download button
                                     br()
                              )) # end fluid row
                   )
                 ) # tabbox
               ) # main panel
               ) # end side bar layout
             ) # end box
             ), # end row 1

# AAlevel1 line plots
#######################################
# CADD
                  fluidRow(
                    box(title = "CADD",
                        solidHeader = T, collapsible = T, collapsed = F,
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
                      solidHeader = T, collapsible = T, collapsed = F,
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
                        solidHeader = T, collapsible = T, collapsed = F,
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
        #)  # end tab panel
      #) # end tabset panel w/ OPTIONS 
    )  # end fluid page 
  ) # end tabItem

########## end #############################
  ) # end tab items
) # end dashboard
ui <- dashboardPage(header, sidebar, body, skin = "blue")
