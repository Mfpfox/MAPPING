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
  title = "CpDAA DB"
)
#######################################
sidebar <- dashboardSidebar( 
  sidebarMenu(id="menu1",
              menuItem("Welcome", icon = icon("user"), tabName="Welcome"),
              menuItem("Analysis", icon = icon("table"), tabName="aalevel1")
              # menuItem("Feedback", icon = icon("comment", lib="glyphicon")
              ) # end sidebar menu
  ) # end dashboard


body <- dashboardBody(
  ### changing theme
  shinyDashboardThemes(
    theme = "grey_light" #blue_gradient grey_light
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
                  width = 6,
                  column(12, 
                         br(),
                         div(img(img(src ="workflow_shiny.png", height = "600", width = "600")), 
                         style="text-align: center"), 
                         br()
                         ) # end col
                  ), # end box
              box(solidHeader = F, collapsible = F, collapsed = F,
                  width = 6,
                  column(12,
                         DT::dataTableOutput("AnnoTable")
                  )
                  )
              ), # end row
            fluidRow(
              box(title = "CpDAA Sources",
                  solidHeader = T, collapsible = F, collapsed = F,
                  width = 6,
                  column(12, 
                         
                         
                         tags$ol(
                           tags$li("Weerapana, E., C. Wang, G. M. Simon, F. Richter, S. Khare, M. B. Dillon, D. A. Bachovchin, K. Mowen, D. Baker, and B. F. Cravatt. 2010. 'Quantitative Reactivity Profiling Predicts Functional Cysteines in Proteomes.' Nature 468 (7325): 790–95."),
                           tags$li("Hacker, S. M., K. M. Backus, M. R. Lazear, S. Forli, B. E. Correia, and B. F. Cravatt. 2017. 'Global Profiling of Lysine Reactivity and Ligandability in the Human Proteome.' Nature Chemistry 9 (12): 1181–90."), 
                           tags$li("Cysteine isoTOP 2019 coming soon"), style="padding-left: 2em"),
                         br(),
                         h4("Summary of CpDAA in 3,840 Proteins"), 
                         DT::dataTableOutput("comboCKmean")
                  ) # end col
              ), # end box
              
              tabBox(
                title = "Database Files", width = 6,
                # The id lets us use input$tabset1 on the server to find the current tab
                id = "databaseTabs", height = "250px",
                tabPanel("AA-level", 
                         DT::dataTableOutput("aa1Table"),
                         fluidRow(
                           column(width=5),
                           column(width=7,
                                  br(),
                                  uiOutput("downloadaa1Table"), #  delay showing download button
                                  br()
                           )) # end fluid row
                         ),
                tabPanel("Missense-level", 
                         DT::dataTableOutput("aa2Table"),
                         fluidRow(
                           column(width=5),
                           column(width=7,
                                  br(),
                                  uiOutput("downloadaa2Table"), #  delay showing download button
                                  br()
                           ))
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
              )
            ), # end row
            fluidRow(
              box(title = "Protein detected frequency and clinvar plot",
                  solidHeader = T, collapsible = F, collapsed = F,
                  width = 6,
                  column(12, 
                         h5("...")
                  ) # end col
              ) # end box
            ), # end row
            fluidRow(
              box(title = "App Info",
                  solidHeader = T, collapsible = F, collapsed = F,
                  width = 6,
                  column(12, 
                         h6("CpDAA is currently in development (Version 0.1)"),
                         h6("Contact email: "),
                         h6("citation hyperlink")
                  ) # end col
              ) # end box
            ) # end row
    ), # end tabItem

############## ANALYSIS #####################################
  tabItem(tabName="aalevel1",
    fluidPage(title = "Detected vs Undetected AA Analysis",
      tabsetPanel(     
        tabPanel(title = h4("Line plots"),  
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
                     h3("aa level df was here")
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
        )  # end tab panel
      ) # end tabset panel w/ OPTIONS 
    )  # end fluid page 
  ) # end tabItem

########## end #############################
  ) # end tab items
) # end dashboard
ui <- dashboardPage(header, sidebar, body, skin = "blue")
