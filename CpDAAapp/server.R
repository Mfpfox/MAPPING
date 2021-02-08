# source("install/prep.R")
## peers certificate error fix- useEnsembl() biomart issue
library(httr) # added set_config(config(ssl_verifypeer = 0L)) to ID mapping sections
library(shiny)
library(shinydashboard)
library(DT)
library(gplots)
library(grDevices)
library(tibble)
library(filesstrings)
library(grid)
library(gridExtra)
library(XML)
library(xlsx)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(plotly)
library(tidyr)


## trim function used for comparision group names' processing
trim <- function (x) gsub("^\\s+|\\s+$", "", x)


## Uploading file's size requirement (<500MB)
options(shiny.maxRequestSize = 500*1024^2)


shinyServer(function(input, output, session) {
# START observing data  ###############
  progress <- reactiveValues(time=shiny::Progress$new(style = "old"))

  # annotation sources table
  annoTable <- read.csv(paste(getwd(), "data/annoTable2.csv", sep="/"))

  # missense level df for CpDAA overlapping pathogenic/likely pathogenic ClinVar missense var
  clinvar <- read.csv(paste(getwd(), 
                            "data/clinvar_CpDAAoverlap_GRCh37_64missense.csv", sep="/"), stringsAsFactors = FALSE)

  # unique counts summary table
  countTable <- read.csv(paste(getwd(), 
                                "data/summary_Det_clinvar_accounting.csv", sep="/"))

  # all aa level annotations for CpDAA and undetected equivalent AA
  dataObs <- reactiveValues(
    aadf = read.csv(paste(getwd(),
                          "data/AAlevel_allLet_4526UKB_16999DET_2690597UNDET_v4.csv",sep="/"))
  )

  # protein level df
  xref <- read.csv(paste(getwd(),
                        "data/proteinLevel_xref_4526IDs.csv",sep="/")) %>% 
    select(matched.UKBID, CpDK.count, CpDC.count, 
           everCpDAA.count,Gene.UKB, Protein.UKB, 
           OMIM, gene.UKBID, Canonical.pro.isoform, CCDS, 
           RefSeq, GO.biologicalProcess,GO.cellularComponent, 
           GO.moleculeFunction, PDB.xref.UKB,
           ENST.85,ENSP.85,
           pLI, oe_lof_upper)

  # drop down menu to select 1 of 4526 genename ukbID (TMPO redundancy)
  HGNCdf <- read.csv(paste(getwd(), 
    "data/geneName_sort_4526.csv", sep="/"))

  # annotation table
  output$AnnoTable <- DT::renderDataTable({
    annoTable
  }, escape = FALSE, 
  options = list(dom = 't', rownames = FALSE, scrollX = TRUE, columnDefs = list(list(targets = "_all")))
  )

  ## unique CpDAA counts table
  output$countTable <- DT::renderDataTable({
    countTable
  },
  options = list(dom = 't', scrollX = TRUE, columnDefs = list(list(targets = "_all")))
  )

  ## protein level table
  output$proTable <- DT::renderDataTable(
    DT::datatable(xref, escape=FALSE, 
                  options = list(
                    pageLength = 5, autoWidth = TRUE,
                    columnDefs = list(list(targets = c(6,12,13,14), 
                                           width = '600px')), scrollX = TRUE
                  )))

  ## DT clinvar missense overlapping AA positions
  output$clinvarTable <- DT::renderDataTable({
    clinvar
  }, filter="none",
  options = list(scrollX = TRUE, columnDefs = list(list(targets = "_all")), pageLength = 5)
  )

  ## AA level DF
  aadf <- reactive({
    df <- dataObs$aadf %>% 
      mutate(AAgroup = factor(AAgroup, levels = c( # TODO dont use col
        "Detected C",
        "Detected K"
      ))) %>%
      mutate(detected = factor(detected, levels = c(
        "Undetected", # 1
        "Detected" # 2
      ))) %>%
      transform(detectedInt = as.integer(detected)) %>% 
      mutate(matched.aapos = as.integer(matched.aapos))
    return(df)
  }) # end aadf()
  
  ## DT table aadf
  output$aa1Table <- DT::renderDataTable({
    df <- aadf() %>% 
      filter(detected == "Detected")
    df <- df %>%
      select(matched.UKBID,AApos,reactivity.2012,react.threshold.2012,reactivity.2019,react.threshold.2019,target.label.2012,CV.AlleleID,overlap.CV.sumIDsPerAA,CADD38.phred.mean,CADD38.phred.max,DANN.mean,DANN.max,fathmmMKL.mean,fathmmMKL.max)
    df
  },
  options = list(scrollX = TRUE, columnDefs = list(list(targets = "_all")),
                 pageLength = 5, rownames= FALSE)
  )

  ## select gene name 1 aa level
  updateSelectizeInput(session, "selectGene1",
                       choices = as.vector(HGNCdf$HGNCname), server = TRUE)
  
  ## following submit select protein id  
  proteinInfo <- eventReactive(input$aasubmit1 , {
    gsym = as.character(trim(input$selectGene1))
    ukb <- xref %>% 
      filter(gene.UKBID == gsym) # updated colname for gene + UKBID
    return(ukb)
  })

  ## DT table showing UKB ID
  output$proteinInfo <- DT::renderDataTable({
    df <- data.frame(t(proteinInfo()))
    names(df) <- "CpDAA Protein Summary"
    df
  },
  options = list(scrollX = TRUE, columnDefs = list(list(targets = "_all")),
                 pageLength = 5)
  )

# select only 1 protin by gene UKB ID column for line plots
  proaadf <- eventReactive(input$aasubmit1 , {
    gsym = as.character(trim(input$selectGene1))
    ukb <- aadf() %>% 
      filter(gene.UKBID == gsym) # updated colname for gene + UKBID
    return(ukb)
  })

# table of lineplot data
  output$aa2Table <- DT::renderDataTable({
    proaadf()
  },
  options = list(scrollX = TRUE, columnDefs = list(list(targets = "_all")),
                 pageLength = 5, rownames= FALSE)
  )

  ## download buttons ##
  output$downloadProTable <- renderUI({
    req(xref) # check if dt = NULL
    downloadButton("DownloadProteins", class="btn-default")
  })
  output$DownloadProteins <- downloadHandler(
    filename = function() {
      paste("CpDAA_proteinLevel_", Sys.Date(),".csv", sep="")
    },
    content = function(file) {
      write.table(xref, file, sep =",", quote = T, row.names = F, col.names = T)
    }, contentType = "text"
  ) # end downloadhandler
  
  output$downloadAaTable <- renderUI({
    req(aadf()) 
    downloadButton("DownloadCpDAA", class="btn-default")
  })
  output$DownloadCpDAA <- downloadHandler(
    filename = function() {
      paste("CpDAA_", Sys.Date(),".csv", sep="")
    },
    content = function(file) {
      write.table(aadf(), file, sep =",", quote = T, row.names = F, col.names = T)
    }, contentType = "text"
  ) # end downloadhandler
  
  output$downloadCVTable <- renderUI({
    req(clinvar)
    downloadButton("DownloadClinVarOverlap", class="btn-default")
  })
  output$DownloadClinVarOverlap <- downloadHandler(
    filename = function() {
      paste("CpDAA_wPathoClinVar_", Sys.Date(),".csv", sep="")
    },
    content = function(file) {
      write.table(clinvar, file, sep =",", quote = T, row.names = F, col.names = T)
    }, contentType = "text"
  ) # end downloadhandler

  
  
  # page 2 downloads #
  
  output$downloadProTable2 <- renderUI({
    req(proteinInfo()) # check if dt = NULL
    downloadButton("DownloadProteins2", class="btn-default")
  })
  output$DownloadProteins2 <- downloadHandler(
    filename = function() {
      paste(input$selectGene1,".csv", sep="")
    },
    content = function(file) {
      write.table(proteinInfo(), file, sep =",", quote = T, row.names = F, col.names = T)
    }, contentType = "text"
  ) # end downloadhandler
  
  output$downloadAaTable2 <- renderUI({
    req(proaadf()) 
    downloadButton("DownloadCpDAA2", class="btn-default")
  })
  output$DownloadCpDAA2 <- downloadHandler(
    filename = function() {
      paste(input$selectGene1, "_CpDAA",".csv", sep="")
    },
    content = function(file) {
      write.table(proaadf(), file, sep =",", quote = T, row.names = F, col.names = T)
    }, contentType = "text"
  ) # end downloadhandler
  
  
  
  ########################### plots ##################
  # ALL LETTERS DF HAS FOLLOWING COLUMNS: 
  # matched.aapos,
  # missingScores (False, True values)
  # CADD38.phred.max, CADD38.phred.mean, 
  # fathmmMKL.max, fathmmMKL.mean, 
  # DANN.max, DANN.mean,
  # AApos, AAgroup

  plotCADD <- function(df, x_var, 
                       y1_var,  y2_var, 
                       ckpos_var, col_var, PROTEINname) {
    qx_var <- enquo(x_var) # aa pos
    qy1_var <- enquo(y1_var) # max
    qy2_var <- enquo(y2_var) # mean
    qckpos_var <- enquo(ckpos_var) 
    qcol_var <- enquo(col_var) # AAgroup
    ### max ###    
    plot1 <- plot_ly(data=df, x = qx_var, colors=c("#d53e4f", "#08519c"))
    plot1 <- plot1 %>% add_trace(y = qy1_var, type = 'scatter',
                                 name="Max",
                                 mode = 'lines',
                                 hoverinfo='text',
                                 text= ~paste0("<b>",AApos," ", CADD38.phred.max),
                                 line = list(color = 'rgb(37, 37, 37)',width = 1), 
                                 showlegend = TRUE,
                                 legendgroup = 1)
    ### mean ##
    plot1 <- plot1 %>% add_trace(y = qy2_var, type = 'scatter',
                                 name="Mean",
                                 mode = 'lines',
                                 line = list(color = 'rgb(189,189,189)', width = 1),
                                 hoverinfo='text',
                                 text= ~paste0("<b>",AApos," ", CADD38.phred.mean),
                                 legendgroup=2,
                                 showlegend=TRUE)
    # colored dots on max38
    plot1<- plot1 %>%  add_markers(size = ~detectedInt,
                                   opacity=0.8,
                                   stroke = I("black"), span = I(1),
                                   inherit=FALSE,
                                   y = qy1_var,
                                   x=qx_var,
                                   color= qcol_var,
                                   hoverinfo = "text",
                                   text = ~paste0("<b>",detected," ", AApos, "</b><br>",
                                                  "<b>Max: ", "</b>",CADD38.phred.max, " ", 
                                                  "<b>Mean: ","</b>", CADD38.phred.mean,
                                                  "<br><b>Reactivity ","</b>", 
                                                  reactivity.2012, " ",react.threshold.2012,
                                                  "/",reactivity.2019," ",react.threshold.2019,
                                                  "<b> Ligandable: ", "</b>",target.label.2012,
                                                  "<br><b>ClinVar Pathogenic: ", "</b>",
                                                  overlap.CV.sumIDsPerAA," ", 
                                                  "<br>", CV.AlleleID),
                                   showlegend=TRUE)
    ### LAYOUT OF ALL PLOTS ###
    TITLEplot = paste(PROTEINname)
    yaxis1 <- list(title="CADD38 Phred Score" , visible=TRUE,
                   showgrid = FALSE,zeroline = FALSE, showline = FALSE, showticklabels = TRUE
    ) #range=c(0,55)
    xaxis <- list(title = "Amino Acid Position",
                  showline = FALSE,
                  showgrid = FALSE,
                  zeroline=FALSE,
                  showticklabels = TRUE,
                  mirror=TRUE,
                  showspikes=TRUE,
                  spikethickness=2,
                  hovermode =  "closest",
                  spikedash="dot",
                  spikemode="across+toaxis+marker", # "across"
                  spikesnap= "cursor") #"cursor"
    plot1 <- plot1 %>% layout(
      yaxis = yaxis1, 
      autosize = FALSE, width=1000, height=400
    ) 
    m <- list(
      l = 100,
      r = 100,
      b = 100,
      t = 100
    )
    plot1 <- plot1 %>% layout(xaxis=xaxis,
                              title=TITLEplot,
                              autosize=FALSE,
                              margin = m,
                              hovermode =  "closest",
                              spikedistance=-1,
                              showlegend = TRUE,
                              legend = list(orientation = 'h', y = 5)) 
    return(plot1)
  }
  
  plotFathmm <- function(df, x_var, 
                         y1_var,  y2_var, 
                         ckpos_var, col_var, PROTEINname) {
    qx_var <- enquo(x_var) # aa pos
    qy1_var <- enquo(y1_var) # max
    qy2_var <- enquo(y2_var) # mean
    qckpos_var <- enquo(ckpos_var) 
    qcol_var <- enquo(col_var) # AAgroup
    ### max ###    
    plot1 <- plot_ly(data=df, x = qx_var, colors=c("#d53e4f", "#08519c"))
    plot1 <- plot1 %>% add_trace(y = qy1_var, type = 'scatter',
                                 name="Max",
                                 mode = 'lines',
                                 hoverinfo='text',
                                 text= ~paste0("<b>",AApos," ", fathmmMKL.max),
                                 line = list(color = 'rgb(37, 37, 37)',width = 1), 
                                 showlegend = TRUE,
                                 legendgroup = 1)
    ### mean ##
    plot1 <- plot1 %>% add_trace(y = qy2_var, type = 'scatter',
                                 name="Mean",
                                 mode = 'lines',
                                 line = list(color = 'rgb(189,189,189)', width = 1),
                                 hoverinfo='text',
                                 text= ~paste0("<b>",AApos," ", fathmmMKL.mean),
                                 legendgroup=2,
                                 showlegend=TRUE)
    # colored dots on max
    plot1<- plot1 %>%  add_markers(size = ~detectedInt,
                                   opacity=0.8,
                                   stroke = I("black"), span = I(1),
                                   inherit=FALSE,
                                   y = qy1_var,
                                   x=qx_var,
                                   color= qcol_var,
                                   hoverinfo = "text",
                                   text = ~paste0("<b>",detected," ", AApos, "</b><br>",
                                                  "<b>Max: ", "</b>",fathmmMKL.max, " ", 
                                                  "<b>Mean: ","</b>", fathmmMKL.mean,
                                                  "<br><b>Reactivity ","</b>", 
                                                  reactivity.2012, " ",react.threshold.2012,
                                                  "/",reactivity.2019," ",react.threshold.2019,
                                                  "<b> Ligandable: ", "</b>",target.label.2012,
                                                  "<br><b>ClinVar Pathogenic: ", "</b>",
                                                  overlap.CV.sumIDsPerAA," ", 
                                                  "<br>", CV.AlleleID),
                                   showlegend=TRUE)
    ### LAYOUT OF ALL PLOTS ###
    TITLEplot = paste(PROTEINname)
    yaxis1 <- list(title="Fathmm-MKL Score" , visible=TRUE,
                   showgrid = FALSE,zeroline = FALSE, showline = FALSE, showticklabels = TRUE
    ) #range=c(0,55)
    xaxis <- list(title = "Amino Acid Position",
                  showline = FALSE,
                  showgrid = FALSE,
                  zeroline=FALSE,
                  showticklabels = TRUE,
                  mirror=TRUE,
                  showspikes=TRUE,
                  spikethickness=2,
                  hovermode =  "closest",
                  spikedash="dot",
                  spikemode="across+toaxis+marker", # "across"
                  spikesnap= "cursor") #"cursor"
    plot1 <- plot1 %>% layout(
      yaxis = yaxis1, 
      autosize = FALSE, width=1000, height=400
    ) 
    m <- list(
      l = 100,
      r = 100,
      b = 100,
      t = 100
    )
    plot1 <- plot1 %>% layout(xaxis=xaxis,
                              title=TITLEplot,
                              autosize=FALSE,
                              margin = m,
                              hovermode =  "closest",
                              spikedistance=-1,
                              showlegend = TRUE,
                              legend = list(orientation = 'h', y = 5)) 
    return(plot1)
  }
  
  plotDANN <- function(df, x_var, 
                       y1_var,  y2_var, 
                       ckpos_var, col_var, PROTEINname) {
    
    qx_var <- enquo(x_var) # aa pos
    qy1_var <- enquo(y1_var) # max
    qy2_var <- enquo(y2_var) # mean
    qckpos_var <- enquo(ckpos_var) 
    qcol_var <- enquo(col_var) # AAgroup
    ### max ###    
    plot1 <- plot_ly(data=df, x = qx_var, colors=c("#d53e4f", "#08519c"))
    plot1 <- plot1 %>% add_trace(y = qy1_var, type = 'scatter',
                                 name="Max",
                                 mode = 'lines',
                                 hoverinfo='text',
                                 text= ~paste0("<b>",AApos," ", DANN.max),
                                 line = list(color = 'rgb(37, 37, 37)',width = 1), 
                                 showlegend = TRUE,
                                 legendgroup = 1)
    ### mean ##
    plot1 <- plot1 %>% add_trace(y = qy2_var, type = 'scatter',
                                 name="Mean",
                                 mode = 'lines',
                                 line = list(color = 'rgb(189,189,189)', width = 1),
                                 hoverinfo='text',
                                 text= ~paste0("<b>",AApos," ", DANN.mean),
                                 legendgroup=2,
                                 showlegend=TRUE)
    # colored dots on max
    plot1<- plot1 %>%  add_markers(size = ~detectedInt,
                                   opacity=0.8,
                                   stroke = I("black"), span = I(1),
                                   inherit=FALSE,
                                   y = qy1_var,
                                   x=qx_var,
                                   color= qcol_var,
                                   hoverinfo = "text",
                                   text = ~paste0("<b>",detected," ", AApos, "</b><br>",
                                                  "<b>Max: ", "</b>",DANN.max, " ", 
                                                  "<b>Mean: ","</b>", DANN.mean,
                                                  "<br><b>Reactivity ","</b>", 
                                                  reactivity.2012, " ",react.threshold.2012,
                                                  "/",reactivity.2019," ",react.threshold.2019,
                                                  "<b> Ligandable: ", "</b>",target.label.2012,
                                                  "<br><b>ClinVar Pathogenic: ", "</b>",
                                                  overlap.CV.sumIDsPerAA," ", 
                                                  "<br>", CV.AlleleID),
                                   showlegend=TRUE)
    ### LAYOUT OF ALL PLOTS ###
    TITLEplot = paste(PROTEINname)
    yaxis1 <- list(title="DANN Score" , visible=TRUE,
                   showgrid = FALSE,zeroline = FALSE, showline = FALSE, showticklabels = TRUE
    ) 
    xaxis <- list(title = "Amino Acid Position",
                  showline = FALSE,
                  showgrid = FALSE,
                  zeroline=FALSE,
                  showticklabels = TRUE,
                  mirror=TRUE,
                  showspikes=TRUE,
                  spikethickness=2,
                  hovermode =  "closest",
                  spikedash="dot",
                  spikemode="across+toaxis+marker", # "across"
                  spikesnap= "cursor") 
    plot1 <- plot1 %>% layout(
      yaxis = yaxis1, 
      autosize = FALSE, width=1000, height=400
    ) 
    m <- list(
      l = 100,
      r = 100,
      b = 100,
      t = 100
    )
    plot1 <- plot1 %>% layout(xaxis=xaxis,
                              title=TITLEplot,
                              autosize=FALSE,
                              margin = m,
                              hovermode =  "closest",
                              spikedistance=-1,
                              showlegend = TRUE,
                              legend = list(orientation = 'h', y = 5)) 
    return(plot1)
  }
  
  

  
  # func calls
  output$caddPlotly1 <- renderPlotly({
    v2 <- plotCADD(proaadf(), matched.aapos,
                         CADD38.phred.max, CADD38.phred.mean, 
                         AApos, AAgroup,
                         isolate(input$selectGene1))
    
    v2})
  
  output$fathmmPlotly1 <- renderPlotly({
    v2 <- plotFathmm(proaadf(), matched.aapos,
                         fathmmMKL.max, fathmmMKL.mean, 
                         AApos, AAgroup,
                         isolate(input$selectGene1))
    
    v2})
  
  output$dannPlotly1 <- renderPlotly({
    v2 <- plotDANN(proaadf(), matched.aapos,
                        DANN.max, DANN.mean,
                        AApos, AAgroup,
                        isolate(input$selectGene1))
    
    v2})
  
#######################################
  
  output$welcomeLink <- renderUI({
      isolate({
        actionButton('switchTab1', label = "GO TO DATA", class="btn-basic")
      })
  })
  
  output$lineplotsLink <- renderUI({
      isolate({
        actionButton('switchTab2', label = "BACK TO WELCOME", class="btn-basic")
      })
  })
  
  observeEvent(input$switchTab1, {
    newtab <- switch(input$menu1, "Welcome" = "aalevel1","aalevel1" = "Welcome")
    updateTabItems(session, "menu1", newtab)
  })
  
  observeEvent(input$switchTab2, {
    newtab <- switch(input$menu1, "aalevel1" = "Welcome","Welcome" = "aalevel1")
    updateTabItems(session, "menu1", newtab)
  })
  
})
### END shiny server ###
