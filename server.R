
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
# library(tidyverse)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(plotly)
library(tidyr)
#---bioconductor packages
# library(DO.db)
#library(clusterProfiler)
#library(enrichplot)
#library(pathview)
#library(DOSE)
#library(Homo.sapiens)
library(enrichR)
#organism_human = "org.Hs.eg.db"
#organism_mouse ='org.Mm.eg.db'
#library(organism_human, character.only = TRUE)
#library(organism_mouse, character.only = TRUE)

## trim function used for comparision group names' processing
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

## Uploading file's size requirement (<30MB)
options(shiny.maxRequestSize = 30*1024^2)

## Start shiny server
shinyServer(function(input, output, session) {


# START observing data  ####
#######################################
  #missdf

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
  aadf <- read.csv(paste(getwd(),
                        "data/AAlevel_CpDAA_undetected_4526UKBID.csv",sep="/")) 
  
  dataObs <- reactiveValues(
    aadf = read.csv(paste(getwd(),
                          "data/AAlevel_CpDAA_undetected_4526UKBID.csv",sep="/"))
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


  ## Reactive expression objects
  aadf <- reactive({
    df <- dataObs$aadf %>% 
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
      mutate(detected = factor(detected, levels = c(
        "Detected",
        "Undetected"
      ))) %>%
      mutate(matched.aapos = as.integer(matched.aapos))
    return(df)
  }) # end aadf()
  
  ## DT table aadf
  output$aa1Table <- DT::renderDataTable({
    df <- aadf() %>% 
      filter(detected == "Detected")
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
  options = list(scrollX = TRUE, columnDefs = list(list(targets = "_all")))
  )

  proaadf <- eventReactive(input$aasubmit1 , {
    gsym = as.character(trim(input$selectGene1))
    ukb <- aadf() %>% 
      filter(gene.UKBID == gsym) # updated colname for gene + UKBID
    return(ukb)
  })
  
  output$aa2Table <- DT::renderDataTable({
    proaadf()
  },
  options = list(scrollX = TRUE, columnDefs = list(list(targets = "_all")),
                 pageLength = 5, rownames= FALSE)
  )

 
  ###### LINE PLOTS ##########################

  # [1] CADD single line plot
  makeCaddLines <- function(df, x_var, y1_var, y2_var, y3_var, y4_var, 
                            ckpos_var, col_var, PROTEINname, SCOREname) {
    qx_var <- enquo(x_var)
    qy1_var <- enquo(y1_var)
    qy2_var <- enquo(y2_var)
    qy3_var <- enquo(y3_var)
    qy4_var <- enquo(y4_var)
    qckpos_var <- enquo(ckpos_var)
    qcol_var <- enquo(col_var)
    TITLEplot = paste("AA-level summary of",
                      SCOREname, "scores for", PROTEINname)
    yaxis <- list(title = SCOREname,
                  showgrid = FALSE,zeroline = FALSE,showline = FALSE,
                  showticklabels = TRUE)
    
    margin <- list(autoexpand = TRUE,t = 110,b = 110) # copy 
    
    xaxis <- list(title = "AA position",showline = FALSE,
                  showgrid = FALSE,zeroline=FALSE,showticklabels = TRUE)
    
    plot <- plot_ly(data=df, x = qx_var, colors=c("#d53e4f", "#fcc5c0", "#08519c", "#9ecae1"))
    #cadd38 max
    plot <- plot %>% add_trace(y = qy1_var, type = 'scatter',
                               name="GRCh38 max",
                               mode = 'lines',
                               line = list(color = 'rgb(37, 37, 37)',
                                           width = 2), legendgroup = 'group1')
    # # cadd38 mean
    plot <- plot %>% add_trace(y = qy2_var, type = 'scatter',
                               name="GRCh38 mean",
                               mode = 'lines',
                               line = list(color = 'rgb(37, 37, 37)',
                                           dash = 'dot', width = 1), legendgroup = 'group1')
    # # cadd37 max
    plot <- plot %>% add_trace(y = qy3_var, type = 'scatter',
                               name="GRCh37 max",
                               mode = 'lines',
                               line = list(color = 'rgb(189,189,189)',
                                           width = 2), legendgroup = 'group2')
    # # cadd37 mean
    plot <- plot %>% add_trace(y = qy4_var, type = 'scatter',
                               name="GRCh37 mean",
                               mode = 'lines',
                               line = list(color = 'rgb(189,189,189)',
                                           dash='dot', width = 1), legendgroup = 'group2')
    # # # colored dots
    plot<- plot %>%  add_markers(size = 5,
                                 opacity=0.95,
                                 stroke = I("black"), span = I(1),
                                 inherit=FALSE,
                                 y = qy1_var,
                                 x=qx_var,
                                 color= qcol_var,
                                 #text = qckpos_var, 
                                 hovertemplate = ~paste0("<b>",CKpos,"</b>"),
                                 showlegend=TRUE)
    
    plot <- plot %>% layout(title = TITLEplot,
                            xaxis = xaxis,
                            yaxis = yaxis,
                            margin = margin,
                            autosize = TRUE,
                            showlegend = TRUE,
                            hovermode = "x unified",
                            legend = list(orientation = 'h', y = -0.3)) 
    return(plot)
  }
  
  
  # [2] for dann and fathmm ck combo plots
  makeLinePlots1assembly <- function(df, x_var, y1_var, y2_var, ckpos_var, col_var, PROTEINname, SCOREname) {
    qx_var <- enquo(x_var)
    qy1_var <- enquo(y1_var)
    qy2_var <- enquo(y2_var)
    qckpos_var <- enquo(ckpos_var)
    qcol_var <- enquo(col_var)
    
    TITLEplot = paste("AA-level summary of",
                      SCOREname, "scores for", PROTEINname)
    
    yaxis <- list(title = SCOREname,
                  showgrid = FALSE,zeroline = FALSE,showline = FALSE,
                  showticklabels = TRUE)
    
    margin <- list(autoexpand = TRUE,t = 110, b=110)
    
    xaxis <- list(title = "AA position",showline = FALSE,
                  showgrid = FALSE,zeroline=FALSE,showticklabels = TRUE)
    
    plot <- plot_ly(data=df, x = qx_var, colors=c("#d53e4f", "#fcc5c0", "#08519c", "#9ecae1"))
    
    plot <- plot %>% add_trace(y = qy1_var, type = 'scatter',
                               name="Max",
                               mode = 'lines',
                               line = list(color = 'rgb(37, 37, 37)',
                                           width = 2))
    
    plot <- plot %>% add_trace(y = qy2_var, type = 'scatter',
                               name="Mean",
                               mode = 'lines',
                               line = list(color = 'rgb(37, 37, 37)',
                                           dash = 'dot', width = 1))
    
    # # # colored dots
    plot<- plot %>%  add_markers(size = 5,
                                 opacity=0.95,
                                 stroke = I("black"), span = I(1),
                                 inherit=FALSE,
                                 y = qy2_var,
                                 x=qx_var,
                                 color= qcol_var,
                                 #text = qckpos_var,
                                 hovertemplate = ~paste0("<b>",CKpos,"</b>"),
                                 showlegend=TRUE)
    
    
    plot <- plot %>% layout(title = TITLEplot,
                            xaxis = xaxis,
                            yaxis = yaxis,
                            margin = margin,
                            autosize = TRUE,
                            showlegend = TRUE,
                            hovermode = "x unified",
                            legend = list(orientation = 'h', y = -0.3))
    return(plot)
  }
  
  
  
  # [3] CADD CK specific line plots
  makeLinePlotsCK <- function(df, x_var, y_var, reactivity,line_var, group_var, aa_var, PROTEINname, SCOREname) {
    reactivity <- enquo(reactivity)  
    xvar <- enquo(x_var)
    yvar <- enquo(y_var)
    linevar <- enquo(line_var)
    groupvar <- enquo(group_var)
    aavar <- enquo(aa_var)
    
    TITLEplot = paste("Cysteine & Lysine mean", SCOREname, "scores for", PROTEINname)
    
    yaxis <- list(title = SCOREname, showgrid = FALSE, zeroline = FALSE, showline = FALSE, 
                  showticklabels = TRUE)
    
    margin <- list(autoexpand = TRUE, t = 110, b=110)
    
    xaxis <- list(title = "AA position",showline = FALSE,
                  showgrid = FALSE,zeroline=FALSE,showticklabels = TRUE)
    
    ck <- plot_ly(data=df, x = xvar, 
                  colors=c("#e31a1c", "#225ea8", # ck
                           "#dd3497","#ae017e","#7a0177", # high
                           "#efedf5","#bcbddc","#756bb1", # med
                           "#deebf7","#9ecae1","#3182bd", # low
                           "#02818a", "#016c59", # tar pan
                           "#ece2f0")) # undetected
    
    ck <- ck %>% 
      add_trace(y=yvar, color=linevar, type='scatter', mode='lines',
                line=list(width=1), hoverinfo="none")
    
    #RAW CADD
    if(SCOREname == "CADD38 RAW"){
      ck <- ck %>% add_markers(size = ~inClinVar.5.10,
                               opacity=0.95,
                               stroke = I("black"), span = I(1),
                               x=xvar,
                               y=yvar,
                               inherit=FALSE,
                               color=groupvar,  
                               hoverinfo = "text",
                               text = ~paste0("<b>",detected.group," ",CKpos,"</b><br><br>",
                                              "Max score: ", CADD38.raw.max,
                                              "<br>Mean score: ", CADD38.raw.mean,
                                              "<br>Reactivity: ",reactivity,
                                              "<br>ClinVar patho&likelypatho AA overlap: ", overlap.CV.sumIDsPerAA,
                                              "<br>ClinVar IDs: ", CV.AlleleID))
    }
    #PHRED CADD
    if(SCOREname == "CADD38 PHRED"){
      ck <- ck %>% add_markers(size = ~inClinVar.5.10,
                               opacity=0.95,
                               stroke = I("black"), span = I(1),
                               x=xvar,
                               y=yvar,
                               inherit=FALSE,
                               color=groupvar,  
                               hoverinfo = "text",
                               text = ~paste0("<b>", detected.group, " ", CKpos, "</b><br><br>",
                                              "Max score: ", CADD38.phred.max,
                                              "<br>Mean score: ", CADD38.phred.mean,
                                              "<br>Reactivity: ",reactivity,
                                              "<br>ClinVar patho&likelypatho AA overlap: ", overlap.CV.sumIDsPerAA,
                                              "<br>ClinVar IDs: ", CV.AlleleID))
    }
    #FATHMM
    if(SCOREname == "FATHMM-mkl"){
      ck <- ck %>% add_markers(size = ~inClinVar.5.10,
                               opacity=0.95,
                               stroke = I("black"), span = I(1),
                               x=xvar,
                               y=yvar,
                               inherit=FALSE,
                               color=groupvar,  
                               hoverinfo = "text",
                               text = ~paste0("<b>", detected.group, " ", CKpos, "</b><br><br>",
                                              "Max score: ", fathmmMKL.max,
                                              "<br>Mean score: ", fathmmMKL.mean,
                                              "<br>Reactivity: ",reactivity,
                                              "<br>ClinVar patho&likelypatho AA overlap: ", overlap.CV.sumIDsPerAA,
                                              "<br>ClinVar IDs: ", CV.AlleleID))
    }
    #DANN
    if(SCOREname == "DANN"){
      ck <- ck %>% add_markers(size = ~inClinVar.5.10,
                               opacity=0.95,
                               stroke = I("black"), span = I(1),
                               x=xvar,
                               y=yvar,
                               inherit=FALSE,
                               color=groupvar, 
                               hoverinfo = "text",
                               text = ~paste0("<b>", detected.group, " ", CKpos, "</b><br><br>",
                                              "Max score: ", DANN.max,
                                              "<br>Mean score: ", DANN.mean,
                                              "<br>Reactivity: ",reactivity,
                                              "<br>ClinVar patho&likelypatho AA overlap: ", overlap.CV.sumIDsPerAA,
                                              "<br>ClinVar IDs: ", CV.AlleleID))
                               
                               
    }
    ck <- ck %>% layout(title = TITLEplot,
                        xaxis = xaxis,
                        yaxis = yaxis,
                        margin = margin,
                        autosize = TRUE,
                        legend = list(orientation = 'h', y = -0.3,
                                      legendgroup="",
                                      font = list(family = "sans-serif", 
                                                  size = 12, color = "#000")))
    
                                     #  bgcolor = "#E2E2E2")), bordercolor = "black",borderwidth = 1))
    return(ck)
  }
  
  

  ####### CADD PLOTLY 1
  output$caddPlotly1 <- renderPlotly({
    if(input$caddTypeInput == "Raw"){
      if(input$lineTypeCADD == "Combo"){
        cadd <- makeCaddLines(proaadf(), matched.aapos, CADD38.raw.max,
                              CADD38.raw.mean, CADD37.raw.max,
                              CADD37.raw.mean, CKpos, AAgroup,
                              input$selectGene1, "CADD RAW")
       }
      if(input$lineTypeCADD == "CK-specific"){
        cadd <- makeLinePlotsCK(proaadf(), matched.aapos, CADD38.raw.mean,
                                reactivity, aaref, detected.group, AAgroup,
                                input$selectGene1, "CADD38 RAW")
      }
    } # end raw
    if(input$caddTypeInput == "PHRED"){
      if(input$lineTypeCADD == "Combo"){
        cadd <- makeCaddLines(proaadf(), matched.aapos, CADD38.phred.max,
                              CADD38.phred.mean, CADD37.phred.max,
                              CADD37.phred.mean, CKpos, AAgroup,
                              input$selectGene1, "CADD PHRED")
      }
      if(input$lineTypeCADD == "CK-specific"){
        cadd <- makeLinePlotsCK(proaadf(), matched.aapos, CADD38.phred.mean,
                                  reactivity, aaref, detected.group, AAgroup,
                                  input$selectGene1, "CADD38 PHRED")
      }
    } # end phred
    cadd
  })

  # ###### FATHMM PLOTLY 1
  output$fathmmPlotly1 <- renderPlotly({
    if(input$lineTypeFathmm == "Combo"){
      fathmm <- makeLinePlots1assembly(proaadf(), matched.aapos, fathmmMKL.max,
                                       fathmmMKL.mean, CKpos, AAgroup,
                                       input$selectGene1, "FATHMM-mkl")
    }
    if(input$lineTypeFathmm == "CK-specific"){
      fathmm <- makeLinePlotsCK(proaadf(), matched.aapos, fathmmMKL.mean,
                                reactivity, aaref, detected.group, AAgroup,
                                input$selectGene1, "FATHMM-mkl")
    }
    fathmm
  })

  # ######## DANN PLOTLY 1
  output$dannPlotly1 <- renderPlotly({
    if(input$lineTypeDann == "Combo"){
      dann <- makeLinePlots1assembly(proaadf(), matched.aapos,
                                     DANN.max, DANN.mean,
                                     CKpos, AAgroup,input$selectGene1, "DANN")
    }
    if(input$lineTypeDann == "CK-specific"){
      dann <- makeLinePlotsCK(proaadf(), matched.aapos, DANN.mean,
                         reactivity, aaref, detected.group, AAgroup,
                         input$selectGene1, "DANN")
    }
    dann
  })

#######################################


})
### END shiny server ###
