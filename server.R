
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
  ## aalevel1 file
  ## aalevel2 file
  ## clinvar missense-level det/undet overlap file
  ## clinvar summary file for AA level groups
  ## missense level all possible score summary for AA level groups
  dataObs <- reactiveValues(
    aadf = read.csv(paste(getwd(),
                          "data/AAlevel6cvIDfix_gnomadScores_detectedCount_189627.csv",sep="/")),
    missdf = read.csv(paste(getwd(),
                            "data/shiny_detectedOnly_missenseLevel_scores_VEPpopulation_104475.csv", sep="/")),
    clinvar = read.csv(paste(getwd(), 
                              "data/clinvarcombo_patho&likelyFiltered_det&undet_allScores_389.csv", sep="/"), stringsAsFactors = FALSE),
    clinvarMean = read.csv(paste(getwd(), 
                                  "data/clinvarSummary_389patho&likely_overlap_det&undet_scoreMeans.csv", sep="/")),
    comboCKmean = read.csv(paste(getwd(), 
                                  "data/missenseScoreSummary_det&undet_Means_1327386.csv", sep="/")),
    enrichdf = read.csv(paste(getwd(),
                              "data/UKB_HGNCdatamapped_sumCpDAA_3840_TMPOabdropped.csv", sep = "/"))
    ) 
  
  
  
  
  # drop down menu to select 1 of 3840 genenames
  HGNCdf <- read.csv(paste(getwd(), "data/HGNCname_sorted_unique_3840_TMPOab_added.csv", sep="/"))
  
  progress <- reactiveValues(time=shiny::Progress$new(style = "old"))

  ## Reactive expression objects
  aadf <- reactive({
    df <- dataObs$aadf
    df <- df %>%
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
      mutate(detected.group = factor(detected.group, levels = c(
        "High","High Target","High panReactive",
        "Medium","Medium Target","Medium panReactive",
        "Low","Low Target","Low panReactive",
        "Target","panReactive",
        "Undetected"
      ))) %>%
      mutate(matched.aapos = as.integer(matched.aapos))
    return(df)

  }) # end aadf()
  
  missdf <-  reactive({
    miss <- dataObs$missdf
    miss <- miss %>%
      mutate(matched.aapos = as.integer(matched.aapos)) %>% 
      mutate(detonly.group = factor(detonly.group, levels = c(
        "High C",
        "High Target C",
        "High panReactive C",
        "Medium C",
        "Medium Target C",
        "Medium panReactive C",
        "Low C",
        "Low Target C",
        "Low panReactive C",
        "Target C",
        "panReactive C",
        "High K",
        "High Target K",
        "High panReactive K",
        "Medium K",
        "Medium Target K",
        "Medium panReactive K",
        "Low K",
        "Low Target K",
        "Low panReactive K",
        "Target K",
        "panReactive K"
      ))) 
    
    miss <- miss[
      with(miss, order(matched.UKBID, matched.aapos)),
    ]
    return(miss)
  }) # end missdf()

  
  ## select gene name 1 aa level
  updateSelectizeInput(session, "selectGene1",
                       choices = as.vector(HGNCdf$HGNCname), server = TRUE)
  
  ## select gene name 2 missense level
  updateSelectizeInput(session, "selectGene2",
                       choices = as.vector(HGNCdf$HGNCname), server = TRUE)
  
  
  ## following submit select protein id  
  proaadf <- eventReactive(input$aasubmit1 , {
    gsym = as.character(trim(input$selectGene1))
    ukb <- aadf() %>% 
      filter(HGNCname == gsym)
    return(ukb)
  })
  
  ## following submit select protein id  
  popdf <- eventReactive(input$aasubmit2 , {
    gsym = as.character(trim(input$selectGene2))
    ukb <- missdf() %>% 
      filter(HGNCname == gsym)
    ukb <- droplevels(ukb)
    return(ukb)
  })
  
  enrichdf <-  reactive({
    df <- dataObs$enrichdf
    return(df)
  })

  ###### enrichr ##################
  ## following submit gene set
  detectedDF <- eventReactive(input$queryENR , {
    df <- enrichdf()
    geneset = input$selectENRsymbol # as.character(trim(input$selectENRsymbol))
    if(geneset == "Cys Detected"){
      df <- df %>% filter(sumCpDC > 0)
      return(df)
    }
    if(geneset == "Lys Detected"){
      df <- df %>% filter(sumCpDK > 0)
      return(df)
    }
    if(geneset == "All Detected"){
      return(df)
    }
  })


  ## DT det and undetect ck mean score summary table
  output$DEToutput <- DT::renderDataTable({
    df <- detectedDF() %>% select(-X)
    df
  },
  options = list(scrollX = TRUE, columnDefs = list(list(targets = "_all")), pageLength = 5)
  )
  
  ## Enrichr query based on user input
  enrichrObj <- eventReactive(input$queryENR , {
    df <- detectedDF()
    # drop duplicate symbol rows
    dfX = df[!duplicated(df[c("HGNC.symbol")]),] 
    dfX <- dfX$HGNC.symbol
    progress$time$set(value = 0.5, detail = "processing 50%")
    dbs<- listEnrichrDbs()
    dbs <- c("GO_Biological_Process_2018", #1
             "GO_Cellular_Component_2018", #2
             "GO_Molecular_Function_2018", #3
             "KEGG_2019_Human", #4
             "ChEA_2016", #5
             "Tissue_Protein_Expression_from_Human_Proteome_Map" #6
    )
    enriched <- tryCatch(
      expr={enrichr(c("Runx1", "Gfi1", "Gfi1b", "Spi1", "Gata1", "Kdr"), dbs)},
      error = function(e){ 
        ERRORMSG("Enrichr - Failure during outbound connection!")
        return(NULL)
      })
    if(is.null(enriched)) return(NULL)
    progress$time$set(value = 0.75, detail = "processing 75%")
    enriched <- tryCatch(
      expr={enrichr(dfX, dbs)},
      error = function(e){ 
        ERRORMSG("Enrichr - Failure during outbound connection!")
        return(NULL)
      })
    if(is.null(enriched)) return(NULL)
    progress$time$set(message = "Querying EnrichR", value = 0.75)
    progress$time$set(value = 1, detail = "processing 100%")
    return(enriched)
  })
  
  
  
  # preview results
  output$ENRoutput <- DT::renderDataTable({
    enrichrDF <- enrichrObj()[[as.numeric(input$selectENRresults)]] # default preview GO mf
    df <- enrichrDF[,c("Term", "Overlap", "P.value", 
                       "Adjusted.P.value", "Odds.Ratio",
                       "Combined.Score")] # Genes
    df['P.value'] <- round(df['P.value'], digits=3)
    df['Adjusted.P.value'] <- round(df['Adjusted.P.value'], digits=3)
    df['Odds.Ratio'] <- round(df['Odds.Ratio'], digits=3)
    df['Combined.Score'] <- round(df['Combined.Score'], digits=3)
    df
  },
  options = list(columnDefs = list(list(className = 'dt-center', targets = c(2,3,4,5,6))), pageLength = 5)
  )
  
  
  output$downloadDataENR <- renderUI({
    req(input$queryENR, enrichrObj())
    downloadButton("enrichrResultDF", label="Download", class="btn-default")
  })
  
  output$enrichrResultDF <- downloadHandler(
    filename = function() {
      paste("enrichR_query_", Sys.Date(), ".xlsx", sep="")
    },
    content = function(file) {
      shname = c("GO_Biological_Process_2018", #1
                 "GO_Cellular_Component_2018", #2
                 "GO_Molecular_Function_2018", #3
                 "KEGG_2019_Human", #4
                 "ChEA_2016", #5
                 "Tissue_Protein_Expression_from_Human_Proteome_Map" #6
      )
      for(i in 1:length(enrichrObj()))
      {write.xlsx(enrichrObj()[[i]], file, sheetName = shname[i], append=T, row.names=F)}
    }
  ) # end downloadhandler
  

  #--enrichr plots------------------
  ENRGObp <- reactive({ #update
    if(is.null(enrichrObj())) return(NULL)
    i <- 1 # update
    namesLS <- names(enrichrObj())
    sortcol <- as.character(trim(input$sortByGObp)) # update
    res <- enrichrObj()[[i]]
    # option 1
    if (sortcol == "Combined.Score") {
      res.cs <- res[order(res$Combined.Score, decreasing = TRUE), ]
      res.cs <- res.cs[1:5,]
      cs <- ggplot(res.cs) + geom_bar(aes(x=Combined.Score, y=reorder(Term, Combined.Score), fill =  Adjusted.P.value), stat = "identity") + scale_x_continuous(expand = c(0,0)) + scale_y_discrete(position="left") + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "white")) + ylab("") + xlab("Combined score") + theme(legend.position = "right") +  scale_fill_gradient(low = "red3", high = "lightblue", name = "P.adjust") + ggtitle(namesLS[i]) +
        theme(plot.title = element_text(hjust = 0.5))
      return(cs)
    } # option 2 default
    else {
      res.pv <- res[order(res$P.value, decreasing = FALSE), ]
      res.pv <- res.pv[1:5,]
      pv <- ggplot(res.pv) + geom_bar(aes(x=P.value, y=reorder(Term, -P.value),
                                          fill =  Adjusted.P.value), stat = "identity") +
        scale_x_continuous(trans = "log10",expand = c(0,0)) + scale_y_discrete(position="right") + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "white")) + ylab("") + xlab("P-value") +
        theme(legend.position = "left") + scale_fill_gradient(low = "red3", high = "lightblue", name = "P.adjust") + ggtitle(namesLS[i]) +
        theme(plot.title = element_text(hjust = 0.5))
      return(pv)
    }
  })
  
  output$ENRGObp <- renderPlot({ # update
    if(is.null(enrichrObj())){
      return(NULL)
    }
    ENRGObp() # update
  })
  
  
  ENRGOcc <- reactive({ #update
    if(is.null(enrichrObj())) return(NULL)
    i <- 2 # update
    namesLS <- names(enrichrObj())
    sortcol <- as.character(trim(input$sortByGOcc)) # update
    res <- enrichrObj()[[i]]
    # option 1
    if (sortcol == "Combined.Score") {
      res.cs <- res[order(res$Combined.Score, decreasing = TRUE), ]
      res.cs <- res.cs[1:5,]
      cs <- ggplot(res.cs) + geom_bar(aes(x=Combined.Score, y=reorder(Term, Combined.Score), fill =  Adjusted.P.value), stat = "identity") + scale_x_continuous(expand = c(0,0)) + scale_y_discrete(position="left") + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "white")) + ylab("") + xlab("Combined score") + theme(legend.position = "right") +  scale_fill_gradient(low = "red3", high = "lightblue", name = "P.adjust") + ggtitle(namesLS[i]) +
        theme(plot.title = element_text(hjust = 0.5))
      return(cs)
    } # option 2 default
    else {
      res.pv <- res[order(res$P.value, decreasing = FALSE), ]
      res.pv <- res.pv[1:5,]
      pv <- ggplot(res.pv) + geom_bar(aes(x=P.value, y=reorder(Term, -P.value),
                                          fill =  Adjusted.P.value), stat = "identity") +
        scale_x_continuous(trans = "log10",expand = c(0,0)) + scale_y_discrete(position="right") + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "white")) + ylab("") + xlab("P-value") +
        theme(legend.position = "left") + scale_fill_gradient(low = "red3", high = "lightblue", name = "P.adjust") + ggtitle(namesLS[i]) +
        theme(plot.title = element_text(hjust = 0.5))
      return(pv)
    }
  })
  
  output$ENRGOcc <- renderPlot({ # update
    if(is.null(enrichrObj())){
      return(NULL)
    }
    ENRGOcc() # update
  })
  
  ENRGOmf <- reactive({ #update
    if(is.null(enrichrObj())) return(NULL)
    i <- 3 # update
    namesLS <- names(enrichrObj())
    sortcol <- as.character(trim(input$sortByGOmf)) # update
    res <- enrichrObj()[[i]]
    # option 1
    if (sortcol == "Combined.Score") {
      res.cs <- res[order(res$Combined.Score, decreasing = TRUE), ]
      res.cs <- res.cs[1:5,]
      cs <- ggplot(res.cs) + geom_bar(aes(x=Combined.Score, y=reorder(Term, Combined.Score), fill =  Adjusted.P.value), stat = "identity") + scale_x_continuous(expand = c(0,0)) + scale_y_discrete(position="left") + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "white")) + ylab("") + xlab("Combined score") + theme(legend.position = "right") +  scale_fill_gradient(low = "red3", high = "lightblue", name = "P.adjust") + ggtitle(namesLS[i]) +
        theme(plot.title = element_text(hjust = 0.5))
      return(cs)
    } # option 2 default
    else {
      res.pv <- res[order(res$P.value, decreasing = FALSE), ]
      res.pv <- res.pv[1:5,]
      pv <- ggplot(res.pv) + geom_bar(aes(x=P.value, y=reorder(Term, -P.value),
                                          fill =  Adjusted.P.value), stat = "identity") +
        scale_x_continuous(trans = "log10",expand = c(0,0)) + scale_y_discrete(position="right") + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "white")) + ylab("") + xlab("P-value") +
        theme(legend.position = "left") + scale_fill_gradient(low = "red3", high = "lightblue", name = "P.adjust") + ggtitle(namesLS[i]) +
        theme(plot.title = element_text(hjust = 0.5))
      return(pv)
    }
  })
  
  output$ENRGOmf <- renderPlot({ # update
    if(is.null(enrichrObj())){
      return(NULL)
    }
    ENRGOmf() # update
  })
  
  
  ENRkegg <- reactive({ #update
    if(is.null(enrichrObj())) return(NULL)
    i <- 4 # update
    namesLS <- names(enrichrObj())
    sortcol <- as.character(trim(input$sortByKEGG)) # update
    res <- enrichrObj()[[i]]
    # option 1
    if (sortcol == "Combined.Score") {
      res.cs <- res[order(res$Combined.Score, decreasing = TRUE), ]
      res.cs <- res.cs[1:5,]
      cs <- ggplot(res.cs) + geom_bar(aes(x=Combined.Score, y=reorder(Term, Combined.Score), fill =  Adjusted.P.value), stat = "identity") + scale_x_continuous(expand = c(0,0)) + scale_y_discrete(position="left") + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "white")) + ylab("") + xlab("Combined score") + theme(legend.position = "right") +  scale_fill_gradient(low = "red3", high = "lightblue", name = "P.adjust") + ggtitle(namesLS[i]) +
        theme(plot.title = element_text(hjust = 0.5))
      return(cs)
    } # option 2 default
    else {
      res.pv <- res[order(res$P.value, decreasing = FALSE), ]
      res.pv <- res.pv[1:5,]
      pv <- ggplot(res.pv) + geom_bar(aes(x=P.value, y=reorder(Term, -P.value),
                                          fill =  Adjusted.P.value), stat = "identity") +
        scale_x_continuous(trans = "log10",expand = c(0,0)) + scale_y_discrete(position="right") + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "white")) + ylab("") + xlab("P-value") +
        theme(legend.position = "left") + scale_fill_gradient(low = "red3", high = "lightblue", name = "P.adjust") + ggtitle(namesLS[i]) +
        theme(plot.title = element_text(hjust = 0.5))
      return(pv)
    }
  })
  
  output$ENRkegg <- renderPlot({ # update
    if(is.null(enrichrObj())){
      return(NULL)
    }
    ENRkegg() # update
  })
  
  ENRchea16 <- reactive({
    if(is.null(enrichrObj())) return(NULL)
    i <- 5
    namesLS <- names(enrichrObj())
    sortcol <- as.character(trim(input$sortByChea))
    res <- enrichrObj()[[i]]
    # option 1
    if (sortcol == "Combined.Score") {
      res.cs <- res[order(res$Combined.Score, decreasing = TRUE), ]
      res.cs <- res.cs[1:5,]
      cs <- ggplot(res.cs) + geom_bar(aes(x=Combined.Score, y=reorder(Term, Combined.Score), fill =  Adjusted.P.value), stat = "identity") + scale_x_continuous(expand = c(0,0)) + scale_y_discrete(position="left") + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "white")) + ylab("") + xlab("Combined score") + theme(legend.position = "right") +  scale_fill_gradient(low = "red3", high = "lightblue", name = "P.adjust") + ggtitle(namesLS[i]) +
        theme(plot.title = element_text(hjust = 0.5))
      return(cs)
    } # option 2 default
    else {
      res.pv <- res[order(res$P.value, decreasing = FALSE), ]
      res.pv <- res.pv[1:5,]
      pv <- ggplot(res.pv) + geom_bar(aes(x=P.value, y=reorder(Term, -P.value),
                                          fill =  Adjusted.P.value), stat = "identity") +
        scale_x_continuous(trans = "log10",expand = c(0,0)) + scale_y_discrete(position="right") + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "white")) + ylab("") + xlab("P-value") +
        theme(legend.position = "left") + scale_fill_gradient(low = "red3", high = "lightblue", name = "P.adjust") + ggtitle(namesLS[i]) +
        theme(plot.title = element_text(hjust = 0.5))
      return(pv)
    }
  })
  
  output$ENRchea16 <- renderPlot({
    if(is.null(enrichrObj())){
      return(NULL)
    }
    ENRchea16()
  })
  
  
  ENRhpm <- reactive({
    if(is.null(enrichrObj())) return(NULL)
    i <- 6 # update
    namesLS <- names(enrichrObj())
    sortcol <- as.character(trim(input$sortByHpm)) # update
    res <- enrichrObj()[[i]]
    # option 1
    if (sortcol == "Combined.Score") {
      res.cs <- res[order(res$Combined.Score, decreasing = TRUE), ]
      res.cs <- res.cs[1:5,]
      cs <- ggplot(res.cs) + geom_bar(aes(x=Combined.Score, y=reorder(Term, Combined.Score), fill =  Adjusted.P.value), stat = "identity") + scale_x_continuous(expand = c(0,0)) + scale_y_discrete(position="left") + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "white")) + ylab("") + xlab("Combined score") + theme(legend.position = "right") +  scale_fill_gradient(low = "red3", high = "lightblue", name = "P.adjust") + ggtitle(namesLS[i]) +
        theme(plot.title = element_text(hjust = 0.5))
      return(cs)
    } # option 2 default
    else {
      res.pv <- res[order(res$P.value, decreasing = FALSE), ]
      res.pv <- res.pv[1:5,]
      pv <- ggplot(res.pv) + geom_bar(aes(x=P.value, y=reorder(Term, -P.value),
                                          fill =  Adjusted.P.value), stat = "identity") +
        scale_x_continuous(trans = "log10",expand = c(0,0)) + scale_y_discrete(position="right") + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "white")) + ylab("") + xlab("P-value") +
        theme(legend.position = "left") + scale_fill_gradient(low = "red3", high = "lightblue", name = "P.adjust") + ggtitle(namesLS[i]) +
        theme(plot.title = element_text(hjust = 0.5))
      return(pv)
    }
  })
  
  output$ENRhpm <- renderPlot({
    if(is.null(enrichrObj())){
      return(NULL)
    }
    ENRhpm()
  })
  
  
  #--------------------
  ###### LINE PLOTS ##########################
  
  ## DT det and undetect ck mean score summary table
  output$comboCKmean <- DT::renderDataTable({
    df <- dataObs$comboCKmean
    df
  },
  options = list(scrollX = TRUE, columnDefs = list(list(targets = "_all")), pageLength = 5)
  )
  
  ## DT clinvar overlap det and undet ck mean score summary table
  output$clinvarMean <- DT::renderDataTable({
    df <- dataObs$clinvarMean
    df
  },
  options = list(scrollX = TRUE, columnDefs = list(list(targets = "_all")), pageLength = 5)
  )
  
  
  ## DT clinvar missense overlapping det or undet CK mean score summary table
  output$clinvarTable <- DT::renderDataTable({
    clinvar <- dataObs$clinvar
    # sort clinvar df
    clinvar <- clinvar[
      with(clinvar, order(matched.UKBID, matched.aapos)),
    ]
    clinvar
  },
  options = list(scrollX = TRUE, columnDefs = list(list(targets = "_all")), pageLength = 5)
  )
  
  
  
  
  
  
  
  
  # TODO add TITLE with GENE and PROTEIN ID information # HGNCname# matched.UKBID
  
  ## DT clinvar overlap 
  output$clinvarTable <- DT::renderDataTable({
    clinvar <- dataObs$clinvar
    # sort clinvar df
    clinvar <- clinvar[
      with(clinvar, order(matched.UKBID, matched.aapos)),
    ]
    clinvar
  },
  options = list(scrollX = TRUE, columnDefs = list(list(targets = "_all")), pageLength = 8)
  )
  
  
  
  ## DT table showing UKB ID
  output$aa1Table <- DT::renderDataTable({
    df <- proaadf() %>% 
      select(pos.ID, HGNCname,detected.group,CKpos,
             CADD.raw.hg38.mean,CADD.raw.hg38.max,CADD.raw.hg19.mean,CADD.raw.hg19.max,
             CADD.phred.hg38.mean,CADD.phred.hg38.max,CADD.phred.hg19.mean,CADD.phred.hg19.max,
             DANN.score.mean,DANN.score.max,
             fathmmMKL.coding.score.mean,fathmmMKL.coding.score.max
      )
    names(df) <- c("AApos ID", "Gene", "Group", "CKpos", "CADDraw38 mean", "CADDraw38 max", 
                "CADDraw37 mean", "CADDraw37 max", "CADD38 mean", "CADD38 max", "CADD37 mean", "CADD37 max", 
                "DANN mean", "DANN max", "FATHMM mean", "FATHMM max")

    df
  },
  options = list(scrollX = TRUE, columnDefs = list(list(targets = "_all")), pageLength = 5)
  )

  # # fix me
  # ## DT table MISSENSE LEVEL
  output$aa2Table <- DT::renderDataTable({
    df <- popdf() %>%
      select(pos.id19,pos.id38,pos.ID,Amino.acids,rs.dbSNP151,Ensembl.transcriptid,cds.strand,refcodon,codonpos,codon.degeneracy,CADD.phred.hg38,CADD.phred.hg19,CADD.phreddiff.38minus19,CADD.raw.hg38,CADD.raw.hg19,CADD.rawdiff.38minus19,CV.AlleleID,CV.Name,CV.ClinicalSignificance,CV.ClinSigSimple,CV.PhenotypeIDS,CV.PhenotypeList,CV.OriginSimple,CV.Protein.position,COSMIC.ACCESSION.NUMBER,COSMIC.ONC.TSG,COSMIC.CGC.TIER,COSMIC.LEGACY.MUTATION.ID,COSMIC.AA.MUT.START,COSMIC.GENOMIC.MUTATION.ID,DISEASE,MUTATION.SIGNIFICANCE.TIER,overlap.CV.sumIDsPerAA,COS.LEGACY.MUTID.list,AAgroup,reactivity,detected.group,react.threshold,target.label,CKpos,HGNCname,SIFT4G.score,SIFT4G.pred,Polyphen2.HDIV.score,Polyphen2.HDIV.pred,Polyphen2.HVAR.score,Polyphen2.HVAR.pred,LRT.score,LRT.pred,LRT.Omega,MutationTaster.score,MutationTaster.pred,MutationTaster.model,MutationTaster.AAE,MutationAssessor.score,MutationAssessor.pred,PROVEAN.score,PROVEAN.pred,VEST4.score,MetaSVM.score,MetaSVM.pred,MetaLR.score,MetaLR.pred,Reliability.index,M.CAP.score,M.CAP.pred,REVEL.score,MutPred.score,MutPred.protID,MutPred.AAchange,MutPred.Top5features,MPC.score,PrimateAI.score,PrimateAI.pred,DANN.score,fathmm.MKL.coding.score,GC.hg19,CpG.hg19,priPhCons.hg19,mamPhCons.hg19,verPhCons.hg19,priPhyloP.hg19,mamPhyloP.hg19,verPhyloP.hg19,GerpRS.hg19,GerpRSpval.hg19,GerpN.hg19,GerpS.hg19,GC.hg38,CpG.hg38,priPhCons.hg38,mamPhCons.hg38,verPhCons.hg38,priPhyloP.hg38,mamPhyloP.hg38,verPhyloP.hg38,GerpRS.hg38,GerpRSpval.hg38,GerpN.hg38,GerpS.hg38,Gene,Codons,Existing.variation,STRAND,CCDS,ENSP,UNIPARC,RefSeq,GENE.PHENO,DOMAINS,gnomAD.AF,gnomAD.AFR.AF,gnomAD.AMR.AF,gnomAD.ASJ.AF,gnomAD.EAS.AF,gnomAD.FIN.AF,gnomAD.NFE.AF,gnomAD.OTH.AF,gnomAD.SAS.AF,MAX.AF,MAX.AF.POPS,CLIN.SIG,SOMATIC,PHENO,BLOSUM62,detonly.group,variation.size,inClinVar.orCLINSIG
      )
    # names(df) <- c("AApos ID", "Gene", "Group", "CKpos", "CADDraw38 mean", "CADDraw38 max",
    #                "CADDraw37 mean", "CADDraw37 max", "CADD38 mean", "CADD38 max", "CADD37 mean", "CADD37 max",
    #                "DANN mean", "DANN max", "FATHMM mean", "FATHMM max")

    df
  },
  options = list(scrollX = TRUE, columnDefs = list(list(targets = "_all")), pageLength = 5)
  )


  

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
                                              "Max score: ", CADD.raw.hg38.max,
                                              "<br>Mean score: ", CADD.raw.hg38.mean,
                                              "<br>Reactivity: ",reactivity,
                                              "<br>ClinVar patho&likelypatho AA overlap: ", overlap.CV.sumIDsPerAA,
                                              "<br>ClinVar IDs: ", CV.AlleleID.list,
                                              "<br>COSMIC IDs: ", COS.LEGACY.MUTID.list))
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
                                              "Max score: ", CADD.phred.hg38.max,
                                              "<br>Mean score: ", CADD.phred.hg38.mean,
                                              "<br>Reactivity: ",reactivity,
                                              "<br>ClinVar patho&likelypatho AA overlap: ", overlap.CV.sumIDsPerAA,
                                              "<br>ClinVar IDs: ", CV.AlleleID.list,
                                              "<br>COSMIC IDs: ", COS.LEGACY.MUTID.list))
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
                                              "Max score: ", fathmmMKL.coding.score.max,
                                              "<br>Mean score: ", fathmmMKL.coding.score.mean,
                                              "<br>Reactivity: ",reactivity,
                                              "<br>ClinVar patho&likelypatho AA overlap: ", overlap.CV.sumIDsPerAA,
                                              "<br>ClinVar IDs: ", CV.AlleleID.list,
                                              "<br>COSMIC IDs: ", COS.LEGACY.MUTID.list))
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
                                              "Max score: ", DANN.score.max,
                                              "<br>Mean score: ", DANN.score.mean,
                                              "<br>Reactivity: ",reactivity,
                                              "<br>ClinVar patho&likelypatho AA overlap: ", overlap.CV.sumIDsPerAA,
                                              "<br>ClinVar IDs: ", CV.AlleleID.list,
                                              "<br>COSMIC IDs: ", COS.LEGACY.MUTID.list))
                               
                               
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
  
  
  
  # [4] missense level BOXPLOTS
  makeBoxplotsMissense <- function(df, xvar, yvar,
                                   boxvar, aavar, sizevar, PROTEINname, SCOREname) {
    xvar <- enquo(xvar) #. detonly.group
    yvar <- enquo(yvar)
    boxvar <- enquo(boxvar) # detected.group
    aavar <- enquo(aavar) # Amino.acids
    sizevar <- enquo(sizevar) # variation.size
    
    
    # after unused factors  dropped
    Cdf <- df %>% filter(aaref == "C")
    Kdf <- df %>% filter(aaref == "K")
    
    
    CKsubs <- subplot(
      plot_ly(Cdf, x = ~as.numeric(detonly.group), y=yvar, color = boxvar, 
              type = "box", boxpoints = FALSE, hoverinfo = 'y', alpha = 0.3, 
              colors=c("#e31a1c", "#225ea8", # ck
                       "#dd3497","#ae017e","#7a0177", # high
                       "#efedf5","#bcbddc","#756bb1", # med
                       "#deebf7","#9ecae1","#3182bd", # low
                       "#02818a", "#016c59", # tar pan
                       "#ece2f0"),
              legendgroup= "g1") %>%            
        add_markers(alpha = 0.8,
                    x=~jitter(as.numeric(detonly.group)),
                    color = ~as.factor(Amino.acids), 
                    size = sizevar,
                    stroke = I("black"), span = I(1),
                    marker = list(sizes = c(6,7), sizemode = 'area'),
                    hoverinfo = "text",
                    text = ~paste0("<b>",Amino.acids," ",CKpos,"</b><br><br>",
                                   "Group: ",detonly.group,
                                   "<br>Blosum62: ", BLOSUM62, 
                                   "<br>ClinVar: ", CV.AlleleID,
                                   "<br>ClinVar Sig.: ", CLIN.SIG, 
                                   "<br>COSMIC ID: ", COSMIC.GENOMIC.MUTATION.ID, 
                                   "<br>GnomAD AF: ", gnomAD.AF,
                                   "<br>Mutation Prediction: ", MutPred.Top5features
                    ),showlegend=FALSE),
      
      plot_ly(Kdf, x = ~as.numeric(detonly.group), y=yvar, color = boxvar, 
               type = "box", boxpoints = FALSE, hoverinfo = 'y', alpha = 0.3, 
              colors=c("#e31a1c", "#225ea8", # ck
                       "#dd3497","#ae017e","#7a0177", # high
                       "#efedf5","#bcbddc","#756bb1", # med
                       "#deebf7","#9ecae1","#3182bd", # low
                       "#02818a", "#016c59", # tar pan
                       "#ece2f0"),
              legendgroup="g2") %>%     
        add_markers(alpha = 0.8, 
                    x = ~jitter(as.numeric(detonly.group)), 
                    color = ~as.factor(Amino.acids), 
                    size = sizevar,
                    stroke = I("black"), span = I(1),
                    marker = list(sizes = c(6,7), sizemode = 'area'),
                    hoverinfo = "text",
                    text = ~paste0("<b>",Amino.acids," ",CKpos,"</b><br><br>",
                                   "Group: ",detonly.group,
                                   "<br>Blosum62: ", BLOSUM62, 
                                   "<br>ClinVar: ", CV.AlleleID,
                                   "<br>ClinVar Sig.: ", CLIN.SIG, 
                                   "<br>COSMIC ID: ", COSMIC.GENOMIC.MUTATION.ID, 
                                   "<br>GnomAD AF: ", gnomAD.AF,
                                   "<br>Mutation Prediction: ", MutPred.Top5features
                    ),showlegend=FALSE),
      
      shareY = TRUE, widths = c(0.4, 0.4), margin = 0, titleX=TRUE
    )
    
    TITLEplot = paste("Distribution of missense-level", SCOREname, "\nscores for", PROTEINname, "CpDAA's")
    yaxis <- list(title = SCOREname, showgrid = FALSE, zeroline = FALSE, showline = FALSE, showticklabels = TRUE)
    margin <- list(autoexpand = TRUE, t = 110, b=110)
    CKsubs <- CKsubs %>%
      layout(title = TITLEplot, autosize = TRUE,
             yaxis = yaxis,
             margin=margin,
             legend = list(orientation = 'h', y = -0.3, x=0, 
                           font = list(family = "sans-serif", size = 12, color = "#000")), 
             xaxis = list(title = "Detected Cysteine", showticklabels = FALSE), 
             xaxis2 = list(title = "Detected Lysine", showticklabels = FALSE))
    return(CKsubs)
  }
  
  ####### CADD PLOTLY 2 equation 4
  output$caddPlotly2 <- renderPlotly({
    if(input$caddTypeInput2 == "Raw"){
      if(input$assemblyCADD2 == "GRCh37"){
        cadd <- makeBoxplotsMissense(popdf(), detonly.group, 
                                     CADD.raw.hg19, detected.group, Amino.acids, 
                                     variation.size, input$selectGene2, "CADD37 RAW")
      }
      if(input$assemblyCADD2 == "GRCh38"){
        cadd <- makeBoxplotsMissense(popdf(), detonly.group, 
                                     CADD.raw.hg38, detected.group, Amino.acids, 
                                     variation.size, input$selectGene2, "CADD38 RAW")
      }
    } # end raw
    if(input$caddTypeInput2 == "PHRED"){
      if(input$assemblyCADD2 == "GRCh37"){
        cadd <- makeBoxplotsMissense(popdf(), detonly.group, 
                                     CADD.phred.hg19, detected.group, Amino.acids, 
                                     variation.size, input$selectGene2, "CADD37 PHRED")
      }
      if(input$assemblyCADD2 == "GRCh38"){
        cadd <- makeBoxplotsMissense(popdf(), detonly.group, 
                                     CADD.phred.hg38, detected.group, Amino.acids, 
                                     variation.size, input$selectGene2, "CADD38 PHRED")
      }
    } # end phred
    cadd
  })
  
  ###### FATHMM PLOTLY 2
  output$fathmmPlotly2 <- renderPlotly({
    fathmm <- makeBoxplotsMissense(popdf(), detonly.group, fathmm.MKL.coding.score, 
                                   detected.group, Amino.acids, variation.size, 
                                   input$selectGene2, "FATHMM-mkl")
    fathmm
  })
  
  ######## DANN PLOTLY 2
  output$dannPlotly2 <- renderPlotly({
    dann <- makeBoxplotsMissense(popdf(), detonly.group, 
                                 DANN.score, 
                                 detected.group, Amino.acids, variation.size, 
                                 input$selectGene2, "DANN")
    dann
  })
  
  
  
  
  
  
  
  
  
  
  
  
  ####### CADD PLOTLY 1
  output$caddPlotly1 <- renderPlotly({
    if(input$caddTypeInput == "Raw"){
      if(input$lineTypeCADD == "Combo"){
        cadd <- makeCaddLines(proaadf(), matched.aapos, CADD.raw.hg38.max, 
                              CADD.raw.hg38.mean, CADD.raw.hg19.max, 
                              CADD.raw.hg19.mean, CKpos, AAgroup,
                              input$selectGene1, "CADD RAW")
        
      }
      if(input$lineTypeCADD == "CK-specific"){
        cadd <- makeLinePlotsCK(proaadf(), matched.aapos, CADD.raw.hg38.mean, 
                                reactivity, aaref, detected.group, AAgroup,  
                                input$selectGene1, "CADD38 RAW")
      }
    } # end raw
    if(input$caddTypeInput == "PHRED"){
      if(input$lineTypeCADD == "Combo"){
        cadd <- makeCaddLines(proaadf(), matched.aapos, CADD.phred.hg38.max, 
                              CADD.phred.hg38.mean, CADD.phred.hg19.max, 
                              CADD.phred.hg19.mean, CKpos, AAgroup,
                              input$selectGene1, "CADD PHRED")
      }
      if(input$lineTypeCADD == "CK-specific"){
        cadd <- makeLinePlotsCK(proaadf(), matched.aapos, CADD.phred.hg38.mean, 
                                  reactivity, aaref, detected.group, AAgroup,  
                                  input$selectGene1, "CADD38 PHRED")
      }
    } # end phred
    cadd
  })
  
  ###### FATHMM PLOTLY 1
  output$fathmmPlotly1 <- renderPlotly({
    if(input$lineTypeFathmm == "Combo"){
      fathmm <- makeLinePlots1assembly(proaadf(), matched.aapos, fathmmMKL.coding.score.max, 
                                       fathmmMKL.coding.score.mean, CKpos, AAgroup,
                                       input$selectGene1, "FATHMM-mkl")
    }
    if(input$lineTypeFathmm == "CK-specific"){
      fathmm <- makeLinePlotsCK(proaadf(), matched.aapos, fathmmMKL.coding.score.mean, 
                                reactivity, aaref, detected.group, AAgroup,  
                                input$selectGene1, "FATHMM-mkl")
    }
    fathmm
  })

  ######## DANN PLOTLY 1
  output$dannPlotly1 <- renderPlotly({
    if(input$lineTypeDann == "Combo"){
      dann <- makeLinePlots1assembly(proaadf(), matched.aapos, 
                                     DANN.score.max, DANN.score.mean, 
                                     CKpos, AAgroup,input$selectGene1, "DANN")
    }
    if(input$lineTypeDann == "CK-specific"){
      dann <- makeLinePlotsCK(proaadf(), matched.aapos, DANN.score.mean, 
                         reactivity, aaref, detected.group, AAgroup,  
                         input$selectGene1, "DANN")
    }
    dann
  })

#######################################


})
### END shiny server ###
