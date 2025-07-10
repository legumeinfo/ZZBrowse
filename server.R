source("gChart.R")
source("gChartMacro.R")
source("pChart.R")
source("pChartMacro.R")
source("zChart.R")
source("gcvMicroservices.R")

library(promises)
library(future)
plan(multisession)

library(DT)

shinyServer(function(input, output, session) {
  #Load any saved datasets
  values <- reactiveValues(urlStatus = logical(nrow(gwas.sources)))
  dataPath <- "./www/config/data/"
  dataFiles <- list.files(dataPath)
  # ignore subdirectories - we can hide unused datasets there
  dataFiles <- dataFiles[sapply(dataFiles, function(f) !dir.exists(paste0(dataPath, f)))]
  # Append names of any data we will create on the fly:
  # Add your organism to legumeInfo.gwas if its GWAS files live on a server instead of locally.
  legumeInfo.gwas <- c("Arabidopsis thaliana GWAS", "Medicago truncatula GWAS")
  # Add your organism to lis.datastore.gwas if its GWAS files live in the LIS data store.
  lis.datastore.gwas <- c("Common Bean GWAS", "Cowpea GWAS", "Mung bean GWAS", "Peanut GWAS", "Soybean GWAS")
  # unique() removes duplicates from lis.datastore.gwas, if already cached
  lis.datastore.qtl <- c("Common Bean QTL", "Cowpea QTL", "Peanut QTL", "Soybean QTL")
  dataFiles <- unique(c(dataFiles, legumeInfo.gwas, lis.datastore.gwas, lis.datastore.qtl))
  for(i in dataFiles){
    dataFile.i <- paste0(dataPath, i)
    if (i %in% legumeInfo.gwas) {
      df.data <- init.gwas(i)
    } else if (i %in% lis.datastore.gwas) {
      if (file.exists(dataFile.i)) {
        df.data <- read.csv(dataFile.i, header = TRUE, stringsAsFactors = FALSE)
      } else {
        # assemble and cache it
        df.data <- build.gwas.from.lis.datastore(i)
        write.csv(df.data, file = dataFile.i, row.names = FALSE)
      }
    } else if (i %in% lis.datastore.qtl) {
      if (file.exists(dataFile.i)) {
        df.data <- read.csv(dataFile.i, header = TRUE, stringsAsFactors = FALSE)
      } else {
        # assemble and cache it
        df.data <- build.qtl.from.lis.datastore(i)
        write.csv(df.data, file = dataFile.i, row.names = FALSE)
      }
    } else {
      df.data <- read.table(dataFile.i,sep=",",stringsAsFactors=FALSE,head=TRUE)
    }
    # Force the chromosome column to be a string
    # (Note: there must be exactly one column whose name starts with "chr", case-insensitive)
    cols <- tolower(names(df.data))
    k <- which(startsWith(cols, "chr"))
    df.data[, k] <- as.character(df.data[, k])

    values[[i]] <- df.data
  }  
  values$datasetlist <- dataFiles
  values$datasetToOrganism <- NULL # map each dataset to an organism
  # TODO: Do we really need legumeInfo.organisms?
  legumeInfo.organisms <- c("Arabidopsis thaliana", "Medicago truncatula", "Soybean", "Cowpea", "Pigeonpea")

  # Call the GCV chromosome service to assign gene families to annotations (on launch)
  for (org in organismNames) {
    # skip if gene families already exist (on reload, or with values specified in URL)
    if ("family" %in% names(org.annotGeneLoc[[org]])) next
    # initialize gene family to "" to avoid NULL values
    org.annotGeneLoc[[org]][, "family"] <<- ""

    for (chr in chrName[[org]]) {
      gcvChr <- sprintf(org.annotChrFormat[[org]], trailingInteger(chr))
      results <- chromosomeMicroservice(org.gcvUrlBase[[org]], gcvChr)
      if (!is.null(results$error)) {
        print(results$error)
      } else {
        gg <- sapply(unlist(results$chromosome$genes), geneNameFromId)
        # workaround for A. thaliana: convert gene names like "arath.Col.AT1G28130" back to "AT1G28130"
        if (org == "Arabidopsis thaliana") {
          gg <- stri_match_first(gg, regex = "^arath.Col.(.+)$")[, 2]
        }
        ff <- unlist(results$chromosome$families)
        bb <- gg %in% org.annotGeneLoc[[org]]$name
        org.annotGeneLoc[[org]][match(gg[bb], org.annotGeneLoc[[org]]$name), "family"] <<- ff[bb]
      }
    }
  }

  # Extract initial values specified in the URL
  isolate({
    values$urlFields <- parseQueryString(session$clientData$url_search)
    # set GCV client URL if necessary
    if (is.null(userConfig$gcv_client_url)) {
      colon_port <- ifelse(is.null(session$clientData$url_port), "", paste0(":", session$clientData$url_port))
      userConfig$gcv_client_url <- sprintf("%s//%s%s/gcv2", session$clientData$url_protocol, session$clientData$url_hostname, colon_port)
    }
  })

  # what to do when the user changes the jth dataset selection
  datasetChanged <- function(j) {
    if (is.null(input[[jth_ref("datasets", j)]])) return()
    nid <- jth_ref("datasets", j)
    showNotification(paste0("Loading dataset ", input[[jth_ref("datasets", j)]], ". Please wait."),
      duration = NULL, id = nid, type = "message")

    # Make sure the organism corresponds to the selected dataset
    # (invoke values$organism first to force a single reaction)
    if (is.null(values[[jth_ref("organism", j)]])) {
      values[[jth_ref("organism", j)]] <- ifelse(j == 2, "Arabidopsis thaliana", "Medicago truncatula")
    } else {
      isolate({
        values[[jth_ref("organism", j)]] <- values$datasetToOrganism[[input[[jth_ref("datasets", j)]]]]
      })
    }
    if (j == 1) {
      chrNames1 <- "All"
      if (!is.null(values$organism)) chrNames1 <- c(chrNames1, chrName[[values$organism]])
      updateSelectInput(session, "macroChromosome", choices = chrNames1)
    }
    shinyjs::click("submitMacroSynteny")

    # Reset GWAS/QTL traits
    updateTextInput(session, inputId = jth_ref("traitFilter", j), value = "")
    trait.choices <- gwas.traits[[values[[jth_ref("organism", j)]]]]
    if (is.null(trait.choices)) trait.choices <- character(0)
    # Initialize remote GWAS trait files specified in the URL (if any)
    gwj <- jth_ref("gwasTraits", j)
    selectedGwasTraits <- isolate(values$urlFields[[gwj]])
    if (!is.null(selectedGwasTraits)) {
      selectedGwasTraits <- unlist(strsplit(selectedGwasTraits, split = ";"))
    }
    # This update triggers a call to remoteTraitsSelected(j) (below), which loads the remote data
    updateSelectizeInput(session, gwj, choices = trait.choices, selected = selectedGwasTraits)
    # Clear loaded remote GWAS traits if the jth organism changes
    values[[gwj]] <- NULL
    # Set 'Trait names from ontology' checkbox value if specified in the URL
    tfoj <- jth_ref("traitsFromOntology", j)
    tfo <- isolate(values$urlFields[[tfoj]])
    updateCheckboxInput(session, tfoj, value = is.null(tfo) || !(tfo %in% c("n", "N")))

    # Check the Append to Current Dataset checkbox and clear any previously selected GWAS files.
    # As there is no updateFileInput() method, send a custom message (to clear the progress bar) and clear values$needsToUploadFiles.
    session$sendCustomMessage(type = "resetFileInputHandler", jth_ref("uploadfile", j))
    values[[jth_ref("needsToUploadFiles", j)]] <- FALSE
    updateCheckboxInput(session, jth_ref("appendSNPs", j), value = TRUE)

    # Clear all genomic linkages if either organism changes
    values$glSelectedGene <- NULL
    clearGenomicLinkages()

    # Clear macrosynteny bands
    clearMacrosyntenyBands()

    removeNotification(nid)
  }
  # This should be the first code block to detect a change in input$datasets
  observe(datasetChanged(1))
  # This should be the first code block to detect a change in input$datasets2
  observe(datasetChanged(2))

  # Load remote GWAS trait files specified in the URL (if any)
  remoteTraitsSelected <- function(j) {
    gwj <- jth_ref("gwasTraits", j)
    if (is.null(values$urlFields[[gwj]])) return()

    selectedGwasTraits <- isolate(input[[gwj]])
    organism <- jth_ref("organism", j)
    for (trait in selectedGwasTraits) {
      trait.url <- gwas.filenames[[values[[organism]]]][which(gwas.traits[[values[[organism]]]] == trait)]
      loadRemoteData(trait, trait.url, j)
    }
    values$urlFields[[gwj]] <- NULL # to reset it
  }
  observeEvent(input$gwasTraits, remoteTraitsSelected(1))
  observeEvent(input$gwasTraits2, remoteTraitsSelected(2))

  # The user selected one or more local GWAS files (not yet loaded)
  gwasFilesSelected <- function(j) {
    input[[jth_ref("uploadfile", j)]]
    values[[jth_ref("needsToUploadFiles", j)]] <- TRUE
  }
  observe(gwasFilesSelected(1))
  observe(gwasFilesSelected(2))

  #handles what displays in the sidebar based on what tab is selected
  manageTabSelected <- "input.datatabs == 'Manage'"
  dataTableTabSelected <- "input.datatabs == 'Table'"
  wholeGenomeTabSelected <- "input.datatabs == 'WhGen'"
  #the second part of the statement [chromosomeTabSelected] is what is allowing the detection of changing panels due to a click event, i [original ZBrowse author] couldn't figure out how to update input$datatabs with javascript
  chromosomeTabSelected <- "input.datatabs == 'Chrom' || $('li.active a').first().html()==='Chromosome View'"
  annotationsTabSelected <- "input.datatabs == 'Annot'"
  wholeGenomeOrChromosomeTabSelected <- paste(wholeGenomeTabSelected, "||", chromosomeTabSelected)
  annotationsOrChromosomeTabSelected <- paste(annotationsTabSelected, "||", chromosomeTabSelected)

  # for one- or two-species mode
  comparing2Species <- "input.compare2species"
  andComparing2Species <- function(condition) {
    paste0("(", condition, ") && ", comparing2Species)
  }

  createManageSidebar <- function(j) {
    list(
      tags$div(id = jth_ref("tour-dataset", j), wellPanel(
        style = paste0("background-color: ", bgColors[j], ";"),
        uiOutput(jth_ref("datasets", j))
      )),
      tags$div(id = jth_ref("tour-loaddata", j), wellPanel(
        style = paste0("background-color: ", bgColors[j], ";"),
        radioButtons(inputId = jth_ref("dataType", j), label = "Load data (Max. 5MB):", c(".csv" = "csv", ".rda" = "rda", "examples" = "examples"), selected = "csv"),
        conditionalPanel(condition = paste0("input.", jth_ref("dataType", j), " != 'examples'"),
          conditionalPanel(condition = paste0("input.", jth_ref("dataType", j), " == 'csv'"),
            checkboxInput(jth_ref('header', j), 'Header', TRUE),
            radioButtons(jth_ref('sep', j), '', c(Comma=',', Semicolon=';', Tab='\t'), ',')
          ),
          selectizeInput(jth_ref("gwasTraits", j), "Remote Trait Files:", choices = NULL, multiple = TRUE),
          fileInput(jth_ref("uploadfile", j), "Local Trait Files:", multiple = TRUE),
          checkboxInput(jth_ref("appendSNPs", j), "Append to current dataset (combine multiple traits)", TRUE),
          actionButton(jth_ref("loadTraits", j), "Load Data")
        ),
        conditionalPanel(condition = paste0("input.", jth_ref("dataType", j), " == 'examples'"),
          actionButton(jth_ref('loadExampleData', j), 'Load examples')
        )
      )) #,
      # wellPanel(
      #   style = paste0("background-color: ", bgColors[j], ";"),
      #   h6("Once your file is finished uploading, press the Save Dataset button below and reload ZBrowse."),
      #   actionButton(jth_ref('saveDatasetButton', j), 'Save Current Dataset'),
      #   conditionalPanel(condition = paste0("input.", jth_ref("saveDatasetButton", j), " > 0"),
      #     h5("Dataset successfully saved!")
      #   )
      # )
    )
  }
  observeEvent(input$appendSNPs, shinyjs::disable("appendSNPs"))
  observeEvent(input$appendSNPs2, shinyjs::disable("appendSNPs2"))
  
  createDataTableSidebar <- function(j) {
    wellPanel(
      style = paste0("background-color: ", bgColors[j], ";"),
      uiOutput(jth_ref("columns", j)),
      tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(), #add some space between selection columns and subset search
      # uiOutput("view_order"), checkboxInput("view_order_desc", "DESC", value = FALSE),
      returnTextInput(jth_ref("dv_select", j), "Subset (e.g., RMIP > 20 & Location == 'FL06')", '')
    )
  }
  
  createTraitsSidebar <- function(j) {
    wellPanel(
      style = paste0("background-color: ", bgColors[j], ";"),
      checkboxInput(jth_ref("traitsFromOntology", j), label = "Trait names from ontology", value = TRUE),
      uiOutput(jth_ref("traitFilter", j)),
      uiOutput(jth_ref("traitColBoxes", j)),
      uiOutput(jth_ref("legend", j)),
      uiOutput(jth_ref("overlaps", j)),
      conditionalPanel(condition = paste0("input.", jth_ref("overlaps", j), "==true"),
        uiOutput(jth_ref("overlapSize", j)),
        uiOutput(jth_ref("numOverlapTraits", j))
      )#end conditional for plotting only overlaps
      #submitButton("Update View"),
    )
  }
  createChromosomeSidebar <- function(j) {
    wellPanel(
      style = paste0("background-color: ", bgColors[j], ";"),
      uiOutput(jth_ref("selectChr", j))
    )
  }
  createAnnotationWindowSidebar <- function(j) {
    wellPanel(
      style = paste0("background-color: ", bgColors[j], ";"),
      h5("Annotation window options:"),
      h6("Click a point or type a basepair value:"),
      uiOutput(jth_ref("selectedOut", j)),
      uiOutput(jth_ref("windowOut", j))
    )
  }
  createDownloadAnnotationsSidebar <- function(j) {
    wellPanel(
      style = paste0("background-color: ", bgColors[j], ";"),
      helpText(h5(p(paste("Download a CSV of the annotations in the selected window.")))),
      downloadButton(jth_ref('downloadAnnot', j), 'Download')
    )
  }
  createMacrosyntenyPanel <- function() {
    val.m <- isolate(values$urlFields$macroMatched)
    if (is.null(val.m)) val.m <- userConfig$macroMatched
    val.i <- isolate(values$urlFields$macroIntermediate)
    if (is.null(val.i)) val.i <- userConfig$macroIntermediate
    val.mask <- isolate(values$urlFields$macroMask)
    if (is.null(val.mask)) val.mask <- userConfig$macroMask
    val.dd <- isolate(values$urlFields$macroDistance)
    if (is.null(val.dd)) {
      val.dd <- unlist(strsplit(userConfig$macroDistance, split = ":"))
    } else {
      val.dd <- unlist(strsplit(val.dd, split = ":"))
    }
    val.dist <- val.dd[1]
    if (val.dist == "jaccard" && length(val.dd) == 3) {
      val.ngram <- as.integer(val.dd[2])
      val.rev <- as.logical(val.dd[3])
    } else {
      val.ngram <- 1
      val.rev <- FALSE
    }
    tags$div(id = "tour-macrosynteny", wellPanel(
      h5("Macro-synteny options"),
      numericInput("macroMatched", "Matched:", min = 1, max = 100, value = val.m),
      numericInput("macroIntermediate", "Intermediate:", min = 0, max = 100, value = val.i),
      numericInput("macroMask", "Mask:", min = 1, max = 100, value = val.mask),
      selectInput("macroDistance", "Distance metric:", choices = c("Jaccard", "Levenshtein"), selected = stri_trans_totitle(val.dist)),
      conditionalPanel("input.macroDistance == 'Jaccard'",
        numericInput("macroNgram", "n-gram size:", min = 1, max = 2, value = val.ngram),
        conditionalPanel("input.macroNgram > 1",
          checkboxInput('macroReversals', 'Reversals?', value = val.rev)
        )
      ),
      actionButton("submitMacroSynteny", "Update"),
      conditionalPanel(condition = wholeGenomeTabSelected,
        selectInput("macroChromosome", "Species 1 Chromosome:", choices = "All")
      ),
      conditionalPanel(condition = wholeGenomeOrChromosomeTabSelected,
        actionLink("clearMacrosynteny", "Clear Macrosynteny Bands")
      )
      # style = paste0("background-color: ", bgColors[1], ";")
    ))
  }
  createGenomicLinkageSidebar <- function() {
    val.bcName <- isolate(values$urlFields$bcName)
    if (is.null(val.bcName)) val.bcName <- userConfig$bcName
    val.n <- isolate(values$urlFields$neighbors)
    if (is.null(val.n)) val.n <- userConfig$neighbors
    val.m <- isolate(values$urlFields$matched)
    if (is.null(val.m)) val.m <- userConfig$matched
    val.i <- isolate(values$urlFields$intermediate)
    if (is.null(val.i)) val.i <- userConfig$intermediate
    tags$div(id = "tour-genLink", wellPanel(
      h5("Genomic Linkage options"),
      textInput("bcName", "Broadcast Channel name:", value = val.bcName),
      checkboxInput('boolBroadcastToBC', 'Broadcast', TRUE),
      checkboxInput('boolListenToBC', 'Listen', TRUE),
      tags$div(id = "tour-microsynteny", wellPanel(
        uiOutput("selectedGene"),
        numericInput("neighbors", "Neighbors:", min = 1, max = 100, value = val.n),
        textInput("matched", "Matched:", value = val.m),
        numericInput("intermediate", "Intermediate:", min = 0, max = 100, value = val.i),
        conditionalPanel("input.boolBroadcastToBC == true",
          actionLink("viewInGCV", "View in GCV")
        ),
        actionButton("clearGenomicLinkage", "Clear"),
        style = paste0("background-color: ", bgColors[1], ";")
      )),
      tags$div(id = "tour-genLink-2", wellPanel(
        uiOutput("relatedRegions"),
        style = paste0("background-color: ", bgColors[2], ";")
      ))
    ))
  }

  observe({
    # test connections every minute
    invalidateLater(60000, session)
    future({
      sapply(gwas.sources$url, url.exists)
    }) %...>% (function(result) {
      values$urlStatus <- result
    })
  })

  displayConnectionStatus <- function() {
    gwas.sources$status <- values$urlStatus
    apply(gwas.sources, FUN = function(ci) {
      # id for each connection's label
      nn <- paste("gwasSource", ci["name"], sep = "_")
      # server part comes first
      output[[nn]] <- renderText({
        if (ci["status"]) {
          cc <- "green"
          ss <- ""
        } else {
          cc <- "red"
          ss <- " unavailable"
        }
        paste0("<span style=\"color:", cc, "\">", ci["name"], ss, "</span>")
      })
      # then the ui part
      htmlOutput(nn)
    }, MARGIN = 1)
  }
  output$ui_All <- renderUI({
    list(
      conditionalPanel(condition = manageTabSelected,
        actionLink("zbrowseTour", "Start Tour"),
        wellPanel(h5("Connection Status"), displayConnectionStatus()),
        checkboxInput("compare2species", "Two Species",
          value = isolate(is.null(values$urlFields$twoSpecies) || !(values$urlFields$twoSpecies %in% c("n", "N")))),
        conditionalPanel(condition = comparing2Species,
          actionButton("swap", "\u21f5 Swap Datasets")
        ),
        createManageSidebar(1),
        conditionalPanel(condition = comparing2Species,
          createManageSidebar(2)
        ),
        tags$script("Shiny.addCustomMessageHandler('resetFileInputHandler', function(x) {
          var id = '#' + x;
          var idProgress = id + '_progress';
          var idBar = id + ' .bar';
          $(idProgress).css('visibility', 'hidden');
          $(idBar).css('width', '0%');
          // TODO: clear the loaded file(s) (none of the following methods work, for Javascript security reasons)
          /*
          $(id).val(''); // method 1
          $(id).replaceWith($(id).val('').clone(true)); // method 2
          $(id).replaceWith($(id) = $(id).clone(true)); // method 3
          // method 4
          $(id).wrap('<form></form>').closest('form').reset();
          $(id).unwrap();
          $(id).stopPropagation();
          $(id).preventDefault();
          */
        });"),
        helpModal('Manage','manage', includeMarkdown("tools/manage.md")),
        HTML(paste('<p style="font-size:10px;">Powered by',
          '<a href="http://www.rstudio.com/shiny/" target=_blank>Shiny</a>,',
          '<a href="https://ramnathv.github.io/rCharts/" target=_blank>rCharts</a>,',
          '<a href="http://www.highcharts.com" target=_blank>Highcharts</a>,',
          'and <a href="https://github.com/carlganz/rintrojs" target=_blank>rintrojs</a>',
          '</p>'))
      ),#end conditional Manage

      conditionalPanel(condition = dataTableTabSelected,
        tags$div(id = "tour-datatableSidebar", createDataTableSidebar(1)),
        conditionalPanel(condition = comparing2Species,
          createDataTableSidebar(2)
        ),
        helpModal('Data Table View','view',includeMarkdown("tools/manage.md"))      
      ),#end conditional Table

      # Sidebar panels for the remaining tabs, in conditional order
      conditionalPanel(condition = wholeGenomeOrChromosomeTabSelected,
        helpText(h5(p("Interactive Graphs for GWAS Data"))),
        createTraitsSidebar(1)
      ),
      conditionalPanel(condition = chromosomeTabSelected,
        createChromosomeSidebar(1)
      ),
      conditionalPanel(condition = annotationsOrChromosomeTabSelected,
        createAnnotationWindowSidebar(1)
      ),
      conditionalPanel(condition = annotationsTabSelected,
        createDownloadAnnotationsSidebar(1)
      ),
      conditionalPanel(condition = andComparing2Species(chromosomeTabSelected),
        createGenomicLinkageSidebar() # goes between the sidebars for organism 1 (above) and organism 2 (below)
      ),
      conditionalPanel(condition = andComparing2Species(wholeGenomeOrChromosomeTabSelected),
        createTraitsSidebar(2)
      ),
      conditionalPanel(condition = andComparing2Species(chromosomeTabSelected),
        createChromosomeSidebar(2)
      ),
      conditionalPanel(condition = andComparing2Species(annotationsOrChromosomeTabSelected),
        createAnnotationWindowSidebar(2)
      ),
      conditionalPanel(condition = andComparing2Species(annotationsTabSelected),
        createDownloadAnnotationsSidebar(2)
      ),
      conditionalPanel(condition = andComparing2Species(wholeGenomeOrChromosomeTabSelected),
        createMacrosyntenyPanel()
      ),
      # end sidebar panels for remaining tabs

      conditionalPanel(condition = wholeGenomeOrChromosomeTabSelected,
        helpModal('Browser Help','browser',includeMarkdown("tools/manage.md"))        
      )#add help button for browser tabs
    )#end list
  }) #end ui_All
  outputOptions(output, "ui_All", suspendWhenHidden=FALSE)

  # find the appropriate UI
  output$ui_finder <- renderUI({
#    if(is.null(input$datatabs)){      
#      get("ui_All")()
#    }else{
#      get(paste0('ui_',input$datatabs))()
#    }
#     if(input$tool == "data") {
#       if(!is.null(input$datatabs)) get(paste0('ui_',input$datatabs))()
#     } else {
#       if(!is.null(input$tool)) get(paste0('ui_',input$tool))()
#     }
  })  
  outputOptions(output, "ui_finder", suspendWhenHidden=FALSE)

  observeEvent(input$swap, {
    dataset1 <- input$datasets
    updateSelectInput(session, "datasets", selected = input$datasets2)
    updateSelectInput(session, "datasets2", selected = dataset1)
  })

  # create the menu for selecting the jth dataset
  renderDatasets <- function(j, defaultDataset) {
    if (is.null(input[[jth_ref("loadTraits", j)]])) return()
    
    # Load any requested GWAS files
    isolate({
      # Local GWAS files
      inFile <- NULL
      if (values[[jth_ref("needsToUploadFiles", j)]]) inFile <- input[[jth_ref("uploadfile", j)]]
      if (!is.null(inFile)) {
        # iterating through the files to upload
        for (i in 1:(dim(inFile)[1])) {
          loadUserData(inFile[i, 'name'], inFile[i, 'datapath'], j)
        }
        values[[jth_ref("needsToUploadFiles", j)]] <- FALSE # since we just loaded them
      }
      # Remote GWAS files
      inTraits <- input[[jth_ref("gwasTraits", j)]]
      if (!is.null(inTraits)) {
        for (trait in inTraits) {
          trait.url <- gwas.filenames[[values[[jth_ref("organism", j)]]]][which(gwas.traits[[values[[jth_ref("organism", j)]]]] == trait)]
          loadRemoteData(trait, trait.url, j)
        }
      }
    })
    
    dat <- isolate(input[[jth_ref("datasets", j)]])
    if (!is.null(dat)) {
      appendSNPs <- isolate(input[[jth_ref("appendSNPs", j)]])
      if (appendSNPs || (is.null(inFile) && is.null(inTraits))) {
        val <- dat
      } else {
        val <- values$datasetlist[1]
      }
    }else{
      val <- defaultDataset
    }
    
    # Drop-down selection of data set
    selectInput(inputId = jth_ref("datasets", j), label = paste0("Dataset ", j, ":"), choices = values$datasetlist, selected = values$datasetlist[values$datasetlist == val], multiple = FALSE, selectize = FALSE)
  }
  initialDataset <- isolate(values$urlFields$datasets)
  if (is.null(initialDataset)) initialDataset <- "Medicago truncatula GWAS"
  output$datasets <- renderUI(renderDatasets(1, initialDataset))
  initialDataset2 <- isolate(values$urlFields$datasets2)
  if (is.null(initialDataset2)) initialDataset2 <- "Arabidopsis thaliana GWAS"
  output$datasets2 <- renderUI(renderDatasets(2, initialDataset2))

  createAnnotTable <- function(j) {
    if (is.null(input[[jth_ref("datasets", j)]])) return()
    centerBP <- as.numeric(input[[jth_ref("selected", j)]])
    winHigh <- centerBP + input[[jth_ref("window", j)]]
    winLow <- centerBP - input[[jth_ref("window", j)]]
    if (winLow < 0) { winLow <- 0 }
    thisChrAnnot <- subset(org.annotGeneLoc[[values[[jth_ref("organism", j)]]]], chromosome == input[[jth_ref("chr", j)]])
    thisAnnot <- thisChrAnnot[thisChrAnnot$transcript_start >= winLow & thisChrAnnot$transcript_end <= winHigh, ]
    thisAnnot
  }
  
  #Returns the nicely formatted preview table
  createHtmlDataExample <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])) return()

    dat <- getdata(j)

    # necessary when deleting a dataset
    if(is.null(dat) || nrow(dat) == 0) return()

    # Show only the first 10 rows
    # (and hide the publication column, which is an HTML link)
    nr <- min(10,nrow(dat))
    dat <- data.frame(dat[1:nr, names(dat) != "publication", drop = FALSE])

    #dat <- date2character_dat(dat) #may be needed to print table if there is a data column

    # Now that the data exist, select the tab specified in the URL (if any)
    initialTab <- isolate(values$urlFields$tab)
    if (!is.null(initialTab)) {
      updateTabsetPanel(session, "datatabs", selected = initialTab)
      values$urlFields$tab <- NULL # to reset it
    }

    html <- print(xtable::xtable(dat), type='html', include.rownames = FALSE, print.results = FALSE)
    html <- sub("<TABLE border=1>","<table class='table table-condensed table-hover'>", html)
    Encoding(html) <- 'UTF-8'
    html
  }
  output$htmlDataExample <- reactive(createHtmlDataExample(1))
  output$htmlDataExample2 <- reactive(createHtmlDataExample(2))
  
  createNrowDataset <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])) return()
    dat <- getdata(j)
    if(is.null(dat) || nrow(dat) == 0) return()
    nr <- nrow(dat)
    
    if(nr>2500){
      paste0('<p>First 10 rows shown of ',nr,' rows. See Data Table tab for details.<br>More than 2500 rows found, only the top 2500 will be plotted.</p>')
    }else{
      paste0('<p>First 10 rows shown of ',nr,' rows. See Data Table tab for details.</p>')
    }
  }
  output$nrowDataset <- reactive(createNrowDataset(1))
  output$nrowDataset2 <- reactive(createNrowDataset(2))

  createDatasetSynopsis <- function(j) {
    if (is.null(input[[jth_ref("datasets", j)]])) return()
    dat <- getdata(j)
    if (is.null(dat) || nrow(dat) == 0 || is.null(dat$publication)) return()

    pubs <- unique(dat$publication)
    ss <- sapply(pubs, function(pub) stri_match_first(pub,
      regex = "title='(.+)' target=_blank>(.+)</a>"))
    pp <- sort(sprintf("<p>%s. %s</p>", ss[3, ], ss[2, ]))
    paste(pp, collapse = "")
  }
  output$datasetSynopsis <- reactive(createDatasetSynopsis(1))
  output$datasetSynopsis2 <- reactive(createDatasetSynopsis(2))

  createColumnSettingsPanel <- function(j) {
    tags$div(id = jth_ref("tour-columnSettings", j), wellPanel(
      style = paste0("background-color: ", bgColors[j], ";"),

      htmlOutput(jth_ref("htmlDataExample", j)),
      htmlOutput(jth_ref("nrowDataset", j)),
      htmlOutput(jth_ref("datasetSynopsis", j)),

      #tags$p(tags$br()),
      row(
        col(3, tags$br()),
        col(7, h4('Select appropriate columns to be used for plotting.'))
        #HTML('<h4>Select appropriate columns to be used for plotting.</h4>'),
      ),
      tags$hr(),
      row(
        #col(2, tags$br()),
        col(2,uiOutput(jth_ref("chrColumn", j)),uiOutput(jth_ref("bpColumn", j))),
        col(2,uiOutput(jth_ref("plotAll", j)),uiOutput(jth_ref("traitColumns", j))),
        col(2,uiOutput(jth_ref("yAxisColumn", j)),uiOutput(jth_ref("logP", j))),
        col(2,uiOutput(jth_ref("axisLimBool", j)),uiOutput(jth_ref("axisLim", j))),
        col(2,actionButton(jth_ref("SubmitColsButton", j),"Submit"))
      ),
      tags$hr(),
      row(
        col(7, uiOutput(jth_ref("supportInterval", j)))#
      ),
      row(
        col(2, uiOutput(jth_ref("SIbpStart", j))),
        col(2, uiOutput(jth_ref("SIyAxisColumn", j))),
        col(2, uiOutput(jth_ref("SIaxisLimBool", j)),uiOutput(jth_ref("SIaxisLim", j)))
      )
    ))
  }
  output$ui_data_tabs <- renderUI({
    tabsetPanel(id = "datatabs",      
      tabPanel(title="Manage",value="Manage",
        createColumnSettingsPanel(1),
        conditionalPanel(condition = comparing2Species,
          createColumnSettingsPanel(2)
        )
      ),
      tabPanel(title="Data Table",value="Table",
        tags$div(id = "tour-datatable",
          wellPanel(DTOutput("dataviewer"), style = paste0("background-color: ", bgColors[1], ";"))
        ),
        conditionalPanel(condition = comparing2Species,
          wellPanel(DTOutput("dataviewer2"), style = paste0("background-color: ", bgColors[2], ";"))
        )
      ),
      tabPanel(title="Whole Genome View",value="WhGen",
        conditionalPanel(condition = comparing2Species,
          wellPanel(
            showOutput("gChartMacro", "highcharts"),
            style = paste0("background-color: ", bgColors[1], ";")
          ),
          tags$div(id = "tour-wholegenome-macrosynteny",
            wellPanel(
              showOutput("gChartMacro2", "highcharts"),
              style = paste0("background-color: ", bgColors[2], ";")
            )
          )
        ),
        tags$div(id = "tour-wholegenome",
          wellPanel(showOutput("gChart", "highcharts"), style = paste0("background-color: ", bgColors[1], ";"))
        ),
        conditionalPanel(condition = comparing2Species,
          wellPanel(showOutput("gChart2", "highcharts"), style = paste0("background-color: ", bgColors[2], ";"))
        ),
        HTML("<style>.rChart { width: 100%; }</style>") # to expand charts to full width
      ),
      tabPanel(title="Chromosome View",value="Chrom",
        conditionalPanel(condition = comparing2Species,
          wellPanel(showOutput("pChartMacro", "highcharts"), style = paste0("background-color: ", bgColors[1], ";")),
          wellPanel(showOutput("pChartMacro2", "highcharts"), style = paste0("background-color: ", bgColors[2], ";"))
        ),
        wellPanel(
          tags$div(id = "tour-pChart", showOutput("pChart", "highcharts")),
          tags$div(id = "tour-zChart", showOutput("zChart", "highcharts")),
          tags$script('Shiny.addCustomMessageHandler("customMsg", function(bandOpts){
            hc = $("#pChartMacro").highcharts();
            if (hc != undefined) {
              chartXAxis = hc.xAxis[0];
              chartXAxis.removePlotBand();
              chartXAxis.addPlotBand(bandOpts);
            }
            hc = $("#pChart").highcharts();
            if (hc != undefined) {
              chartXAxis = hc.xAxis[0];
              chartXAxis.removePlotBand();
              chartXAxis.addPlotBand(bandOpts);
            }
          })'),
          style = paste0("background-color: ", bgColors[1], ";")
        ),
        conditionalPanel(condition = comparing2Species,
          wellPanel(showOutput("zChart2", "highcharts"), showOutput("pChart2", "highcharts"),
          tags$script('Shiny.addCustomMessageHandler("customMsg2", function(bandOpts){
            hc = $("#pChartMacro2").highcharts();
            if (hc != undefined) {
              chartXAxis = hc.xAxis[0];
              chartXAxis.removePlotBand();
              chartXAxis.addPlotBand(bandOpts);
            }
            hc = $("#pChart2").highcharts();
            if (hc != undefined) {
              chartXAxis = hc.xAxis[0];
              chartXAxis.removePlotBand();
              chartXAxis.addPlotBand(bandOpts);
            }
          })'),
          style = paste0("background-color: ", bgColors[2], ";")
        )),
        HTML("<style>.rChart { width: 100%; }</style>") # to expand charts to full width
      ),
      tabPanel(title="Annotations Table",value="Annot",
        tags$div(id = "tour-annotations",
          wellPanel(DTOutput("annotViewer"), style = paste0("background-color: ", bgColors[1], ";"))
        ),
        conditionalPanel(condition = comparing2Species,
          wellPanel(DTOutput("annotViewer2"), style = paste0("background-color: ", bgColors[2], ";"))
        )
      )
    )#end tabsetPanel
  })#end data tabs
  outputOptions(output, "ui_data_tabs", suspendWhenHidden=FALSE)

  annotViewer.options <- list(orderClasses = TRUE, bCaseInsensitive = TRUE,
    lengthMenu = c(15, 50, 100, 200, 500), pageLength = 15,
    # make odd rows white
    initComplete = JS(
      "function(settings, json) {",
        "$('table.dataTable.display tbody tr.odd').css('background-color', 'white');",
      "}"
    ),
    "dom" = 'T<"clear">lfrtip',
    "oTableTools" = list(
      "sSwfPath" = "/tabletools/swf/copy_csv_xls_pdf.swf",
      "aButtons" = list(
        "copy",
        "print",
        list("sExtends" = "collection",
          "sButtonText" = "Save",
          "aButtons" = c("csv","xls","pdf")
        )
      )
    )                    
  )
  output$annotViewer <- renderDT({
    createAnnotTable(1)
  }, options = annotViewer.options, escape = FALSE)
  output$annotViewer2 <- renderDT({
    createAnnotTable(2)
  }, options = annotViewer.options, escape = FALSE)

  dataViewer.options <- annotViewer.options # they happen to be identical, but we could make them different
  createDataViewer <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]]) || is.null(input[[jth_ref("columns", j)]])) return()
    
    dat <- getdata(j)
    #dat <- date2character()
    
    if(!all(input[[jth_ref("columns", j)]] %in% colnames(dat))) return()
    
    if(input[[jth_ref("dv_select", j)]] != '') {
      selcom <- input[[jth_ref("dv_select", j)]]
      selcom <- gsub(" ", "", selcom)
      
      seldat <- try(do.call(subset, list(dat,parse(text = selcom))), silent = TRUE)
      
      if(!is(seldat, 'try-error')) {
        if(is.data.frame(seldat)) {
          dat <- seldat
          seldat <- NULL
        }
      }
    }
    
    dat <- data.frame(dat[, input[[jth_ref("columns", j)]], drop = FALSE])
    dat
    
    # html <- print(xtable::xtable(dat), type='html', include.rownames = FALSE, print.results = FALSE)
    # html <- sub("<TABLE border=1>","<table class='table table-condensed table-hover'>", html)
    # html
  }
  output$dataviewer <-renderDT(createDataViewer(1), options = dataViewer.options, escape = FALSE)
  output$dataviewer2 <-renderDT(createDataViewer(2), options = dataViewer.options, escape = FALSE)
  
  createDownloadAnnot <- function(j) {
    downloadHandler(
      filename = function() {paste0("AnnotationsAround.chr",input[[jth_ref("chr", j)]],".",input[[jth_ref("selected", j)]],"bp.",values[[jth_ref("organism", j)]],".csv")},
      content = function(file) {write.csv(createAnnotTable(j),file,row.names=F)}
    )
  }
  output$downloadAnnot <- createDownloadAnnot(1)
  output$downloadAnnot2 <- createDownloadAnnot(2)

  createColumns <- function(j) {
    cols <- varnames(j)    
    selectInput(jth_ref("columns", j), "Select columns to show:", choices = as.list(cols), selected = cols, multiple = TRUE)
  }
  output$columns <- renderUI(createColumns(1))
  output$columns2 <- renderUI(createColumns(2))
  
  createAxisLimBool <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    if(input[[jth_ref("datasets", j)]] %in% datasetProp()$dataset){
      val = datasetProp()$axisLim[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
    } else if (values[[jth_ref("organism", j)]] %in% legumeInfo.organisms) {
      val <- FALSE
    }else{
      val = TRUE}
    checkboxInput(jth_ref('axisLimBool', j), 'Set Y-axis Limits?', val)
  }
  output$axisLimBool <- renderUI(createAxisLimBool(1))
  output$axisLimBool2 <- renderUI(createAxisLimBool(2))

  createLogP <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    if(input[[jth_ref("datasets", j)]] %in% datasetProp()$dataset){
      val = datasetProp()$logP[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
    } else if (values[[jth_ref("organism", j)]] %in% legumeInfo.organisms) {
      val <- TRUE
    }else{
      val = FALSE}
    checkboxInput(jth_ref('logP', j), 'Take -log10 of column?', val)
  }
  output$logP <- renderUI(createLogP(1))
  output$logP2 <- renderUI(createLogP(2))

  createChrColumn <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    cols <- varnames(j)    
    if(input[[jth_ref("datasets", j)]] %in% datasetProp()$dataset){
      selected = datasetProp()$chrColumn[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
    }else{
      selected = names(cols[1])
    }    
    selected <- stri_match(selected, regex = ".*(?=\\ \\{)")[, 1]
    selectizeInput(jth_ref("chrColumn", j), "Chromosome Column:", choices = as.list(cols), selected = selected, multiple = FALSE, options = list(dropdownParent="body"))
  }
  output$chrColumn <- renderUI(createChrColumn(1))
  output$chrColumn2 <- renderUI(createChrColumn(2))

  createBpColumn <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    cols <- varnames(j)
    if(input[[jth_ref("datasets", j)]] %in% datasetProp()$dataset){
      selected = datasetProp()$bpColumn[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
    }else{
      selected = names(cols[2])
    }
    selected <- stri_match(selected, regex = ".*(?=\\ \\{)")[, 1]
    selectizeInput(jth_ref("bpColumn", j), "Base Pair Column:", choices = as.list(cols), selected = selected, multiple = FALSE, options = list(dropdownParent="body"))
  }
  output$bpColumn <- renderUI(createBpColumn(1))
  output$bpColumn2 <- renderUI(createBpColumn(2))

  createTraitColumns <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    cols <- varnames(j)
    if(input[[jth_ref("datasets", j)]] %in% datasetProp()$dataset){
      selected = unlist(strsplit(datasetProp()$traitCol[datasetProp()$dataset == input[[jth_ref("datasets", j)]]],";"))
    } else if (values[[jth_ref("organism", j)]] %in% legumeInfo.organisms) {
      selected <- names(cols[3])
    }else{
      selected = names(cols[3:4])
    }
    selected <- sapply(selected, FUN = function(s) stri_match(selected, regex = ".*(?=\\ \\{)")[, 1], USE.NAMES = FALSE)
    conditionalPanel(condition = paste0("input.", jth_ref("plotAll", j), "==false"),
      selectizeInput(jth_ref("traitColumns", j), "Group by these trait column(s):", choices = as.list(cols), selected = selected, multiple = TRUE, options = list(dropdownParent="body"))
    )        
  }
  output$traitColumns <- renderUI(createTraitColumns(1))
  output$traitColumns2 <- renderUI(createTraitColumns(2))

  createPlotAll <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    if(input[[jth_ref("datasets", j)]] %in% datasetProp()$dataset){      
      val = datasetProp()$plotAll[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
    }else{
      val = FALSE
    }    
    checkboxInput(jth_ref('plotAll', j), 'All data is the same trait', val)
  }
  output$plotAll <- renderUI(createPlotAll(1))
  output$plotAll2 <- renderUI(createPlotAll(2))

  createYAxisColumn <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    cols <- varnames(j)
    if(input[[jth_ref("datasets", j)]] %in% datasetProp()$dataset){      
      selected = datasetProp()$yAxisColumn[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
    } else if (values[[jth_ref("organism", j)]] %in% legumeInfo.organisms) {
      selected <- names(cols[4])
    }else{
      #selected = names(cols[10])
      selected = as.character(cols[10])
    }    
    selected <- stri_match(selected, regex = ".*(?=\\ \\{)")[, 1]
    selectizeInput(jth_ref("yAxisColumn", j), "Y-axis column:", choices = as.list(cols), selected = selected, multiple = FALSE, options = list(dropdownParent="body"))
  }
  output$yAxisColumn <- renderUI(createYAxisColumn(1))
  outputOptions(output, "yAxisColumn", suspendWhenHidden=FALSE)
  output$yAxisColumn2 <- renderUI(createYAxisColumn(2))
  outputOptions(output, "yAxisColumn2", suspendWhenHidden=FALSE)
  
  createAxisLim <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    if(input[[jth_ref("datasets", j)]] %in% datasetProp()$dataset){
      min = datasetProp()$axisMin[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
      max = datasetProp()$axisMax[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
    }else{
      min = 0
      max = 1
    }
    conditionalPanel(condition = paste0("input.", jth_ref("axisLimBool", j), "==true"),
      numericInput(jth_ref("axisMin", j), "Min:", value=min),
      numericInput(jth_ref("axisMax", j), "Max:", value=max)
    )
  }
  output$axisLim <- renderUI(createAxisLim(1))
  output$axisLim2 <- renderUI(createAxisLim(2))

  createSupportInterval <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    if(input[[jth_ref("datasets", j)]] %in% datasetProp()$dataset){
      val = datasetProp()$supportInterval[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
    }else{
      val = FALSE}
    checkboxInput(jth_ref('supportInterval', j), 'Plot base pair intervals (e.g., Joint linkage support intervals)?', val)
  }
  output$supportInterval <- renderUI(createSupportInterval(1))
  output$supportInterval2 <- renderUI(createSupportInterval(2))

  createSIbpStart <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    cols <- varnames(j)
    if(input[[jth_ref("datasets", j)]] %in% datasetProp()$dataset){
      selected = datasetProp()$SIbpStart[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
      selectedEnd = datasetProp()$SIbpEnd[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
    }else{
      selected = names(cols[2])
      selectedEnd = names(cols[2])
    }
    conditionalPanel(condition = paste0("input.", jth_ref("supportInterval", j), "==true"),
      selectizeInput(jth_ref("SIbpStart", j), "Interval Base Pair Start:", choices = as.list(cols), selected = selected, multiple = FALSE, options = list(dropdownParent="body")),
      selectizeInput(jth_ref("SIbpEnd", j), "Interval Base Pair End:", choices = as.list(cols), selected = selectedEnd, multiple = FALSE, options = list(dropdownParent="body"))
    )
  }
  output$SIbpStart <- renderUI(createSIbpStart(1))
  output$SIbpStart2 <- renderUI(createSIbpStart(2))

  createSIyAxisColumn <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    cols <- varnames(j)
    if(input[[jth_ref("datasets", j)]] %in% datasetProp()$dataset){  
      selected = datasetProp()$SIyAxisColumn[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
    }else{
      #selected = names(cols[10])
      selected = as.character(cols[10])
    }        
    conditionalPanel(condition = paste0("input.", jth_ref("supportInterval", j), "==true"),
      selectizeInput(jth_ref("SIyAxisColumn", j), "Support Interval Y-axis column:", choices = as.list(cols), selected = selected, multiple = FALSE, options = list(dropdownParent="body"))
    )
  }
  output$SIyAxisColumn <- renderUI(createSIyAxisColumn(1))
  outputOptions(output, "SIyAxisColumn", suspendWhenHidden=FALSE)  
  output$SIyAxisColumn2 <- renderUI(createSIyAxisColumn(2))
  outputOptions(output, "SIyAxisColumn2", suspendWhenHidden=FALSE)
  
  createSIaxisLimBool <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    if(input[[jth_ref("datasets", j)]] %in% datasetProp()$dataset){
      val = datasetProp()$SIaxisLimBool[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
    } else if (values[[jth_ref("organism", j)]] %in% legumeInfo.organisms) {
      val <- FALSE
    }else{
      val = TRUE
    }
    conditionalPanel(condition = paste0("input.", jth_ref("supportInterval", j), "==true"),
      checkboxInput(jth_ref('SIaxisLimBool', j), 'Set Support Interval Y-axis Limits?', val)
    )
  }
  output$SIaxisLimBool <- renderUI(createSIaxisLimBool(1))
  output$SIaxisLimBool2 <- renderUI(createSIaxisLimBool(2))

  createSIaxisLim <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    if(input[[jth_ref("datasets", j)]] %in% datasetProp()$dataset){
      min = datasetProp()$axisMin[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
      max = datasetProp()$axisMax[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
    }else{
      min = 0
      max = 1
    }
    conditionalPanel(condition = paste0("input.", jth_ref("supportInterval", j), "==true && input.", jth_ref("SIaxisLimBool", j), "==true"),
      numericInput(jth_ref("SIaxisMin", j), "Min:", value=min),
      numericInput(jth_ref("SIaxisMax", j), "Max:", value=max)
    )
  }
  output$SIaxisLim <- renderUI(createSIaxisLim(1))
  output$SIaxisLim2 <- renderUI(createSIaxisLim(2))

  # Text input that filters traits if their name matches
  createTraitFilter <- function(j) {
    textInput(inputId = jth_ref("traitFilter", j), label = "Filter traits:", placeholder = "e.g. Seed")
  }
  output$traitFilter <- renderUI(createTraitFilter(1))
  output$traitFilter2 <- renderUI(createTraitFilter(2))

  #builds list of multiple selection boxes for traits that have multiple columns in dataset
  createTraitColBoxes <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    if(input[[jth_ref("plotAll", j)]] == TRUE){return()}
    df.data <- values[[input[[jth_ref("datasets", j)]]]]
    hasGwasTraits <- any(is.na(df.data[, input[[jth_ref("SIbpStart", j)]]]))
    hasQtlTraits <- any(!is.na(df.data[, input[[jth_ref("SIbpStart", j)]]]))
    lapply(input[[jth_ref("traitColumns", j)]], function(i) {
      traits <- c("Select All", "Deselect All", sort(unique(values[[input[[jth_ref("datasets", j)]]]][,i])))
      # trait index is 1 for Select All, 2 for Deselect All, 3 for the first trait, etc
      if (hasGwasTraits && hasQtlTraits) {
        traits <- c(traits[1], "Select only GWAS traits", "Select only QTL traits", traits[-1])
      }
      selectedTraits <- traits[1] # default is Select All
      # Select traits specified in the URL (if any)
      traitsFromUrl <- isolate(values$urlFields[[jth_ref("traits", j)]])
      if (!is.null(traitsFromUrl)) {
        selectedTraits <- unlist(strsplit(traitsFromUrl, split = ";"))
      }
      selectizeInput(inputId = jth_ref(i, j), label = paste0("Select ", i),
        choices = traits, selected = selectedTraits, multiple = TRUE,
        options = list(dropdownParent = "body", plugins = list("remove_button")))
    })
  }
  output$traitColBoxes <- renderUI(createTraitColBoxes(1))
  output$traitColBoxes2 <- renderUI(createTraitColBoxes(2))

  updateTraitsMenu <- function(j) {
    lapply(input[[jth_ref("traitColumns", j)]], function(i){
      df.data <- values[[input[[jth_ref("datasets", j)]]]]
      tf <- input[[jth_ref("traitFilter", j)]]
      if ("Select All" %in% input[[jth_ref(i, j)]]) {
        selected_choices <- sort(unique(df.data[, i]))
        selected_choices <- stringsThatMatchPattern(selected_choices, tf)
        updateSelectizeInput(session, jth_ref(i, j), selected = selected_choices)
      } else if ("Select only GWAS traits" %in% input[[jth_ref(i, j)]]) {
        df.data <- df.data[is.na(df.data[, input[[jth_ref("SIbpStart", j)]]]), ]
        selected_choices <- sort(unique(df.data[, i]))
        selected_choices <- stringsThatMatchPattern(selected_choices, tf)
        updateSelectizeInput(session, jth_ref(i, j), selected = selected_choices)
      } else if ("Select only QTL traits" %in% input[[jth_ref(i, j)]]) {
        df.data <- df.data[!is.na(df.data[, input[[jth_ref("SIbpStart", j)]]]), ]
        selected_choices <- sort(unique(df.data[, i]))
        selected_choices <- stringsThatMatchPattern(selected_choices, tf)
        updateSelectizeInput(session, jth_ref(i, j), selected = selected_choices)
      } else if ("Deselect All" %in% input[[jth_ref(i, j)]]) {
        updateSelectizeInput(session, jth_ref(i, j), selected = character(0))
      }
    })
  }
  observe(updateTraitsMenu(1))
  observe(updateTraitsMenu(2))

  #checkbox to suppress plot legend
  createLegend <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    checkboxInput(jth_ref('legend', j), 'Suppress Legend', FALSE)
  }
  output$legend <- renderUI(createLegend(1))
  outputOptions(output, "legend", suspendWhenHidden=FALSE)
  output$legend2 <- renderUI(createLegend(2))
  outputOptions(output, "legend2", suspendWhenHidden=FALSE)
  
  #checkbox for whether to filter for only overlapping SNPs
  createOverlaps <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
#    if(input[[jth_ref("plotAll", j)]] == TRUE){return()}
    checkboxInput(jth_ref('overlaps', j), 'Show only overlapping SNPs', FALSE)
  }
  output$overlaps <- renderUI(createOverlaps(1))
  outputOptions(output, "overlaps", suspendWhenHidden=FALSE)
  output$overlaps2 <- renderUI(createOverlaps(2))
  outputOptions(output, "overlaps2", suspendWhenHidden=FALSE)
  
  #how many traits must overlap to be included in output, 1 means traits that overlap with themselves will be included
  createNumOverlapTraits <- function(j) {
    numericInput(jth_ref("numOverlaps", j), "Minimum number of overlaps?", value=2,min=1)
  }
  output$numOverlapTraits <- renderUI(createNumOverlapTraits(1))
  output$numOverlapTraits2 <- renderUI(createNumOverlapTraits(2))

  #how big is the window when calculating whether two snps overlap
  createOverlapSize <- function(j) {
    numericInput(inputId=jth_ref("overlapSize", j), label="Overlap size around each point:",min=1,max=.5e6,value=10000)
  }
  output$overlapSize <- renderUI(createOverlapSize(1))
  output$overlapSize2 <- renderUI(createOverlapSize(2))

  createSelectChr <- function(j) {
    organism <- values[[jth_ref("organism", j)]]
    if (is.null(organism)) { return() }
    chrChoices <- chrName[[organism]]
    chr <- jth_ref("chr", j)
    # Select the chromosome specified in the URL (if any)
    initialChromosome <- isolate(values$urlFields[[chr]])
    if (is.null(initialChromosome)) { initialChromosome <- NULL }
    selectInput(chr, "Chromosome:", chrChoices, selected = initialChromosome, selectize = FALSE)
  }
  output$selectChr <- renderUI(createSelectChr(1))
  outputOptions(output, "selectChr", suspendWhenHidden=FALSE)
  output$selectChr2 <- renderUI(createSelectChr(2))
  outputOptions(output, "selectChr2", suspendWhenHidden=FALSE)
  
  createSelectedOut <- function(j) {
    selected <- jth_ref("selected", j)
    # Select the center position specified in the URL (if any)
    initialCenter <- isolate(values$urlFields[[selected]])
    if (is.null(initialCenter)) { initialCenter <- defaultCenter }
    numericInput(selected, "", min = 1000, value = initialCenter, step = 1000)
  }
  output$selectedOut <- renderUI(createSelectedOut(1))
  outputOptions(output, "selectedOut", suspendWhenHidden=FALSE)
  output$selectedOut2 <- renderUI(createSelectedOut(2))
  outputOptions(output, "selectedOut2", suspendWhenHidden=FALSE)

  observeEvent(input$chr, {
    chr <- trailingInteger(input$chr)
    chr_length <- chrSize[[values$organism]][chr]
    # prevent exceeding the chromosome length
    updateNumericInput(session, "selected", max = chr_length)
    if (input$selected > chr_length - input$window) {
      updateNumericInput(session, "selected", value = chr_length - input$window)
    }
  })
  observeEvent(input$chr2, {
    chr2 <- trailingInteger(input$chr2)
    chr2_length <- chrSize[[values$organism2]][chr2]
    # prevent exceeding the chromosome length
    updateNumericInput(session, "selected2", max = chr2_length)
    if (input$selected2 > chr2_length - input$window2) {
      updateNumericInput(session, "selected2", value = chr2_length - input$window2)
    }
  })
  observeEvent(input$selected, {
    validNumericInput("selected", defaultCenter, 1000)
  })
  observeEvent(input$selected2, {
    validNumericInput("selected2", defaultCenter, 1000)
  })

  createWindowOut <- function(j) {
    window <- jth_ref("window", j)
    # Select the window size (half-width) specified in the URL (if any)
    initialWindowSize <- isolate(values$urlFields[[window]])
    if (is.null(initialWindowSize)) { initialWindowSize <- defaultWindowSize }
    sliderInput(inputId = window, label="Window size around selected point:",
      min = 1000, max = 2000000, value = initialWindowSize, step = 1000)
  }
  output$windowOut <- renderUI(createWindowOut(1))
  outputOptions(output, "windowOut", suspendWhenHidden=FALSE)
  output$windowOut2 <- renderUI(createWindowOut(2))
  outputOptions(output, "windowOut2", suspendWhenHidden=FALSE)
  
  #returns datasets from uploaded file
  getdata <- function(j) {
    if (is.null(input[[jth_ref("datasets", j)]])) { return() }
    values[[input[[jth_ref("datasets", j)]]]]
  }

  #builds list of column names and type in dataset
  varnames <- function(j) {
    if (is.null(input[[jth_ref("datasets", j)]])) return()
    dat <- getdata_class(j)
    vars <- names(dat)
    names(vars) <- paste(vars, " {", dat, "}", sep = "")
    vars
  }

  getdata_class <- function(j) {
    if (is.null(input[[jth_ref("datasets", j)]])) return()
    cls <- sapply(getdata(j), function(x) class(x)[1])
    gsub("ordered","factor", cls)
  }

  # Function to handle loading of GWAS data from a file or rObject, for the jth dataset
  loadUserData <- function(filename, uFile, j) {
    ext <- file_ext(filename)
    objname <- sub(paste(".", ext, sep = ""), "", basename(filename))
    ext <- tolower(ext)

    if (ext == 'rda' || ext == 'rdata') {
      # objname will hold the name of the object(s) inside the R datafile
      robjname <- load(uFile)
      if (length(robjname) > 1) {
        loaded.values <- data.frame(get(robjname[1]))
      } else {
        loaded.values <- data.frame(get(robjname))  # only work with data.frames
      }
    }

    appendSNPs <- isolate(input[[jth_ref("appendSNPs", j)]])
    if (!appendSNPs) {
      # add new datasets to the datasetToOrganism map
      values$datasetToOrganism[[objname]] <- values[[jth_ref("organism", j)]]
    }
    if (length(values[['datasetlist']]) == 0 || values[['datasetlist']][1] == '') {
      values[['datasetlist']] <- c(objname)
    } else if (!appendSNPs) {
      values[['datasetlist']] <- unique(c(objname, values[['datasetlist']]))
    }

    if (ext == 'sav') {
      loaded.values <- as.data.frame(as.data.set(spss.system.file(uFile)))
    } else if (ext == 'dta') {
      loaded.values <- read.dta(uFile)
    } else {
      loaded.values <- read.csv(uFile, header = input[[jth_ref("header", j)]], sep = input[[jth_ref("sep", j)]], stringsAsFactors = FALSE)
    }
    # Force the chromosome column to be a string
    # (Note: there must be exactly one column whose name starts with "chr", case-insensitive)
    cols <- tolower(names(loaded.values))
    k <- which(startsWith(cols, "chr"))
    loaded.values[, k] <- as.character(loaded.values[, k])

    if (appendSNPs) {
      dsj <- jth_ref("datasets", j)
      values[[input[[dsj]]]]$totalBP <- NULL
      names(loaded.values) <- names(values[[input[[dsj]]]])
      values[[input[[dsj]]]] <- unique(rbind(values[[input[[dsj]]]], loaded.values))
    } else {
      values[[objname]] <- loaded.values
    }
  }
  
  # Load GWAS data from .csv files at a remote URL, for the jth dataset
  loadRemoteData <- function(trait, traitUrl, j) {
    ext <- file_ext(traitUrl)
    objname <- sub(paste(".", ext, sep = ""), "", basename(traitUrl))

    appendSNPs <- isolate(input[[jth_ref("appendSNPs", j)]])
    if (!appendSNPs) {
      # add new datasets to the datasetToOrganism map
      values$datasetToOrganism[[objname]] <- values[[jth_ref("organism", j)]]
    }
    if (length(values[['datasetlist']]) == 0 || values[['datasetlist']][1] == '') {
      values[['datasetlist']] <- c(objname)
    } else if (!appendSNPs) {
      values[['datasetlist']] <- unique(c(objname, values[['datasetlist']]))
    }

    # Note: there must be exactly one column whose name starts with "chr", case-insensitive
    loaded.values <- load.gwas.remote(values[[jth_ref("organism", j)]], traitUrl, trait)

    if (appendSNPs) {
      dsj <- jth_ref("datasets", j)
      values[[input[[dsj]]]]$totalBP <- NULL
      names(loaded.values) <- names(values[[input[[dsj]]]])
      values[[input[[dsj]]]] <- unique(rbind(values[[input[[dsj]]]], loaded.values))
    } else {
      values[[objname]] <- loaded.values
    }

    # Keep track of which remote GWAS traits are loaded, for use in the URL
    # (do not use unique(values[[input[[dsj]]]]$Trait) which may include local traits)
    gwj <- jth_ref("gwasTraits", j)
    values[[gwj]] <- union(values[[gwj]], trait)
  }

  #this function makes the chromsomeview chart  
  #subset whole chart based on selection
  output$pChart <- renderChart(create_pChart(1, input, values))
  output$pChart2 <- renderChart(create_pChart(2, input, values))
  output$pChartMacro <- renderChart(create_pChartMacro(1, input, values))
  output$pChartMacro2 <- renderChart(create_pChartMacro(2, input, values))
  
  #Genome wide chart
  output$gChart <- renderChart(create_gChart(1, input, values))
  output$gChart2 <- renderChart(create_gChart(2, input, values))
  output$gChartMacro <- renderChart(create_gChartMacro(1, input, values))
  output$gChartMacro2 <- renderChart(create_gChartMacro(2, input, values))

  output$zChart <- renderChart(create_zChart(1, input, values))
  output$zChart2 <- renderChart(create_zChart(2, input, values))

  #highcharts test chart
  output$testChart <- renderChart({
    h1 <- hPlot(x = "Wr.Hnd", y = "NW.Hnd", data = MASS::survey, type = c("line", 
                                                                          "bubble", "scatter"), group = "Clap", size = "Age")
    h1$set(dom = 'testChart')
    return(h1)     
  })

  handleTraitFilter <- function(j) {
    if (is.null(input[[jth_ref("traitFilter", j)]])) return()
    isolate({
      tc <- input[[jth_ref("traitColumns", j)]]
      tc.id <- jth_ref(tc, j)
      df.data <- values[[input[[jth_ref("datasets", j)]]]]
      all.choices <- sort(unique(df.data[, tc]))
      tf <- input[[jth_ref("traitFilter", j)]]
      # Select traits specified in the URL (if any)
      traitsFromUrl <- isolate(values$urlFields[[jth_ref("traits", j)]])
      if (!is.null(traitsFromUrl)) {
        filtered.choices <- unlist(strsplit(traitsFromUrl, split = ";"))
      } else if (nchar(tf) == 0) {
        filtered.choices <- "Select All"
      } else {
        filtered.choices <- all.choices[grep(tf, all.choices, ignore.case = TRUE)]
      }
      all.choices <- c("Select All", "Deselect All", all.choices)
      hasGwasTraits <- any(is.na(df.data[, input[[jth_ref("SIbpStart", j)]]]))
      hasQtlTraits <- any(!is.na(df.data[, input[[jth_ref("SIbpStart", j)]]]))
      if (hasGwasTraits && hasQtlTraits) {
        all.choices <- c(all.choices[1], "Select only GWAS traits", "Select only QTL traits", all.choices[-1])
      }
      updateSelectizeInput(session, tc.id, choices = all.choices, selected = filtered.choices)
    })
  }
  observe(handleTraitFilter(1))
  observe(handleTraitFilter(2))

  handleSubmitColsButton <- function(j) {
    if(is.null(input[[jth_ref("SubmitColsButton", j)]]) || input[[jth_ref("SubmitColsButton", j)]] == 0){return()}
    isolate({
      currDatasetProp <- datasetProp()
      #print("before")
      #print(currDatasetProp)
      if(as.character(input[[jth_ref("datasets", j)]]) %in% currDatasetProp$dataset){
        currDatasetProp <- currDatasetProp[currDatasetProp$dataset != as.character(input[[jth_ref("datasets", j)]]),]
      }
      #print("after")
      #print(currDatasetProp)
      cols <- varnames(j)
      #print("data.frame")
      #print(data.frame(dataset=input$datasets,chrColumn=names(cols[cols==input$chrColumn]),bpColumn=names(cols[cols==input$bpColumn]),
      #                 traitCol=paste(names(cols[cols %in% input$traitColumns]),collapse=";"),yAxisColumn=names(cols[cols==input$yAxisColumn]),axisLim=input$axisLimBool,axisMin=input$axisMin,axisMax=input$axisMax,stringsAsFactors=FALSE))
      currDatasetProp <-  rbind(currDatasetProp,data.frame(dataset=input[[jth_ref("datasets", j)]],chrColumn=names(cols[cols==input[[jth_ref("chrColumn", j)]]]),bpColumn=names(cols[cols==input[[jth_ref("bpColumn", j)]]]),
        traitCol=paste(names(cols[cols %in% input[[jth_ref("traitColumns", j)]]]),collapse=";"),yAxisColumn=names(cols[cols==input[[jth_ref("yAxisColumn", j)]]]),
        logP=input[[jth_ref("logP", j)]],axisLim=input[[jth_ref("axisLimBool", j)]],axisMin=input[[jth_ref("axisMin", j)]],axisMax=input[[jth_ref("axisMax", j)]],organism=values[[jth_ref("organism", j)]],plotAll=input[[jth_ref("plotAll", j)]],
        supportInterval=input[[jth_ref("supportInterval", j)]],SIyAxisColumn=input[[jth_ref("SIyAxisColumn", j)]],SIbpStart=input[[jth_ref("SIbpStart", j)]],SIbpEnd=input[[jth_ref("SIbpEnd", j)]],
        SIaxisLimBool=input[[jth_ref("SIaxisLimBool", j)]],SIaxisMin=input[[jth_ref("SIaxisMin", j)]],SIaxisMax=input[[jth_ref("SIaxisMax", j)]],stringsAsFactors=FALSE))      
      #print("rbind")
      #print(currDatasetProp)
      # write.table(file="./www/config/datasetProperties.csv",x=currDatasetProp,col.names=TRUE,row.names=FALSE,sep=",")
      updateTabsetPanel(session, "datatabs", selected = "WhGen")
    })
#    if(input$selected != 1e5){
#      updateTabsetPanel(session, "datatabs", selected = "WhGen")  
#    }
  }
  observe(handleSubmitColsButton(1))
  observe(handleSubmitColsButton(2))

  # saveDataset <- function(j) {
  #   if (is.null(input[[jth_ref("saveDatasetButton", j)]]) || input[[jth_ref("saveDatasetButton", j)]] == 0) { return() }
  #   isolate({
  #     if(!file.exists(paste0("./www/config/data/",input[[jth_ref("datasets", j)]]))){
  #       write.table(getdata(j),paste0("./www/config/data/",input[[jth_ref("datasets", j)]]),sep=",",col.names=TRUE,row.names=FALSE)
  #     }
  #   })
  # }
  # observe(saveDataset(1))
  # observe(saveDataset(2))
  
   datasetProp <- function(){
     dp <- read.table("./www/config/datasetProperties.csv",sep=",",head=TRUE,stringsAsFactors=FALSE)
     if (is.null(values$datasetToOrganism)) {
       # populate the datasetToOrganism map from stored values
       values$datasetToOrganism <- list()
       for (i in 1:nrow(dp)) {
         values$datasetToOrganism[[dp$dataset[i]]] <- dp$organism[i]
       }
     }
     # replace hyphens in strings like "p-value" as they are not allowed in column names
     # (note that the above read.table() converts them to a period)
     dp$yAxisColumn <- sapply(dp$yAxisColumn, function(s) gsub("-", ".", s))
     dp$SIyAxisColumn <- sapply(dp$SIyAxisColumn, function(s) gsub("-", ".", s))
     return(dp)
   }
   #this provides the functionality to update the plotband window after a user clicks without rerendering the whole plot
   #From:
   #http://stackoverflow.com/questions/20247759/add-highcharts-plotband-after-render-in-r-shiny-rcharts/20249933?noredirect=1#20249933
   updatePlotbandWindow <- function(j) {
     center <- as.numeric(input[[jth_ref("selected", j)]])
     winHigh <- center + input[[jth_ref("window", j)]]
     winLow <- center - input[[jth_ref("window", j)]]
     #eventually I would use winLow/winHigh to change the plotband range
     band = list(from = winLow, to = winHigh, color = windowPlotBandColor)
     #print(band)
     session$sendCustomMessage(type = jth_ref("customMsg", j), band)
   }
   observe(updatePlotbandWindow(1))
   observe(updatePlotbandWindow(2))

  output$selectedGene <- renderUI({
    # break lines at '.' to ensure visibility
    h5(HTML(gsub("\\.", "\\.<br>", values$glSelectedGene)))
  })
  output$relatedRegions <- renderUI({
    selectInput("relatedRegions", "Related Regions:", choices = NULL, selectize = FALSE)
  })
  observe({
    if (is.null(input$relatedRegions) || length(input$relatedRegions) == 0) {
      updateNumericInput(session, "selected2", value = input$selected2)
    } else {
      # parse from the format "[Chr] [minBP]-[maxBP] Mbp"
      # (or "chr[Chr]" if [Chr] is a number)
      ss <- unlist(strsplit(input$relatedRegions, split = " "))
      chr <- ifelse(hasNumericChromosomeNames(values$organism2), stri_sub(ss[1], 4), ss[1])
      ss2 <- unlist(strsplit(ss[2], split = "-"))
      centerBP <- as.integer(1.0e6*mean(as.numeric(ss2)))
      updateSelectInput(session, "chr2", selected = chr)
      updateNumericInput(session, "selected2", value = centerBP)
    }
  })

  observe({
    # Change the Broadcast Channel name
    if (!is.null(input$bcName)) {
      runjs(paste(
        # Close any existing channel
        "try {",
          "bc.close();",
          # "console.log('Closed GCV broadcast channel \"' + bc.name + '\"');",
        "} catch (ex) {}",
        # Open the new one
        "try {",
          sprintf("bc = new BroadcastChannel('%s');", input$bcName),
          "bc.onmessage = function(e) {",
            "Shiny.onInputChange('bc_gcv', e.data);",
            # "console.log(e.data);",
          "};",
          # "console.log('Opened GCV broadcast channel \"' + bc.name + '\"');",
        "} catch (ex) {",
          # "console.log('Could not open GCV broadcast channel \"' + bc.name + '\"');",
        "}"
      ))
    }
  })

  # Recompute macro-synteny blocks
  observeEvent(input$submitMacroSynteny, {
    # Clear any existing macro-synteny blocks before proceeding
    values$pairwiseBlocks <- list(NULL, NULL)

    org1 <- values$organism
    org2 <- values$organism2
    if (is.null(org1) || is.null(org2)) return()
    if (!("family" %in% names(org.annotGeneLoc[[org1]]) && "family" %in% names(org.annotGeneLoc[[org2]]))) return()

    # Consider all pairs of chromosomes between species 1 and 2
    for (chr1 in chrName[[org1]]) {
      # drop = TRUE to avoid contaminating the JSON in the macroSyntenyBlocksMicroservice() call
      ff <- subset(org.annotGeneLoc[[org1]], chromosome == chr1, select = family, drop = TRUE)
      tt <- sapply(1:length(chrName[[org2]]), function(chr2) {
        sprintf(org.annotChrFormat[[org2]], chr2)
      })
      distanceMetric <- tolower(input$macroDistance)
      if (distanceMetric == "jaccard") distanceMetric <- paste(distanceMetric, input$macroNgram, tolower(input$macroReversals), sep = ":")
      blocks <- macroSyntenyBlocksMicroservice(org.gcvUrlBase[[org1]], ff, input$macroMatched, input$macroIntermediate, input$macroMask, tt, distanceMetric, chr1, org1, org2)
      if (!is.null(blocks$error)) {
        print(blocks$error)
      } else {
        for (j in 1:2) {
          if (is.null(blocks[[j]])) next
          if (is.null(values$pairwiseBlocks[[j]])) {
            values$pairwiseBlocks[[j]] <- blocks[[j]]
          } else {
            values$pairwiseBlocks[[j]] <- rbind(values$pairwiseBlocks[[j]], blocks[[j]])
          }
        }
      }
    }
    if (!is.null(values$pairwiseBlocks[[1]])) values$pairwiseBlocks[[1]] <- unique(values$pairwiseBlocks[[1]])
    updateURL()
  })

  clearGenomicLinkages <- function() {
    values$glGenes <- values$glGenes2 <- values$glColors <- NULL
    values$highlightGenes <- NULL
    updateSelectInput(session, "relatedRegions", choices = character(0))
  }

  # Clear any existing genomic linkage query results
  observeEvent(input$clearGenomicLinkage, {
    values$glSelectedGene <- NULL
    clearGenomicLinkages()
  })

  # Determine genomic linkages using the GCV microservices
  observeEvent(input$selectedGene, {
    # User selected a gene in the organism 1 zChart, or in the URL
    if (is.null(input$selectedGene)) return()

    org1 <- values$organism
    values$glSelectedGene <- input$selectedGene
    df.genes <- subset(org.annotGeneLoc[[org1]], chromosome == input$chr,
      select = c(id, chromosome, transcript_start, transcript_end, strand, family))
    nt <- nrow(df.genes)
    n0 <- which(df.genes$id == values$glSelectedGene)
    n <- input$neighbors
    nn <- n0 + (-n:n)
    # if at either end of the chromosome, adjust to retain (2*n + 1) genes
    if (nn[1] < 1) {
      nn <- nn - (nn[1] - 1)
    } else if (tail(nn, 1) > nt) {
      nn <- nn - (tail(nn, 1) - nt)
    }
    values$glGenes <- df.genes[nn, ]
    values$glGenes2 <- NULL # to disable redrawing the zCharts until we have the micro-synteny-search results

    # Send to the micro-synteny-search service
    results <- microSyntenySearchMicroservice(org.gcvUrlBase[[org1]], values$glGenes$family, input$matched, input$intermediate)
    if (!is.null(results$error)) {
      print(results$error)
    } else {
      if (length(results$tracks) == 0) {
        clearGenomicLinkages()
        return()
      }
      # Parse related genes from species 2
      for (i in 1:length(results$tracks)) {
        results$tracks[[i]]$id <- i
      }
      org2 <- values$organism2
      df.annot <- subset(org.annotGeneLoc[[org2]], select = c(name, chromosome, transcript_start, transcript_end, strand, family))
      values$glGenes2 <- do.call(rbind, lapply(results$tracks, FUN = function(tr) {
        sp <- unlist(strsplit(tr$species, split = ":"))[1]
        if (paste(substr(tr$genus, 1, 1), sp, sep = ".") == org.G.species[[org2]]) {
          df.genes <- data.frame(name = sapply(unlist(tr$genes), geneNameFromId), trackId = tr$id, stringsAsFactors = FALSE)
          # workaround for A. thaliana: convert gene names like "arath.Col.AT1G28130" back to "AT1G28130"
          if (org2 == "Arabidopsis thaliana") {
            df.genes$name <- stri_match_first(df.genes$name, regex = "^arath.Col.(.+)$")[, 2]
          }
          df.genes <- merge(df.genes, df.annot)
          if (nrow(df.genes) == 0) return()
          df.genes <- df.genes[nchar(df.genes$family) > 0, ]
          df.genes
        }
      }))

      # Highlight families common to both genomes
      families <- intersect(values$glGenes$family, values$glGenes2$family)
      nf <- length(families)
      if (nf == 0) {
        clearGenomicLinkages()
        return()
      }

      # Create nf colors
      fc <- getRainbowColors(nf)
      familyColors <- list()
      lapply(1:nf, FUN = function(i) familyColors[[families[i]]] <<- stri_sub(fc[i], 1, 7))

      # Construct the chart data
      values$glGenes <- values$glGenes[values$glGenes$family %in% families, ]
      values$glGenes$color <- familyColors[values$glGenes$family]
      values$glGenes2 <- values$glGenes2[values$glGenes2$family %in% families, ]
      values$glGenes2$color <- familyColors[values$glGenes2$family]
      values$glColors <- familyColors

      # Construct the related regions (each corresponds to a track from results$tracks)
      glTrackIds <- unique(values$glGenes2$trackId)
      glRelatedRegions <- do.call(rbind.data.frame, lapply(glTrackIds, FUN = function(tr.id) {
        tr.genes <- subset(values$glGenes2, trackId == tr.id)
        tr.chr <- tr.genes$chr[1]
        # Prepend "chr" if chromosome name is a number
        chrd <- ifelse(hasNumericChromosomeNames(org2), paste0("chr", tr.chr), tr.chr)
        tr.minBP <- min(tr.genes$transcript_start)
        tr.maxBP <- max(tr.genes$transcript_end)
        list(region = sprintf("%s %3.2f-%3.2f Mbp", chrd, tr.minBP*1.0e-6, tr.maxBP*1.0e-6),
          chr = tr.chr, minBP = tr.minBP, maxBP = tr.maxBP)
      }))
      # Sort the related regions
      glRelatedRegions <- glRelatedRegions[with(glRelatedRegions, order(chr, minBP)), ]

      # Recenter the window around the selected gene
      centerBP1 <- (as.integer(values$glGenes$transcript_start[1]) + as.integer(values$glGenes$transcript_end[nrow(values$glGenes)])) %/% 2
      updateNumericInput(session, "selected", value = centerBP1)
      # Select the related region specified in the URL (if any)
      selectedRegion <- values$urlFields$relatedRegion
      if (!(is.null(selectedRegion) || selectedRegion %in% glRelatedRegions$region)) {
        # If the user-specified region does not exist, select the first one
        # (TODO: select the best match instead)
        selectedRegion <- NULL
      }
      values$urlFields$relatedRegion <- NULL # to reset it
      updateSelectInput(session, "relatedRegions", choices = glRelatedRegions$region, selected = selectedRegion)
    }
  })

  # Handle Broadcast Channel messages from the Genome Context Viewer
  # (move to a separate file?...)
  observeEvent(input$bc_gcv, {
    if (!input$boolListenToBC) return()
    if (input$bc_gcv$type == 'select') {
      # Parse the message
      flags = 0
# or handle each one separately...
# TODO: add a SNPs field
      if (!is.null(input$bc_gcv$targets$organism)) flags <- flags + 1
      if (!is.null(input$bc_gcv$targets$chromosome)) flags <- flags + 2
      if (!is.null(input$bc_gcv$targets$genes)) flags <- flags + 4
      if (!is.null(input$bc_gcv$targets$family)) flags <- flags + 8
      if (!is.null(input$bc_gcv$targets$extent)) flags <- flags + 16
      if (!is.null(input$bc_gcv$targets$block)) flags <- flags + 32
      isMacroSyntenyRow <- (flags == 3)
      isMacroSyntenyBlock <- (flags == 33)
      isMacroSyntenyOrganism <- (flags == 1)
      isMicroSyntenyRow <- (flags == 23)
      isMicroSyntenyGene <- (flags == 12)
      isMicroSyntenyFamily <- (flags == 8)
      # TODO:
      # isDotPlot <- isMicroSyntenyGene # they have the same fields
      # Multi-view fields:
      # isMacroSyntenyCircosInnerBlock <- isMacroSyntenyRow # they have the same fields
      # isMacroSyntenyCircosOuterBlock <- (flags == 32)

      # Determine the selected organism and its index
      j <- 0
      if (!is.null(input$bc_gcv$targets$organism)) {
        org <- input$bc_gcv$targets$organism # always in "Genus species" format
        if (org == org.Genus_species[values$organism]) {
          j <- 1
        } else if (org == org.Genus_species[values$organism2]) {
          j <- 2
        }
      }

      # Do something
      if (isMacroSyntenyRow) {
        if (j > 0) {
          chr <- trailingChromosomeName(input$bc_gcv$targets$chromosome, values[[jth_ref("organism", j)]])
          # Adjust the Chromosome window to match the selection
          updateTabsetPanel(session, "datatabs", selected = "Chrom")
          clearGenomicLinkages()
          updateSelectInput(session, jth_ref("chr", j), selected = chr)
        }

      } else if (isMacroSyntenyBlock) {
        # org <- input$bc_gcv$targets$organism
        blk <- input$bc_gcv$targets$block
        ref <- blk$reference
        srx <- blk$source
        genspRef <- stri_trans_totitle(substr(ref$chromosome, 1, 5))
        genspSrc <- stri_trans_totitle(substr(srx$chromosome, 1, 5))
        if (genspRef == org.Gensp[values$organism] && genspSrc == org.Gensp[values$organism2]) {
          updateTabsetPanel(session, "datatabs", selected = "Chrom")
          clearGenomicLinkages()
          # Adjust the organism 1 Chromosome window to match the block reference
          chrRef <- trailingChromosomeName(ref$chromosome, values$organism)
          bpStartRef <- ref$locus[[1]]
          bpEndRef <- ref$locus[[2]]
          bpCenterRef <- (bpStartRef + bpEndRef) %/% 2
          bpWidthRef <- (bpEndRef - bpStartRef) %/% 2
          updateSelectInput(session, "chr", selected = chrRef)
          updateNumericInput(session, "selected", value = bpCenterRef)
          # updateSliderInput(session, "window", value = bpWidthRef)
          values$chrNumber <- trailingInteger(chrRef)
          values$blockStart <- bpCenterRef - bpWidthRef
          values$blockEnd <- bpCenterRef + bpWidthRef
          cumBp <- c(0, cumsum(as.numeric(chrSize[[values$organism]])))
          values$blockCumStart <- values$blockStart + cumBp[values$chrNumber]
          values$blockCumEnd <- values$blockEnd + cumBp[values$chrNumber]
          # Adjust the organism 2 Chromosome window to match the block source
          chrSrc <- trailingChromosomeName(srx$chromosome, values$organism2)
          bpStartSrc <- srx$locus[[1]]
          bpEndSrc <- srx$locus[[2]]
          bpCenterSrc <- (bpStartSrc + bpEndSrc) %/% 2
          bpWidthSrc <- (bpEndSrc - bpStartSrc) %/% 2
          updateSelectInput(session, "chr2", selected = chrSrc)
          updateNumericInput(session, "selected2", value = bpCenterSrc)
          # updateSliderInput(session, "window2", value = bpWidthSrc)
          values$chrNumber2 <- trailingInteger(chrSrc)
          values$blockStart2 <- bpCenterSrc - bpWidthSrc
          values$blockEnd2 <- bpCenterSrc + bpWidthSrc
          cumBp2 <- c(0, cumsum(as.numeric(chrSize[[values$organism2]])))
          values$blockCumStart2 <- values$blockStart2 + cumBp2[values$chrNumber2]
          values$blockCumEnd2 <- values$blockEnd2 + cumBp2[values$chrNumber2]
        }

      } else if (isMacroSyntenyOrganism) {
        org <- input$bc_gcv$targets$organism
        # TODO: ...

      } else if (isMicroSyntenyRow) {
        if (j > 0) {
          chr <- trailingChromosomeName(input$bc_gcv$targets$chromosome, values[[jth_ref("organism", j)]])
          # Range of base pairs
          bpMin <- input$bc_gcv$targets$extent[[1]]
          bpMax <- input$bc_gcv$targets$extent[[2]]
          centerBP <- (bpMax + bpMin) %/% 2
          widthBP <- (bpMax - bpMin) %/% 2 + 50000 # give it 50k BPs on either side to ensure visibility
          # Adjust the Chromosome window to match the selection
          updateTabsetPanel(session, jth_ref("datatabs", j), selected = "Chrom")
          clearGenomicLinkages()
          updateSelectInput(session, jth_ref("chr", j), selected = chr)
          updateNumericInput(session, jth_ref("selected", j), value = centerBP)
          updateSliderInput(session, jth_ref("window", j), value = widthBP)

          # TODO: Do something with the genes?
          # genes <- input$bc_gcv$targets$genes
        }

      } else if (isMicroSyntenyGene) {
        # Highlight the selected gene
        values$highlightGenes <- input$bc_gcv$targets$genes[[1]]
        fam <- input$bc_gcv$targets$family

      } else if (isMicroSyntenyFamily) {
        fam <- input$bc_gcv$targets$family
        # Check for singleton and orphan genes
        isSingletons <- startsWith(fam, "singleton")
        # Parse "singleton,<family 1>,<family 2>,..."
        if (isSingletons) fam <- unlist(strsplit(fam, split = ","))[-1]
        # isOrphans <- (fam == "")

        # Determine and highlight the selected genes
        values$highlightGenes <- subset(org.annotGeneLoc[[values$organism]], family %in% fam, select = id)
        values$highlightGenes <- unlist(union(values$highlightGenes, subset(org.annotGeneLoc[[values$organism2]], family %in% fam, select = id)))
      }
    }
  })
  # The messages GCV currently implements conform to the following schema:
    # {
    #   type: "select" | "deselect",
    #   targets: {
    #     organism?: String,  // an organism identifier of the form "<genus> <species>"
    #     chromosome?: String,  // a chromosome identifier
    #     genes?: String[],  // an array of gene identifiers
    #     family?: String,  // a gene family identifier
    #     extent?: Integer[2],  // a length 2 array representing a genomic interval
    #     block?: {  // an object representing a pairwise synteny block
    #       source: {
    #         chromosome: String,
    #         locus: Integer[2]
    #       },
    #       reference: {
    #         chromosome: String,
    #         locus: Integer[2]
    #       },
    #       orientation: "+" | "-"
    #     }
    #   }
    # }

  # Turn off (deselect) any old highlighted targets
  clearGCVHighlights <- function() {
    if (!is.null(values$gcvTargets)) {
      robj <- list(type = "deselect", targets = values$gcvTargets)
      msg <- sprintf("try { bc.postMessage(JSON.parse('%s')); } catch (ex) {}", toJSON(robj, auto_unbox = TRUE))
      runjs(msg)
    }
  }

  # View the current genomic linkage query in the Genome Context Viewer
  # (this does not involve Broadcast Channel)
  observeEvent(input$viewInGCV, {
    if (is.null(values$glSelectedGene)) {
      alert("No gene selected.")
    } else {
      # filter by organism 2 genome prefix to avoid displaying similar tracks from different genomes
      gene2 <- org.annotGeneLoc[[values$organism2]]$id[1]
      regexp2 <- stri_match_first(gene2, regex = "(^.+)\\.ann\\d+\\.")[, 2]
      gcvQuery <- sprintf("window.open('%s/gene;lis=%s?neighbors=%d&matched=%s&intermediate=%d&regexp=%s', 'gcv');",
        userConfig$gcv_client_url, values$glSelectedGene, input$neighbors, input$matched, input$intermediate, regexp2)
      runjs(gcvQuery)
    }
  })

  # Post a Broadcast Channel message to the Genome Context Viewer
  # (which handles the same kind of messages it sends)
  observe({
    if (is.null(input$gcvGeneFamily)) return()
    isolate(clearGCVHighlights())

    # User selected a family (in the Chromosome View) to highlight in the Genome Context Viewer
    gcvTargets = list(family = input$gcvGeneFamily)
    robj <- list(type = "select", targets = gcvTargets)
    msg <- sprintf("try { bc.postMessage(JSON.parse('%s')); } catch (ex) {}", toJSON(robj, auto_unbox = TRUE))
    values$gcvTargets <- gcvTargets
    runjs(msg)

    # Then reset it to enable selecting the same family in the future
    runjs("Shiny.onInputChange('gcvGeneFamily', null);")
  })

  # Tour based on rintrojs
  df.tour <- read.csv(file = "tour.csv", header = TRUE)
  observeEvent(input$zbrowseTour, {
    introjs(session,
      options = list(steps = df.tour, showBullets = FALSE, showStepNumbers = FALSE, skipLabel = "End Tour"),
      # allow switching tabs
      events = list(onbeforechange = readCallback("switchTabs"))
    )
  })

  # Update application state to URL on the fly
  updateURL <- function() {
    # required for query string: selected tab; inputs for each organism
    url.q <- paste0("?tab=", isolate(input$datatabs))
    # number of species to display
    cmp2sp <- isolate(input$compare2species)
    twoSpecies <- (is.null(cmp2sp) || cmp2sp)
    ts <- ifelse(twoSpecies, "y", "n")
    ns <- ifelse(twoSpecies, 2, 1)
    url.q <- paste0(url.q, "&twoSpecies=", ts)
    ss <- outer(c("datasets", "chr", "selected", "window"), 1:ns, FUN = jth_ref)
    for (s in ss) {
      x <- isolate(input[[s]])
      if (!is.null(x)) {
        url.q <- paste0(url.q, "&", s, "=", x)
      }
    }
    # remote GWAS traits to load
    for (j in 1:ns) {
      tfoj <- jth_ref("traitsFromOntology", j)
      url.q <- paste0(url.q, "&", tfoj, "=", ifelse(input[[tfoj]], "y", "n"))
      gwj <- jth_ref("gwasTraits", j)
      if (!is.null(values[[gwj]])) {
        url.q <- paste0(url.q, "&", gwj, "=", paste(values[[gwj]], collapse = ";"))
      }
    }
    # selected traits
    for (j in 1:ns) {
      numTraits <- 0
      allTraits <- isolate(input[[jth_ref("datasets", j)]])
      traits <- ""
      for (i in input[[jth_ref("traitColumns", j)]]) {
        try({
          numTraits <- length(unique(values[[allTraits]][,i]))
        }, silent = TRUE)
        x <- input[[jth_ref(i, j)]]
        if (!(is.null(x) || length(x) == 0 || length(x) == numTraits)) {
          traits <- paste0(x, collapse=";")
        }
      }
      trj <- jth_ref("traits", j)
      if (traits != "") url.q <- paste0(url.q, "&", trj, "=", traits)
    }
    # optional: genomic linkage information
    if (twoSpecies) {
      if (isolate(!is.null(values$glSelectedGene))) {
        ss.gl <- c("neighbors", "matched", "intermediate")
        for (s in ss.gl) {
          x <- isolate(input[[s]])
          if (!is.null(x)) {
            url.q <- paste0(url.q, "&", s, "=", x)
          }
        }
        # optional: selected gene and related region
        if (!is.null(values$glSelectedGene)) {
          url.q <- paste0(url.q, "&selectedGene=", values$glSelectedGene)
        }
        if (!is.null(isolate(input$relatedRegions))) {
          url.q <- paste0(url.q, "&relatedRegion=", isolate(input$relatedRegions))
        }
      }
      # optional: macro-synteny information
      ss.gl <- c("macroMatched", "macroIntermediate", "macroMask", "macroDistance")
      for (s in ss.gl) {
        x <- isolate(input[[s]])
        if (!is.null(x)) {
          if (s == "macroDistance") {
            x <- tolower(x)
            if (x == "jaccard") {
              x <- paste(x, isolate(input$macroNgram), tolower(isolate(input$macroReversals)), sep = ":")
            }
          }
          url.q <- paste0(url.q, "&", s, "=", x)
        }
      }
    }
    updateQueryString(url.q)
    # updateQueryString(url.q, mode = "push") # to save history
  }
  observeEvent(
    eventExpr = {
      # any input that can trigger the handlerExpr
      input$datasets
      input$datasets2
      values$gwasTraits
      values$gwasTraits2
      input$traitsFromOntology
      input$traitsFromOntology2
      sapply(input$traitColumns, function(i) input[[jth_ref(i, 1)]])
      sapply(input$traitColumns2, function(i) input[[jth_ref(i, 2)]])
      input$chr
      input$selected
      input$window
      input$chr2
      input$selected2
      input$window2
      values$glSelectedGene
      input$relatedRegions
      input$datatabs
      input$compare2species
    },
    handlerExpr = updateURL(), ignoreInit = TRUE)

  # Get genomic linkages for the selected gene specified in the URL (if any)
  observeEvent(input$datatabs, {
    if (input$datatabs != "Chrom") return()

    if (is.null(values$urlFields$selectedGene)) return()

    # Query the Genome Context Viewer for genomic linkages
    runjs(sprintf("Shiny.onInputChange('selectedGene', '%s');", values$urlFields$selectedGene))

    values$urlFields$selectedGene <- NULL # to reset it
  }, ignoreInit = TRUE)

  # Validate numericInputs
  validNumericInput <- function(inputId, default, min, max = Inf, FUN = updateNumericInput) {
    ii <- input[[inputId]]
    if (length(ii) == 0) return(FALSE)
    x <- suppressWarnings(as.numeric(ii))
    if (is.na(x)) {
      FUN(session, inputId, value = default)
      return(FALSE)
    } else if (x < min) {
      FUN(session, inputId, value = min)
      return(FALSE)
    } else if (max < Inf && x > max) {
      FUN(session, inputId, value = max)
      return(FALSE)
    }
    TRUE
  }

  observeEvent(input$macroMatched, {
    validNumericInput("macroMatched", userConfig$macroMatched, 1, 100)
  })
  observeEvent(input$macroIntermediate, {
    validNumericInput("macroIntermediate", userConfig$macroIntermediate, 0, 100)
  })
  observeEvent(input$macroMask, {
    validNumericInput("macroMask", userConfig$macroMask, 1, 100)
  })
  observeEvent(input$macroNgram, {
    validNumericInput("macroNgram", 1, 1, 2)
  })
  observeEvent(input$neighbors, {
    validNumericInput("neighbors", userConfig$neighbors, 1, 100)
  })
  observeEvent(input$matched, {
    validNumericInput("matched", userConfig$matched, 0.01, 100, updateTextInput)
  })
  observeEvent(input$intermediate, {
    validNumericInput("intermediate", userConfig$intermediate, 0, 100)
  })

  observeEvent(input$set_gChartMacro, {
    j <- input$set_gChartMacro$j
    values[[jth_ref("blockCumStart", j)]] <- unlist(input$set_gChartMacro$min)
    values[[jth_ref("blockCumEnd", j)]] <- unlist(input$set_gChartMacro$max)
    cumBp_j <- c(0, cumsum(as.numeric(chrSize[[values[[jth_ref("organism", j)]]]])))
    values[[jth_ref("chrNumber", j)]] <- chr_j <- which(values[[jth_ref("blockCumStart", j)]] <= cumBp_j)[1] - 1
    values[[jth_ref("blockStart", j)]] <- values[[jth_ref("blockCumStart", j)]] - cumBp_j[chr_j]
    values[[jth_ref("blockEnd", j)]] <- values[[jth_ref("blockCumEnd", j)]] - cumBp_j[chr_j]
    if (j == 1) {
      values$blockCumStart2 <- unlist(input$set_gChartMacro$minSrc)
      values$blockCumEnd2 <- unlist(input$set_gChartMacro$maxSrc)
      cumBp2 <- c(0, cumsum(as.numeric(chrSize[[values$organism2]])))
      # there may be multiple species 2 blocks in the pCharts
      values$chrNumber2 <- sapply(values$blockCumStart2, function(bp) which(bp <= cumBp2)[1] - 1)
      values$blockStart2 <- values$blockCumStart2 - cumBp2[values$chrNumber2]
      values$blockEnd2 <- values$blockCumEnd2 - cumBp2[values$chrNumber2]
    } else if (j == 2) {
      values$blockCumStart <- unlist(input$set_gChartMacro$minRef)
      values$blockCumEnd <- unlist(input$set_gChartMacro$maxRef)
      cumBp <- c(0, cumsum(as.numeric(chrSize[[values$organism]])))
      # there is only one species 1 block, otherwise we could handle multiple blocks as we do above for species 2
      values$chrNumber <- which(values$blockCumStart <= cumBp)[1] - 1
      values$blockStart <- values$blockCumStart - cumBp[values$chrNumber]
      values$blockEnd <- values$blockCumEnd - cumBp[values$chrNumber]
    }
    updateSelectInput(session, "chr", selected = chrName[[values$organism]][values$chrNumber])
    updateNumericInput(session, "selected", value = (values$blockStart + values$blockEnd) %/% 2)
    updateSelectInput(session, "chr2", selected = chrName[[values$organism2]][values$chrNumber2[1]])
    updateNumericInput(session, "selected2", value = (values$blockStart2[1] + values$blockEnd2[1]) %/% 2)
  })

  observeEvent(input$set_pChartMacro, {
    j <- input$set_pChartMacro$j
    values$chrNumber <- unlist(input$set_pChartMacro$chrNumber)
    values$chrNumber2 <- unlist(input$set_pChartMacro$chrNumber2)
    values[[jth_ref("blockStart", j)]] <- unlist(input$set_pChartMacro$min)
    values[[jth_ref("blockEnd", j)]] <- unlist(input$set_pChartMacro$max)
    cumBp_j <- c(0, cumsum(as.numeric(chrSize[[values[[jth_ref("organism", j)]]]])))
    chr_j <- values[[jth_ref("chrNumber", j)]]
    values[[jth_ref("blockCumStart", j)]] <- values[[jth_ref("blockStart", j)]] + cumBp_j[chr_j]
    values[[jth_ref("blockCumEnd", j)]] <- values[[jth_ref("blockEnd", j)]] + cumBp_j[chr_j]
    if (j == 1) {
      values$blockStart2 <- unlist(input$set_pChartMacro$minSrc)
      values$blockEnd2 <- unlist(input$set_pChartMacro$maxSrc)
      cumBp2 <- c(0, cumsum(as.numeric(chrSize[[values$organism2]])))
      values$blockCumStart2 <- values$blockStart2 + cumBp2[values$chrNumber2]
      values$blockCumEnd2 <- values$blockEnd2 + cumBp2[values$chrNumber2]
      # also update species 2 chromosome in Chromosome view
      updateSelectInput(session, "chr2", selected = chrName[[values$organism2]][values$chrNumber2[1]])
      updateNumericInput(session, "selected2", value = (values$blockStart2 + values$blockEnd2) %/% 2)
    } else if (j == 2) {
      values$blockStart <- unlist(input$set_pChartMacro$minRef)
      values$blockEnd <- unlist(input$set_pChartMacro$maxRef)
      cumBp <- c(0, cumsum(as.numeric(chrSize[[values$organism]])))
      values$blockCumStart <- values$blockStart + cumBp[values$chrNumber]
      values$blockCumEnd <- values$blockEnd + cumBp[values$chrNumber]
      # also update species 1 chromosome in Chromosome view
      updateSelectInput(session, "chr", selected = chrName[[values$organism]][values$chrNumber])
      updateNumericInput(session, "selected", value = (values$blockStart + values$blockEnd) %/% 2)
    }
  })

  clearMacrosyntenyBands <- function() {
    values$blockCumStart <- values$blockCumStart2 <- NULL
    values$blockCumEnd <- values$blockCumEnd2 <- NULL
    values$chrNumber <- values$chrNumber2 <- NULL
    values$blockStart <- values$blockStart2 <- NULL
    values$blockEnd <- values$blockEnd2 <- NULL
  }
  observeEvent(input$clearMacrosynteny, {
    clearMacrosyntenyBands()
  })

  observeEvent(input$traitsFromOntology, {
    if (!("trait_id" %in% varnames(1))) return()
    updateSelectizeInput(session, "traitColumns", selected = ifelse(input$traitsFromOntology, "trait", "trait_id"))
  })
  observeEvent(input$traitsFromOntology2, {
    if (!("trait_id" %in% varnames(2))) return()
    updateSelectizeInput(session, "traitColumns2", selected = ifelse(input$traitsFromOntology2, "trait", "trait_id"))
  })

})#end server
