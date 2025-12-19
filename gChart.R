source("common.R")

#Genome wide chart
create_gChart <- function(j, input, values) {
  nid <- jth_ref("notify.create_gChart", j)
  showNotification(paste0("Creating Whole Genome chart for ", values[[jth_ref("organism", j)]], ". Please wait."),
    duration = NULL, id = nid, type = "message")

  calculateTotalBP(j, input, values)
  
  #subset whole chart based on selection
  genomeChart <- values[[input[[jth_ref("datasets", j)]]]]
  selcom <- trimws(input[[jth_ref("dv_select", j)]])
  if (selcom != "") genomeChart <- applyGlobalFilter(genomeChart, selcom)
  if(input[[jth_ref("plotAll", j)]] == FALSE){
    for(i in input[[jth_ref("traitColumns", j)]]){
      genomeChart <- genomeChart[genomeChart[,i] %in% input[[jth_ref(i, j)]],]
    }    
    if(length(input[[jth_ref("traitColumns", j)]]) > 1){
      genomeChart$trait <- do.call(paste,c(genomeChart[,input[[jth_ref("traitColumns", j)]]],sep="_"))    
    }else{
      genomeChart$trait <- genomeChart[,input[[jth_ref("traitColumns", j)]]]
    }
  }else{
    genomeChart$trait <- input[[jth_ref("datasets", j)]]
  }
  
  #Separate Support Interval data from GWAS data, if support, GWAS data is assumed to be anything that has an NA in the SIbpStart column
  if(input[[jth_ref("supportInterval", j)]] == TRUE){
    SIchart <- genomeChart[!(is.na(genomeChart[,input[[jth_ref("SIbpStart", j)]]])),]
    genomeChart <- genomeChart[is.na(genomeChart[,input[[jth_ref("SIbpStart", j)]]]),]
  }
  
  #filter genomeChart for only rows that have a base pair and yaxis value
  genomeChart <- genomeChart[!(is.na(genomeChart[,input[[jth_ref("bpColumn", j)]]])),]
  genomeChart <- genomeChart[!(is.na(genomeChart[,input[[jth_ref("yAxisColumn", j)]]])),]
  
  #if checked, filter for only overlapping SNPs
  if(!is.null(input[[jth_ref("overlaps", j)]]) & input[[jth_ref("overlaps", j)]] == TRUE){
    genomeChart <- findGWASOverlaps(genomeChart, j, input)
  }
  
  #check if there is any data for the selected traits
  gwas.data.exist <- TRUE
  if(nrow(genomeChart)==0){ #nothing is in the window, but lets still make a data.frame
    gwas.data.exist <- FALSE
    genomeChart <- values[[input[[jth_ref("datasets", j)]]]][1,]
    genomeChart[,input[[jth_ref("yAxisColumn", j)]]] <- -1    
    genomeChart$totalBP <- 1
    if(length(input[[jth_ref("traitColumns", j)]]) > 1){
      genomeChart$trait <- do.call(paste,c(genomeChart[,input[[jth_ref("traitColumns", j)]]],sep="_"))    
    }else{
      genomeChart$trait <- genomeChart[,input[[jth_ref("traitColumns", j)]]]
    }
    # workaround to display x axis when no data exist
    genomeChart[is.na(genomeChart)] <- -1
  }
  
  #take -log10 of y-axis column if requested
  if (input[[jth_ref("logP", j)]] == TRUE && gwas.data.exist) {
    genomeChart[,input[[jth_ref("yAxisColumn", j)]]] <- -log(genomeChart[,input[[jth_ref("yAxisColumn", j)]]],10)
  }
  
  #check if there is too much data (>2500 data points), trim to 2500
  if(nrow(genomeChart)>2500){
    cutVal <- sort(genomeChart[,input[[jth_ref("yAxisColumn", j)]]],decreasing = T)[2500]
    genomeChart <- genomeChart[genomeChart[,input[[jth_ref("yAxisColumn", j)]]] >= cutVal,]
  }
  
  colorTable <- getColorTable(j, input)
  publication <- genomeChart$publication
  publication[is.null(publication)] <- ""
  genomeTable <- data.frame(
    x = genomeChart$totalBP,
    y = genomeChart[, input[[jth_ref("yAxisColumn", j)]]],
    trait = genomeChart$trait,
    name = sprintf("Base Pair: %1$s<br/>Chromosome: %2$s<br/>",
      prettyNum(genomeChart[, input[[jth_ref("bpColumn", j)]]], big.mark = ","),
      genomeChart[, input[[jth_ref("chrColumn", j)]]]
    ),
    chr = genomeChart[,input[[jth_ref("chrColumn", j)]]],
    bp = genomeChart[,input[[jth_ref("bpColumn", j)]]],
    pub = publication,
    stringsAsFactors = FALSE
  )
  genomeSeries <- lapply(split(genomeTable, genomeTable$trait), function(x) {
    res <- lapply(split(x, rownames(x)), as.list)
    names(res) <- NULL
    res <- res[order(sapply(res, function(x) x$x))]
    return(res)
  })

  #build JL series
  qtl.data.exist <- FALSE
  if(input[[jth_ref("supportInterval", j)]]==TRUE){
    qtl.data.exist <- TRUE
    if(nrow(SIchart)==0){ #nothing is in the window, but lets still make a data.frame
      qtl.data.exist <- FALSE
      SIchart <- values[[input[[jth_ref("datasets", j)]]]][1,]
      SIchart[,input[[jth_ref("SIyAxisColumn", j)]]] <- -1    
      if(length(input[[jth_ref("traitColumns", j)]]) > 1){
        SIchart$trait <- do.call(paste,c(SIchart[,input[[jth_ref("traitColumns", j)]]],sep="_"))    
      }else{
        SIchart$trait <- SIchart[,input[[jth_ref("traitColumns", j)]]]
      }             
    }
    SIchart$loc_el <- SIchart$trait
    if (dynamic.interval.height) {
      SIchart$h <- computeBlockHeights(SIchart, c("SIbpStartTotal", "SIbpEndTotal"), chart.max.height, interval.bar.height)
    } else {
      SIchart$h <- SIchart[[input[[jth_ref("SIyAxisColumn", j)]]]]
    }
    SIchart$publication[is.null(SIchart$publication)] <- ""
    SIchart <- SIchart[order(SIchart$SIbpStartTotal),]
    jlTable <- adply(SIchart, 1, function(x) {
      data.frame(
        x = c(x$SIbpStartTotal, x$SIbpEndTotal, x$SIbpEndTotal),
        y = c(x$h, x$h, NA),
        trait = x$trait,
        name = sprintf("<table cellpadding='4' style='line-height:1.5'><tr><td align='left'>Interval: %1$s-%2$s<br>Chromosome: %3$s</td></tr></table>",
          prettyNum(x[[input[[jth_ref("SIbpStart", j)]]]], big.mark = ","),
          prettyNum(x[[input[[jth_ref("SIbpEnd", j)]]]], big.mark = ","),
          x[[input[[jth_ref("chrColumn", j)]]]]
        ),
        loc_el = x$loc_el,
        bp = x[[input[[jth_ref("bpColumn", j)]]]],
        chr = x[[input[[jth_ref("chrColumn", j)]]]],
        pub = x$publication,
        stringsAsFactors = FALSE
      ) #end jlTable
    }) #end adply
    jlTable <- jlTable[, c("x", "y", "trait", "name", "loc_el", "bp", "chr", "pub")]
  } #end build jlTable if support intervals
  
  #build list for where to put plotbands for this organism
  bigList <- list()
  title.margin <- 15
  cumBP<-c(0,cumsum(as.numeric(chrSize[[values[[jth_ref("organism", j)]]]])))
  for(i in 1:(length(cumBP)-1)){
    chrLabel <- chrName[[values[[jth_ref("organism", j)]]]][i]
    # split character names like "Arahy.01"
    k.dot <- unlist(gregexpr("\\.", chrLabel))
    if (k.dot > 1) {
      chrLabel <- paste(substring(chrLabel, 1, k.dot - 1), substring(chrLabel, k.dot), sep = "<br>")
      title.margin <- 30
    }
    bigList[[i]] <- list(
      from = cumBP[i] + 1, to = cumBP[i + 1],
      label = list(
        text = chrLabel,
        style = list(color = "#6D869F", fontSize = "10px"),
        verticalAlign = "bottom"
      )
    )
    # different color for odd plot bands
    if (i %% 2 != 0) bigList[[i]]$color <- chromosomePlotBandColor
  }    
  
  c <- rCharts::Highcharts$new()
  c$LIB$url <- 'highcharts/'

  if(input[[jth_ref("axisLimBool", j)]] == TRUE){       
    c$yAxis(title=list(text=input[[jth_ref("yAxisColumn", j)]]),min=input[[jth_ref("axisMin", j)]],max=input[[jth_ref("axisMax", j)]],startOnTick=FALSE)
  }else{
    c$yAxis(title=list(text=input[[jth_ref("yAxisColumn", j)]]),min=0,startOnTick=FALSE)
  }
  
  if(input[[jth_ref("supportInterval", j)]]==TRUE){
    if(input[[jth_ref("SIaxisLimBool", j)]] == TRUE){
      c$yAxis(visible=FALSE,title=list(text=input[[jth_ref("SIyAxisColumn", j)]]),min=input[[jth_ref("SIaxisMin", j)]],max=input[[jth_ref("SIaxisMax", j)]],gridLineWidth=0,minorGridLineWidth=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)
    }else{
      c$yAxis(visible=FALSE,title=list(text=input[[jth_ref("SIyAxisColumn", j)]]),min=0,max=chart.max.height,gridLineWidth=0,minorGridLineWidth=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)
    }
    
    if (qtl.data.exist) {
      d_ply(jlTable,.(trait),function(x){
        c$series(
          data = toJSONArray2(x,json=F,names=T),
          type = "line",
          name = unique(x$trait),
          dashStyle = 'Solid',
          marker = list(enabled=F),
          yAxis = 1,
          tooltip = list(
            pointFormat = '<span style="color:{point.color}">\u25a0</span> {series.name}<br>{point.pub}',
            followPointer = TRUE
          ),
          color = colorTable$color[colorTable$trait == as.character(unique(x$loc_el))])})            
    }
  }

  if (gwas.data.exist) {
    invisible(sapply(genomeSeries, function(x) {
      c$series(
        data = x,
        turboThreshold = 5000,
        type = "scatter",
        color = colorTable$color[colorTable$trait == as.character(x[[1]]$trait)],
        name = x[[1]]$trait
      )
    }))
  }

  if (!(gwas.data.exist || qtl.data.exist)) {
    # force it to show the x axis when no data exist
    c$series(data = list(x = 0, y = -1), type = "scatter")
  }

  c$chart(zoomType="x",alignTicks=FALSE,events=list(click = "#!function(event) {this.tooltip.hide();}!#"))
  c$title(text=paste(input[[jth_ref("datasets", j)]]," Results",sep=" "))
  c$subtitle(text="Rollover for more info. Drag chart area to zoom. Click point to switch to chromosome and annotation view.")
  
  js <- jth_ref("", j)
  doClickOnPoint <- sprintf(
    paste(
      "#! function(){$('select#chr%s').val(this.options.chr); $('select#chr%s').trigger('change'); $('input#selected%s').val(this.options.bp);",
      "$('input#selected%s').trigger('change'); $('ul#datatabs li').eq(0).removeClass('active');",
      "$('ul#datatabs li').eq(1).removeClass('active'); $('ul#datatabs li').eq(2).removeClass('active');",
      "$('ul#datatabs li').eq(4).removeClass('active');",
      "$('ul#datatabs li').eq(3).addClass('active');",
      "$('#pChart%s').trigger('change');$('#pChart%s').trigger('shown');",
      "$('.tab-content div').toggleClass(function(){if(this.getAttribute('data-value')=='Chrom' || this.getAttribute('data-value')=='WhGen'){return 'active';}else{return '';}});",
      "$('.tab-content div').trigger('change');$('ul#datatabs li').trigger('change');} !#"
    ), js, js, js, js, js, js)
  doClickOnLine <- sprintf(paste(
    "#! function() {",
      "$('select#chr%s').val(this.options.chr); $('select#chr%s').trigger('change'); $('input#selected%s').val(this.options.bp);",
      "$('input#selected%s').trigger('change'); $('ul#datatabs li').eq(0).removeClass('active');",
      "$('ul#datatabs li').eq(1).removeClass('active'); $('ul#datatabs li').eq(2).removeClass('active');",
      "$('ul#datatabs li').eq(4).removeClass('active');",
      "$('ul#datatabs li').eq(3).addClass('active');",
      "$('#pChart%s').trigger('change');$('#pChart%s').trigger('shown');",
      "$('.tab-content div').toggleClass(function(){if(this.getAttribute('data-value')=='Chrom' || this.getAttribute('data-value')=='WhGen'){return 'active';}else{return '';}});",
      "$('.tab-content div').trigger('change');$('ul#datatabs li').trigger('change');",
    "} !#"
  ), js, js, js, js, js, js)
  c$plotOptions(
    scatter = list(
      cursor = "pointer",
      point = list(
        events = list(
          #click = "#! function() { window.open(this.options.url); } !#")), #open webpage
          #click = "#! function(event) {alert(this.name);} !#")), #display popup
          #click = "#! function(event) {console.log(this);} !#")), #write object to log
          #click = "#! function(){$('input#selected').val(134); $('input#selected').trigger('change');} !#")),
          click = doClickOnPoint
        )
      ),
      marker = list(
        symbol = "circle",
        radius = 5
      ),
      tooltip = list(
        headerFormat = "<b>{series.name}</b><br/>{point.key}",
        pointFormat = "Y-value: {point.y}<br/>{point.pub}",
        valueDecimals = 2,
        followPointer = TRUE
      )
    ),
    line = list(
      lineWidth = interval.bar.height,
      cursor = "pointer",
      point = list(
        events = list(
          #click = "#! function() { window.open(this.url); } !#")), #open webpage
          #click = "#! function(event) {alert(this.url);} !#")), #display popup
          #click = "#! function(event) {console.log(this);} !#")), #write object to log
          #click = "#! function(){$('input#selected').val(134); $('input#selected').trigger('change');} !#")),
          #click = "#! function(){$('input#selected').val(this.options.bp); $('input#selected').trigger('change');} !#")),
          #click = "#! function(){$('select#chr').val(this.options.chr); $('select#chr').trigger('change'); $('input#selected').val(this.options.bp); $('input#selected').trigger('change');  $('ul#datatabs li').eq(2).removeClass('active'); $('ul#datatabs li').eq(3).addClass('active'); $('.tab-content div').toggleClass('active'); $('#pChart').trigger('shown')} !#")),
          click = doClickOnLine
        )
      ),
      marker = list(
        enabled = FALSE,
        states = list(hover = list(enabled=FALSE)),
        symbol = "square"
      )
    )            
  )#end plotOptions        
  #c$tooltip(useHTML = T, formatter = "#! function() { return this.point.name; } !#")
  #c$tooltip(formatter = "#! function() { return this.point.name; } !#")
  c$exporting(enabled=TRUE,filename='genomeChart',sourceWidth=2000)
  hideLegend <- input[[jth_ref("legend", j)]]
  if (is.null(hideLegend)) hideLegend <- FALSE
  hideLegend <- hideLegend || !(gwas.data.exist || qtl.data.exist)
  c$legend(enabled = !hideLegend)

  # set c$xAxis here to add macrosynteny plot bands, if any
  bcs <- values[[jth_ref("blockCumStart", j)]]
  bce <- values[[jth_ref("blockCumEnd", j)]]
  nBlocks <- length(bcs)
  if (!is.null(bcs) && nBlocks > 0) {
    for (i in 1:nBlocks) {
      bigList[[nBlocks + i]] <- list(id = paste0("band", i), from = bcs[i], to = bce[i], color = macrosyntenyPlotBandColor)
    }
  }
  c$xAxis(title = list(text = "Chromosome", margin = title.margin), startOnTick = TRUE,
    min = 0, max = sum(as.numeric(chrSize[[values[[jth_ref("organism", j)]]]])),
    endOnTick = FALSE, labels = list(enabled = FALSE), tickWidth = 0, plotBands = bigList)

  removeNotification(nid)

  c$credits(enabled=TRUE)
  c$set(dom = jth_ref('gChart', j))
  return(c)
}
