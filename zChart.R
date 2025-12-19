source("common.R")

provideMultipleURLs <- function(includeGenomicLinkage) {
  paste(
    # From the JSON at this.url, extract the URLs related to this gene.
    # Note that this.url = 'https://services.lis.ncgr.org/gene_linkouts?genes=' + geneId
    # And for now, add the gene family phylogram URL by hand.
    "var geneFamily = this.family;",
    "$.getJSON(this.url, function(data) {",
      "const geneIdRegex = /\\?genes=(.+)$/;",
      "var geneId = geneIdRegex.exec(this.url)[1];",
      "var content = '';",
      "if (data.length == 0) {",
        "content = '<p>No ' + geneId + ' links found.</p>';",
      "} else {",
        "$.each(data, function(i, obj) {",
          "content = content + '<p><a href=' + obj.href + ' target=_blank>' + obj.text + '</a></p>';",
          "if (i == 0) {",
            "var urlPhylogram = 'https://funnotate.legumeinfo.org/?family=' + geneFamily + '&gene_name=' + geneId;",
            "var textPhylogram = 'View LIS gene family phylogram page for : ' + geneId;",
            "content = content + '<p><a href=' + urlPhylogram + ' target=_blank>' + textPhylogram + '</a></p>';",
          "}",
        "});",
      "}",
      # Add the Genomic Linkage button if requested
      ifelse(!includeGenomicLinkage, "", paste(
        "content = content + '<p><button",
          # button actions: first close the dialog, then handle genomic linkages
          "onclick=\"$(this).closest(&quot;.ui-dialog-content&quot;).dialog(&quot;close&quot;);",
          "Shiny.onInputChange(&quot;selectedGene&quot;, &quot;' + geneId + '&quot;);\"",
        ">Genomic Linkage</button></p>';"
      )),
      "var $div = $('<div></div>');",
      "$div.html(content);",
      "$div.dialog({",
        "title: geneId + ' Links',",
        "width: 512,",
        "height: 'auto',",
        "modal: true",
      "});",
    "});"
  )
}

create_zChart <- function(j, input, values) {
  if (!is.numeric(input[[jth_ref("selected", j)]])) return(rCharts::Highcharts$new())

  org.j <- values[[jth_ref("organism", j)]]
  nid <- jth_ref("notify.create_zChart", j)
  showNotification(paste0("Creating Annotation chart for ", org.j, ". Please wait."),
    duration = NULL, id = nid, type = "message")

  centerBP <- as.numeric(input[[jth_ref("selected", j)]])
  winHigh <- centerBP + input[[jth_ref("window", j)]]
  winLow <- centerBP - input[[jth_ref("window", j)]]
  if (winLow < 0) {winLow <- 0}
  
  zoomChart <- values[[input[[jth_ref("datasets", j)]]]]
  zoomChart <- zoomChart[zoomChart[,input[[jth_ref("chrColumn", j)]]]==input[[jth_ref("chr", j)]],]    
  selcom <- trimws(input[[jth_ref("dv_select", j)]])
  if (selcom != "") zoomChart <- applyGlobalFilter(zoomChart, selcom)
  
  if(input[[jth_ref("plotAll", j)]] == FALSE){
    for(i in input[[jth_ref("traitColumns", j)]]){
      zoomChart <- zoomChart[zoomChart[,i] %in% input[[jth_ref(i, j)]],]
    }
    
    if(length(input[[jth_ref("traitColumns", j)]]) > 1){
      zoomChart$trait <- do.call(paste,c(zoomChart[,input[[jth_ref("traitColumns", j)]]],sep="_"))    
    }else{
      zoomChart$trait <- zoomChart[,input[[jth_ref("traitColumns", j)]]]
    }
  }else{
    zoomChart$trait <- input[[jth_ref("datasets", j)]]
  }
  
  #Separate Support Interval data from GWAS data, if support, GWAS data is assumed to be anything that has an NA in the SIbpStart column
  if(input[[jth_ref("supportInterval", j)]] == TRUE){
    SIchart <- zoomChart[!(is.na(zoomChart[,input[[jth_ref("SIbpStart", j)]]])),]
    zoomChart <- zoomChart[is.na(zoomChart[,input[[jth_ref("SIbpStart", j)]]]),]
    if (dynamic.interval.height) {
      # Retain only visible intervals, in order to dynamically assign their height
      SIchart <- SIchart[((SIchart[,input[[jth_ref("SIbpStart", j)]]] <= winHigh & SIchart[,input[[jth_ref("SIbpEnd", j)]]] >= winLow)), ]
    }
  }    
  
  zoomChart <- zoomChart[(zoomChart[,input[[jth_ref("bpColumn", j)]]] <= winHigh) & (zoomChart[,input[[jth_ref("bpColumn", j)]]] >= winLow),]    
  
  #filter for only rows that have a base pair value
  zoomChart <- zoomChart[!(is.na(zoomChart[,input[[jth_ref("bpColumn", j)]]])),]
  zoomChart <- zoomChart[!(is.na(zoomChart[,input[[jth_ref("yAxisColumn", j)]]])),]
  
  #if checked, filter for only overlapping SNPs
  if(!is.null(input[[jth_ref("overlaps", j)]]) & input[[jth_ref("overlaps", j)]] == TRUE){
    zoomChart <- findGWASOverlaps(zoomChart, j, input)
  }                    
  
  if(nrow(zoomChart)==0){ #nothing is in the window, but lets still make a data.frame
    zoomChart <- values[[input[[jth_ref("datasets", j)]]]][1,]
    zoomChart[,input[[jth_ref("yAxisColumn", j)]]] <- -1    
    if(length(input[[jth_ref("traitColumns", j)]]) > 1){
      zoomChart$trait <- do.call(paste,c(zoomChart[,input[[jth_ref("traitColumns", j)]]],sep="_"))
    }else{
      zoomChart$trait <- zoomChart[,input[[jth_ref("traitColumns", j)]]]
    }                   
  }
  colorTable <- getColorTable(j, input)
  
  #take -log10 of y-axis column if requested
  if(input[[jth_ref("logP", j)]] == TRUE && zoomChart[1,input[[jth_ref("yAxisColumn", j)]]] != -1){
    zoomChart[,input[[jth_ref("yAxisColumn", j)]]] <- -log(zoomChart[,input[[jth_ref("yAxisColumn", j)]]],10)
  }                
  
  #check if there is too much data (>2500 data points), trim to 2500
  if(nrow(zoomChart)>2500){
    cutVal <- sort(zoomChart[,input[[jth_ref("yAxisColumn", j)]]],decreasing = T)[2500]
    zoomChart <- zoomChart[zoomChart[,input[[jth_ref("yAxisColumn", j)]]] >= cutVal,]
  }                
  
  zoomChart$publication[is.null(zoomChart$publication)] <- ""
  zoomTable <- data.frame(
    x = zoomChart[, input[[jth_ref("bpColumn", j)]]],
    y = zoomChart[, input[[jth_ref("yAxisColumn", j)]]],
    trait = zoomChart$trait,
    name = sprintf("<table cellpadding='4' style='line-height:1.5'><tr><td align='left'><b>%1$s</b><br>Base Pair: %2$s<br>Chromosome: %3$s<br>Y-value: %4$.2f<br>%5$s</td></tr></table>",
      zoomChart$trait,
      prettyNum(zoomChart[, input[[jth_ref("bpColumn", j)]]], big.mark = ","),
      zoomChart[, input[[jth_ref("chrColumn", j)]]],
      zoomChart[, input[[jth_ref("yAxisColumn", j)]]],
      zoomChart$publication
    ),
    chr = zoomChart[, input[[jth_ref("chrColumn", j)]]],
    bp = zoomChart[, input[[jth_ref("bpColumn", j)]]],
    stringsAsFactors = FALSE
  )
  zoomSeries <- lapply(split(zoomTable, zoomTable$trait), function(x) {
    res <- lapply(split(x, rownames(x)), as.list)
    names(res) <- NULL
    res <- res[order(sapply(res, function(x) x$x))]
    return(res)
  })
  
  #     #build the sliding window GWAS summary line
  #     if(input$SlidingWinCheckbox==TRUE){
  #       winTable <- fullWinTable[fullWinTable$chr==input$chr & fullWinTable$el %in% els & fullWinTable$loc %in% input$locs]
  #     }else{
  #       winTable <- data.frame()
  #     }      
  #     
  #     if(nrow(winTable)==0){ #make a dummy point if there are none in this plot
  #       winTable <- fullWinTable[1,]
  #       winTable$chr <- input$chr
  #       winTable$env <- input$locs[1]
  #       winTable$el <- els[1]
  #       winTable$sumRMIP <- -5
  #     }
  #     
  #     winTable$loc_el <- paste(winTable$loc,winTable$el,"Window",sep="_")
  #     winTableSeries <- adply(winTable,1,function(x) {data.table(x=x$bp,y=x$sumRMIP,element=x$loc_el,
  #                                                                name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>sumRMIP: %2$s<br>Num Points: %3$s<br>Loc: %4$s<br>Base Pair: %5$s<br>cM: %6$s<br>Chromosome: %7$s</td></tr></table>",
  #                                                                             x$el,
  #                                                                             x$sumRMIP,
  #                                                                             x$numPoints,
  #                                                                             x$loc,
  #                                                                             prettyNum(x$bp, big.mark = ","),
  #                                                                             round(x$cM,digits=3),
  #                                                                             x$chr
  #                                                                ),
  #                                                                bp=x$bp,
  #                                                                chr=x$chr                                                    
  #     )})
  #     #this removes columns not needed (using data.table which is why its complicated)
  #     winTableSeries[,colnames(winTableSeries)[which(!colnames(winTableSeries) %in% c("x","y","name","element","bp","chr"))] := NULL,with=FALSE]
  #     #winTableSeries <- winTableSeries[order(winTableSeries$bp,winTableSeries$x),] #should be ordered from the file
  #     #build JL series
  if(input[[jth_ref("supportInterval", j)]]==TRUE){
    if(nrow(SIchart) == 0){ #make a dummy table, but we won't plot the series anyways
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
      SIchart$h <- computeBlockHeights(SIchart, c(input[[jth_ref("SIbpStart", j)]], input[[jth_ref("SIbpEnd", j)]]), chart.max.height, interval.bar.height)
    } else {
      SIchart$h <- SIchart[[input[[jth_ref("SIyAxisColumn", j)]]]]
    }
    SIchart$publication[is.null(SIchart$publication)] <- ""
    SIchart <- SIchart[order(SIchart[[input[[jth_ref("SIbpStart", j)]]]]),]
    jlTable <- adply(SIchart, 1, function(x) {
      data.frame(
        x = c(x[[input[[jth_ref("SIbpStart", j)]]]], x[[input[[jth_ref("SIbpEnd", j)]]]], x[[input[[jth_ref("SIbpEnd", j)]]]]),
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
  } #end if support interval

  #     jl$loc_el <- paste(jl$env,jl$el,"JL",sep="_")      
  #     jlTable <- adply(jl,1,function(x) {data.frame(x=c(x$lowerCIbp,x$upperCIbp,x$upperCIbp),y=c(x$F,x$F,NA),element=x$loc_el,url="http://danforthcenter.org",
  #                                                   name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>F: %2$.2f  -logP: %3$.2f<br>Location: %4$s<br>Base Pair: %5$s<br>SI: %6$s-%7$s<br>SNP: %8$s<br>Chromosome: %9$s</td></tr></table>",       
  #                                                                x$el,
  #                                                                x$F,
  #                                                                x$negLogP,
  #                                                                x$env,
  #                                                                prettyNum(x$bp, big.mark = ","),
  #                                                                prettyNum(x$lowerCIbp, big.mark = ","),
  #                                                                prettyNum(x$upperCIbp, big.mark = ","),
  #                                                                x$Name,
  #                                                                x$Chromosome
  #                                                   ),
  #                                                   bp=x$bp,
  #                                                   chr=x$Chromosome,stringsAsFactors=FALSE
  #     )}
  #     )
  #     jlTable <- jlTable[,c("x","y","name","element","bp","chr","url")]
  #     
  #     jlTable <- jlTable[order(jlTable$bp,jlTable$x),]
  #     #     jlSeries <- lapply(split(jlTable, jlTable$element), function(x) {
  #     #       res <- lapply(split(x, rownames(x)), as.list)
  #     #       names(res) <- NULL
  #     #       #res <- res[order(sapply(res, function(x) x$x))] #
  #     #       return(res)
  #     #     })
  
  #build annotation series
  #thisChrAnnot <- subset(org.annotGeneLoc,chromosome==input[[jth_ref("chr", j)]])
  thisChrAnnot <- subset(org.annotGeneLoc[[org.j]], chromosome == input[[jth_ref("chr", j)]])
  thisAnnot <- thisChrAnnot[thisChrAnnot$transcript_start >= winLow & thisChrAnnot$transcript_end <= winHigh,]
  if(nrow(thisAnnot)==0){ #nothing is in the window, but lets still make a data.frame (actually make it big just to hopefully pick up one row from each strand...)
    thisAnnot <- thisChrAnnot[1:100,]
  }
  thisAnnot <- thisAnnot[order(thisAnnot$transcript_start),]

  annotYvalReverse <- 0.02
  #if(input[[jth_ref("axisLimBool", j)]] == TRUE){annotYvalReverse <- input[[jth_ref("axisMin", j)]] + 0.01}
  annotYvalForward <- annotYvalReverse + 0.04
  annotTable <- adply(thisAnnot[thisAnnot[[org.tag_strand[[org.j]]]] == org.strand_fwd[[org.j]], ], 1, function(x) {
    data.frame(
      x = c(x[[org.tag_start[[org.j]]]], x[[org.tag_end[[org.j]]]], x[[org.tag_end[[org.j]]]]),
      y = c(annotYvalForward, annotYvalForward, NA),
      url = sprintf(org.urlFormat[[org.j]], x[[org.tag_id[[org.j]]]]),
      name = sprintf("<table cellpadding='4' style='line-height:1.5'><tr><td align='left'><b>%1$s</b><br>Location: %2$s-%3$s<br>Chromosome: %4$s, Strand: %5$s<br>%6$s</td></tr></table>",
        x[[org.tag_id[[org.j]]]],
        prettyNum(x[[org.tag_start[[org.j]]]], big.mark = ","),
        prettyNum(x[[org.tag_end[[org.j]]]], big.mark = ","),
        x[[org.tag_chr[[org.j]]]],
        x[[org.tag_strand[[org.j]]]],
        brAt(x[[org.tag_desc[[org.j]]]])
      ),
      gene = x[[org.tag_id[[org.j]]]],
      family = x[["family"]],
      marker = c(NA, "Arrow", NA),
      stringsAsFactors = FALSE
    )
  })
  annotTableReverse <- adply(thisAnnot[thisAnnot[[org.tag_strand[[org.j]]]] == org.strand_rev[[org.j]], ], 1, function(x) {
    data.frame(
      x = c(x[[org.tag_start[[org.j]]]], x[[org.tag_end[[org.j]]]], x[[org.tag_end[[org.j]]]]),
      y = c(annotYvalReverse, annotYvalReverse, NA),
      url = sprintf(org.urlFormat[[org.j]], x[[org.tag_id[[org.j]]]]),
      name = sprintf("<table cellpadding='4' style='line-height:1.5'><tr><td align='left'><b>%1$s</b><br>Location: %2$s-%3$s<br>Chromosome: %4$s, Strand: %5$s<br>%6$s</td></tr></table>",
        x[[org.tag_id[[org.j]]]],
        prettyNum(x[[org.tag_start[[org.j]]]], big.mark = ","),
        prettyNum(x[[org.tag_end[[org.j]]]], big.mark = ","),
        x[[org.tag_chr[[org.j]]]],
        x[[org.tag_strand[[org.j]]]],
        brAt(x[[org.tag_desc[[org.j]]]])
      ),
      gene = x[[org.tag_id[[org.j]]]],
      family = x[["family"]],
      marker = c(NA, "Arrow", NA),
      stringsAsFactors = FALSE
    )
  })
  #annotTable <- annotTable[,c("x","y","name","url","marker")]
  if (nrow(annotTable) == 0) {
    annotTable <- data.frame(x = character(0), y = character(0), name = character(0),
      url = character(0), gene = character(0), family = character(0), stringsAsFactors = FALSE)
  }
  highlight <- (annotTable$gene %in% values$highlightGenes)
  hasHighlightsForward <- (sum(highlight) > 0)
  annotTable <- annotTable[, c("x", "y", "name", "url", "gene", "family")]
  if (hasHighlightsForward) {
    annotTableH <- annotTable[highlight, ]
    annotTable <- annotTable[!highlight, ]
  }
  #annotTable <- annotTable[order(annotTable$x),]
  
  #annotTableReverse <- annotTableReverse[,c("x","y","name","url","marker")]
  if (nrow(annotTableReverse) == 0) {
    annotTableReverse <- data.frame(x = character(0), y = character(0), name = character(0),
      url = character(0), gene = character(0), family = character(0), stringsAsFactors = FALSE)
  }
  highlight <- (annotTableReverse$gene %in% values$highlightGenes)
  hasHighlightsReverse <- (sum(highlight) > 0)
  annotTableReverse <- annotTableReverse[, c("x", "y", "name", "url", "gene", "family")]
  if (hasHighlightsReverse) {
    annotTableReverseH <- annotTableReverse[highlight, ]
    annotTableReverse <- annotTableReverse[!highlight, ]
  }
  #annotTableReverse <- annotTableReverse[order(annotTableReverse$x),]
  
  annotArray <- toJSONArray2(annotTable, json = F, names = T)
  if (hasHighlightsForward) annotArrayH <- toJSONArray2(annotTableH, json = F, names = T)
  #     for(i in 1:length(annotArray)){ #use this to add a symbol before or after the gene track
  #       if(is.na(annotArray[[i]]$marker)){
  #         annotArray[[i]]$marker <- NULL
  #       }else{
  #         annotArray[[i]]$marker <- NULL
  #         annotArray[[i]]$marker$symbol <- "url(./forwardArrow.svg)"
  #       }
  #     }
  
  annotArrayReverse <- toJSONArray2(annotTableReverse, json = F, names = T)
  if (hasHighlightsReverse) annotArrayReverseH <- toJSONArray2(annotTableReverseH, json = F, names = T)
  #     for(i in 1:length(annotArrayReverse)){
  #       if(is.na(annotArrayReverse[[i]]$marker)){
  #         annotArrayReverse[[i]]$marker <- NULL
  #       }else{
  #         annotArrayReverse[[i]]$marker <- NULL
  #         annotArrayReverse[[i]]$marker$symbol <- "url(./reverseArrow.svg)"
  #       }
  #     }
  
  
  b <- rCharts::Highcharts$new()
  b$LIB$url <- 'highcharts/'
  #b$xAxis(title=list(text=zoomTitle), min= -0, max=1,startOnTick = FALSE,reversed=FALSE)
  #b$yAxis(title = list(text = "Base Pairs"),min=winLow,max=winHigh,startOnTick=TRUE,endOnTick=TRUE)
  #b$xAxis(list(title=list(text="RMIP"),min=-0,max=1,startOnTick=FALSE,reversed=FALSE),list(title=list(text="F-score"),min=0,max=1,startOnTick=FALSE,opposite=TRUE,reversed=FALSE))      
  b$chart(zoomType="x",alignTicks=FALSE,events=list(click = "#!function(event) {this.tooltip.hide();}!#"))
  b$xAxis(title = list(text = "Base Pairs"),startOnTick=FALSE,min=winLow,max=winHigh,endOnTick=FALSE)      
  #     if(winTable$sumRMIP[1] != -5){
  #       b$yAxis(list(title=list(text="RMIP"),min=-0,max=1,startOnTick=FALSE),list(title=list(text="F-score"),min=0,startOnTick=FALSE,opposite=TRUE,gridLineWidth=0),
  #               list(title=list(text="sumRMIP"),min=-0,startOnTick=FALSE,opposite=TRUE,gridLineWidth=0))          
  #     }else{
  #      b$yAxis(list(title=list(text="RMIP"),min=-0,max=1,startOnTick=FALSE),list(title=list(text="F-score"),min=0,startOnTick=FALSE,opposite=TRUE,gridLineWidth=0))
  #    }
  if(input[[jth_ref("axisLimBool", j)]] == TRUE){
    b$yAxis(title=list(text=input[[jth_ref("yAxisColumn", j)]]),min=input[[jth_ref("axisMin", j)]],max=input[[jth_ref("axisMax", j)]],startOnTick=FALSE)
    #create a hidden axis to put the gene track on, all the options are setting to hide everything from the axis 
    b$yAxis(labels=list(enabled=FALSE),title=list(text=NULL),min=0,max=1,lineWidth=0,gridLineWidth=0,minorGridLineWidth=0,lineColor="transparent",minorTickLength=0,tickLength=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)
  }else{      
    b$yAxis(title=list(text=input[[jth_ref("yAxisColumn", j)]]),min=0,startOnTick=FALSE)
    #create a hidden axis to put the gene track on, all the options are setting to hide everything from the axis
    b$yAxis(labels=list(enabled=FALSE),title=list(text=NULL),min=0,max=1,lineWidth=0,gridLineWidth=0,minorGridLineWidth=0,lineColor="transparent",minorTickLength=0,tickLength=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)
  }
  
  if(input[[jth_ref("supportInterval", j)]]==TRUE){
    if(input[[jth_ref("SIaxisLimBool", j)]] == TRUE){
      b$yAxis(visible=FALSE,title=list(text=input[[jth_ref("SIyAxisColumn", j)]]),min=input[[jth_ref("SIaxisMin", j)]],max=input[[jth_ref("SIaxisMax", j)]],gridLineWidth=0,minorGridLineWidth=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)
    }else{
      b$yAxis(visible=FALSE,title=list(text=input[[jth_ref("SIyAxisColumn", j)]]),min=0,max=chart.max.height,gridLineWidth=0,minorGridLineWidth=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)
    }
    
    if (SIchart[1, input[[jth_ref("SIyAxisColumn", j)]]] != -1) {
      d_ply(jlTable, .(trait), function(x) {
        b$series(
          data = toJSONArray2(x, json = FALSE, names = TRUE),
          type = "line",
          showInLegend = FALSE,
          dashStyle = 'Solid',
          lineWidth = interval.bar.height,
          name = unique(x$trait),
          yAxis = 2,
          tooltip = list(
            pointFormat = '<span style="color:{point.color}">\u25a0</span> {series.name}<br>{point.pub}',
            followPointer = TRUE
          ),
          color = colorTable$color[colorTable$trait == as.character(unique(x$loc_el))]
        )
      })
    }
  }
  
  #b$series(
  #  data = toJSONArray2(annotTable, json = F, names = F),
  #  type = "columnrange",
  #  pointWidth=15,
  #  pointPadding=0,
  #  pointPlacement=0,
  #  name="Genes"    
  #)    
  
  #invisible(sapply(jlSeries, function(x) {b$series(data = x, type = "line", yAxis=1,name = x[[1]]$element)}))
  #     if(jl$F[1] != -5){
  #       d_ply(jlTable,.(element),function(x){
  #         b$series(
  #           data = toJSONArray2(x,json=F,names=T),
  #           type = "line",
  #           name = unique(x$element),
  #           yAxis=1,
  #           visible = if(input$JLCheckbox==FALSE) FALSE else TRUE,
  #           color = colorTable$color[colorTable$loc_el == as.character(unique(x$element))])})
  #     }
  #     
  #     if(winTable$sumRMIP[1] != -5){      
  #       d_ply(winTableSeries,.(element),function(x){
  #         b$series(
  #           data = toJSONArray2(x,json=F,names=T),
  #           type = "line",
  #           name = unique(x$element),
  #           yAxis = 2,
  #           marker = list(enabled = FALSE),
  #           lineWidth = 3,
  #           turboThreshold=5500,
  #           color = colorTable$color[colorTable$loc_el == as.character(sub("_Window","",unique(x$element)))]
  #         )
  #       })
  #     }

  if(zoomChart[1,input[[jth_ref("yAxisColumn", j)]]] != -1){
    invisible(sapply(zoomSeries, function(x) {
      if (length(x) == 0) {
        return()
      }
      b$series(
        data = x,
        type = "scatter",
        showInLegend = FALSE,
        color = colorTable$color[colorTable$trait == as.character(x[[1]]$trait)],
        name = x[[1]]$trait
      )
    }))
  }

  geneColor <- "#53377A"
  geneColorH <- "#C00080" # between magenta and red, should not conflict with rainbow colors from red to magenta
  b$series(
    data = annotArray,
    type = "line",
    showInLegend = FALSE,
    name = "Forward Genes",
    id = "forward-genes",
    zIndex = 1,
    color = geneColor,
    marker = list(symbol = "circle", enabled = FALSE),
    tooltip = list(pointFormat = '', followPointer = TRUE),
    yAxis = 1
  )
  if (hasHighlightsForward) {
    b$series(
      data = annotArrayH,
      type = "line",
      showInLegend = FALSE,
      name = "Forward Genes H",
      id = "forward-genes-H",
      zIndex = 1,
      color = geneColorH,
      marker = list(symbol = "circle", enabled = FALSE),
      tooltip = list(pointFormat = '', followPointer = TRUE),
      yAxis = 1
    )
  }

  b$series(
    data = annotArrayReverse,
    type = "line",
    showInLegend = FALSE,
    name = "Reverse Genes",
    id = "reverse-genes",
    zIndex = 1,
    color = geneColor,
    marker = list(symbol = "circle", enabled = FALSE),
    tooltip = list(pointFormat = '', followPointer = TRUE),
    yAxis = 1
  )
  if (hasHighlightsReverse) {
    b$series(
      data = annotArrayReverseH,
      type = "line",
      showInLegend = FALSE,
      name = "Reverse Genes H",
      id = "reverse-genes-H",
      zIndex = 1,
      color = geneColorH,
      marker = list(symbol = "circle", enabled = FALSE),
      tooltip = list(pointFormat = '', followPointer = TRUE),
      yAxis = 1
    )
  }

  # Wait for values$glGenes2 to be populated before drawing the zChart
  if (!is.null(values$glGenes2)) {
    apply(values[[jth_ref("glGenes", j)]], 1, FUN = function(g) {
      g <- data.frame(as.list(g), stringsAsFactors = FALSE) # to avoid "$ operator is invalid for atomic vectors" warning
      yh <- -1
      if (g$strand == org.strand_fwd[[org.j]]) {
        yh <- annotYvalForward
        sid <- "forward-genes"
      } else if (g$strand == org.strand_rev[[org.j]]) {
        yh <- annotYvalReverse
        sid <- "reverse-genes"
      }
      if (yh > 0 && !(is.na(g$chromosome)) && g$chromosome == input[[jth_ref("chr", j)]]) {
        g.data <- vector("list", 2)
        g.data[[1]]$x <- as.integer(g$transcript_start)
        g.data[[2]]$x <- as.integer(g$transcript_end)
        g.data[[1]]$y <- g.data[[2]]$y <- yh
        b$series(
          type = "line",
          data = g.data,
          color = g$color,
          linkedTo = sid,
          lineWidth = 12,
          yAxis = 1,
          tooltip = list(
            headerFormat = sprintf("<b>%s</b><br>%s", g$id, g$family),
            pointFormat = '',
            followPointer = TRUE
          ),
          zIndex = 0
        )
      }
    })
  }
  
  #b$series(data = annotSeries[[1]], type = "columnrange", pointWidth=15,pointPadding=0,pointPlacement=0,name = "Genes")    
  b$chart(zoomType="x",alignTicks=FALSE,events=list(click = "#!function(event) {this.tooltip.hide();}!#"))
  #b$title(text=paste("NAM GWAS Results",sep=" "))
  #b$subtitle(text="Rollover for more info. Drag chart area to zoom. Click point for zoomed annotated plot.")
  
  # User clicked on a point -> display trait in popup
  #doClickOnPoint <- "#! function(event) { alert(this.trait); } !#"
  # User clicked on a line -> various possible responses:
  bGenomicLinkage <- ifelse(j == 1 && input$compare2species, 1, 0)
  doClickOnLine <- sprintf(paste(
    "#! function() {",
      "if (this.url.includes('services.lis.ncgr.org')) {",
        provideMultipleURLs(bGenomicLinkage),
      "} else {",
        # for all other species
        "window.open(this.url);", #open webpage
      "}",
    "} !#"
  ))
  
  # Create the gene family legend. For organism 2, show only the gene families on the currently selected chromosome.
  glFamilies <- names(values$glColors)
  if (j == 2 && !is.null(values$glGenes2) && !(is.null(input$relatedRegions) || length(input$relatedRegions) == 0)) {
    glGenes2 <- values$glGenes2
    # parse from the format "chr[Chr] [minBP]-[maxBP] Mbp"
    ss <- unlist(strsplit(input$relatedRegions, split = " "))
    if (startsWith(ss[1], "Gm") || startsWith(ss[1], "Vu")) {
      chr <- ss[1]
    } else {
      chr <- stri_sub(ss[1], 4)
    }
    glFamilies <- intersect(glFamilies, glGenes2$family[glGenes2$chromosome == chr])
  }
  doClickOnColumn <- paste(
    "#! function() {",
      "if ($('input#boolBroadcastToBC').prop('checked')) {",
        # Broadcast the gene family to the Genome Context Viewer
        "Shiny.onInputChange('gcvGeneFamily', this.name);",
      "} else {",
        # Go to the gene family's web page (in Funnotate)
        "window.open('https://funnotate.legumeinfo.org/?family=' + this.name);",
      "}",
      "return false;", # and disable toggling the legend item
    "} !#"
  )
  sapply(glFamilies, FUN = function(f) {
    b$series(
      data = list(),
      type = "column",
      name = f,
      color = values$glColors[[f]],
      events = list(legendItemClick = doClickOnColumn)
    )
  })
  
  b$plotOptions(
    scatter = list(
      cursor = "pointer",
      # Disable clicking on point, as the tooltip is sufficient to identify the trait
      # point = list(
      #   events = list(
      #     #click = "#! function() { window.open(this.options.url); } !#")), #open webpage
      #     click = doClickOnPoint
      #   )
      # ),
      #click = "#! function(event) {console.log(this);} !#")), #write object to log
      #click = "#! function(){$('input#selected').val(134); $('input#selected').trigger('change');} !#")),
      #click = "#! function(){$('input#selected').val(this.options.bp); $('input#selected').trigger('change');} !#")),
      tooltip = list(headerFormat = '', pointFormat = '{point.name}', followPointer = TRUE),
      marker = list(
        symbol = "circle",
        radius = 5
      )
    ),
    line = list(
      lineWidth = 6,
      cursor = "pointer",
      #stickyTracking=FALSE,
      point = list(
        events = list(
          click = doClickOnLine
        )
      ),
      marker = list(symbol = "square")
    ),
    #click = "#! function(event) {alert(this.url);} !#")), #display popup
    #click = "#! function(event) {console.log(this);} !#")), #write object to log
    #click = "#! function(){$('input#selected').val(134); $('input#selected').trigger('change');} !#")),
    #click = "#! function(){$('input#selected').val(this.options.bp); $('input#selected').trigger('change');} !#")),
    marker = list(
      enabled = FALSE,
      radius = 1,
      #symbol = "url(./sun.png)",
      states = list(hover = list(enabled=FALSE))
    )
  )

  removeNotification(nid)

  #it seems almost impossible to get the tooltip to hover along the chart with this version of highcharts (4.0.1), perhaps a question to stackoverflow could solve it.
  #see an example of the problem here: http://jsfiddle.net/N5ymb/
  #one hack/fix would be to add dummy points to the middle of the line that show up when moused over
  #b$tooltip(snap=5, useHTML = T, formatter = "#! function() { return this.point.name; } !#") #followTouchMove = T, shared=T, followPointer = T
  b$exporting(enabled=TRUE,filename='zoomChart',sourceWidth=2000)
  b$credits(enabled=TRUE)
  b$set(dom = jth_ref('zChart', j))
  return(b)
}
