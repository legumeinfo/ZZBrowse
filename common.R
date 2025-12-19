# --------------------------------------------------------------
library(rtracklayer)
library(stringi)
library(yaml)
# --------------------------------------------------------------

# Construct the full name of a named list item for species j (j = 1 or 2).
# for example, input$datasets == input[[jth_ref("datasets", 1)]]
# and input$datasets2 == input[[jth_ref("datasets", 2)]]
jth_ref <- function(name, j) {
  ifelse(j == 1, name, paste0(name, j))
}

#add a totalBP column to an input dataset if not already present
calculateTotalBP <- function(j, input, values) {
  if(nrow(values[[input[[jth_ref("datasets", j)]]]]) == 0 || "totalBP" %in% colnames(values[[input[[jth_ref("datasets", j)]]]])){
    
  }else{
    #      progress <- Progress$new(session, min=1, max=1)
    #      on.exit(progress$close())
    
    #      progress$set(message = 'Preparing data for plotting.',
    #                   detail = 'This may take a few seconds...')
    #       withProgress(session, min=1, max=400, expr={
    #         for(i in 1:400) {
    #           setProgress(message = 'Preparing Data for Plotting',
    #                       detail = 'This may take a few seconds...',
    #                       value=i)
    #           #print(i)
    #           Sys.sleep(0.1)
    #         }
    #       })     
    cumBP<-c(0,cumsum(as.numeric(chrSize[[values[[jth_ref("organism", j)]]]])))
    #to order by desired chromosome add factor levels in the desired order to the chrColumn, any chr names that differ in gwas file compared
    #to organism file will turn into NA
    ff <- factor(values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("chrColumn", j)]]],levels=chrName[[values[[jth_ref("organism", j)]]]])
    values[[input[[jth_ref("datasets", j)]]]] <- values[[input[[jth_ref("datasets", j)]]]][order(ff,values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("bpColumn", j)]]]),]
    numeachchr<-aggregate(values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("bpColumn", j)]]],list(ff),length)
    #      adjust<-rep(cumBP[1],numeachchr$x[numeachchr$Group.1==1])            
    adjust <- numeric()
    for (i in 1:(length(cumBP)-1)){#max(unique(ff))){
      if(length(numeachchr$x[numeachchr$Group.1==chrName[[values[[jth_ref("organism", j)]]]][i]])==0){next;}
      adjust<-c(adjust,rep(cumBP[i],numeachchr$x[numeachchr$Group.1==chrName[[values[[jth_ref("organism", j)]]]][i]]))
    }
    #newval <- values[[input[[jth_ref("datasets", j)]]]][600,input[[jth_ref("bpColumn", j)]]]+adjust[600]      
    values[[input[[jth_ref("datasets", j)]]]]$totalBP <- values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("bpColumn", j)]]]+adjust
    
    #values[[input[[jth_ref("datasets", j)]]]] <- adply(values[[input[[jth_ref("datasets", j)]]]],1,function(x){data.frame(totalBP=sum(x[[input[[jth_ref("bpColumn", j)]]]],chrSize$bp[chrSize$chr %in% if(x[[input[[jth_ref("chrColumn", j)]]]]==1) 0 else c(1:(x[[input[[jth_ref("chrColumn", j)]]]]-1))]))})
  }
  if(input[[jth_ref("supportInterval", j)]] == TRUE){
    if("SIbpStartTotal" %in% colnames(values[[input[[jth_ref("datasets", j)]]]])){
      
    }else{
      
      cumBP<-c(0,cumsum(as.numeric(chrSize[[values[[jth_ref("organism", j)]]]])))
      ff <- factor(values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("chrColumn", j)]]],levels=chrName[[values[[jth_ref("organism", j)]]]])
      values[[input[[jth_ref("datasets", j)]]]] <- values[[input[[jth_ref("datasets", j)]]]][order(ff,values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("SIbpStart", j)]]]),]
      numeachchr<-aggregate(values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("SIbpStart", j)]]],list(ff),length)
      adjust <- numeric()
      for (i in 1:(length(cumBP)-1)){#max(unique(ff))){
        if(length(numeachchr$x[numeachchr$Group.1==chrName[[values[[jth_ref("organism", j)]]]][i]])==0){next;}
        adjust<-c(adjust,rep(cumBP[i],numeachchr$x[numeachchr$Group.1==chrName[[values[[jth_ref("organism", j)]]]][i]]))
      }
      values[[input[[jth_ref("datasets", j)]]]]$SIbpStartTotal <- values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("SIbpStart", j)]]]+adjust    
    }
    
    if("SIbpEndTotal" %in% colnames(values[[input[[jth_ref("datasets", j)]]]])){
      
    }else{
      
      cumBP<-c(0,cumsum(as.numeric(chrSize[[values[[jth_ref("organism", j)]]]])))
      ff <- factor(values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("chrColumn", j)]]],levels=chrName[[values[[jth_ref("organism", j)]]]])
      values[[input[[jth_ref("datasets", j)]]]] <- values[[input[[jth_ref("datasets", j)]]]][order(ff,values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("SIbpEnd", j)]]]),]
      numeachchr<-aggregate(values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("SIbpEnd", j)]]],list(ff),length)
      adjust <- numeric()
      for (i in 1:(length(cumBP)-1)){#max(unique(values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("chrColumn", j)]]]))){
        if(length(numeachchr$x[numeachchr$Group.1==chrName[[values[[jth_ref("organism", j)]]]][i]])==0){next;}
        adjust<-c(adjust,rep(cumBP[i],numeachchr$x[numeachchr$Group.1==chrName[[values[[jth_ref("organism", j)]]]][i]]))
      }
      values[[input[[jth_ref("datasets", j)]]]]$SIbpEndTotal <- values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("SIbpEnd", j)]]]+adjust    
    }
  } #end SI total bp calculation
}#end calculateTotalBP

#return color table (which sets colors for series in charts) up to date with all combinations of traits
#this uses allColors set up in global, colors will repeat after 30series
getColorTable <- function(j, input) {
  traitVals <- list()
  if(input[[jth_ref("plotAll", j)]] == FALSE){
    # Consider traits from both datasets,
    # so as to assign any common ones the same color in both charts
    for (jk in 1:2) {
      for (i in input[[jth_ref("traitColumns", jk)]]) {
        traitVals[[i]] <- c(traitVals[[i]], input[[jth_ref(i, jk)]])
      }
    }
    
    traits <- unique(unlist(traitVals))
    if(length(traits)==0){return(NULL)}
    
    colorTable <- data.frame(trait=traits,color=rep(allColors,ceiling(length(traits)/30))[1:length(traits)])
  }else{
    colorTable <- data.frame(trait=input[[jth_ref("datasets", j)]],color=allColors[1])
  }
  colorTable
}

findGWASOverlaps <- function(genomeChart, j, input) {
  if(is.null(input[[jth_ref("overlapSize", j)]])){return(genomeChart[1,])}
  tableIn <- genomeChart
  tableIn$winStart <- tableIn[, input[[jth_ref("bpColumn", j)]]] - input[[jth_ref("overlapSize", j)]]
  tableIn$winStop <- tableIn[, input[[jth_ref("bpColumn", j)]]] + input[[jth_ref("overlapSize", j)]]
  
  allGr <- GRanges(tableIn[,input[[jth_ref("chrColumn", j)]]], IRanges(start=tableIn$winStart,end=tableIn$winStop))
  
  tableIn$group <- subjectHits(findOverlaps(allGr, reduce(allGr)))
  
  #just groups that have more than one unique SNP
  gwasDataOverlap <- tableIn[tableIn$group %in% as.data.frame(table(tableIn$group))[as.data.frame(table(tableIn$group))$Freq>1,"Var1"],]
  
  #just groups that have more than one unique phenotype
  gwasDataOverlapDiffPheno <- ddply(gwasDataOverlap,.(group),function(x){if(nrow(unique(as.data.frame(x[,"trait"])))>=input[[jth_ref("numOverlaps", j)]]){x}else{x[0,]}})
  
  return(gwasDataOverlapDiffPheno)    
}

# Test whether chromosome name is just its number
hasNumericChromosomeNames <- function(organism) {
  is.na(org.chrFormat[organism])
}

# Return the trailing integer value of a string like "phavu.Chr02",
# or NA if it does not exist
trailingInteger <- function(s) {
  n <- as.integer(regmatches(s, regexpr("\\d+$", s)))
  if (length(n) == 0) return(NA)
  n
}

# Return the trailing part ("Chr02") of a chromosome name like "phavu.Chr02"
# (for peanut, extract trailing "Arahy.02")
trailingChromosomeName <- function(s, organism) {
  if (hasNumericChromosomeNames(organism)) {
    chrName <- as.character(trailingInteger(s))
  } else {
    chrName <- stri_match(s, regex = ".+\\.([A-Za-z]+\\.?\\d+)\\s*.*$")[, 2]
  }
  chrName
}

getRainbowColors <- function(n) {
  # Return n rainbow colors (from red to magenta)
  # TODO: a more clearly distinguishable set of colors
  rainbow(n, end = 5/6)
}

# Break the string s into lines at most n characters long,
# replacing the breaks (last whitespace in each line) by <br>
brAt <- function(s, n = 90) {
  gsub(sprintf('(.{1,%d})(\\s|$)', n), '\\1<br>', s)
}

# Return the subset of string vector ss that matches pattern p (case-insensitive)
stringsThatMatchPattern <- function(ss, p) {
  if (is.null(p) || p == "") return(ss)
  ss[!is.na(unlist(
    sapply(ss, function(s) {
      stri_match_all(s, regex = paste0("(?i)", p))
    }, USE.NAMES = FALSE)
  ))]
}

# Uses interval scheduling algorithm to assign the minimum number of rows to the macrosynteny blocks and QTLs.
#   blocks = data frame representing the macrosynteny blocks or QTLs
#   columns = column names (of blocks) containing the start and end base pairs
#   ymax = top of chart
#   h = vertical width of each block
# Output is the chart height for each block = (ymax - j*h) for the jth block.
computeBlockHeights <- function(blocks, columns, ymax, h) {
  bb <- blocks[, columns]
  cmin <- columns[1]
  cmax <- columns[2]
  n <- nrow(bb)
  bb$i0 <- 1:n
  yh <- numeric(n)

  # order by nondecreasing fmax
  bb <- bb[!is.na(bb[, cmax]), ]
  bb <- bb[order(bb[, cmax]), ]
  j <- 1
  while (nrow(bb) > 0) {
    f <- 0
    ii <- c()
    for (i in 1:nrow(bb)) {
      if (f < bb[i, cmin]) {
        ii <- c(ii, i)
        f <- bb[i, cmax]
      }
    }
    yh[bb$i0[ii]] <- ymax - j*h
    bb <- bb[-ii, ]
    j <- j + 1
  }

  yh
}

# Determine best contrast color (black or white) for color, according to its luminance
contrastColor <- function(color, threshold = 0.50) {
  rgb <- col2rgb(color)/255
  rgb <- ifelse(rgb <= 0.03928, rgb/12.92, ((rgb + 0.055)/1.055)^2.4)
  luminance <- as.vector(c(0.2126, 0.7152, 0.0722) %*% rgb) # flatten 1 x 1 matrix
  ifelse(luminance <= threshold, "white", "black")
}

geneNameFromId <- function(geneId) {
  paste0(stri_match_first(geneId, regex = "(\\w+)\\..+ann\\d+(.+)")[, -1], collapse = "")
}

# --------------------------------------------------------------

# LIS Data Store methods

read.gff3.lis.datastore <- function(fin) {
  # download a temporary file as R cannot handle https
  tmp <- tempfile()
  download.file(fin, tmp, method = "wget", quiet = TRUE)
  df.gff <- readGFF(tmp, columns = c("seqid", "start"), tags = "Name", filter = list(type = "genetic_marker"))
  unlink(tmp)
  df.gff
}

scrub.gff <- function(df.gff.in, chrRegex) {
  # Extract chromosome
  chromosome <- sapply(df.gff.in$seqid, function(s) stri_match_first(s, regex = chrRegex)[, 2])
  # Restore internal commas (if any) to marker name
  marker <- sapply(df.gff.in$Name, function(name) paste(name, collapse = ","))
  # output
  df.gff <- data.frame(marker, chromosome, position = df.gff.in$start, row.names = NULL, stringsAsFactors = FALSE)
  df.gff[!is.na(df.gff$chromosome), ]
}

# Read metadata from README file in same directory as input file
read.metadata <- function(fin) {
  directory <- dirname(fin)
  collection <- basename(directory)
  fReadme <- sprintf("%s/README.%s.yml", directory, collection)
  tmp <- tempfile()
  download.file(fReadme, tmp, method = "wget", quiet = TRUE)
  yy <- read_yaml(tmp)
  unlink(tmp)

  metadata <- list()
  metadata$publication <- sprintf("<a href='https://doi.org/%s' title='%s' target=_blank>%s</a>",
    yy$publication_doi, yy$synopsis, yy$identifier)
  metadata
}

# Apply fltr (a valid R logical expression) to dataset
applyGlobalFilter <- function(dataset, fltr) {
  if (fltr != '') {
    seldat <- try(do.call(subset, list(dataset, parse(text = fltr))), silent = TRUE)
    if (!is(seldat, 'try-error')) {
      if (is.data.frame(seldat)) {
        dataset <- seldat
        seldat <- NULL
      }
    }
  }
  dataset
}

# --------------------------------------------------------------
