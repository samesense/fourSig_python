#########################################################
##
## GLOBAL CONSTANTS
##
#########################################################

TAB.CHR        <- 1
TAB.POS        <- 2
TAB.READ.COUNT <- 3
TAB.MAP        <- 4
TAB.4BP        <- 5
TAB.SHORT      <- 6
TAB.6BP        <- 7
TAB.LENGTH     <- 8

DATA.CHR        <- 1
DATA.POS        <- 2
DATA.READ.COUNT <- 3
DATA.MAP        <- 4
DATA.4BP        <- 5
DATA.SHORT      <- 6
DATA.6BP        <- 7
DATA.LENGTH     <- 8

#########################################################
##
## FUNCTIONS
##
#########################################################

calc.windows <- function(data, window.size) {
  sum.v <- c(0, cumsum(data))

  ## diff, in this case is a vector of
  ## sum.v[(1+n):length(sum.v)] - sum.v[1:(length(sum.v)-n)]
  ## i.e.
  ## sum.v = 0 0 1 1 1 2 2 2 3 3 3
  ##         ^       ^               1 - 0 = 1
  ## sum.v = 0 0 1 1 1 2 2 2 3 3 3
  ##           ^       ^             2 - 0 = 2
  ## sum.v = 0 0 1 1 1 2 2 2 3 3 3
  ##             ^       ^           2 - 1 = 1 etc...
  diff(sum.v, window.size)
}	

saveWindowsToWig <- function(window.data=NULL, filename=NULL, trackname=NULL) {
  if(is.null(window.data)) {
    return("Error: You must provide a value for window.data... Exiting")
  }
  
  if(is.null(filename)) {
    return("Error: You must specify a value for filename... Exiting")
  }
  
  window.data$x <- window.data$x+1
  
  trackline.string <- "track type=wiggle_0 name=\"fourSig\" graphType=bar\n"
  if (!is.null(trackname)) {
    trackline.string <- paste("track type=wiggle_0 name=\"", trackname, "\" graphType=bar\n", sep="")
  }
  variable.string <- "variableStep chrom=chr"
  
  cat(trackline.string, file=filename)
  
  for(chrNum in levels(as.ordered(window.data$chr))) {
    chr.variable.string <- paste(variable.string, chrNum, "\n", sep="")
    cat(chr.variable.string, file=filename, append=TRUE)
    write.table(window.data[window.data$chr == chrNum,2:3], 
                file=filename, append=TRUE, 
                col.names=FALSE, row.names=FALSE,
                sep="\t", quote=FALSE)
  }
}

readsToWindows <- function(filename=NULL, chr, window.size=100,
                           flatten=FALSE, chr.coords=TRUE, mask.start=NULL,
                           mask.end=NULL, only.mappable=TRUE, bed.output=NULL) {

  print(paste("Reading in", filename))
  data <- readTab(filename, onlyMappable=only.mappable) # read in the data...
  print("Done reading the file")

  chromosomes <- NULL
  if (chr == "all") {
    chromosomes <- c(1:19, "X")
  } else {
    chromosomes <- chr
  }

  output.df <- NULL
  for(chr in chromosomes) {
    print(paste("Working on chromosme:", chr))
    
    chr.data <- data[data[DATA.CHR] == chr,] 
    
    if((!is.null(mask.start)) & (!is.null(mask.end))) {
      print(paste("...", mask.start, "to", mask.end, "on chromosome",
                  chr, "will be masked out"))
      chr.data[((chr.data[DATA.POS] >= mask.start) &
                (chr.data[DATA.POS] <= mask.end)), DATA.READ.COUNT] <- 0
    }
    
    chr.reads <- chr.data[, DATA.READ.COUNT]
    if(length(chr.reads) < window.size) {
      print("ERROR: The window size is larger than the number of positions!")
      print("...exiting")
      return()
    }
    
    if (flatten) {
      print ("Flattening the reads...")
      chr.reads <- chr.reads > 0
    }
    
    windows <- calc.windows(data=chr.reads,
                            window.size=window.size)
    
    window.num <- length(windows)
    
    if (chr.coords) {
      window.df <- data.frame(chr=chr,
                              x=chr.data[1:window.num, DATA.POS], y=windows)
      ##return(window.df)
      output.df = rbind(output.df, window.df, deparse.level=0)
    } else {
      window.df <- data.frame(chr=chr, x=c(1:window.num), y=windows)
      ##return(window.df)
      output.df = rbind(output.df, window.df, deparse.level=0)
    }
  }
  return(output.df)
}

readTab <- function( fileName, onlyMappable=TRUE){
  tab.out <- read.table(file=fileName, sep="\t", header=TRUE)
  
  ## now remove all of the fragments that can't be mapped...
  if (onlyMappable) {
    tab.out <- tab.out[tab.out[TAB.MAP] == TRUE,] ## in R TRUE == 1
    gUnMappable <<- tab.out[tab.out[TAB.MAP] == FALSE,] ## in R FALSE == 0
  }
  
  tab.out
}
