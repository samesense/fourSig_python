saveSigAsBed <- function(sig.data=NULL, filename=NULL, trackname=NULL, 
                         mask.start=NULL,
                         mask.end=NULL) {
  
  if(is.null(sig.data)) {
    return("Error: You must provide a value for sig.data... Exiting")
  }
  
  if (!is.data.frame(sig.data)) {
    if (is.list(sig.data)) {
      sig.data <- rbind(sig.data$cis, sig.data$trans) 
    } else {
      return("Error: sig.data must either be a list or a data.frame... Exiting")
    }
  }
  
  if(is.null(filename)) {
    return("Error: You must specify a value for filename... Exiting")
  }
  
  sig.data <- cbind(sig.data, score=rep(166, nrow(sig.data)))

  if((!is.null(mask.start)) & (!is.null(mask.end))) {
    print(paste(mask.start, "to", mask.end, "will be masked out"))
    sig.data <- sig.data[((sig.data$end < mask.start) |
                            (sig.data$start > mask.end)),]
  }
  
  if (nrow(sig.data[sig.data$category==3,]) > 0) {
    sig.data[sig.data$category==3,]$score <- 167
    sig.data[sig.data$category==3,]$category <- "254,232,200"
  }
  if (nrow(sig.data[sig.data$category==2,]) > 0) {
    sig.data[sig.data$category==2,]$score <- 500
    sig.data[sig.data$category==2,]$category <- "253,187,132"
  }
  if (nrow(sig.data[sig.data$category==1,]) > 0) {
    sig.data[sig.data$category==1,]$score <- 1000
    sig.data[sig.data$category==1,]$category <- "227,74,51"
  }
    
  output.df <- data.frame(chr=paste("chr", sig.data$chr, sep=""),
                          start=sig.data$start+1,
                          end=sig.data$end+1,
                          id=paste("fourSig_", c(1:nrow(sig.data)), sep=""),
                          score=sig.data$score,
                          blank2=rep(".", nrow(sig.data)),
                          start2=sig.data$start+1,
                          end=sig.data$end+1,
                          color=sig.data$category)
  output.string <- "track type=bed name=\"fourSig\" nextItemButton=on itemRgb=on\n"
  if (!is.null(trackname)) {
    output.string <- paste("track type=bed name=\"", trackname, "\" nextItemButton=on itemRgb=on\n", sep="")
  }
  
  cat(output.string, file=filename)
  write.table(output.df, file=filename, append=TRUE, 
              col.names=FALSE, row.names=FALSE,
              sep="\t", quote=FALSE)
}

plotReads <- function(window.data, sig.data, use.log=FALSE, x.start=NULL, x.end=NULL, y.min=NULL, y.max=NULL, connect.points=FALSE, mask.start=NULL, mask.end=NULL, min.sig.width=1000, draw.reads.first=TRUE, show.fragments=TRUE) {
  
  if (is.null(x.start)) {
    x.start <- 0;
  }
  if (is.null(x.end)) {
    x.end <- max(window.data$x)
  }
  
  window.data <- window.data[window.data$x >= x.start & window.data$x <= x.end,]
  sig.data    <- sig.data[sig.data$end >= x.start & sig.data$start <= x.end,]
  
  if((!is.null(mask.start)) & (!is.null(mask.end))) {
    print(paste(mask.start, "to", mask.end, "will be masked out"))
    window.data[((window.data$x >= mask.start) &
                   (window.data$x <= mask.end)),]$y <- 0
    sig.data <- sig.data[((sig.data$end < mask.start) |
                            (sig.data$start > mask.end)),]
    
    print("NOTE: Significant regions that overlap the masked area will not be displayed")
  }
  
  ##print(window.data)
  ##print(sig.data)
  
  layout(matrix(1:2, ncol = 1), widths = 1, heights = c(2,0.6), respect = FALSE)
  par(mar = c(0, 4.1, 4.1, 2.1))
  
  if (is.null(y.min)) {
    y.min <- 0
  }
  if (is.null(y.max)) {
    y.max <- max(window.data$y)
  }
  
  sig.y.max <- max(window.data$y)
  
  if (use.log) {
    print("Using log10....")
    
    y.max <- max(log(window.data$y))
    
    plot(window.data$x, log10(window.data$y), xlab="", ylab="log10(Reads Per Window)", type="n", xaxt="none", xlim=c(x.start, x.end), ylim=c(y.min, y.max))
  } else {
    plot(window.data$x, window.data$y, xlab="", ylab="Reads Per Window", type="n", xaxt="none", xlim=c(x.start, x.end), ylim=c(y.min, y.max))
  }
  
  if (draw.reads.first) {
    print("Drawing the read counts first, then the category shading")
    if (use.log) {
      if (connect.points) {
        lines(window.data$x, log10(window.data$y), xaxt="n")
      } else {
        lines(window.data$x, log10(window.data$y), type="h", xaxt="n")
      }
    } else {
      if (connect.points) {
        lines(window.data$x, window.data$y, xaxt="n")
      } else {
        lines(window.data$x, window.data$y, type="h", xaxt="n")
      }
    }
  } else {
    print("Drawing the category shading first, then the read counts")
  }
  
  box(col="#F0F0F0")
  
  if(nrow(sig.data[sig.data$category == 1,]) > 0) {
    pad.indices <- which((sig.data[sig.data$category == 1,]$end -
                            sig.data[sig.data$category == 1,]$start)
                         < min.sig.width)
    
    cat.0 <- sig.data[sig.data$category == 1,]
    
    sig.pads <- rep(0, nrow(cat.0))
    sig.pads[pad.indices] <- min.sig.width -
      (sig.data[pad.indices,]$end - sig.data[pad.indices,]$start)
    
    
    max.sig.y <- rep(0, nrow(cat.0))
    
    for (i in 1:nrow(cat.0)) {
      some.windows <- subset(window.data,
                             x >= cat.0[i,]$start & x <= cat.0[i,]$end)
      max.sig.y[i] <- max(some.windows$y)
    }
    
    rect(xleft=cat.0$start,
         xright=(cat.0$end + sig.pads),
         ytop=max.sig.y, ybottom=0, col="#E34A33BB", border=NA)
  }
  
  if(nrow(sig.data[sig.data$category == 2,]) > 0) {
    pad.indices <- which((sig.data[sig.data$category == 2,]$end -
                            sig.data[sig.data$category == 2,]$start)
                         < min.sig.width)
    
    cat.1 <- sig.data[sig.data$category == 2,]
    
    sig.pads <- rep(0, nrow(cat.1))
    sig.pads[pad.indices] <- min.sig.width -
      (sig.data[pad.indices,]$end - sig.data[pad.indices,]$start)
    
    
    max.sig.y <- rep(0, nrow(cat.1))
    
    for (i in 1:nrow(cat.1)) {
      some.windows <- subset(window.data,
                             x >= cat.1[i,]$start & x <= cat.1[i,]$end)
      max.sig.y[i] <- max(some.windows$y)
    }
    
    rect(xleft=cat.1$start,
         xright=(cat.1$end + sig.pads),
         ytop=max.sig.y, ybottom=0, col="#FDBB84BB", border=NA)
  }
  
  if (nrow(sig.data[sig.data$category == 3,]) > 0) {
    pad.indices <- which((sig.data[sig.data$category == 3,]$end -
                            sig.data[sig.data$category == 3,]$start)
                         < min.sig.width)
    
    cat.2 <- sig.data[sig.data$category == 3,]
    
    sig.pads <- rep(0, nrow(cat.2))
    sig.pads[pad.indices] <- min.sig.width -
      (sig.data[pad.indices,]$end - sig.data[pad.indices,]$start)
    
    
    max.sig.y <- rep(0, nrow(cat.2))
    
    for (i in 1:nrow(cat.2)) {
      some.windows <- subset(window.data,
                             x >= cat.2[i,]$start & x <= cat.2[i,]$end)
      max.sig.y[i] <- max(some.windows$y)
    }
    
    rect(xleft=cat.2$start,
         xright=(cat.2$end + sig.pads),
         ytop=max.sig.y, ybottom=0, col="#FEE8C8BB", border=NA)
  }
  
  
  if (!draw.reads.first) {
    if (use.log) {
      if (connect.points) {
        lines(window.data$x, log10(window.data$y), xaxt="n")
      } else {
        lines(window.data$x, log10(window.data$y), type="h", xaxt="n")
      }
    } else {
      if (connect.points) {
        lines(window.data$x, window.data$y, xaxt="n")
      } else {
        lines(window.data$x, window.data$y, type="h", xaxt="n")
      }
    }
  }
  
  if (!show.fragments) {
    legend(x.start, y.max,
           c("Category 1", "Category 2", "Category 3"),
           fill=c("#E34A33", "#FDBB84", "#FEE8C8"), cex=0.75)
    
  } else {
    legend(x.start, y.max,
           c("Category 1", "Category 2", "Category 3", "3C Fragments"),
           fill=c("#E34A33", "#FDBB84", "#FEE8C8", "#888888"), cex=0.75)
    
    par(mar = c(4.1, 4.1, 0, 2.1))
    
    plot(c(1,1), c(-1,1), xlab="", ylab="", type="n", yaxt="none", xaxt="none", xlim=c(x.start, x.end), ylim=c(-1, 2))
    box(col="#F0F0F0")
    
    window.start.pos <- window.data[seq(from=1, to=nrow(window.data), by=2),]$x
    window.start.neg <- window.data[seq(from=2, to=nrow(window.data), by=2),]$x
    
    if (nrow(window.data) %% 2) {
      ## even number of rows in window.data
      
      rect(xleft=window.start.pos[1:(length(window.start.pos)-1)],
           xright=window.start.neg,
           ytop=1, ybottom=0, col="#888888")
      
      rect(xleft=window.start.neg,
           xright=window.start.pos[-1],
           ytop=0, ybottom=-1, col="#888888")
      
    } else {
      ## odd number of rows in window.data
      
      rect(xleft=window.start.pos,
           xright=window.start.neg,
           ytop=1, ybottom=0, col="#888888")
      
      rect(xleft=window.start.neg[1:(length(window.start.neg)-1)],
           xright=window.start.pos[-1],
           ytop=0, ybottom=-1, col="#888888")
      
    }
  }
  
  tick.locs <- round(seq(x.start, x.end, length.out=4))
  axis(side=1, at=tick.locs, labels=tick.locs)
  
}
