#Ruping Sun, regularhand@gmail.com

library(KernSmooth)
library(fields)

.smoothScatterCalcDensity1 <- function(x, nbin, bandwidth, range.x) {
  
  if (length(nbin) == 1)
    nbin <- c(nbin, nbin)
  if (!is.numeric(nbin) || (length(nbin)!=2))
    stop("'nbin' must be numeric of length 1 or 2")

  if (missing(bandwidth)) {
    bandwidth <- diff(apply(x, 2, quantile, probs=c(0.05, 0.95), na.rm=TRUE)) / 25
  } else {
    if(!is.numeric(bandwidth))
      stop("'bandwidth' must be numeric")
  }
  ## create density map
  if(missing(range.x))
     rv <- bkde2D(x, gridsize=nbin, bandwidth=bandwidth)
  else
     rv <- bkde2D(x, gridsize=nbin, bandwidth=bandwidth, range.x=range.x) 
  rv$bandwidth <- bandwidth
  return(rv)
}

image.plot2 = function (..., add = FALSE, nlevel = 64, legend.shrink = 0.9, 
    legend.width = 1.2, legend.mar = NULL, graphics.reset = FALSE, 
    horizontal = FALSE, bigplot = NULL, smallplot = NULL, legend.only = FALSE, 
    col = tim.colors(nlevel)) 
{
    old.par <- par(no.readonly = TRUE)
    info <- image.plot.info(...)
    if (add) {
        big.plot <- old.par$plt
    }
    if (legend.only) {
        graphics.reset <- TRUE
    }
    if (is.null(legend.mar)) {
        legend.mar <- ifelse(horizontal, 3.1, 5.1)
    }
    temp <- image.plot.plt(add = add, legend.shrink = legend.shrink, 
        legend.width = legend.width, legend.mar = legend.mar, 
        horizontal = horizontal, bigplot = bigplot, smallplot = smallplot)
    smallplot <- temp$smallplot
    bigplot <- temp$bigplot
    if (!legend.only) {
        if (!add) {
            par(plt = bigplot)
        }
	##par(bty = 'n')
        image(..., add = add, col = col)
        box()
	
        big.par <- par(no.readonly = TRUE)	
    }
    if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
        par(old.par)
        stop("plot region too small to add legend\n")
    }
    ix <- 1
    minz <- info$zlim[1]
    maxz <- info$zlim[2]
    binwidth <- (maxz - minz)/nlevel
    midpoints <- seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
    iy <- midpoints
    iz <- matrix(iy, nrow = 1, ncol = length(iy))
    breaks <- list(...)$breaks
    par(new = TRUE, pty = "m", plt = smallplot, err = -1)
    if (!horizontal) {                         #for the scale bar
        if (is.null(breaks)) {
            image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
                ylab = "", col = col)
            axis(4, mgp = c(3, 1, 0), las = 2)
        }
        else {
            image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
                ylab = "", col = col, breaks = breaks)
            axis(4, at = breaks, labels = format(breaks), mgp = c(3, 
                1, 0), las = 2)
        }
    }
    else {
        if (is.null(breaks)) {
            image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
                ylab = "", col = col)
            axis(1, mgp = c(3, 1, 0))
        }
        else {
            image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
                ylab = "", col = col, breaks = breaks)
            axis(1, at = breaks, labels = format(breaks), mgp = c(3, 
                1, 0))
        }
    }
    box()
    mfg.save <- par()$mfg
    if (graphics.reset | add) {
        par(old.par)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
    else {
        par(big.par)
        par(plt = big.par$plt, xpd = TRUE)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
}


smkey <- function(x, y=NULL, 
                          nbin=128,
                          bandwidth,
                          colramp=colorRampPalette(c("white", brewer.pal(9, "Blues"))),
                          nrpoints=100,
                          transformation=function(x) x^.25,
                          xlab=NULL, ylab=NULL, postPlotHook=box,
                          pch=".", cex=1,
                          xlim, ylim, col="black",
                          hline=NULL, vline=NULL,
                          xaxs=par("xaxs"), yaxs=par("yaxs"), ...) {
  
  if (!is.numeric(nrpoints) | (nrpoints<0) | (length(nrpoints)!=1) )
    stop("'nrpoints' should be numeric scalar with value >= 0.")

  ## similar as in plot.default
  xlabel <- if (!missing(x)) 
    deparse(substitute(x))
  ylabel <- if (!missing(y)) 
    deparse(substitute(y))
  xy <- xy.coords(x, y, xlabel, ylabel)
  xlab <- if (is.null(xlab)) 
    xy$xlab
  else xlab
  ylab <- if (is.null(ylab)) 
    xy$ylab
  else ylab


  ## eliminate NA
  x <- cbind(xy$x, xy$y)[!(is.na(xy$x)|is.na(xy$y)), ]

  ## xlim and ylim
  if(!missing(xlim)) {
    stopifnot(is.numeric(xlim), length(xlim)==2, !any(is.na(xlim)))
    x <- x[ (x[,1]>=xlim[1]) & (x[,1]<=xlim[2]), ]
  } else {
    xlim <- range(x[,1], na.rm=TRUE)
  }
  if(!missing(ylim)) {
    stopifnot(is.numeric(ylim), length(ylim)==2, !any(is.na(ylim)))
    x <- x[ (x[,2]>=ylim[1]) & (x[,2]<=ylim[2]), ]
  } else {
    ylim <- range(x[,2], na.rm=TRUE)
  }
  
  ## create density map
  map  <- .smoothScatterCalcDensity1(x, nbin, bandwidth)
  xm   <- map$x1
  ym   <- map$x2
  dens <- map$fhat
  dens <- array(transformation(dens), dim=dim(dens))
  	  
  ## plot color image
  image.plot2(xm, ym, z=dens, legend.shrink = 1.0,
        xlab = xlab, ylab = ylab, nlevel = 256,...)
  if(!is.null(postPlotHook)) postPlotHook()
  if(!is.null(hline)){segments(x0=xlim[1],y0=hline,x1=xlim[2],y1=hline,lty=2)}
  if(!is.null(vline)){segments(x0=vline,y0=ylim[1],x1=vline,y1=ylim[2],lty=2)}
  
  ## plot selection of dots
  if (nrpoints!=0){
    ## we assume that map$x1 and map$x2 go linearly from
    ## their first to their last value in nbin steps
    stopifnot(length(xm)==nrow(dens), length(ym)==ncol(dens))
    ixm <- round((x[,1]-xm[1])/(xm[length(xm)]-xm[1])*(length(xm)-1))
    iym <- round((x[,2]-ym[1])/(ym[length(ym)]-ym[1])*(length(ym)-1))
    idens <- dens[1 + iym*length(xm) + ixm]
    nrpoints <- min(nrow(x), ceiling(nrpoints))
    sel <- order(idens, decreasing=FALSE)[1:nrpoints]
    points(x[sel,1:2], pch=pch, cex=cex, col=col)
  }
}


scatterDensityPlot <- function(x, y, xlim=c(0,1), ylim=c(0,1), div=0.02, xlab="x", ylab="y", main="Scatter Density", cex=1.5, cex.axis=1.5, cex.lab=1.5, cex.main=1.7, abline=TRUE, pch=19,
                               drx=c(), dry=c(), drlabels=c(), denscolor=vector(), groups=list(), groupColors=list(), colScaleLabel="# sSNV", xaxisat=vector(), xaxislb=vector(), yaxisat=vector(), yaxislb=vector(),
                               legend=c("Public","Pvt-Shared","Pvt-Rgn Specific"), legendCol=c(rgb(0,0,0,1/4),rgb(178/255,223/255,138/255,1),rgb(31/255,120/255,180/255,1)), layout=TRUE, alpha=1, box=TRUE) {

    if (length(groups) > 0) {
        colLegends = list()
        if (layout == TRUE) {
            layout(matrix(rep(1:2, c(7,1)), nrow=1))
            par(mar=c(4.5,5,5.5,0))
        }
        else {
            par(mar=c(4.5,5,5,3))
        }
        for (i in 1:length(groups)) {                                  #plot dense scatter for each group
            colpanel = groupColors[[i]]
            indexes = groups[[i]]
            nbinx = round((max(x[indexes]) - min(x[indexes]))/div)
            nbiny = round((max(y[indexes]) - min(y[indexes]))/div)
            #if (nbinx < 20) {
            #    nbinx = 20
            #}
            #if (nbiny < 20) {
            #    nbiny = 20
            #}
            message(nbinx)
            message(nbiny)
            denscolor = densCols(x[indexes], y[indexes], colramp = colorRampPalette(colpanel), nbin=c(nbinx,nbiny))
            denscolor = add.alpha(denscolor, alpha)
            
            dd = grDevices:::.smoothScatterCalcDensity(cbind(x[indexes], y[indexes]), nbin=c(nbinx,nbiny))
            dens <- as.numeric(dd$fhat)
            dens <- dens[dens>0]
            ktotal = length(x[indexes])
            dens = ktotal/sum(dens) * dens

            colLegend <- data.frame(density=ceiling(seq(min(dens), max(dens), len=5)+0.8),0,
                                    color=I(colorRampPalette(colpanel)(5)))
            colLegend = subset(colLegend,!duplicated(colLegend$density))
            colLegends[[i]] = colLegend
            if (i == 1) {
                if (length(xaxisat) > 0) {
                    plot(x[indexes],y[indexes],col=denscolor, bg=denscolor, pch=pch, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, main=main, cex=cex, cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main, axes = F)
                    axis(side=1, at=xaxisat, labels=xaxislb, cex.axis=cex.axis)
                    axis(side=2, at=yaxisat, labels=yaxislb, cex.axis=cex.axis)
                    if (box){
                        box("plot")
                    }
                } else {
                    message("b")
                    plot(x[indexes],y[indexes],col=denscolor, bg=denscolor, pch=pch, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, main=main, cex=cex, cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main)
                }
            } else {
                points(x[indexes],y[indexes],col=denscolor, bg=denscolor, pch=pch, cex=cex)
            }
        }
        if (abline == TRUE) {
            abline(0,1,lty=3,col=rgb(0,0,0,2/3))
        }
        if (length(drlabels) > 0) {
            points(drx,dry,col=rgb(238/255,59/255,59/255,4/5), pch=1, cex=cex+.3)
            if (layout == TRUE){
                pointLabel(drx, dry, labels=drlabels, col=rgb(238/255,59/255,59/255,1), cex=cex.axis)
            }
        }
        if (layout == TRUE) {
            nbottomright = length(which(x > 0.6 & y < 0.2))
            ntopleft = length(which(x < 0.4 & y > 0.8))
            lposition = "topleft"
            if (ntopleft > nbottomright) {
                lposition = "bottomright"
            }
            legend(lposition, legend=legend, col=legendCol, pch=19, bty="n", cex=cex.axis-.3)
        }
        if (layout == TRUE) {
            par(mar=c(4,1,4,0))
            for (i in rev(1:length(colLegends))) {
                si = 5*(length(colLegends)-i)+1
                ri = dim(colLegends[[i]])[1]
                if (i == length(colLegends)) {
                    plot(NA, xlim=c(0,5), ylim=c(0,5*(length(colLegends)+1)), type="n", ann=FALSE, axes=FALSE)
                    rect(0, si:(si+ri-1), 1, (si+1):(si+ri), border=NA, col=colLegends[[i]]$col)
                    text(2, si:(si+ri-1)+0.5, signif(colLegends[[i]]$density, 2), adj=0, cex=cex.axis)
                } else {                
                    rect(0, si:(si+ri-1), 1, (si+1):(si+ri), border=NA, col=colLegends[[i]]$col)
                    text(2, si:(si+ri-1)+0.5, signif(colLegends[[i]]$density, 2), adj=0, cex=cex.axis)
                }
            }
            text(1, 5*(length(colLegends)+1)-2, label=colScaleLabel, srt=90, cex=cex.axis)
        }
    } else {

        colpanel = rev(brewer.pal(11, "RdYlBu"))[2:11]    #RdYlBu as default
        
        #determine dense color
        if (length(denscolor) == 0){
            denscolor = densCols(x, y, colramp = colorRampPalette(colpanel))
        }

        
        #determine density range
        dd = grDevices:::.smoothScatterCalcDensity(cbind(x, y), nbin=128)
        dens <- as.numeric(dd$fhat)
        dens <- dens[dens>0]
        ktotal = length(x)
        dens = ktotal/sum(dens) * dens

        colLegend <- data.frame(density=ceiling(seq(min(dens), max(dens), len=7)),0,
                                color=I(colorRampPalette(colpanel)(7)))
        layout(matrix(rep(1:2, c(7,1)), nrow=1))
        par(mar=c(4.5,5,5.5,0))
        plot(x,y,col=denscolor, pch=19, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, main=main, cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main)
        abline(0,1,lty=3)
        if (length(drlabels) > 0) {
            pointLabel(drx, dry, labels=drlabels, col=rgb(1,0,0,1/4), cex=cex.axis)
        }
        
        par(mar=c(4,1,4,0))
        plot(NA, xlim=c(0,7), ylim=c(0,8), type="n", ann=FALSE, axes=FALSE)
        rect(0, 0:6, 1, 1:7, border=NA, col=colLegend$col)
        text(2, (0:6)+0.5, signif(colLegend$density, 2), adj=0, cex=cex.axis)
        #text(1, -.1, label="Density", srt=90, cex=cex.axis)
    }
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}




## Add an alpha value to a colour
add.alpha <- function(col, alpha=1) {
    if(missing(col))
        stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2, 
          function(x)
              rgb(x[1], x[2], x[3], alpha=alpha))  
}
