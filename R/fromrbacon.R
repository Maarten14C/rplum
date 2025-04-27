# tmp file to repair functions with bugs in rbacon

# temporarily adding PlotSuppPost to do unlist(set$ps)

# temporarily addind PlotLogPost to deal with axis limits (removing Inf values)

# temporarily adding draw.pbmodelled, owing to changes in how images are plotted (implemented in rbacon version 3.4.1 but this has problems when BCAD=TRUE)

# temporarily adding rbacon's agedepth, to accommodate set <- draw.pbmodelled() (a post-3.4.1 change) (adding rbacon::: where needed)

# Time series of the log of the posterior
PlotLogPost <- function(set, from=0, to=set$Tr, xaxs="i", yaxs="i", panel.size=.9, col=grey(0.4)) {
  y <- -set$Us[(from+1):to]
  y[is.infinite(y)] <- NA
  plot(from:(to-1), y, type="l",
    ylab="Log of Objective", xlab="Iteration", main="", xaxs=xaxs, yaxs=yaxs, col=col, cex.axis=panel.size)
}



# plot the Supported data (for plum)
PlotSuppPost_repaired <- function(set=get('info'), xaxs="i", yaxs="i", legend=TRUE, supp.xlim=c(), supp.ylim=c(), yaxt="n", prior.size=.9, panel.size=.9, line.col=3, line.width=2, text.col=2, hist.col=grey(0.8), hist.border=grey(0.4), data.col=rgb(.5,0,.5,.5)) {
  lab <- ifelse(set$Bqkg, "Bq/kg", "dpm/g")

  if(length(supp.xlim) == 0)
    supp.xlim <- c(min( set$detsPlum[,4]), max(set$detsPlum[,4]))

  post.mn <- mean(unlist(set$ps)) # added unlist April 2025
  post.shape <- post.mn^2 / var(unlist(set$ps))

  if(set$nPs > 1) {
    rng <- array(NA, dim=c(set$nPs, 22)) # 22 is the number of segments to draw. Always???
    for(i in 1:set$nPs) {
      rng[i,1:21] <- quantile( set$ps[,i] , seq(0,2,0.1)/2)
      rng[i,22] <- mean( set$ps[,i] )
    }

    if(length(supp.ylim) == 0)
      supp.ylim <- c(min( rng[,1]), max(rng[,21]))

    plot(0, type="n", ylim=supp.ylim, xlim=supp.xlim, main="", xlab="Depth (cm)", ylab=lab, yaxt=yaxt, cex.axis=panel.size)
    n = 21
    colorby = 1.0 / (n/2)
    for(i in 1:(n/2)) {
      segments(set$detsPlum[,4], rng[,i], set$detsPlum[,4], rng[,(i+1)], grey(1.0-colorby*i), lwd=3)
      segments(set$detsPlum[,4], rng[,n-i], set$detsPlum[,4], rng[,n-(i-1)], grey(1.0-colorby*i), lwd=3)
    }
    lines(set$detsPlum[,4], rng[,22], col="red", lty=12) # mean

  } else {
    post <- density(set$ps)
    suppdata <- set$supportedData[,1:2]
    rng <- range(post$x, suppdata[,1]-suppdata[,2], suppdata[,1]+suppdata[,2])

    if(length(supp.ylim) == 0)
      supp.ylim <- c(0, max(post$y))

    plot(post, type="n", xlab=lab, xlim=rng, main="", ylab="", yaxt=yaxt, cex.axis=panel.size)
    polygon(post, col=hist.col, border=hist.border)
    lines(seq(min(set$ps),max(set$ps),.05), dgamma(seq(min(set$ps), max(set$ps), .05), shape=set$s.shape, scale=set$s.mean/set$s.shape), col=line.col, lwd=line.width)

    # plot the measurements of supported
    # extend the plot's horizontal axis, range(...)
    y <- seq(0, max(post$y), length=nrow(suppdata)+50)[1+(1:nrow(suppdata))] # at bottom graph
    segments(suppdata[,1]-suppdata[,2], y, suppdata[,1]+suppdata[,2], y, col=data.col)
    points(suppdata[,1], y, pch=20, col=data.col)
  }
  txt <- paste0("supported", "\ns.shape: ", set$s.shape, "\ns.mean: ", set$s.mean)  

  if(legend)
    legend("topright", txt, bty="n", cex=prior.size, text.col=text.col, adj=c(0,.2))
  
  invisible(c(post.mn, post.shape))
}



#' @name draw.pbmodelled
#' @title Plot the 210Pb data
#' @description Produce a plot of the 210Pb data and their depths
#' @details This function is generally called internally to produce the age-depth graph.
#' It can be used to produce custom-built graphs.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param BCAD The calendar scale of graphs is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param rotate.axes The default of plotting age on the horizontal axis and event probability on the vertical one can be changed with \code{rotate.axes=TRUE}.
#' @param rev.d The direction of the depth axis can be reversed from the default (\code{rev.d=TRUE}).
#' @param rev.age The direction of the calendar age axis can be reversed from the default (\code{rev.age=TRUE})
#' @param pb.lim Minimum and maximum of the 210Pb axis ranges, calculated automatically by default (\code{pb.lim=c()}).
#' @param d.lim Minimum and maximum depths to plot; calculated automatically by default (\code{d.lim=c()}).
#' @param d.lab The labels for the depth axis. Default \code{d.lab="Depth (cm)"}.
#' @param pb.lab The label for the 210Pb axis (default \code{pb.lab="210Pb (Bq/kg)"} or \code{"210Pb (dpm/g)"}).
#' @param supp.col Colour of the supported 210Pb data. Defaults to red: \code{supp.col="red"}.
#' @param pbmodelled.col Colour of the modelled 210Pb values. Defaults to scales of blue: \code{pbmodelled.col=function(x) rgb(0,0,1,x)}.
#' @param pbmeasured.col Colour of the measured 210Pb values. Defaults to blue.
#' @param plot.measured Plot the measured 210Pb values (default \code{plot.measured=TRUE}).
#' @param age.lim values of the age axis. Used to calculate where to plot the pb values on the secondary axis
#' @param mgp Axis text margins (where should titles, labels and tick marks be plotted). Defaults to \code{mgp=c(1.7, .7, .0)}.
#' @param pb.lty Line type of measured Pb-210 data.
#' @author Maarten Blaauw, J. Andres Christen, Marco Aquino-Lopez
#' @return A plot of the modelled (and optionally the measured) 210Pb values
#' @export
draw.pbmodelled <- function(set=get('info'), BCAD=set$BCAD, rotate.axes=FALSE, rev.d=FALSE, rev.age=FALSE, pb.lim=c(), d.lim=c(), d.lab=c(), pb.lab=c(), pbmodelled.col=function(x) rgb(0,0,1,x), pbmeasured.col="blue", supp.col="purple", plot.measured=TRUE, age.lim=c(), mgp=mgp, pb.lty=1) {

  pb <- set$dets[set$dets[,9] == 5,]
  depths <- pb[,4] # set$detsOrig[,2]
  dns <- pb[,6] # set$detsOrig[,3]
  Pb <- pb[,2] # set$detsOrig[,4]
  err <- pb[,3] # set$detsOrig[,5]
  thickness <- pb[,5] # set$detsOrig[,6]
  n <- nrow(pb)

  if(ncol(pb) > 6) {
    supp <- pb[,7] # set$detsOrig[,7]
    supperr <- pb[,8] # set$detsOrig[,8]
  } else {
    supp <- set$supportedData[,1]
    supperr <- set$supportedData[,2]
    suppd <- set$supportedData[,3]
    suppthick <- set$supportedData[,4]
  }

  if(length(d.lab) == 0)
    d.lab <- paste0("depth (", set$depth.unit, ")")
  if(length(pb.lab) == 0)
    pb.lab <- ifelse(set$Bqkg,
      expression(""^210*"Pb (Bq/kg)"),
        expression(""^210*"Pb (dpm/g)"))

  if(length(d.lim) == 0)
    d.lim <- range(depths)
  if(rev.d)
    d.lim <- d.lim[2:1]

  if(length(set$phi) > 0) {
    Ai <- list(x=NULL, y=NULL)
    hght <- 0; pbmin <- c(); pbmax <- 0
    A.rng <- array(0, dim=c(n,2))
    for(i in 1:length(depths)) {
      A <- A.modelled(depths[i]-thickness[i], depths[i], dns[i], set)
      tmp <- density(A)
      Ai$x[[i]] <- tmp$x
      Ai$y[[i]] <- tmp$y
      hght <- max(hght, Ai$y[[i]])
      pbmin <- min(pbmin, Ai$x[[i]])
      pbmax <- max(pbmax, Ai$x[[i]])
      A.rng[i,] <- quantile(A, c((1-set$prob)/2, 1-(1-set$prob)/2))
    }

    if(length(pb.lim) == 0)
      pb.lim <- extendrange(c(0, Pb-2*err, Pb+2*err, pbmax), f=c(0, 0.05))

    pb2bp <- function(pb, pb.min=pb.lim[1], pb.max=pb.lim[2], agemin=min(age.lim), agemax=max(age.lim), AD=BCAD) {
      if(AD) {
        ex <- (agemin - agemax) / (pb.max - pb.min)
        return(agemax + ex*pb)
      } else {
          ex <- (agemax - agemin) / (pb.max - pb.min)
          return(agemin + ex*pb)
      }
    }

    # save the values for later
    set$Ai <- Ai
    set$A.rng <- A.rng

    this <- ifelse(rotate.axes, 3, 4)
    pretty.pb <- pretty(pb.lim)
    onbp <- pb2bp(pretty.pb)
    if(BCAD)
      axis(this, rev(onbp), rev(pretty.pb), col=pbmeasured.col, col.axis=pbmeasured.col, col.lab=pbmeasured.col) else
        axis(this, onbp, pretty.pb, col=pbmeasured.col, col.axis=pbmeasured.col, col.lab=pbmeasured.col)

    mtext(pb.lab, this, 2.5, col=pbmeasured.col, cex=.8)

    for(i in 1:length(depths)) {
      ages <- pb2bp(Ai$x[[i]], AD=BCAD)
      if(BCAD) 
		ages <- rev(ages)
      z <- matrix(Ai$y[[i]]/hght, nrow=1)
      d_slice <- c(depths[i]-thickness[i], depths[i])
      modelled_col <- pbmodelled.col(seq(0, 1-max(z), length=50))
	  
	  # we're not using ghost.mirror, because of BC/AD axis reversal issues
      if(rotate.axes) {
        if(BCAD)	  
          image(ages, d_slice, z, add=TRUE, col=modelled_col, useRaster=TRUE) else
            image(ages, d_slice, z, add=TRUE, col=modelled_col, useRaster=TRUE)
        } else {
		d_slice <<- d_slice; ages <<- ages; z <<- z	
	     if(BCAD)  
  	       image(d_slice, ages, z, add=TRUE, col=modelled_col, useRaster=TRUE) else
	         image(d_slice, ages, z, add=TRUE, col=modelled_col, useRaster=TRUE)
        }	  
    }
  }

  if(plot.measured)
    draw.pbmeasured(set=set, newplot=FALSE, rotate.axes=rotate.axes, BCAD=BCAD, on.agescale=TRUE, pb.lim=pb.lim, age.lim=age.lim, supp.col=supp.col)	
  
  invisible(set)
}



#' @title Plot an age-depth model
#' @description Plot the age-depth model of a core.
#' @details After loading a previous run, or after running either the scissors or thinner command, plot the age-model
#' again using the command \code{agedepth()}.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param d.lab The labels for the depth axis. Default \code{d.lab="Depth (cm)"}. See also \code{depth.unit}.
#' @param age.lab The labels for the calendar axis (default \code{age.lab="cal BP"} or \code{"BC/AD"} if \code{BCAD=TRUE}).
#' @param kcal Use kcal BP. Default is \code{kcal=FALSE}.
#' @param yr.lab Deprecated - use age.lab instead
#' @param acc.lab The labels for the accumulation rate plot (top middle). Default \code{d.lab="Acc. rate (yr/cm)"} (or whatever units you're using).
#' @param mem.lab The labels for the memory plot (top right). Default \code{d.lab="Memory"}.
#' @param d.min Minimum depth of age-depth model (use this to extrapolate to depths higher than the top dated depth). 
#' @param d.max Maximum depth of age-depth model (use this to extrapolate to depths below the bottom dated depth).
#' @param d.by Depth intervals at which ages are calculated. Default 1. Alternative depth intervals can be provided using, e.g., d.\code{d.by=0.5}.
#' @param depths By default, Bacon will calculate the ages for the depths \code{d.min} to \code{d.max} in steps of \code{d.by}. Alternative depths can be provided as, e.g., \code{depths=seq(0, 100, length=500)} or as a file, e.g., \code{depths=read.table("CoreDepths.txt"}. See also \code{depths.file}.
#' @param depths.file By default, Bacon will calculate the ages for the depths \code{d.min} to \code{d.max} in steps of \code{d.by}.
#' If \code{depths.file=TRUE}, Bacon will read a file containing the depths for which you require ages.
#' This file, containing the depths in a single column without a header, should be stored within \code{coredir},
#' and its name should start with the core's name and end with '_depths.txt'. Then specify \code{depths.file=TRUE} (default \code{FALSE}). See also \code{depths}.
#' @param accordion An experimental option to squeeze and stretch depths below a boundary. Needs 2 parameters: boundary depth and compression ratio, e.g., \code{accordion=c(25,20)}. Defaults to inactive, \code{accordion=c()}.
#' @param plotatthesedepths An option to plot ages at different depths (e.g., if a core has been compressed during a run). Use with extreme caution!
#' @param age.min Minimum age of the age-depth plot.
#' @param yr.min Deprecated - use age.min instead.
#' @param age.max Maximum age of the age-depth plot.
#' @param yr.max Deprecated - use age.min instead.
#' @param hiatus.option How to calculate accumulation rates and ages for sections with hiatuses. Either extrapolate from surrounding sections (default, \code{hiatus.option=1}), use a w-weighted mix between the prior and posterior values for depths below the hiatus and prior information only for above the hiatus (\code{hiatus.option=2}), or use the originally calculated slopes (\code{hiatus.option=0}).
#' @param dark Darkness of the greyscale age-depth model. By default, the darkest grey value is calculated as 10 times the height of the lowest-precision age estimate \code{dark=c()}. Lower values will result in lighter grey but values >1 are not allowed.
#' @param prob Confidence interval to report (between 0 and 1, default 0.95 or 95\%).
#' @param rounded Rounding of years. Default is to round to single years (1 digit for plum models).
#' @param d.res Resolution or amount of greyscale pixels to cover the depth scale of the age-model plot. Default \code{d.res=200}.
#' @param age.res Resolution or amount of greyscale pixels to cover the age scale of the age-model plot. Default \code{yr.res=200}.
#' @param yr.res Deprecated - use age.res instead.
#' @param date.res Date distributions are plotted using \code{date.res=100} points by default.
#' @param rotate.axes By default, the age-depth model is plotted with the depths on the horizontal axis and ages on the vertical axis. This can be changed with \code{rotate.axes=TRUE}.
#' @param rev.age The direction of the age axis, which can be reversed using \code{rev.age=TRUE}.
#' @param rev.yr Deprecated - use rev.age instead.
#' @param rev.d The direction of the depth axis, which can be reversed using \code{rev.d=TRUE}.
#' @param depth.unit Units of the depths. Defaults to the one provided in the Bacon() command, \code{depth.unit=set$depth.unit}.
#' @param age.unit Units of the ages. Defaults to \code{age.unit="yr"}.
#' @param unit Deprecated and replaced by \code{depth.unit}.
#' @param maxcalc Number of depths to calculate ages for. If this is more than \code{maxcalc=500}, a warning will be shown that calculations will take time.
#' @param height The maximum heights of the distributions of the dates on the plot. See also \code{normalise.dists}.
#' @param calheight Multiplier for the heights of the distributions of dates on the calendar scale. Defaults to \code{calheight=1}.
#' @param ex As an alternative to 'height' and 'calheight', the distribution heights can be set as either a single value (e.g., \code{ex=1}) or for each date (for example, \code{myex=rep(1, nrow(info$dets)); myex[1:5] <- 10; agedepth(ex=myex)}).
#' @param mirror Plot the dates as 'blobs'. Set to \code{mirror=FALSE} to plot simple distributions.
#' @param up Directions of distributions if they are plotted non-mirrored. Default \code{up=TRUE}.
#' @param cutoff Avoid plotting very low probabilities of date distributions (default \code{cutoff=0.1}).
#' @param plot.range Whether or not to plot the curves showing the confidence ranges of the age-model. Defaults to (\code{plot.range=TRUE}).
#' @param range.col The colour of the curves showing the confidence ranges of the age-model. Defaults to medium grey (\code{range.col=grey(0.5)}).
#' @param range.lty The line type of the curves showing the confidence ranges of the age-model. Defaults to \code{range.lty=12}.
#' @param range.lwd Widths of the lines of the ranges of the age-depth model. Default \code{range.lwd=1}.
#' @param mn.col The colour of the mean age-depth model: default \code{mn.col="red"}.
#' @param mn.lty The line type of the mean age-depth model. Default \code{mn.lty=12}.
#' @param mn.lwd Width of the line of the mean age-depth model. Default \code{mn.lwd=1}.
#' @param med.col The colour of the median age-depth model: not drawn by default \code{med.col=NA}.
#' @param med.lty The line type of the median age-depth model. Default \code{med.lty=12}.
#' @param med.lwd Width of the line of the median age-depth model. Default \code{med.lwd=1}.
#' @param C14.col The colour of the calibrated ranges of the dates. Default is semi-transparent blue: \code{C14.col=rgb(0,0,1,.35)}.
#' @param C14.border The colours of the borders of calibrated 14C dates. Default is semi-transparent dark blue: \code{C14.border=rgb(0, 0, 1, 0.5)}.
#' @param cal.col The colour of the non-14C dates. Default is semi-transparent blue-green: \code{cal.col=rgb(0,.5,.5,.35)}.
#' @param cal.border The colour of the border of non-14C dates in the age-depth plot: default semi-transparent dark blue-green: \code{cal.border=rgb(0,.5,.5,.5)}. Not used by default.
#' @param dates.col As an alternative to colouring dates based on whether they are 14C or not, sets of dates can be coloured as, e.g., \code{dates.col=colours()[2:100]}.
#' @param pb.background Probability at which total Pb values are considered to have reached background values, or in other words, that their modelled values are at or below supported + detection limit (Al)). Setting this at 0.5 means that any depth with a Pb measurement, where at least half of the iterations model Pb values reaching background values, is flagged as having reached background. The age-model is not extended to any Pb measurements that have reached background.
#' @param pbmodelled.col Colour of the modelled 210Pb values. Defaults to shades of blue: \code{pbmodelled.col=function(x) rgb(0,0,1,x)}.
#' @param pbmeasured.col Colour of the measured 210Pb values (default \code{pbmeasured.col="blue"}). Draws rectangles of the upper and lower depths as well as the Pb values with 95 percent error ranges. 
#' @param pb.lim Axis limits for the Pb-210 data. Calculated automatically by default (\code{pblim=c()}).
#' @param supp.col Colour of supported Pb-210. Defaults to semi-transparent purple, because why not.
#' @param remove.tail Whether or not to remove the tail measurements when plotting. Sometimes automated removal might go wrong, so then this option can be used to avoid removing the tail measurements.
#' @param MCMC.resample After the MCMC run, if there are more MCMC iterations than requested, only the last 'ssize' iterations will be retained. Defaults to TRUE.
#' @param hiatus.col The colour of the depths of any hiatuses. Default \code{hiatus.col=grey(0.5)}.
#' @param hiatus.lty The line type of the depths of any hiatuses. Default \code{hiatus.lty=12}.
#' @param rgb.scale The function to produce a coloured representation of all age-models. Needs 3 values for the intensity of red, green and blue. Defaults to grey-scales: \code{rgb.scale=c(0,0,0)}, but could also be, say, scales of red (\code{rgb.scale=c(1,0,0)}). 
#' @param rgb.res Resolution of the colour spectrum depicting the age-depth model. Default \code{rgb.res=100}.
#' @param slump.col Colour of slumps. Defaults to \code{slump.col=grey(0.8)}.
#' @param normalise.dists By default, the distributions of more precise dates will cover less time and will thus peak higher than less precise dates. This can be avoided by specifying \code{normalise.dists=FALSE}.
#' @param same.heights Plot the distributions of the dates all at the same maximum height (default \code{same.height=FALSE}).
#' @param cc Calibration curve for 14C dates: \code{cc=1} for IntCal20 (northern hemisphere terrestrial), \code{cc=2} for Marine20 (marine), \code{cc=3} for SHCal20 (southern hemisphere terrestrial). For dates that are already on the cal BP scale use \code{cc=0}.
#' @param title The title of the age-depth model is plotted on the main panel. By default this is the core's name. To leave empty: \code{title=""}.
#' @param title.location Location of the title. Default \code{title.location='topleft'}.
#' @param title.size Size of the title font. Defaults to \code{title.size=1.5}.
#' @param plot.labels Whether or not to plot labels next to the dated depths. Defaults to \code{FALSE}.
#' @param labels Add labels to the dates (as given by the first column of the .csv file). \code{FALSE} by default.
#' @param label.age Position on the age axis of the date labels. By default draws them before the youngest age (1), but can also draw them after the oldest age (2), or above its mean (3).
#' @param label.offset Offsets of the positions of the labels, giving the x and y offsets. Defaults to c(0,0).
#' @param label.size Size of labels.
#' @param label.col Colour of the labels. Defaults to the colour given to the borders of the dates.
#' @param label.adj  Justification of the labels. Follows R's adj option: A value of '0' produces left-justified text, '0.5' (the default) centered text and '1' right-justified text.
#' @param label.rot Rotation of the label. 0 by default (horizontal).
#' @param after Sets a short section above and below hiatus.depths within which to calculate ages. For internal calculations - do not change.
#' @param bty Type of box to be drawn around plots (\code{"n"} for none, and \code{"l"} (default), \code{"7"}, \code{"c"}, \code{"u"}, or \code{"o"} for correspondingly shaped boxes).
#' @param mar.left Plot margins for the topleft panel (amount of white space along edges of axes 1-4). Default \code{mar.left=c(3,3,1,1)}.
#' @param mar.middle Plot margins for the middle panel(s) at the top (amount of white space along edges of axes 1-4). Default \code{mar.middle=c(3,3,1,1)}.
#' @param mar.right Plot margins for the topright panel (amount of white space along edges of axes 1-4). Default \code{mar.right=c(3,3,1,1)}.
#' @param mar.main Plot margins for the main panel (amount of white space along edges of axes 1-4). Default \code{mar.main=c(3,3,1,1)}.
#' @param righthand Adapt the righthand margins by a certain amount (default 2) to allow a righthand axis to be plotted (for plum)
#' @param mgp Axis text margins (where should titles, labels and tick marks be plotted). Defaults to \code{mgp=c(1.7, 0.7, 0.0)}.
#' @param xaxs Extension of x-axis. By default, add some extra white-space at both extremes (\code{xaxs="r"}). See ?par for other options.
#' @param yaxs Extension of y-axis. By default, add no extra white-space at both extremes (\code{yaxs="i"}). See ?par for other options.
#' @param MCMC.col Colour of the MCMC output. Defaults to \code{post.col=grey(0.4))}.
#' @param post.col Colour of the posterior histogram. Defaults to \code{post.col=grey(0.8))}.
#' @param post.border Colour of the posterior border. Defaults to \code{post.border=grey(0.4))}.
#' @param prior.col Colour of the prior curve. Defaults to light green, \code{prior.col=3)}.
#' @param prior.lwd Line width of the prior curve. Defaults to \code{prior.lwd=2)}.
#' @param prior.fontcol Colour of the font accompanying the posterior histograms. Defaults to red, \code{prior.fontcol=2)}.
#' @param prior.ticks Plot tickmarks and values on the vertical axes for the prior and posterior distributions. Defaults to no tick marks (\code{prior.ticks="n"}). Set to \code{prior.ticks="s"} to plot the tick marks. Note that these values are of little practical use, as they correspond poorly to, e.g., the mean and strength values. All that matters is that the areas of both the prior and the posterior distributions sum to 1; wider distributions tend to give lower peaks, and narrower distributions higher peaks. 
#' @param prior.fontsize Font size of the prior, relative to R's standard size. Defaults to \code{prior.fontsize=0.9}.
#' @param toppanel.fontsize Font size of the top panels, relative to R's standard size. Defaults to \code{prior.fontsize=0.9}.
#' @param mainpanel.tickfontsize Font size of values at the tick marks in the main panel, relative to R's standard size. Defaults to \code{mainpanel.tickfontsize=1}.
#' @param mainpanel.labelfontsize Font size of axis labels in the main panel, relative to R's standard size. Defaults to \code{mainpanel.labelsize=1}.
#' @param acc.xlim Horizontal axis limits of the accumulation rate panel. Calculated automatically by default.
#' @param acc.ylim Vertical axis limits of the accumulation rate panel. Calculated automatically by default.
#' @param mem.xlim Horizontal axis limits of the memory panel. Calculated automatically by default.
#' @param mem.ylim Vertical axis limits of the memory panel. Calculated automatically by default.
#' @param hiatus.xlim Horizontal axis limits of the hiatus size panel. Calculated automatically by default.
#' @param hiatus.ylim Vertical axis limits of the hiatus size panel. Calculated automatically by default.
#' @param phi.xlim Horizontal axis limits of the phi panel. Calculated automatically by default.
#' @param phi.ylim Vertical axis limits of the phi panel. Calculated automatically by default.
#' @param supp.xlim Horizontal axis limits of the supported-Pb panel. Calculated automatically by default.
#' @param supp.ylim Vertical axis limits of the supported-Pb panel. Calculated automatically by default.
#' @param xaxt Whether or not to plot the x-axis. Can be used to adapt axes after a plot. See ?par for other options.
#' @param yaxt Whether or not to plot the y-axis. Can be used to adapt axes after a plot. See ?par for other options.
#' @param plot.pb Plot the 210Pb data (if present). Defaults to \code{plot.pb=TRUE}.
#' @param pb.lty Line type of measured Pb-210 data.
#' @param plot.pdf Produce a pdf file of the age-depth plot.
#' @param model.only By default, panels showing the MCMC iterations and the priors and posteriors for accumulation rate and memory are plotted above the main age-depth model panel. This can be avoided by supplying \code{model.only=TRUE}. Note however that this removes relevant information to evaluate the age-depth model, so we do recommend to present age-models together with these upper panels.
#' @param dates.only By default, the age-depth model is plotted on top of the dates. This can be avoided by supplying \code{dates.only=TRUE}.
#' @param verbose Provide a summary of the age ranges after producing the age-depth model graph; default \code{verbose=FALSE}.
#' @param roundby Rounding of the values reported at the end of the function. Defaults to 2 decimals.
#' @param save.info By default, a variable called `info' with relevant information about the run (e.g., core name, priors, settings, ages, output) is saved into the working directory. Note that this will overwrite any existing variable with the same name - as an alternative, one could run, e.g., \code{myvar <- Bacon()}, followed by supplying the variable \code{myvar} in any subsequent commands.
#' @author Maarten Blaauw, J. Andres Christen
#' @return A plot of the age-depth model, and estimated ages incl. confidence ranges for each depth.
#' @examples
#' \donttest{
#'   rbacon::Bacon(ask=FALSE, coredir=tempfile())
#'   agedepth()
#' }
#' @export
agedepth <- function(set=get('info'), BCAD=set$BCAD, depth.unit=set$depth.unit, age.unit="yr", unit=depth.unit, d.lab=c(), age.lab=c(), yr.lab=age.lab, kcal=FALSE, acc.lab=c(), mem.lab=c(), d.min=c(), d.max=c(), d.by=c(), depths=set$depths, depths.file=FALSE, accordion=c(), plotatthesedepths=c(), age.min=c(), yr.min=age.min, age.max=c(), yr.max=age.max, hiatus.option=1, dark=c(), prob=set$prob, rounded=c(), d.res=400, age.res=400, yr.res=age.res, date.res=100, rotate.axes=FALSE, rev.age=FALSE, rev.yr=rev.age, rev.d=FALSE, maxcalc=500, height=1, calheight=1, ex=1, mirror=TRUE, up=TRUE, cutoff=.1, plot.range=TRUE, range.col=grey(.5), range.lty="12", range.lwd=1, mn.col="red", mn.lty="12", mn.lwd=1, med.col=NA, med.lty="12", med.lwd=1, C14.col=rgb(0,0,1,.35), C14.border=rgb(0,0,1,.5), cal.col=rgb(0,.5,.5,.35), cal.border=rgb(0,.5,.5,.5), dates.col=c(), pb.background=.5, pbmodelled.col=function(x) rgb(0,0,1,.5*x), pbmeasured.col="blue", pb.lim=c(), supp.col=rgb(.5,0,.5,.5), remove.tail=TRUE, MCMC.resample=TRUE, hiatus.col=grey(0.5), hiatus.lty="12", rgb.scale=c(0,0,0), rgb.res=100, slump.col=grey(0.8), normalise.dists=TRUE, same.heights=FALSE, cc=set$cc, title=set$core, title.location="topleft", title.size=1.5, plot.labels=FALSE, labels=c(), label.age=1, label.size=0.8, label.col="black", label.offset=c(0,0), label.adj=c(0.5,0), label.rot=0, after=set$after, bty="l", mar.left=c(3,3,1,.5), mar.middle=c(3,0,1,.5), mar.right=c(3,0,1,.5), mar.main=c(3,3,1,1), righthand=3, mgp=c(1.7,.7,.0), xaxs="r", yaxs="i", MCMC.col=grey(.4), post.col=grey(.8), post.border=grey(.4), prior.col=3, prior.lwd=2, prior.fontcol=2, prior.ticks="n", prior.fontsize=0.9, toppanel.fontsize=0.9, mainpanel.tickfontsize=1, mainpanel.labelfontsize=1, acc.xlim=c(), acc.ylim=c(), mem.xlim=c(), mem.ylim=c(), hiatus.xlim=c(), hiatus.ylim=c(), phi.xlim=c(), phi.ylim=c(), supp.xlim=c(), supp.ylim=c(), xaxt="s", yaxt="s", plot.pb=TRUE, pb.lty=1, plot.pdf=FALSE, dates.only=FALSE, model.only=FALSE, verbose=TRUE, roundby=2, save.info=set$save.info) {
# Load the output, if it exists
  outp <- paste0(set$prefix, ".out")

  if(file.exists(outp))
    set <- rbacon:::Bacon.AnaOut(outp, set, MCMC.resample)

  # Plum-specific
  if(set$isplum) {
    outPlum <- paste0(set$prefix, "_plum.out")
    if(file.exists(outPlum))
      set <- rbacon:::Plum.AnaOut(outPlum, set, MCMC.resample)
    if(save.info)
      assign_to_global("info", set)

    # sometimes runs don't go well, with the age-model totally lost. This is indicated by a very peaked posterior for memory, very close to 1
    if(mn <- mean(set$output[,set$K+2]) > 0.95)
      if(mn > set$mem.mean)
        message("\nWarning, this run has a very high posterior memory (", round(mn, 2), ") and probably didn't go very well. Please run again\n")
    bg <- background(set) # calculate which Pb data have likely reached background
    set$background <- bg
  }

  set$BCAD <- BCAD # tmp May 2021
  # Adapt ages of sections which contain hiatuses
  if(!is.na(set$hiatus.depths[1]))
    set <- hiatus.slopes(set, hiatus.option)
  if(save.info)
    assign_to_global("info", set)

  oldpar <- par(no.readonly = TRUE)
  #on.exit(par(oldpar))
  newpar <- par(mar=mar.main, mgp=mgp)

  if(!model.only) {
    newpar <- par(mar=mar.left, bty=bty, mgp=mgp, xaxs=xaxs, yaxs=yaxs)
    ifelse(set$isplum, pn <- c(1:5, rep(6,5)), pn <- c(1:3, rep(4,3)))
    if(!is.na(set$hiatus.depths[1]))
      if(is.na(set$boundary[1]))
        ifelse(set$isplum, pn <- c(1:6, rep(7,6)), pn <- c(1:4, rep(5,4)))
    layout(matrix(pn, nrow=2, byrow=TRUE), heights=c(.3,.7))
    PlotLogPost(set, 0, set$Tr, xaxs=xaxs, yaxs=yaxs, panel.size=toppanel.fontsize, col=MCMC.col) # convergence information
    newpar <- par(mar=mar.middle) # reduce white space
    #on.exit(par(oldpar))
    set$post.acc <- PlotAccPost(set, depth.unit=depth.unit, age.unit=age.unit, xaxs=xaxs, yaxs=yaxs, yaxt=prior.ticks, prior.size=prior.fontsize, panel.size=toppanel.fontsize, acc.xlim=acc.xlim, acc.ylim=acc.ylim, acc.lab=acc.lab, line.col=prior.col, line.width=prior.lwd, text.col=prior.fontcol, hist.col=post.col, hist.border=post.border)
    set$post.mem <- PlotMemPost(set, set$core, set$K, "", set$mem.strength, set$mem.mean, ds=1, thick=set$thick, xaxs=xaxs, yaxs=yaxs, yaxt=prior.ticks, prior.size=prior.fontsize, panel.size=toppanel.fontsize, mem.xlim=mem.xlim, mem.ylim=mem.ylim, mem.lab=mem.lab, line.col=prior.col, line.width=prior.lwd, text.col=prior.fontcol, hist.col=post.col, hist.border=post.border)
    if(!is.na(set$hiatus.depths[1]))
      if(is.na(set$boundary[1]))
         set$post.hiatus <- PlotHiatusPost(set, set$hiatus.max, xaxs=xaxs, yaxs=yaxs, yaxt=prior.ticks, prior.size=prior.fontsize, panel.size=toppanel.fontsize, hiatus.xlim=mem.xlim, hiatus.ylim=mem.ylim, line.col=prior.col, line.width=prior.lwd, text.col=prior.fontcol, hist.col=post.col, hist.border=post.border)
    if(set$isplum) {
       set$post.phi <- PlotPhiPost(set, xaxs=xaxs, yaxs=yaxs, yaxt=prior.ticks, prior.size=prior.fontsize, panel.size=toppanel.fontsize, phi.xlim=phi.xlim, phi.ylim=phi.ylim, line.col=prior.col, line.width=prior.lwd, text.col=prior.fontcol, hist.col=post.col, hist.border=post.border)
       if(set$nPs > 1)
         prior.ticks <- "s" # because with varying supported Pb, the y-axis is important
       par(mar=mar.right)
       set$post.supp <- PlotSuppPost_repaired(set, xaxs=xaxs, yaxs=yaxs, yaxt=prior.ticks, prior.size=prior.fontsize, panel.size=toppanel.fontsize, supp.xlim=supp.xlim, supp.ylim=supp.ylim, line.col=prior.col, line.width=prior.lwd, text.col=prior.fontcol, hist.col=post.col, hist.border=post.border, data.col=supp.col)
       mar.main[4] <- mar.main[4] + righthand # to enable space for righthand axis
    }
    par(mar=mar.main) # new May 2021
  }
  
  # calculate and plot the ranges and 'best' estimates for each required depth
  if(length(d.min) == 0)
    d.min <- set$d.min
  if(length(d.max) == 0)
    if(length(accordion) == 2)
      d.max <- stretch(set$d.max, accordion[1], accordion[2]) else
        d.max <- set$d.max
  if(length(d.by) == 0)
    d.by <- set$d.by
  if(depths.file) {
    dfile <- paste0(set$coredir, set$core, "/", set$core, "_depths.txt")
    if(!file.exists(dfile))
      stop("I cannot find the file ", paste0(set$coredir, set$core, "/", set$core, "_depths.txt"), call.=FALSE)
    depths <- rbacon:::fastread(dfile, header=FALSE)[,1]
    if(!is.numeric(depths[1]))
      stop("File should contain numbers only, no headers", call.=FALSE)
  }
  if(length(depths) > 0)
    d <- sort(depths) else
      d <- seq(set$d.min, set$d.max, by=d.by) # not d.min itself as depths < set$d.min cannot be calculated. Same for d.max, best not extrapolate here
  # but if squeezing/stretching has to be done, then :
  if(length(accordion) == 2) {
    d <- seq(set$d.min, stretch(set$d.max, accordion[1], accordion[2]), by=d.by)
    d <- squeeze(d, accordion[1], accordion[2])
  }

  if(length(d) > maxcalc)
    message("Warning, this will take quite some time to calculate. I suggest increasing d.by to, e.g., ", 10*d.by, "\n") # was set$d.by

  for(i in set$hiatus.depths)
    d <- sort(unique(c(i+after, i, d)))
  for(i in set$slump)
    d <- sort(unique(c(i+after, i, d)))
  
  if(set$isplum) # new May 2021
    if(set$ra.case == 0) # renamed from incorrectly named radon.case
      if(!set$hasBaconData) # May 2021
        d <- d[which(d <= max(set$detsOrig[,2]))]

  if(verbose)
    message("Calculating age ranges...")
  modelranges <- c()
  ranges <- Bacon.rng(d, set, BCAD=BCAD, prob=prob)
  d.rng <- d
  # calculate calendar axis limits
  modelranges <- range(ranges[!is.na(ranges)])
  if(length(set$calib$probs) > 0) {
    dates <- set$calib$probs
  dateranges <- c()
  for(i in 1:length(dates))
    if(BCAD)
      dateranges <- range(dateranges, 1950-dates[[i]][,1], na.rm=TRUE) else
        dateranges <- range(dateranges, dates[[i]][,1], na.rm=TRUE)
  } else dateranges <- modelranges # assuming plum with no additional dates

 if(length(age.min) == 0)
    age.min <- min(modelranges, dateranges)
  if(length(age.max) == 0)
    age.max <- max(modelranges, dateranges)
  age.lim <- extendrange(c(age.min, age.max), f=0.01)

  if(BCAD)
    age.lim <- rev(age.lim)
  if(rev.age)
    age.lim <- rev(age.lim)
  d.lim <- rev(extendrange(c(d.max, d.min), f=0.01))
  if(rev.d)
    d.lim <- d.lim[2:1]

  if(length(d.lab) == 0)
    d.lab <- paste0("Depth (", depth.unit, ")")
  if(length(age.lab) == 0)
    age.lab <- ifelse(BCAD, "BC/AD", ifelse(kcal, "kcal BP", paste("cal", age.unit, "BP")))

  if(kcal)
    ifelse(rotate.axes, xaxt <- "n", yaxt <- "n")

  if(set$isplum)
    if(remove.tail) {
      # has to check for non-210Pb dates below the tail
      above <- which(set$background < pb.background)
      if(set$hasBaconData) 
        d.lim <- range(set$dets[above,4], set$detsBacon[,4])[2:1] else
          d.lim[which(d.lim==max(d.lim))] <- max(set$dets[above,4])
    }

  if(rotate.axes)
    plot(0, type="n", ylim=d.lim, xlim=age.lim, ylab=d.lab, xlab=age.lab, bty="n", xaxt=xaxt, yaxt=yaxt, mar=mar.main, cex.axis=mainpanel.tickfontsize, cex.lab=mainpanel.labelfontsize) else
      plot(0, type="n", xlim=d.lim[2:1], ylim=age.lim, xlab=d.lab, ylab=age.lab, bty="n", xaxt=xaxt, yaxt=yaxt, mar=mar.main, cex.axis=mainpanel.tickfontsize, cex.lab=mainpanel.labelfontsize)
  if(kcal)
    axis(ifelse(rotate.axes, 1, 2), pretty(age.lim), pretty(age.lim/1e3), cex.axis=mainpanel.labelfontsize)

  if(set$isplum) { # new May 2021
    top <- min(set$dets[,4]-set$dets[,5])
    if(d.min > top)
      d.min <- top # because the top 210Pb-dated depth should be at or below d.min
    if(!set$hasBaconData) { # needs more work
      above <- which(set$background < pb.background)
      dmax <- max(set$dets[above,4]) # cut depths that have reached background
      if(remove.tail)
        if(!is.infinite(dmax)) {
          ranges <- ranges[which(d.rng <= dmax),]
          d <- d[which(d <= dmax)]
          d.max <- dmax
        }
    }
  }

  if(!dates.only) {
    if(verbose)
      message("\nPreparing ghost graph... ")
    agedepth.ghost(set, rotate.axes=rotate.axes, accordion=accordion, BCAD=BCAD, d.res=d.res, age.res=age.res, rev.d=rev.d, rev.age=rev.age, rgb.res=rgb.res, dark=dark, rgb.scale=rgb.scale, age.lim=age.lim)
  }

  if(length(set$slump) > 0 )
    for(i in 1:nrow(set$slump))
      if(rotate.axes)
        rect(min(age.lim)-1e3, set$slump[i,1], max(age.lim)+1e3, set$slump[i,2], col=slump.col, border=slump.col) else
          rect(set$slump[i,1], min(age.lim)-1e3, set$slump[i,2], max(age.lim)+1e3, col=slump.col, border=slump.col)

  if(length(plotatthesedepths) == 0) 
    thesedateddepths <- c() else {
      dseq <- seq(d.min, d.max, length=length(plotatthesedepths)) # this OK? May 2023
      if(set$isplum)
        thesedateddepths <- approx(dseq, plotatthesedepths, set$detsBacon[,4])$y else
          thesedateddepths <- approx(dseq, plotatthesedepths, set$dets[,4])$y
    }

  if(set$isplum) {
    if(set$hasBaconData)
      calib.plot(set, dets=set$detsBacon, accordion=accordion, BCAD=BCAD, cc=cc, rotate.axes=rotate.axes, height=height, calheight=calheight, ex=ex, mirror=mirror, up=up, date.res=date.res, cutoff=cutoff, C14.col=C14.col, C14.border=C14.border, cal.col=cal.col, cal.border=cal.border, dates.col=dates.col, new.plot=FALSE, same.heights=same.heights)
	
    if(plot.pb)		
      set <- draw.pbmodelled(set, BCAD=BCAD, rotate.axes=rotate.axes, age.lim=age.lim, d.lim=c(d.min, d.max), pbmodelled.col=pbmodelled.col, pbmeasured.col=pbmeasured.col, pb.lim=pb.lim, supp.col=supp.col, mgp=mgp, pb.lty=pb.lty) # pointing to set since April 2025 
  } else
    calib.plot(set, dets=set$dets, accordion=accordion, BCAD=BCAD, cc=cc, rotate.axes=rotate.axes, height=height, calheight=calheight, ex=ex, mirror=mirror, up=up, date.res=date.res, cutoff=cutoff, C14.col=C14.col, C14.border=C14.border, cal.col=cal.col, cal.border=cal.border, dates.col=dates.col, new.plot=FALSE, same.heights=same.heights)

  if(plot.labels) {
    if(length(labels) == 0)
      labels <- set$dets[,1]
    for(i in 1:nrow(set$dets)) {
      age <- set$calib$probs[[i]][,1]
      if(label.age == 1)
        age.pos <- max(age) else
          ages.pos <- min(age)
      d <- set$calib$d[i]
      if(rotate.axes)
         text(age.pos+label.offset[1], d+label.offset[2], labels[i], cex=label.size, col=label.col, adj=label.adj, srt=label.rot) else
           text(d+label.offset[1], age.pos+label.offset[2], labels[i], cex=label.size, col=label.col, adj=label.adj, srt=label.rot)
    }
  }
  legend(title.location, title, bty="n", cex=title.size)
  box(bty=bty)

  hiatus.depths <- set$hiatus.depths
  if(!is.na(set$boundary[1]))
    hiatus.depths <- set$boundary
  if(length(hiatus.depths) > 0) # draw locations of boundaries/hiatuses
    if(rotate.axes)
      abline(h=hiatus.depths, col=hiatus.col, lty=hiatus.lty) else
        abline(v=hiatus.depths, col=hiatus.col, lty=hiatus.lty)

  th <- rbind(1, nrow(ranges))
  if(!is.na(set$hiatus.depths[1])) {
    hi.d <- c()
    for(i in set$hiatus.depths)
      hi.d <- c(hi.d, max(which(d <= i)))
    th <- array(sort(c(1, nrow(ranges), hi.d-1, hi.d)), dim=c(2,length(hi.d)+1))
  }

  if(length(accordion) == 2)
    d <- stretch(d, accordion[1], accordion[2])

  if(!dates.only)
    for(i in 1:ncol(th)) {
      h <- th[1,i] : th[2,i]
      if(rotate.axes) {
        lines(ranges[h,1], d[h], col=range.col, lty=range.lty, lwd=range.lwd)
        lines(ranges[h,2], d[h], col=range.col, lty=range.lty, lwd=range.lwd)
        lines(ranges[h,3], d[h], col=med.col, lty=med.lty, lwd=med.lwd) # median
        lines(ranges[h,4], d[h], col=mn.col, lty=mn.lty, lwd=mn.lwd) # mean
      } else {
          lines(d[h], ranges[h,1], col=range.col, lty=range.lty, lwd=range.lwd)
          lines(d[h], ranges[h,2], col=range.col, lty=range.lty, lwd=range.lwd)
          lines(d[h], ranges[h,3], col=med.col, lty=med.lty, lwd=med.lwd) # median
          lines(d[h], ranges[h,4], col=mn.col, lty=mn.lty, lwd=mn.lwd) # mean
        }
    }

  if(length(rounded) == 0)
    rounded <- ifelse(set$isplum, 1, 0)
  set$ranges <- cbind(d, round(ranges, rounded))
  colnames(set$ranges) <- c("depth", "min", "max", "median", "mean")
  if(save.info)
    assign_to_global("info", set)

  if(plot.pdf)
    if(names(dev.cur()) != "null device")
      rbacon:::export.pdf(paste0(set$prefix, ".pdf"))

  rbacon:::fastwrite(set$ranges, paste0(set$prefix, "_ages.txt"), quote=FALSE, row.names=FALSE, sep="\t") # was write.table
  rng <- abs(round(set$ranges[,3]-set$ranges[,2], rounded))
  min.rng <- d[which(rng==min(rng, na.rm=TRUE))]
  max.rng <- d[which(rng==max(rng, na.rm=TRUE))]

  if(length(min.rng)==1)
    min.rng <- paste(age.unit, "at", min.rng, noquote(depth.unit)) else
      min.rng <- paste(age.unit, "between", min(min.rng), "and", max(min.rng), noquote(depth.unit))
  if(length(max.rng)==1)
    max.rng <- paste(age.unit, "at", max.rng, noquote(depth.unit)) else
      max.rng <- paste(age.unit, "between", min(max.rng), "and", max(max.rng), noquote(depth.unit))

  if(verbose) {
    if(!dates.only)
      message("\nMean ", 100*prob, "% confidence ranges ", round(mean(rng), rounded), " ", age.unit, ", min. ",
        min(rng), " ", min.rng, ", max. ", max(rng), " ", max.rng)

    if(set$isplum) {
      below <- which(bg > pb.background)
      if(length(below) == 0)
        message("\nIt seems that the Pb-measurements haven't reached background") else {
        message("Pb-210 likely at or below detection limit from ", min(set$dets[below,4]), " ", set$depth.unit, " depth onward: ")
        for(i in seq_along(below))
          cat(set$dets[below[i],4], " ", set$depth.unit,
            " (", round(100*bg[below[i]]), "%)",
            if(i < length(below)) ", " else "\n", sep = "")
      }
    } else
        rbacon:::overlap(set)
	
	# report summaries of posteriors	
    posteriors <- c(set$post.acc, set$post.mem)	
    if(!is.na(hiatus.depths[1]))
      posteriors <- c(posteriors, set$post.hiatus)
	posteriors <- round(posteriors, roundby)
	if(!is.na(hiatus.depths[1]))
	  message("Posteriors: accrate mean ", posteriors[1], ", shape ", posteriors[2], 
	    ", memory mean ", posteriors[3], ", strength ", posteriors[4], 
		", hiatus mean ", posteriors[5], ", shape ", posteriors[6]) else
  	  message("Posteriors: accrate mean ", posteriors[1], ", shape ", posteriors[2], 
  	    ", memory mean ", posteriors[3], ", strength ", posteriors[4])
		
    if(set$isplum) {
      plumpost <- round(c(set$post.phi, set$post.supp), roundby)
      message("phi mean ", plumpost[1], ", shape ", plumpost[2], 
	    ", supported mean ", plumpost[3], ", shape ", plumpost[4])
	}	
  }
  
  #par(oldpar) # tmp Jan 2021
  invisible(set)
}
