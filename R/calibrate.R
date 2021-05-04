# these functions are commented because they are already defined in rbacon

# ### for running Plum, but is looked for by generic agedepth() function, so is included in the rbacon code
# #' @name calib.plumbacon.plot
# #' @title Plot the dates
# #' @description Produce a plot of the dated depths and their dates
# #' @details This function is generally called internally to produce the age-depth graph.
# #' It can be used to produce custom-built graphs.
# #' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
# #' @param BCAD The calendar scale of graphs is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
# #' @param cc Calibration curve to be used (defaults to info$cc)
# #' @param rotate.axes The default of plotting age on the horizontal axis and event probability on the vertical one can be changed with \code{rotate.axes=TRUE}.
# #' @param rev.d The direction of the depth axis can be reversed from the default (\code{rev.d=TRUE}).
# #' @param rev.age The direction of the calendar age axis can be reversed from the default (\code{rev.age=TRUE})
# #' @param rev.yr Deprecated - use rev.age instead
# #' @param age.lim Minimum and maximum calendar age ranges, calculated automatically by default (\code{age.lim=c()}).
# #' @param yr.lim Deprecated - use age.lim instead
# #' @param d.lab The labels for the depth axis. Default \code{d.lab="Depth (cm)"}.
# #' @param age.lab The labels for the calendar axis (default \code{yr.lab="cal BP"} or \code{"BC/AD"} if \code{BCAD=TRUE}).
# #' @param yr.lab Deprecated - use age.lab instead
# #' @param height The heights of the distributions of the dates. See also \code{normalise.dists}.
# #' @param calheight Multiplier for the heights of the distributions of dates on the calendar scale. Defaults to \code{calheight=1}.
# #' @param mirror Plot the dates as 'blobs'. Set to \code{mirror=FALSE} to plot simple distributions.
# #' @param up Directions of distributions if they are plotted non-mirrored. Default \code{up=TRUE}.
# #' @param cutoff Avoid plotting very low probabilities of date distributions (default \code{cutoff=0.001}).
# #' @param date.res Date distributions are plotted using \code{date.res=100} points by default.
# #' @param C14.col Colour of the calibrated distributions of the dates. Default is semi-transparent blue: \code{rgb(0,0,1,.35)}.
# #' @param C14.border Colours of the borders of calibrated 14C dates. Default is transparent dark blue: cal.col
# #' @param cal.col Colour of the non-14C dates in the age-depth plot: default semi-transparent blue-green: \code{rgb(0,.5,.5,.35)}.
# #' @param cal.border Colour of the of the border of non-14C dates in the age-depth plot: default semi-transparent dark blue-green: \code{rgb(0,.5,.5,.5)}.
# #' @param dates.col As an alternative to colouring dates based on whether they are 14C or not, sets of dates can be coloured as, e.g., \code{dates.col=colours()[2:100]}.
# #' @param slump.col Colour of slumps. Defaults to \code{slump.col=grey(0.8)}.
# #' @param new.plot Start a new plot (\code{new.plot=TRUE}) or plot over an existing plot (\code{new.plot=FALSE}).
# #' @param plot.dists Plot the distributions of the dates (default \code{plot.dists=TRUE}).
# #' @param same.heights Plot the distributions of the dates all at the same maximum height (default \code{same.height=FALSE}).
# #' @param normalise.dists By default, the distributions of more precise dates will cover less time and will thus peak higher than less precise dates. This can be avoided by specifying \code{normalise.dists=FALSE}.
# #' @author Maarten Blaauw, J. Andres Christen
# #' @return NA
# #' @export
# ### produce plots of the calibrated distributions
# calib.plumbacon.plot <- function(set=get('info'), BCAD=set$BCAD, cc=set$cc, rotate.axes=FALSE, rev.d=FALSE, rev.age=FALSE, rev.yr=rev.age, age.lim=c(), yr.lim=age.lim, date.res=100, d.lab=c(), age.lab=c(), yr.lab=age.lab, height=1, calheight=1, mirror=TRUE, up=TRUE, cutoff=.001, C14.col=rgb(0,0,1,.5), C14.border=rgb(0,0,1,.75), cal.col=rgb(0,.5,.5,.5), cal.border=rgb(0,.5,.5,.75), dates.col=c(), slump.col=grey(0.8), new.plot=TRUE, plot.dists=TRUE, same.heights=FALSE, normalise.dists=TRUE) {
#   #height <- length(set$d.min:set$d.max) * height/50
#   if(length(age.lim) == 0)
#     lims <- c()
#   for(i in 1:length(set$calib$probs))
#     lims <- c(lims, set$calib$probs[[i]][,1])
#   age.min <- min(lims)
#   age.max <- max(lims)
#   if(BCAD) {
#     age.min <- 1950 - age.min
#     age.max <- 1950 - age.max
#     }
#   if(length(age.lab) == 0)
#     age.lab <- ifelse(set$BCAD, "BC/AD", paste("cal", set$age.unit, " BP"))
#   age.lim <- extendrange(c(age.min, age.max), f=0.01)
#   if(rev.age)
#     age.lim <- age.lim[2:1]
#   dlim <- extendrange(set$elbows, f=0.05)
#   if(rev.d)
#     dlim <- dlim[2:1]
#   if(length(d.lab) == 0)
#     d.lab <- paste0("depth (", set$depth.unit, ")")
#
#   if(new.plot)
#     if(rotate.axes)
#       plot(0, type="n", xlim=age.lim, ylim=dlim[2:1], xlab=age.lab, ylab=d.lab, main="") else
#         plot(0, type="n", xlim=dlim, ylim=age.lim, xlab=d.lab, ylab=age.lab, main="")
#
#   if(length(set$slump) > 0)
#     if(rotate.axes)
#       abline(h=set$slump, lty=2, col=slump.col) else
#         abline(v=set$slump, lty=2, col=slump.col)
#
#   if(plot.dists)
#     for(i in 1:length(set$calib$probs)) {
#         cal <- cbind(set$calib$probs[[i]])
#         d <- set$calib$d[[i]]
#         cc <- set$calib$cc[[i]]
#         if(BCAD)
#           cal[,1] <- 1950-cal[,1]
#         o <- order(cal[,1])
#         cal <- cbind(cal[o,1], cal[o,2])
#         if(same.heights)
#           cal[,2] <- cal[,2]/max(cal[,2])
#         if(normalise.dists)
#           cal[,2] <- cal[,2]/sum(cal[,2])
#         cal[,2] <- (height * cal[,2])/1 * (max(dlim) - min(dlim)) # 1 is just a scaling factor that looks nice in the examples I tested. Adapt by changing parameter height
#         cal <- cal[cal[,2] >= cutoff*max(cal[,2]),]
#         cal[,2] <- height*cal[,2]
#         if(ncol(set$dets) > 4 && set$dets[i,9] == 0) # cal BP date
#           cal[,2] <- calheight*cal[,2]
#
#         x = cal[,1]
#         y = cal[,2]
#
#         y = y[!duplicated(x)]
#         x = x[!duplicated(x)]
#         #seq(min(cal[,1]), max(cal[,1]), length= length(cal[,1]) )
#         cal <- approx(x, y, seq(min(x), max(x), length= 100 ) ) # tmp
#
#         if(mirror)
#           pol <- cbind(c(d-cal$y, d+rev(cal$y)), c(cal$x, rev(cal$x))) else
#             if(up)
#               pol <- cbind(d-c(0, cal$y, 0), c(min(cal$x), cal$x, max(cal$x))) else
#                 pol <- cbind(d+c(0, cal$y, 0), c(min(cal$x), cal$x, max(cal$x)))
#         if(rotate.axes)
#           pol <- cbind(pol[,2], pol[,1])
#         if(cc > 0) {
#           col <- C14.col
#           border <- C14.border
#         } else {
#            col <- cal.col
#            border <- cal.border
#         }
#         if(length(dates.col) > 0) {
#           col <- dates.col[i]
#           border <- dates.col[i]
#         }
#         polygon(pol, col=col, border=border)
#       }
# }



# ### for running Plum, but is looked for by generic agedepth() function, so is included in the rbacon code
# #' @name draw.pbmodelled
# #' @title Plot the 210Pb data
# #' @description Produce a plot of the 210Pb data and their depths
# #' @details This function is generally called internally to produce the age-depth graph.
# #' It can be used to produce custom-built graphs.
# #' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
# #' @param BCAD The calendar scale of graphs is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
# #' @param rotate.axes The default of plotting age on the horizontal axis and event probability on the vertical one can be changed with \code{rotate.axes=TRUE}.
# #' @param rev.d The direction of the depth axis can be reversed from the default (\code{rev.d=TRUE}).
# #' @param rev.age The direction of the calendar age axis can be reversed from the default (\code{rev.age=TRUE})
# #' @param pb.lim Minimum and maximum of the 210Pb axis ranges, calculated automatically by default (\code{pb.lim=c()}).
# #' @param d.lim Minimum and maximum depths to plot; calculated automatically by default (\code{d.lim=c()}).
# #' @param d.lab The labels for the depth axis. Default \code{d.lab="Depth (cm)"}.
# #' @param pb.lab The label for the 210Pb axis (default \code{pb.lab="210Pb (Bq/kg)"} or \code{"210Pb (dpm/g)"}).
# #' @param supp.col Colour of the supported 210Pb data. Defaults to red: \code{supp.col="red"}.
# #' @param pbmodelled.col Colour of the modelled 210Pb values. Defaults to scales of blue: \code{pbmodelled.col=function(x) rgb(0,0,1,x)}.
# #' @param pbmeasured.col Colour of the measured 210Pb values. Defaults to blue.
# #' @param plot.measured Plot the measured 210Pb values (default \code{plot.measured=TRUE}).
# #' @param age.lim values of the age axis. Used to calculate where to plot the pb values on the secondary axis
# #' @param mgp Axis text margins (where should titles, labels and tick marks be plotted). Defaults to \code{mgp=c(1.7, .7, .0)}.
# #' @author Maarten Blaauw, J. Andres Christen, Marco Aquino-Lopez
# #' @return A plot of the modelled (and optionally the measured) 210Pb values
# #' @export
# draw.pbmodelled <- function(set=get('info'), BCAD=set$BCAD, rotate.axes=FALSE, rev.d=FALSE, rev.age=FALSE, pb.lim=c(), d.lim=c(), d.lab=c(), pb.lab=c(), pbmodelled.col=function(x) rgb(0,0,1,x), pbmeasured.col="blue", supp.col="purple", plot.measured=TRUE, age.lim=c(), mgp=mgp) {
#   depths <- set$detsOrig[,2]
#   dns <- set$detsOrig[,3]
#   Pb <- set$detsOrig[,4]
#   err <- set$detsOrig[,5]
#   thickness <- set$detsOrig[,6]
#   n <- nrow(set$detsOrig)
#
#   if(ncol(set$detsPlum) > 6) {
#     supp <- set$detsOrig[,7]
#     supperr <- set$detsOrig[,8]
#   } else {
#     supp <- set$supportedData[,1]
#     supperr <- set$supportedData[,2]
#     suppd <- set$supportedData[,3]
#     suppthick <- set$supportedData[,4]
#   }
#
#   if(length(d.lab) == 0)
#     d.lab <- paste0("depth (", set$depth.unit, ")")
#   if(length(pb.lab) == 0)
#     pb.lab <- ifelse(set$Bqkg, "210Pb (Bq/kg)", "210Pb (dpm/g)")
#
#   if(length(d.lim) == 0)
#     d.lim <- range(depths)
#   if(rev.d)
#     d.lim <- d.lim[2:1]
#
#   if(length(set$phi) > 0) {
#     Ai <- list(x=NULL, y=NULL)
#     hght <- 0; pbmin <- c(); pbmax <- 0
#     A.rng <- array(0, dim=c(n,2))
#     for(i in 1:length(depths)) {
#       A <- A.modelled(depths[i]-thickness[i], depths[i], dns[i], set)
#       tmp <- density(A)
#       Ai$x[[i]] <- tmp$x
#       Ai$y[[i]] <- tmp$y
#       hght <- max(hght, Ai$y[[i]])
#       pbmin <- min(pbmin, Ai$y[[i]])
#       pbmax <- max(pbmax, Ai$x[[i]])
#       A.rng[i,] <- quantile(A, c((1-set$prob)/2, 1-(1-set$prob)/2))
#     }
#
#   if(length(pb.lim) == 0)
#     pb.lim <- extendrange(c(0, Pb-2*err, Pb+2*err, pbmax), f=c(0,0.05))
#
#     # translate pb values to cal BP values for plotting on the age axis
#     pb2bp <- function(pb, pb.min=pb.lim[1], pb.max=pb.lim[2], agemin=min(age.lim), agemax=max(age.lim)) {
#       ex <- (agemax-agemin) / (pb.max - pb.min)
#       agemin + ex*pb
#     }
#
#     pb2ad <- function(pb, pb.min=pb.lim[1], pb.max=pb.lim[2], agemin=max(age.lim), agemax=min(age.lim)) {
#       ex <- (agemin-agemax) / (pb.max - pb.min)
#       agemin - ex*pb
#     }
#
#     # save the values for later
#     set$Ai <- Ai
#     set$A.rng <- A.rng
#     assign_to_global("info", set, .GlobalEnv)
#
#     this <- ifelse(rotate.axes, 3, 4)
#     pretty.pb <- pretty(c(pbmin, pbmax))
#     if(BCAD) {
#       onad <- pb2ad(pretty.pb)
#       axis(this, rev(onad), rev(pretty.pb), col=pbmeasured.col, col.axis=pbmeasured.col, col.lab=pbmeasured.col)
#     } else
#       {
#         onbp <- pb2bp(pretty.pb)
#         axis(this, onbp, pretty.pb, col=pbmeasured.col, col.axis=pbmeasured.col, col.lab=pbmeasured.col)
#       }
#     pb.lab <- ifelse(set$Bqkg, "Bq/kg", "dpm/g")
#     mtext(pb.lab, this, 1.4, col=pbmeasured.col, cex=.8)
#
#     for(i in 1:length(depths)) {
#       if(BCAD)
#         ages <- pb2ad(rev(Ai$x[[i]])) else
#           ages <- pb2bp(Ai$x[[i]])
#
#     if(BCAD)
#       z <- t(rev(Ai$y[[i]]))/hght else
#         z <- t(Ai$y[[i]])/hght
#
#     if(rotate.axes)
#       image(ages, c(depths[i]-thickness[i], depths[i]), t(z), col=pbmodelled.col(seq(0, 1-max(z), length=50)), add=TRUE) else
#         image(c(depths[i]-thickness[i], depths[i]), ages, z, col=pbmodelled.col(seq(0, 1-max(z), length=50)), add=TRUE)
#     }
#   }
#
#   if(BCAD)
#     pb2bp <- pb2ad
#
#   if(plot.measured)
#     draw.pbmeasured(newplot=FALSE, rotate.axes=rotate.axes, BCAD=BCAD, on.agescale=TRUE, pb.lim=pb.lim, age.lim=age.lim, supp.col=supp.col)
# }



# ### for running Plum, but is looked for by generic agedepth() function (through draw.pbmodelled()), so is included in the rbacon code
# #' @name A.modelled
# #' @title Calculate modelled 210Pb
# #' @description Calculate modelled 210Pb values of a sample slice, based on the parameters of the age-model (i.e., time passed since deposition of the bottom and top of the slice), supported and influx
# #' @param d.top top depth of the slice
# #' @param d.bottom bottom depth of the slice
# #' @param dens Density of the slice (in g/cm3)
# #' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
# #' @param phi The modelled values of the 210Pb influx
# #' @param sup The modelled values of the supported 210Pb
# #' @author Maarten Blaauw
# #' @return a list of modelled values of A
# #' @export
# A.modelled <- function(d.top, d.bottom, dens, set=get('info'), phi=set$phi, sup=set$ps) {
#   if(d.top >= d.bottom)
#     stop("\n d.top should be above d.bottom", call.=FALSE)
#   dd <- 1
#   if(ncol(cbind(sup)) > 1) { # then multiple, varying estimates of supported, find the one belonging to the specified depth interval
#     dd <- set$supportedData[,3] # bottom depths
#     dd <- max(1, which(dd <= d.bottom))
#     sup <- sup[,dd]
#   }
#
#   t.top <- Bacon.Age.d(d.top, BCAD=F) - set$theta0
#   t.bottom <- Bacon.Age.d(d.bottom, BCAD=F) - set$theta0
#   multiply <- ifelse(set$Bqkg, 10, 500)
#   return(sup + ((phi / (.03114*multiply*dens) ) * (exp( -.03114*t.top) - exp(-.03114*t.bottom)) ) )
# }


