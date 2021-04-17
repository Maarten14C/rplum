#' @name draw.pbmeasured
#' @title Plot the 210Pb data
#' @description Produce a plot of the 210Pb data and their depths
#' @details This function is generally called internally to produce the age-depth graph.
#' It can be used to produce custom-built graphs.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param rotate.axes The default of plotting age on the horizontal axis and event probability on the vertical one can be changed with \code{rotate.axes=TRUE}.
#' @param rev.d The direction of the depth axis can be reversed from the default (\code{rev.d=TRUE}).
#' @param rev.age The direction of the calendar age axis can be reversed from the default (\code{rev.age=TRUE})
#' @param BCAD The calendar scale of graphs and age output-files is in cal BP (calendar or calibrated years before the present, where the present is AD 1950) by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param pb.lim Minimum and maximum of the 210Pb axis ranges, calculated automatically by default (\code{pb.lim=c()}).
#' @param age.lim Minimum and maximum of the age ranges to be used to plot 210Pb values. Calculated automatically by default (\code{age.lim=c()}).
#' @param d.lim Minimum and maximum depths to plot; calculated automatically by default (\code{d.lim=c()}).
#' @param d.lab The labels for the depth axis. Default \code{d.lab="Depth (cm)"}.
#' @param pb.lab The label for the 210Pb axis (default \code{pb.lab="210Pb (Bq/kg)"} or \code{"210Pb (dpm/g)"}).
#' @param pbmeasured.col The label for the measured 210Pb data. \code{pbmeasured.col="blue"}.
#' @param pb.log Use a log scale for the 210Pb-axis (default \code{pb.log=FALSE}).
#' @param supp.col Colour of the supported 210Pb data. Defaults to red: \code{supp.col="red"}.
#' @param newplot make new plot (default TRUE)
#' @param on.agescale Plot the Pb-210 on the cal BP scale. Defaults to FALSE.
#' @author Maarten Blaauw, J. Andres Christen, Marco Aquino-Lopez
#' @return A plot of the measured 210Pb values
#' @export
draw.pbmeasured <- function(set=get('info'), rotate.axes=FALSE, rev.d=FALSE, rev.age=FALSE, BCAD=set$BCAD, pb.lim=c(), age.lim=c(), d.lim=c(), d.lab=c(), pb.lab=c(), pbmeasured.col="blue", pb.log=FALSE, supp.col="purple", newplot=TRUE, on.agescale=FALSE) {
  depths <- set$detsOrig[,2]
  dns <- set$detsOrig[,3]
  Pb <- set$detsOrig[,4]
  err <- set$detsOrig[,5]
  thickness <- set$detsOrig[,6]
  n <- nrow(set$detsOrig)

  if(length(pb.lim) == 0) 
    pb.lim <- extendrange(c(0, Pb+2*err), f=c(0,0.05))
 
  # translate pb values to cal BP/AD values for plotting on the age axis
  pb2bp <- function(pb, pb.min=pb.lim[1], pb.max=pb.lim[2], agemin=min(age.lim), agemax=max(age.lim), AD=BCAD) {
    if(on.agescale) {  
        if(AD) {
          ex <- (agemin - agemax) / (pb.max - pb.min)
          return(agemax + ex*pb) 	
        } else {
            ex <- (agemax - agemin) / (pb.max - pb.min)
            return(agemin + ex*pb)
        }
      } else
        return(pb)
  }

  if(newplot) {
    if(length(d.lab) == 0)
      d.lab <- paste0("depth (", set$depth.unit, ")")
    if(length(pb.lab) == 0)
      pb.lab <- ifelse(set$Bqkg, "210Pb (Bq/kg)", "210Pb (dpm/g)")

    if(length(d.lim) == 0)
      d.lim <- range(depths, set$supportedData[,3])
    if(rev.d)
      d.lim <- d.lim[2:1]
    if(rotate.axes)
      plot(0, type="n", ylim=d.lim, ylab=d.lab, xlim=pb2bp(pb.lim), xlab=pb.lab) else
        plot(0, type="n", xlim=d.lim, xlab=d.lab, ylim=pb2bp(pb.lim), ylab=pb.lab)
  }

  if(rotate.axes)
    rect(pb2bp(Pb-err), depths-thickness, pb2bp(Pb+err), depths, border=pbmeasured.col, lty=3) else
      rect(depths-thickness, pb2bp(Pb-err), depths, pb2bp(Pb+err), lty=3, border=pbmeasured.col)
    
  if(length(set$supportedData) > 0) {
    supp <- set$supportedData[,1]
    supperr <- set$supportedData[,2]
    suppd <- set$supportedData[,3]
    suppthick <- set$supportedData[,4]

    if(rotate.axes)
      rect(pb2bp(supp-supperr), suppd-suppthick, pb2bp(supp+supperr), suppd,
        border=supp.col, lty=3) else
        rect(suppd-suppthick, pb2bp(supp-supperr), suppd, pb2bp(supp+supperr),
          border=supp.col, lty=3) 
  }
}





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
#     stop("\n d.top should be higher than d.bottom", call.=FALSE)
#   t.top <- Bacon.Age.d(d.top, BCAD=F) - set$theta0
#   t.bottom <- Bacon.Age.d(d.bottom, BCAD=F) - set$theta0
#   multiply <- 500
#   if(set$Bqkg)
# 	multiply <- 10  
#   return(sup + ((phi / (.03114*multiply*dens) ) * (exp( -.03114*t.top) - exp(-.03114*t.bottom)) ) )
# } 





.plum.calib <- function(dat, set=get('info'), date.res=100, normal=set$normal, t.a=set$t.a, t.b=set$t.b,
                        delta.R=set$delta.R, delta.STD=set$delta.STD, ccdir="") {
  # read in the curves
  if(set$cc1=="IntCal20" || set$cc1=="\"IntCal20\"")
    cc1 <- read.table(paste0(ccdir, "3Col_intcal20.14C")) else
      cc1 <- read.csv(paste0(ccdir, set$cc1, ".14C"), header=FALSE, skip=11)[,1:3]
  if(set$cc2=="Marine20" || set$cc2=="\"Marine20\"")
    cc2 <- read.table(paste0(ccdir, "3Col_marine20.14C")) else
      cc2 <- read.csv(paste0(ccdir, set$cc2, ".14C"), header=FALSE, skip=11)[,1:3]
  if(set$cc3=="SHCal20" || set$cc3=="\"SHCal20\"")
    cc3 <- read.table(paste0(ccdir, "3Col_shcal20.14C")) else
      cc3 <- read.csv(paste0(ccdir, set$cc3, ".14C"), header=FALSE, skip=11)[,1:3]
  if(set$cc4=="ConstCal" || set$cc4=="\"ConstCal\"") cc4 <- NA else
    cc4 <- read.table(paste0(ccdir, set$cc4))[,1:3]

  if(set$postbomb != 0) {
    if(set$postbomb==1) bomb <- read.table(paste0(ccdir,"postbomb_NH1.14C"))[,1:3] else
      if(set$postbomb==2) bomb <- read.table(paste0(ccdir,"postbomb_NH2.14C"))[,1:3] else
        if(set$postbomb==3) bomb <- read.table(paste0(ccdir,"postbomb_NH3.14C"))[,1:3] else
          if(set$postbomb==4) bomb <- read.table(paste0(ccdir,"postbomb_SH1-2.14C"))[,1:3] else
            if(set$postbomb==5) bomb <- read.table(paste0(ccdir,"postbomb_SH3.14C"))[,1:3] else
              stop("cannot find postbomb curve #", set$postbomb, " (use values of 1 to 5 only)", call.=FALSE)
      bomb.x <- seq(max(bomb[,1]), min(bomb[,1]), by=-.1) # interpolate
      bomb.y <- approx(bomb[,1], bomb[,2], bomb.x)$y
      bomb.z <- approx(bomb[,1], bomb[,3], bomb.x)$y
      bomb <- cbind(bomb.x, bomb.y, bomb.z, deparse.level=0)
      if(set$postbomb < 4)
        cc1 <- rbind(bomb, cc1, deparse.level=0) else
          cc3 <- rbind(bomb, cc3, deparse.level=0)
  }

  ## use Gaussian or t (Christen and Perez Radiocarbon 2009) calibration
  if(round(set$t.b-set$t.a) !=1)
    stop("t.b - t.a should always be 1, check the manual", call.=FALSE)

  d.cal <- function(cc, rcmean, w2, t.a, t.b) {
    if(set$normal)
      cal <- cbind(cc[,1], dnorm(cc[,2], rcmean, sqrt(cc[,3]^2+w2))) else
        cal <- cbind(cc[,1], (t.b+ ((rcmean-cc[,2])^2) / (2*(cc[,3]^2 + w2))) ^ (-1*(t.a+0.5))) # student-t
    cal[,2] <- cal[,2]/sum(cal[,2])
    if(length(which(cal[,2]>set$cutoff)) > 5) # ensure that also very precise dates get a range of probabilities
      cal[which(cal[,2]>set$cutoff),] else {
        calx <- seq(min(cal[,1]), max(cal[,1]), length=100)
        caly <- approx(cal[,1], cal[,2], calx)$y
        cbind(calx, caly/sum(caly), deparse.level = 0)
      }
  }

  # now calibrate all dates
  calib <- list(d=dat[,4])

  for(i in 1:nrow(dat)) {
    dets <- c(NA, as.numeric(dat[i,-1])) # first entry is often not numeric    
    if(dets[9]==0 || dets[9] == 5) { # cal BP or 210Pb data
      x <- seq(dets[2]-(5*dets[3]), dets[2]+(5*dets[3]), by=5) # simplify, May 2019
      if(length(x) < 5 || length(x) > 100) # if too many resulting years, make 100 vals
        x <- seq(dets[2]-(5*dets[3]), dets[2]+(5*dets[3]), length=100)
      ccurve <- cbind(x, x, rep(0,length(x))) # dummy 1:1 curve
    } else
        if(dets[9]==1) ccurve <- cc1 else if(dets[9]==2) ccurve <- cc2 else
          if(dets[9]==3) ccurve <- cc3 else ccurve <- cc4

    delta.R <- set$delta.R; delta.STD <- set$delta.STD; t.a <- set$t.a; t.b <- set$t.b

    if(dets[9] > 0 && dets[9] < 5) { # only for C14 dates
      delta.R <- dets[5]
      delta.STD <- dets[6]
    }

    t.a <- dets[7]
    t.b <- dets[8]

    calib$probs[[i]] <- d.cal(ccurve, dets[2]-delta.R, dets[3]^2+delta.STD^2, t.a, t.b)

  }

  calib
}


.plumbacon.calib <- function(datPlum, dat, set=get('info'), date.res=100, normal=set$normal, t.a=set$t.a, t.b=set$t.b, delta.R=set$delta.R, delta.STD=set$delta.STD, ccdir="") {
  # read in the curves
  if(set$cc1=="IntCal20" || set$cc1=="\"IntCal20\"")
    cc1 <- read.table(paste0(ccdir, "3Col_intcal20.14C")) else
      cc1 <- read.csv(paste0(ccdir, set$cc1, ".14C"), header=FALSE, skip=11)[,1:3]
  if(set$cc2=="Marine20" || set$cc2=="\"Marine20\"")
    cc2 <- read.table(paste0(ccdir, "3Col_marine20.14C",sep="")) else
      cc2 <- read.csv(paste0(ccdir, set$cc2, ".14C"), header=FALSE, skip=11)[,1:3]
  if(set$cc3=="SHCal20" || set$cc3=="\"SHCal20\"")
    cc3 <- read.table(paste0(ccdir, "3Col_shcal20.14C")) else
      cc3 <- read.csv(paste0(ccdir, set$cc3, ".14C"), header=FALSE, skip=11)[,1:3]
  if(set$cc4=="ConstCal" || set$cc4=="\"ConstCal\"") cc4 <- NA else
    cc4 <- read.table(paste0(ccdir, set$cc4))[,1:3]

  if(set$postbomb != 0) {
    if(set$postbomb==1) bomb <- read.table(paste0(ccdir,"postbomb_NH1.14C"))[,1:3] else
      if(set$postbomb==2) bomb <- read.table(paste0(ccdir,"postbomb_NH2.14C"))[,1:3] else
        if(set$postbomb==3) bomb <- read.table(paste0(ccdir,"postbomb_NH3.14C"))[,1:3] else
          if(set$postbomb==4) bomb <- read.table(paste0(ccdir,"postbomb_SH1-2.14C"))[,1:3] else
            if(set$postbomb==5) bomb <- read.table(paste0(ccdir,"postbomb_SH3.14C"))[,1:3] else
              stop("cannot find postbomb curve #", set$postbomb, " (use values of 1 to 5 only)", call.=FALSE)
      bomb.x <- seq(max(bomb[,1]), min(bomb[,1]), by=-.1) # interpolate
      bomb.y <- approx(bomb[,1], bomb[,2], bomb.x)$y
      bomb.z <- approx(bomb[,1], bomb[,3], bomb.x)$y
      bomb <- cbind(bomb.x, bomb.y, bomb.z, deparse.level=0)
      if(set$postbomb < 4)
        cc1 <- rbind(bomb, cc1, deparse.level=0) else
          cc3 <- rbind(bomb, cc3, deparse.level=0)
  }

  ## use Gaussian or t (Christen and Perez Radiocarbon 2009) calibration
  if(round(set$t.b-set$t.a) !=1)
    stop("t.b - t.a should always be 1, check the manual", call.=FALSE)

  d.cal <- function(cc, rcmean, w2, t.a, t.b) {
    if(set$normal)
      cal <- cbind(cc[,1], dnorm(cc[,2], rcmean, sqrt(cc[,3]^2+w2))) else
        cal <- cbind(cc[,1], (t.b+ ((rcmean-cc[,2])^2) / (2*(cc[,3]^2 + w2))) ^ (-1*(t.a+0.5))) # student-t
    cal[,2] <- cal[,2]/sum(cal[,2])
    if(length(which(cal[,2]>set$cutoff)) > 5) # ensure that also very precise dates get a range of probabilities
      cal[which(cal[,2]>set$cutoff),] else {
        calx <- seq(min(cal[,1]), max(cal[,1]), length=100)
        caly <- approx(cal[,1], cal[,2], calx)$y
        cbind(calx, caly/sum(caly), deparse.level = 0)
      }
  }

  # now calibrate all dates
  calib <- list(d=dat[,4])
  if(ncol(dat)==4) { # only one type of dates (e.g., calBP, or all IntCal20 C14 dates)
    if(set$cc==0) {
      x <- seq(min(dat[,2])-(5*max(dat[,3])), max(dat[,2])+(5*max(dat[,3])), by=5) # simplify, May 2019
      if(length(x) > 100) # if too many resulting years, make 100 vals
        x <- seq(min(dat[,2])-(5*max(dat[,3])), max(dat[,2])+(5*max(dat[,3])), length=100)
      ccurve <- cbind(x, x, rep(0,length(x))) # dummy 1:1 curve
    } else {
        if(set$cc==1) ccurve <- cc1 else
          if(set$cc==2) ccurve <- cc2 else
            if(set$cc==3) ccurve <- cc3 else
              ccurve <- cc4
      }
    for(i in 1:nrow(dat))
      calib$probs[[i]] <- d.cal(ccurve, dat[i,2]-delta.R, dat[i,3]^2+delta.STD^2, set$t.a, set$t.b)
  } else
      for(i in 1:nrow(dat)) {
        dets <- c(NA, as.numeric(dat[i,-1])) # first entry is often not numeric
        if(dets[5]==0) {
          x <- seq(dets[2]-(5*dets[3]), dets[2]+(5*dets[3]), by=5) # simplify, May 2019
          if(length(x) < 5 || length(x) > 100) # if too many resulting years, make 100 vals
            x <- seq(dets[2]-(5*dets[3]), dets[2]+(5*dets[3]), length=100)
          ccurve <- cbind(x, x, rep(0,length(x))) # dummy 1:1 curve
        } else
            if(dets[5]==1) ccurve <- cc1 else if(dets[5]==2) ccurve <- cc2 else
              if(dets[5]==3) ccurve <- cc3 else ccurve <- cc4

        delta.R <- set$delta.R; delta.STD <- set$delta.STD; t.a <- set$t.a; t.b <- set$t.b
        if(length(dets) >= 7 && dets[5] > 0) { # the user provided age offsets; only for C14 dates
          delta.R <- dets[6]
          delta.STD <- dets[7]
        }

        if(length(dets) >= 9) { # the user provided t.a and t.b values for each date
          t.a <- dets[8]
          t.b <- dets[9]
          if(round(t.b-t.a) != 1)
            stop("t.b - t.a should always be 1, check the manual", call.=FALSE)
        }
        calib$probs[[i]] <- d.cal(ccurve, dets[2]-delta.R, dets[3]^2+delta.STD^2, t.a, t.b)
      }
  calib
}
