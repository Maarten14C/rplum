

plum.calib <- function(dat, set=get('info'), date.res=100, normal=set$normal, t.a=set$t.a, t.b=set$t.b, delta.R=set$delta.R, delta.STD=set$delta.STD, ccdir="") {
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
      bomb <- bomb[order(bomb[,1], decreasing=FALSE),]
      bomb.x <- seq(max(bomb[,1]), min(bomb[,1]), by=-.1) # interpolate
      bomb.y <- approx(bomb[,1], bomb[,2], bomb.x)$y
      bomb.z <- approx(bomb[,1], bomb[,3], bomb.x)$y
      bomb <- cbind(bomb.x, bomb.y, bomb.z, deparse.level=0)
      if(set$postbomb < 4)
        cc1 <- rbind(bomb, cc1, deparse.level=0) else
          cc3 <- rbind(bomb, cc3, deparse.level=0)
  }

  ## use Gaussian or t (Christen and Perez Radiocarbon 2009) calibration
  if(round(set$t.b-set$t.a) != 1)
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
  calib <- list(d=dat[,4], cc=dat[,9])

  for(i in 1:nrow(dat)) {
    dets <- c(NA, as.numeric(dat[i,-1])) # first entry is often not numeric
    if(dets[9] == 0 || dets[9] == 5) { # cal BP or 210Pb data
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
