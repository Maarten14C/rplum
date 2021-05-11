

### for running Plum, but is looked for by generic agedepth() function (through draw.pbmodelled()), so is included in the rbacon code
#' @name tmpbackground
#' @title calculate probabilities that Pb-210 data have reached background levels
#' @description Checks which of the Pb-210 data most likely have reached background levels and thus are below the detection limit Al (probabilities between 0 and 1)
#' @author Maarten Blaauw
#' @return a list of probabilities for each Pb-210 data point
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param Al The detection limit. Default \code{Al=0.1}.
#' @export
tmpbackground <- function(set=get('info'), Al=set$Al) {
  if(set$isplum) { # works with Pb-210 data only
    pb <- 0
    its <- nrow(set$output)
#    dets <- set$detsOrig[,c(2,6,3)] # we need maxdepth, mindepth, density
    dets <- set$dets[which(set$dets[,9] == 5),4:6] # should leave out any non-Pb data
    ps <- cbind(set$ps)
    for(i in 1:nrow(dets)) {
      As <- A.modelled(dets[i,1]-dets[i,2], dets[i,1], dets[i,3])
      if(set$ra.case == 2)
        ps <- set$ps[,i] else
          ps <- set$ps
      bg <- which((As - ps) <= Al) # which modelled data are at or below the detection limit?
      pb[i] <- length(bg) / its
    }
    return(pb)
  }
}

# function to read plum output files into memory
tmpPlum.AnaOut <- function(fnam, set=get('info')) {
  out <- read.table(fnam)
  n <- ncol(out)-1
  set$nPs  <- n
  set$TrPs <- nrow(out)
  set$phi  <- out[,1]
  set$ps   <- out[,2:(n+1)]
  set
}



# function to read output file into memory
tmpBacon.AnaOut <- function(fnam, set=get('info')) {
  out <- read.table(fnam)
  n <- ncol(out)-1
  set$n <- n
  set$Tr <- nrow(out)
  set$Us <- out[,n+1]
  set$output <- out[,1:n]
  set
}



#' @name Plum_runs
#' @title List the folders present in the current core directory.
#' @description Lists all folders located within the core's directory.
#' @details The directory is either "Plum_runs", "Cores" or a custom-named one.
#' @author Maarten Blaauw, J. Andres Christen
#' @return A list of folders
#' @param coredir The directory where the Bacon runs reside. Defaults to \code{coredir="Plum_runs"}.
#' @export
Plum_runs <- function(coredir=get('info')$coredir)
  list.files(coredir)



# Do a regression to determine the optimal number of supported data to use; the minimum is 3
check.equi <- function(dets, suggest=TRUE) {
  rawdata <- dets[,4] #210Pb
  rawsd <- dets[,5] #sd(210pb)
  # deps    <- dets[,2] #depth

  lendat <- length(rawdata)
  numdat <- as.integer(.5*lendat) # why 0.5?
  usedat <- rawdata[(lendat-3):lendat]
  usesd  <- rawsd[(lendat-3):lendat]
  usex   <- 1:length(usedat)
  usereg <- lm(usedat ~ usex, weights=1/(usesd^2))
  reg    <- coef(summary(usereg))[2,4]
  est    <- coef(summary(usereg))[1,1]
  coe    <- 3

  for(i in 1:numdat) {
    usedat <- rawdata[(lendat-3-i):lendat]
    usesd  <- rawsd[(lendat-3-i):lendat]
    usex   <- 1:length((lendat-3-i):lendat)
    usereg <- lm(usedat ~ as.numeric(scale(usex)), weights=1/(usesd^2))
    reg1   <- coef(summary(usereg))[2,4]
    est1   <- mean(usedat) #coef(summary(usereg))[1,1]
    if(reg1 > reg) {
      reg  <- reg1
      coe  <- (3+i)
      est  <- est1
    }
  }
  
  # from https://www.mail-archive.com/r-help@r-project.org/msg117700.html
  stopQuietly <- function() {
    blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "));
    stop(simpleError(blankMsg));
  }

  if(suggest) {
    ans <- readline(message("The regression process proposes using the last ", as.integer(coe), " data points as estimates of the supported activity, with a p-value of ", round((reg),3), ", OK? (Y/n) "))
    if(!(ans=="y" || ans=="")) {
      message("  OK. Please provide correct n.supp (as Plum option or in the .csv file).")
      stopQuietly()
    }
  }
#  return(c(coe, reg1))
  return(as.integer(coe))
}



# read the 210Pb dets file
read.dets.plum <- function(core, coredir, n.supp=c(), date.sample, sep=",", dec=".", cc=1, Bqkg=TRUE, ra.case=c(), suggest=TRUE) {

  # read the file. Removing the option to read and convert dat files because this is moot for Pb210
  csv.file <- paste0(coredir,  core, "/", core, ".csv")
  changed <- FALSE
  if(file.exists(csv.file)) {
    dets <- read.table(csv.file, header=TRUE, sep=sep)
    message("Reading ", csv.file)
  } else {
    if(file.exists(paste0(csv.file, ".txt"))) {
      file.rename(paste0(csv.file, ".txt"), csv.file)
      message("Removing .txt extension from .csv file")
      dets <- read.table(csv.file, header=TRUE, sep=sep)
      message("Reading", csv.file, "\n")
    } else
        message("No .csv file found. Please check if the name and/or location is correct")
  }

  # some general housekeeping: remove NAs and/or superfluous commas
  commas <- grep(",,", readLines(csv.file)) # check if there are too many commas
  if(length(!is.na(commas)) > 0) # often an artefact of spreadsheet programs
    stop("check the .csv file in a plain-text editor for 'orphan' commas", call.=FALSE)
  nas <- which(is.na(dets[,ncol(dets)]))
  if(length(nas) > 0)
    dets[nas,ncol(dets)] <- "" # replace NAs with empty cells
  # relations between the names of columns and their positions in the .csv file
  idColumn       <- 1 # ID
  depthColumn    <- 2 # depth
  rhoColumn      <- 3 # density
  plumdataColumn <- 4 # means of measurements
  stdColumn      <- 5 # their errors
  deltaColumn    <- 6 # sample thickness
  raColumn    <- 7 # if present
  sdRaColumn  <- 8 # if present

  #check that depths are in ascending order
  if(min(diff(dets[,depthColumn])) < 0) {
    message("Warning, the depths are not in ascending order, I will correct this.")
    dets <- dets[ order(dets[,depthColumn]),]
    changed <- TRUE
  }
  date.infile <- NA; nsupp.infile <- NA; racase.infile <- NA #; Bqkg.infile <- NA
  if(ncol(dets) == 6 || ncol(dets) == 8) # no additional information in file
    detsOrig <- dets else
      if(ncol(dets) == 7 || ncol(dets) == 9) { # additional information in file
        n <- ifelse(ncol(dets) == 7, 7, 9)
        detsOrig <- dets[,-n]
        if(length(dets[1,n]) > 0) # first entry is sampling date (in AD, e.g., 2019.4)
          if(!is.na(dets[1,n]))
            if(dets[1,n] != "")
              date.infile <- dets[1,n]
        if(length(dets[2,n]) > 0) # 2nd, how many 'tail-dates' to use to estimate supported
          if(!is.na(dets[2,n]))
            if(dets[2,n] != "")
              nsupp.infile <- dets[2,n]
        if(length(dets[3,n]) > 0) # 3rd, which radium case to use to estimate supported
          if(!is.na(dets[3,n]))
            if(dets[3,n] != "")
              racase.infile <- dets[3,n]
      } else
        stop(paste(csv.file, "should have between 6 and 9 columns. Please check."), call.=TRUE)

  date.asoption <- date.sample
  nsupp.asoption <- n.supp
  racase.asoption <- ra.case
  Bqkg.asoption <- Bqkg

  # now decide which options to use
  choice <- function(infile, asoption, string1, string2, testnumeric=TRUE, test=c()) {
    if(length(infile) == 0 || is.na(infile) || infile == "") {
      if(length(asoption) == 0) {
        if(length(test) > 0)
          ans <- test else
            ans <- readline(string2)
        if(testnumeric)
          if(grepl("^[0-9.][0-9]*[.]?[0-9]*[0-9.]$",ans[1]) == FALSE)
            if(grepl("^[0-9]+$", ans[1]) == FALSE)
              stop(cat(string1, "should be a numeric value "), call.=FALSE)
        chosen <- as.numeric(ans)
      } else
          chosen <- as.numeric(asoption)
    } else
        if(length(asoption) == 0)
          chosen <- as.numeric(infile) else
            chosen <- as.numeric(asoption)
    message("Using for ", string1, ": ", chosen)
    return(chosen)
  }

  if(core=="HP1C") {
    if(length(date.asoption) == 0) {
      message("Core HP1C was sampled in summer 2018, setting date.sample to 2018.5")
      date.sample <- 2018.5
    }
  } else
    date.sample <- choice(date.infile, date.asoption, "sampling date", "Please provide a date (in AD) for when the Pb210 samples were measured: ")

  # Now test different scenarios. Only check n.supp if ra.case < 2?
  if(ncol(dets) == 6) { # no radium, no information provided
    if(length(racase.asoption) == 0)
      message("No radium-226 data, setting ra.case to 0, using tail data to estimate supported Pb-210") else {
        if(racase.asoption != 0)
          message("Setting ra.case to 0 as no radium-226 data provided, using tail data to estimate supported Pb-210") else
            message("No radium-226 data, setting ra.case to 0")
        }
    ra.case <- 0
  }

  if(ncol(dets) == 7) { # no radium, information provided (thanks!)
    if(length(racase.asoption) == 0) {
      if(length(racase.infile) == 0 || is.na(racase.infile) || racase.infile > 0)
        message("Setting ra.case to 0 as no radium-226 data provided, using tail data to estimate supported Pb-210 (ra.case 0)") else
          message("No radium-226 data, using tail data to estimate supported Pb-210 (ra.case 0)")
      } else {
          if(racase.asoption != 0)
            message("Setting ra.case to 0 as no radium-226 data provided, using tail data to estimate supported Pb-210 (ra.case 0)") else
              message("No radium-226 data, using tail data to estimate supported Pb-210 (ra.case 0)")
        }
    ra.case <- 0
  }

  if(ncol(dets) == 8) { # radium, no information provided
    if(length(racase.asoption) == 0) {
      message("Radium-226 data provided. Should I assume constant (ra.case 1) or varying (ra.case 2) supported Pb-210? Note that using ra.case 2 will greatly increase the computing time and should only be used when clear patterns are observed in the radium-226 data.")
      ans <- readline(" Use ra.case (1 or 2):")
      if(ans == 1)
        ra.case <- 1 else
        if(ans == 2)
          ra.case <- 2 else
          stop("I do not understand this value for ra.case (should be 1 or 2). Please adapt the settings", call.=TRUE)
    } else {
        if(racase.asoption == 1)
          ra.case <- 1 else
          if(racase.asoption == 2)
            ra.case <- 2 else {
              message("Radium-226 data provided. Should I assume constant (ra.case 1) or varying (ra.case 2) supported Pb-210? Note that using ra.case 2 will greatly increase the computing time and should only be used when clear patterns are observed in the radium data.")
              ans <- readline("Use ra.case (1 or 2):")
              if(ans == 1)
                ra.case <- 1 else
                if(ans == 2)
                  ra.case <- 2 else
                    stop("I do not understand this value for ra.case (should be 1 or 2). Please adapt the settings", call.=TRUE)
            }
      }
  }

  if(ncol(dets) == 9) { # radium, information provided (fun!)
    if(length(racase.asoption) == 0) {
      if(length(racase.infile) == 0 || is.na(racase.infile) || racase.infile == 0 || racase.infile == "") {
        message("Radium data provided so ra.case cannot be 0 or empty. Should I assume constant (ra.case 1) or varying (ra.case 2) supported Pb-210? Note that using ra.case 2 will greatly increase the computing time and should only be used when clear patterns are observed in the radium-226 data. ")
        ans <- readline("Use ra.case (1 or 2):")
        if(ans == 1)
          ra.case <- 1 else
          if(ans == 2)
            ra.case <- 2 else
              stop("I do not understand this answer. Please adapt the settings", call.=TRUE)
      } else
        if(racase.infile == 1)
          ra.case <- 1 else
          if(racase.infile == 2)
            ra.case <- 2 else
              stop("I do not understand the radium-226 case value in the .csv file. Please adapt", call.=TRUE)
    } else {
      if(racase.asoption == 0) {
        message("Radium-226 data provided so ra.case cannot be 0. Should I assume constant (ra.case 1) or varying (ra.case 2) supported Pb-210? Note that using ra.case 2 will greatly increase the computing time and should only be used when clear patterns are observed in the radium-226 data.")
        ans <- readline("Use ra.case (1 or 2):")
        if(ans == 1)
          ra.case <- 1 else
          if(ans == 2)
            ra.case <- 2 else
              stop("I do not understand this value for ra.case (should be 1 or 2). Please adapt the settings", call.=TRUE)
      } else
        if(racase.asoption == 1) {
          message("ra.case 1")
          ra.case <- 1
        } else
          if(racase.asoption == 2) {
            message("ra.case 2")
            ra.case <- 2
          } else
            stop("I do not understand this value for ra.case (should be 1 or 2). Please adapt the settings", call.=TRUE)
    }
  }

  if(ra.case < 2) { # n.supp cannot be used if ra.case 2 (check with Marco that this is correct!)
    if(ra.case == 1)
      message("Besides using the radium data, the tail Pb-210 data can also be used to estimate supported Pb-210. ")
    n.supp <- choice(nsupp.infile, nsupp.asoption, "number of supported data", "",, check.equi(dets))
  }
  if(ra.case == 2)
    n.supp <- 0

  if(length(Bqkg) == 0 || !(Bqkg %in% c(0, 1))) {
    message("Assuming that the Pb units are in Bq/kg, Bqkg=1")
    Bqkg <- 1
  }

  # now put the chosen options into the .csv file if they differ from what's in there already
  choices <- c(date.sample, n.supp, ra.case, rep("", nrow(dets)-3)) # empty after line 4
  suggested.names <- c("labID", "depth(cm)","density(g/cm^3)","210Pb(Bq/kg)","sd(210Pb)","thickness(cm)", "226Ra(Bq/kg)", "sd(226Ra)", "settings")
  if(ra.case == 0) # then no radium columns
   suggested.names <- suggested.names[-(7:8)]

  if(ncol(dets) %in% c(6,8)) { # no data provided in the .csv file
    changed <- TRUE
    dets <- cbind(dets, choices)
  } else {
      current <- c(date.infile, nsupp.infile, racase.infile)
      if(length(is.na(current)) > 0 || length(choices[1:4] == current) < 4) # then update .csv file
        changed <- TRUE
      dets[,ncol(dets)] <- choices
    }

  if(changed) {
    message("Writing changes to ", csv.file)
    write.table(dets, csv.file, sep=paste0(sep, "\t"), dec=dec, row.names=FALSE, col.names=suggested.names, quote=FALSE)
  }

  # now identify supported Pb-210 data
  if(ncol(detsOrig) == 6) { # then repeat the Pb210 columns again
    raColumn <- 4
    sdRaColumn <- 5
    supportedData <- detsOrig[(nrow(detsOrig)-n.supp+1):nrow(detsOrig),c(raColumn, sdRaColumn, depthColumn, deltaColumn)]
    detsOrig <- detsOrig[1:(nrow(detsOrig)-n.supp),]
  } else
    if(ncol(detsOrig) == 8) { # columns 7 and 8 are supported data
      raColumn <- 7
      sdRaColumn <- 8
      supportedData <- detsOrig[,c(raColumn, sdRaColumn, depthColumn, deltaColumn)]
      detsOrig <- detsOrig[,-c(raColumn,sdRaColumn)]

      if(length(supportedData[is.na(supportedData)]) > 0) {
        message("Missing values are detected; the radium case is set to 1.")
        ra.case <- 1
        elim <- c() # get rid of data with NAs
        for(i in 1:nrow(supportedData))
          if(length(is.na(supportedData[i,])) > 0)
            elim <- c(elim, i)
        supportedData <- supportedData[-elim,]
      }

      if(length(n.supp) > 0)
        if(n.supp > 0) {
          raColumn <- 4
          sdRaColumn <- 5
          tmp <- detsOrig[(nrow(detsOrig)-n.supp+1):(nrow(detsOrig)),c(raColumn, sdRaColumn, depthColumn, deltaColumn)]
          names(tmp) <- colnames(supportedData)
          supportedData <- rbind(supportedData, tmp) # combine Ra-226 data with tail of Pb-210 data
          detsOrig <- detsOrig[1:(nrow(detsOrig)-n.supp),]
        }

    } else
      if( n.supp > 0 ) { # 9 columns, so radium data and information provided, and n.supp to be added
        raColumn <- 4
        sdRaColumn <- 5
        tmp <- detsOrig[(nrow(detsOrig)-n.supp+1):(nrow(detsOrig)),c(raColumn, sdRaColumn, depthColumn, deltaColumn)]
        names(tmp) <- colnames(supportedData)
        supportedData <- rbind(supportedData, tmp)
        detsOrig <- detsOrig[1:(nrow(detsOrig)-n.supp),]
      } else
          stop("Unexpected column names, order or values in dets file. \nPlease check the manual for how to produce a correct dets file.", call.=FALSE)

  # more sanity checks
  if(!is.numeric(dets[,plumdataColumn]) || !is.numeric(dets[,stdColumn]) || !is.numeric(dets[,depthColumn]))
    stop("unexpected values in dets file, I expected numbers. Check the manual.", call.=FALSE)
  if(!is.numeric(dets[,deltaColumn]) || !is.numeric(dets[,rhoColumn]) )
    stop("unexpected values in dets file, I expected numbers. Check the manual.", call.=FALSE)

  dets <- dets[,c(idColumn, plumdataColumn, stdColumn, depthColumn, deltaColumn, rhoColumn)]

  # find the plot limits
  if(ncol(detsOrig) == 6) {
    age.min <- min( c(detsOrig[,2]-(detsOrig[,6]/2)), c( detsOrig[,4]-detsOrig[,5]) )
    age.max <- max( c(detsOrig[,2]-(detsOrig[,6]/2)), c( detsOrig[,4]+detsOrig[,5]) )
  } else {
    age.min <- min( c(detsOrig[,2]-(detsOrig[,6]/2),detsOrig[,2]-(detsOrig[,6]/2)), c( detsOrig[,4]-detsOrig[,5],detsOrig[,7]-detsOrig[,8]) )
    age.max <- max( c(detsOrig[,2]-(detsOrig[,6]/2),detsOrig[,2]-(detsOrig[,6]/2)), c( detsOrig[,4]+detsOrig[,5],detsOrig[,7]+detsOrig[,8]) )
  }

  # plot the data
  layout(1)
  oldpar <- par(mar=c(3,3,1,1), mgp=c(1.5,.7,.0), bty="l")
  on.exit(par(oldpar))

  age.lim <- extendrange(c(age.min, age.max), f=0.01)
  dlim <- c(0, max(detsOrig[,depthColumn]))
  ylab <- ifelse(Bqkg, '210Pb (Bq/kg)', '210Pb (dpm/g)')
  plot(0, type='n', pch=16,col=c(rep('red',nrow(detsOrig)),rep('red',nrow(detsOrig))),
    cex=.3, ylab=ylab, xlab='depth(cm)', xlim=dlim, ylim=age.lim )

  rect(detsOrig[,2]-detsOrig[,6], detsOrig[,4]-detsOrig[,5],
    detsOrig[,2], detsOrig[,4]+detsOrig[,5],
    lty=3, border=4)
  if(ncol(detsOrig) > 6)
    rect(detsOrig[,2], detsOrig[,7]-detsOrig[,8],
      detsOrig[,2]-detsOrig[,6], detsOrig[,7]+detsOrig[,8],
      lty=3, border=2)

  return(list(dets, supportedData, ra.case, date.sample, detsOrig, n.supp, Bqkg))
}



#' @name Plum.cleanup
#' @title Remove files made to produce the current core's age-depth model.
#' @description Remove files .bacon, .out, .pdf, _ages.txt, and _settings.txt of current core.
#' @details If cores behave badly, you can try cleaning up previous runs and settings, by
#' removing files .bacon, .out, .pdf, _ages.txt, and _settings.txt of current core.
#' @return A message stating that the files and settings of this run have been deleted.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @author Maarten Blaauw, J. Andres Christen
#' @export
Plum.cleanup <- function(set=get('info')) {
  files <- c(paste0(set$prefix, ".bacon"), paste0(set$prefix, ".out"),
    paste0(set$prefix, ".pdf"), paste0(set$prefix, "_ages.txt"),
    paste0(set$coredir,set$core, "/", set$core, "_settings.txt"))
  for(i in files)
    if(file.exists(i))
      tmp <- file.remove(i)
  if(exists("tmp"))
    rm(tmp)
  message("Previous Plum runs of core ", set$core, " with thick=", set$thick, " deleted. Now try running the core again\n")
}



# read in default values, values from previous run, any specified values, and report the desired one. Internal function. Differs from Bacon.settings because it requires Pb-specific settings.
.plum.settings <- function(core, coredir, dets, thick, remember=TRUE, d.min, d.max, d.by, depths.file,
  slump, acc.mean, acc.shape, mem.mean, mem.strength, boundary, hiatus.depths, hiatus.max, hiatus.shape,
  BCAD, cc, postbomb, cc1, cc2, cc3, cc4, depth.unit, normal, t.a, t.b, delta.R, delta.STD, prob,
  defaults, runname, ssize, dark, MinAge, MaxAge, cutoff, age.res, after, age.unit,
  supportedData, date.sample, Al, phi.shape, phi.mean, s.shape, s.mean, ra.case, Bqkg, n.supp) {

  vals <- list(d.min, d.max, d.by, depths.file, slump, acc.mean, acc.shape, mem.mean, mem.strength, boundary, hiatus.depths, hiatus.max, BCAD, cc, postbomb, cc1, cc2, cc3, cc4, depth.unit, normal, t.a, t.b, delta.R, delta.STD, prob, age.unit)
  valnames <- c("d.min", "d.max", "d.by", "depths.file", "slump", "acc.mean", "acc.shape", "mem.mean", "mem.strength", "boundary", "hiatus.depths", "hiatus.max", "BCAD", "cc", "postbomb", "cc1", "cc2", "cc3", "cc4", "depth.unit", "normal", "t.a", "t.b", "delta.R", "delta.STD", "prob", "age.unit")
  #TODO: modificar para que acepte el vector de soportado y los valores propios de plum, como es el nombre de archivo de soportado y los parametros como "Al"
  extr <- function(i, def=deffile, pre=prevfile, exists.pre=prevf, rem=remember, sep=" ", isnum=TRUE) {
    if(length(vals[[i]]) > 0) # tmp
      if(any(is.na(vals[[i]]))) {
        ext.def <- strsplit(def[i], sep)[[1]]
        ext.def <- ext.def[-length(ext.def)] # remove description
        if(exists.pre) {
          ext.pre <- strsplit(pre[i], sep)[[1]]
          ext.pre <- ext.pre[-length(ext.pre)] # remove description
          if(def[i] == pre[i]) # values for dev and pre similar, no worries
            ext <- ext.pre else
              if(rem) {
                if(i==13) ifelse(ext.pre, "using BC/AD", "using cal BP") else
                if(i>2) message(" using previous run's value for ", valnames[i], ", ", ext.pre)
                ext <- ext.pre
              } else {
                  if(i==13) ifelse(ext.def, "using BC/AD", "using cal BP") else
                  if(i>2) message(" using default value for ", valnames[i], ", ", ext.def)
                  ext <- ext.def
                }
        } else ext <- ext.def

        if(any(ext=="NA") || any(is.na(ext))) NA else
          if(isnum) as.numeric(ext) else noquote(ext)
      } else
        if(isnum) as.numeric(vals[[i]]) else vals[[i]]
  }

  # read in default values and those of previous run if available
  deffile <- readLines(defaults, n=-1)
  prevfile <- paste0(coredir, core, "/", core, "_settings.txt")
  prevf <- FALSE
  if(file.exists(prevfile)) {
    prevfile <- readLines(prevfile, n=-1)
    if(length(prevfile) > 0) prevf <- TRUE
  }

  if(is.na(d.min) || d.min=="NA")
    d.min <- min(dets[,4])
  if(is.na(d.max) || d.max=="NA")
    #d.max <- max(dets[,4])
    d.max <- max(dets[,4]) # tmp
  if(length(acc.shape) < length(acc.mean))
    acc.shape <- rep(acc.shape, length(acc.mean)) else
      if(length(acc.shape) > length(acc.mean))
        acc.mean <- rep(acc.mean, length(acc.shape))
  if(length(mem.strength) < length(mem.mean))
    mem.strength <- rep(mem.strength, length(mem.mean)) else
      if(length(mem.strength) > length(mem.mean))
        mem.mean <- rep(mem.mean, length(mem.strength))

  ## produce/update settings file, and return the values
  prevfile <- file(paste0(coredir, core, "/", core, "_settings.txt"), "w")
  scat <- function(m, n="") cat(m, n, sep="", file=prevfile)
  cat(d.min, " #d.min\n", d.max, " #d.max\n", d.by, " #d.by\n",
    depths.file, " #depths.file\n", slump, " #slump\n", sep="", file=prevfile)
  for(i in acc.mean) scat(i, " "); scat("#acc.mean\n")
  for(i in acc.shape) scat(i, " "); scat("#acc.shape\n", "")
  for(i in mem.mean) scat(i, " "); scat("#mem.mean\n", "")
  for(i in mem.strength) scat(i, " "); scat("#mem.strength\n", "")
  for(i in boundary) scat(i, " "); scat("#boundary\n", "")
  for(i in hiatus.depths) scat(i, " "); scat("#hiatus.depths\n", "")
  for(i in hiatus.max) scat(i, " "); scat("#hiatus.max\n", "")
  #for(i in hiatus.max) scat(i, " "); scat("#hiatus.max\n", "") # redundant
  cat(BCAD, " #BCAD\n", cc, " #cc\n", postbomb, " #postbomb\n",
    cc1, " #cc1\n", cc2, " #cc2\n", cc3, " #cc3\n", cc4, " #cc4\n",
    depth.unit, " #depth.unit\n", normal, " #normal\n", t.a, " #t.a\n", t.b, " #t.b\n",
    delta.R, " #delta.R\n", delta.STD, " #d.STD\n", prob, " #prob\n", age.unit, "#age.unit\n", sep="", file=prevfile)

  cat(date.sample, " #date.sample\n", Al, " #Al\n", phi.shape, " #phi.shape\n", phi.mean, " #phi.mean\n",
    s.shape, " #s.shape\n", s.mean, " #s.mean\n", ra.case, " #ra.case\n", Bqkg, " #Bqkg\n", sep="", file=prevfile)

  cat(n.supp, " #n.supp\n", sep="", file=prevfile);

  close(prevfile)

  if(length(MinAge) == 0)
    MinAge <- min(1950 - as.integer(format(Sys.time(), "%Y")), round(dets[,2] - (5*dets[,3])))
  if(length(MaxAge) == 0)
    MaxAge <- max(1e6, round(dets[,2] + (5*dets[,3])))

  theta0 <- 1950 - date.sample

  list(core=core, thick=thick, dets=dets, d.min=d.min, d.max=d.max, coredir=coredir, # was coredir=core
    d.by=d.by, depths.file=depths.file, slump=slump,
    acc.mean=acc.mean, acc.shape=acc.shape, mem.mean=mem.mean,
    mem.strength=mem.strength, boundary=boundary,
    hiatus.depths=hiatus.depths, hiatus.max=hiatus.max,
    BCAD=BCAD, cc=cc, postbomb=postbomb,
    cc1=cc1, cc2=cc2, cc3=cc3, cc4=cc4, depth.unit=noquote(depth.unit), unit=depth.unit, age.unit=noquote(age.unit), normal=normal,
    t.a=t.a, t.b=t.b, delta.R=delta.R, delta.STD=delta.STD, prob=prob, date=date(),
    runname=runname, ssize=ssize, dark=dark, MinAge=MinAge, MaxAge=MaxAge,
    cutoff=cutoff, age.res=age.res, after=after,
    supportedData=supportedData, theta0=theta0, Al=Al, phi.shape=phi.shape, phi.mean=phi.mean, s.shape=s.shape, s.mean=s.mean, ra.case=ra.case, Bqkg=Bqkg)
}



#function to merge dets of plum and bacon data
merge.dets <- function(detsPlum, detsBacon, delta.R, delta.STD, t.a, t.b, cc) {
  if(ncol(detsBacon) >= 5) {
    cc <- detsBacon[,5]
    detsBacon <- detsBacon[,-5]
  } else
      cc <- array(cc, dim=c(nrow(detsBacon),1))

  if(ncol(detsBacon) < 9 ) {
    for(i in(ncol(detsBacon)+1):9) {
      if(i==5) {
        col <- array(delta.R, dim=c(nrow(detsBacon),1))
      } else if(i==6) {
        col <- array(delta.STD, dim=c(nrow(detsBacon),1))
      } else if(i==7) {
        col <- array(t.a, dim=c(nrow(detsBacon),1))
      } else if(i==8) {
        col <- array(t.b, dim=c(nrow(detsBacon),1))
      } else if(i==9) {
        col <- cc
      }
      detsBacon <- cbind(detsBacon, col)
    }
    colnames(detsBacon) <- c("labID", "X210Pb.Bq.kg.", "sd.210Pb.", "depth.cm.", "thickness.cm.", "density.g.cm.3.",  "t.a", "t.b", "cc")
  }

  if(ncol(detsPlum) < 9) {
    for(i in (ncol(detsPlum)+1):9){
      if(i==5) {
        col <- array(delta.R, dim=c(nrow(detsPlum),1))
      } else if(i==6) {
        col <- array(delta.STD, dim=c(nrow(detsPlum),1))
      } else if(i==7) {
        col <- array(t.a, dim=c(nrow(detsPlum),1))
      } else if(i==8) {
        col <- array(t.b, dim=c(nrow(detsPlum),1))
      } else if(i==9) {
        col <- array(5, dim=c(nrow(detsPlum),1))
      }
      detsPlum <- cbind(detsPlum, col)
    }
    colnames(detsPlum) <- c("labID", "X210Pb.Bq.kg.", "sd.210Pb.", "depth.cm.", "thickness.cm.", "density.g.cm.3.",  "t.a", "t.b", "cc")
  }

  dets <- rbind(detsPlum, detsBacon, make.row.names=FALSE)
  dets <- dets[order(dets[,4]),]
}



# write files to be read by the main Bacon age-depth modelling function. Has plum-specific settings so differs from write.Bacon.file
write.plum.file <- function(set=get('info')) {

  if(length(set$slump) > 0) {
    dets <- set$slumpdets
    hiatus.depths <- set$slumphiatus
    boundary <- set$slumpboundary
  } else {
    dets <- set$dets
    hiatus.depths <- set$hiatus.depths
    boundary <- set$boundary
  }

  depthColumn <- 4 # column of the depths
  if(is.na(set$d.min) || set$d.min < min(dets[,depthColumn])) { # repeat relevant row, change error and depth
    dets <- rbind(dets[which(dets[,depthColumn] == min(dets[,depthColumn]))[1],], dets, make.row.names=FALSE)
    dets[1,1] <- NA # calling this "d.min" causes issues
    dets[1,3] <- max(1e5, 1e3*dets[,4], 1e3*dets[,3])
    dets[1,depthColumn] <- set$d.min
  }

  if(is.na(set$d.max) || set$d.max > max(dets[,depthColumn])) { # repeat relevant row, change error and depth
    dets <- rbind(dets, dets[which(dets[,depthColumn] == max(dets[,depthColumn]))[1],], make.row.names=FALSE)
    dets[nrow(dets),1] <- NA # calling this "d.max" causes issues
    dets[nrow(dets),3] <- max(1e5, 1e3*dets[,4], 1e3*dets[,3])
    dets[nrow(dets),depthColumn] <- set$d.max
  }

  supportedData <- set$supportedData

  fl <- file(set$bacon.file, "w")
  cat("## Ran on", set$date, "\n\n", file=fl)
  cat("Cal 0 : ConstCal;\nCal 1 : ",
  if(set$cc1=="IntCal20" || set$cc1=="\"IntCal20\"") "IntCal20"
    else noquote(set$cc1), ", ", set$postbomb, ";\nCal 2 : ",
  if(set$cc2=="Marine20" || set$cc2=="\"Marine20\"") "Marine20"
    else noquote(set$cc2), ";\nCal 3 : ",
  if(set$cc3=="SHCal20" || set$cc3=="\"SHCal20\"") "SHCal20"
    else noquote(set$cc3), ", ", set$postbomb, ";",
  if(set$cc4=="ConstCal" || set$cc4=="\"ConstCal\"") set$cc4 <- c()
    else
      paste0("\nCal 4 : GenericCal, ", set$cc4, ";"), sep="", file=fl)
  cat("\nCal 4 : ConstCal;", sep="", file=fl)
  cat("\n##          alPhi mPhi  alS  mS     Al   theta0  Radon_case  supported_data_file", file=fl)
  cat("\nCal 5 : Plum, ", set$phi.shape, ", ",  set$phi.mean, ", ",  set$s.shape, ", ", set$s.mean, ", ", set$Al, ", ", set$theta0, ", ",
        set$ra.case, ", ", set$plum.file,";", sep="", file=fl)
  cat("\n##    ", colnames(dets), " ... Plum: 210Pb data",sep=", ", file=fl)

  # we need to send the dets with all columns so pre-processing is needed
  for( i in 1:nrow(dets) ) {
    cat( "\nDet ", i-1, " : ", as.character(dets[i,1]),
        " , ", dets[i,2],
        ", ", dets[i,3],
        ", ", dets[i,4],
        ", ", dets[i,5],
        ", ", dets[i,6],
        ", ", dets[i,7],
        ", ", dets[i,8],
        ", ", dets[i,9],
        ";", sep="", file=fl)
  }

  if(!is.na(hiatus.depths[1])) {
    if(is.null(boundary[1]))
      message("\n  Hiatus set at depth(s) ", hiatus.depths, "\n") else
        message("\n  Boundary set at depth(s) ", boundary, "\n")
    if(length(set$acc.shape)==1)
      set$acc.shape <- rep(set$acc.shape, length(hiatus.depths)+1)
    if(length(set$acc.mean)==1)
      set$acc.mean <- rep(set$acc.mean, length(hiatus.depths)+1)
    if(length(set$hiatus.max)==1)
      set$hiatus.max <- rep(set$hiatus.max, length(hiatus.depths))
    assign_to_global("info", set)
    cat("\n\n### Depths and priors for fixed hiatuses, in descending order",
      "\n##### cm  alpha beta      ha     hb", file=fl)
    for(i in length(hiatus.depths):1)
      cat("\nHiatus ", i-1, ":  ", hiatus.depths[i], ",  ", set$acc.shape[i+1],
        ",  ", set$acc.shape[i+1]/set$acc.mean[i+1], ",  ", .1, # last value (h.a) was NA but this conflicts with setting initial values for hiatus length
        ",  ", set$hiatus.max[i], ";", sep="", file=fl)
  }

  cK <- set$d.min+(set$thick*set$K)
  ### final parameters - dmax now calculated as dmin+(dC*K)
  if( is.na(set$seed) ) {
  wrapup <- paste0("\n\n##\t\t K   MinAge   MaxAge   th0   th0p   w.a   w.b   alpha  beta  dmin  dmax",
    "\nBacon 0: ", ifelse(set$normal, "FixNor", "FixT"), ", ", set$K,
    ",  ", set$theta0-.02, ",  ", 26500, ",  ", set$theta0-0.01, ",  ", set$theta0+0.01,
    ",  ", set$mem.strength*set$mem.mean, ",  ", set$mem.strength*(1-set$mem.mean),
    ",  ", set$acc.shape[1], ",  ", set$acc.shape[1]/set$acc.mean[1], ", ", set$d.min,
    ", ", cK, ";\n")
  } else {
    wrapup <- paste0("\n\n##\t\t K   MinAge   MaxAge   th0   th0p   w.a   w.b   alpha  beta  dmin  dmax  seed",
      "\nBacon 0: ", ifelse(set$normal, "FixNor", "FixT"), ", ", set$K,
      ",  ", set$theta0-.02, ",  ", 26500, ",  ", set$theta0-0.01, ",  ", set$theta0+0.01,
      ",  ", set$mem.strength*set$mem.mean, ",  ", set$mem.strength*(1-set$mem.mean),
      ",  ", set$acc.shape[1], ",  ", set$acc.shape[1]/set$acc.mean[1], ", ", set$d.min,
      ", ", cK, ", ", set$seed, ";\n")
  }
  cat(wrapup, file=fl)
  close(fl)

  fl <- file(set$plum.file, "w")
  if(length(supportedData) > 0) # if the .plum file has NA NA, then the output files have no lines
    for(i in 1:nrow(supportedData)) {
      for(j in 1:2) 
        cat(supportedData[i,j], " ", file=fl)
      cat("\n", file=fl)
    }	
  close(fl)
  # we have to check that there are no NAs in the .plum file
}
