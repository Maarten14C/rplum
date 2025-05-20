# updated to rbacon 3.5.2, removed temporary commands within fromrbacon.R

# set.initvals from rbacon doesn't work as expected in rplum. The function makes the bottom-left panel active and initial age-depth points can be selected, but the selected initial values do not run as expected. Probably because additional initvals are also required for Pb-210 pars?

# write an R package to download and plot climate data (grip, ngrip, gisp2, hulu, cariaco, EPICA, ...) working name icecream, pickles, or cream. check pangaear package, also check what rioja provides

# do: add more guidance on acc.mean - what type of site is it? option to enter supported data as file (instead of in parent .csv file), change column order in .csv file??? Adapt default value of dark? .01 works well if a Pb core also has C14 dates. check par righthand toppanel as too much space, A.rng and Ai in calibrate.plum.plot cannot be saved to info (needed to provide post-run info on fit 210Pb data), is it OK that d.min is set to 0 by default?

#' @name Plum
#' @title Main 210Pb age-depth modelling function
#' @description This is the main age-depth modelling function of the rplum package for 210Pb age-modelling.
#' @details Plum is an approach to age-depth modelling that uses Bayesian statistics in order to reconstruct
#' accumulation histories for 210Pb-dated deposits by taking into account prior information, and can combine 210Pb, radiocarbon and other dates (Aquino et al. 2018).
#'
#' Plum handles 210Pb and other dated depths within in a core, by dividing a core into many thin vertical sections (by default of \code{thick=1} cm thickness),
#' and through millions of Markov Chain Monte Carlo (MCMC) iterations estimates the flux of 210Pb and supported 210Pb, as well as 
#' the accumulation rate (in years/cm; so more correctly, sedimentation times) for each of these sections.
#' Combined with an estimated starting date for the first section, these accumulation rates and values for 210Pb then form the age-depth and 210Pb model.
#' The accumulation rates are constrained by prior information on the accumulation rate (\code{acc.mean, acc.shape)} and its
#' variability between neighbouring depths, or "memory" (\code{mem.mean, mem.strength}). Hiatuses can be introduced as well,
#' also constrained by prior information (\code{hiatus.max}). The 210Pb flux (phi) and supported 210Pb (s) are constrained by priors \code{phi.mean, phi.shape, s.mean and s.shape}.  
#'
#' Although Plum was developed for 210Pb dates, it can also include absolute dates (e.g., 14C, OSL, tephra or other dates on a calendar scale).
#' Radiocarbon dates should be calibrated using either IntCal20
#' (for terrestrial northern hemisphere material; Reimer et al., 2020), Marine20 (for marine dates; Hughen et al., 2020),
#' SHCal20 (for southern hemisphere dates; Hogg et al., 2020) or any other calibration curve (see below), while modern 14C
#' dates are calibrated using one of the post-bomb calibration curves (NH1, NH2 or NH3 for the northern hemisphere,
#' SH1-2 or SH3 for the southern hemisphere; Hua et al., 2022). See \url{http://calib.org/CALIBomb/} if you are unsure which
#' postbomb curve you need. If Plum finds postbomb dates (negative 14C ages) and you haven't specified a postbomb curve,
#' you will be prompted. Provide postbomb curves as, e.g., \code{postbomb=1} for the NH1 postbomb curve (2 for NH2, 3 for NH3,
#' 4 for SH1-2, 5 for SH3).
#'
#' For calendar dates, i.e. dates that are already on the calendar scale and thus should not be calibrated, set\code{cc=0}. 
#' Plum also needs the date of sampling, in AD (\code{date.sample}). 
#'
#' rplum works by calling the rbacon package. Since version 3.1.0, Bacon can also handle younger-than and older-than ages, with the model aiming to either go 'above'
#' or 'below' such dates as requested. If the resulting combination of parameters becomes problematic (e.g., no initial
#' combination of parameters can be found that obeys the priors or is in chronological order), then the output will often be wrong.
#' If so, using the function set.initvals could help.
#'
#' By default, the initial MCMC values of the Bacon age-depth model (upper ages and accumulation rate for each model section)
#' are estimated randomly. Since version 3.1.0, these starting values can also be provided in a file with extension _bacon.init,
#' placed within the core's folder. This file will need to have two rows, each for one of the two initial sets of parameters required
#' (the t-walk requires two starting estimates for all MCMC parameters).
#' If such a file is found (and correctly formatted), Bacon will use the values within this file
#' as starting points for the MCMC run. See function set.initvals for more information.
#'
#' @param core Name of the core, given using quotes. Defaults to one of the cores provided with rplum, \code{core="HP1C"} also reported by Aquino-Lopez et al. (2018). Also available is LL14, a core kindly provided by Dr Lysanna Anderson (USGS). LL14 has ra-226 data (so can be run with \code{ra.case=1} or \code{ra.case=2}, see below), and also has additional C-14 and cal BP data (these can be added using \code{otherdates="LL14_14C.csv"}). The original LL14 core has more 14C data than provided here (for reasons of brevity).
#' To run your own core, produce a .csv file with the dates as outlined in the manual, add a folder with the core's name to the default directory for cores (see \code{coredir}), and save the .csv file there. For example, the file's location and name could be \code{Plum_runs/MyCore/MyCore.csv}. Then run Plum as follows: \code{Plum("MyCore")}.
#' Note that for Pb-210 data, the depth in the .csv should be the bottom of the slice, not the mid-point. (For any non-Pb data, depths are the midpoints of their slices). Also make sure that the thickness and density are given correctly for each Pb-210 data point.
#' @param thick Plum will divide the core into sections of equal thickness specified by thick (default \code{thick=1}).
#' @param otherdates Name of (optional) file with radiocarbon dates. This file should have the same format as the one used for rbacon. For example, \code{Bacon("LL14", otherdates="LL14_14C.csv")}.
#' @param coredir Folder where the core's files \code{core} are and/or will be located. This will be a folder with the core's name, within either the folder \code{coredir='Plum_runs/'}, or the folder Cores/ if it already exists within R's working directory, or a custom-built folder. For example, use \code{coredir="."} to place the core's folder within the current working directory, or \code{coredir="F:"} if you want to put the core's folder and files on a USB drive loaded under F:.
#' Thinner (and thus more) sections will result in smoother age-models, but too many sections can cause `run-away' models.
#' @param phi.shape Shape parameter of the prior gamma distribution used for the influx of Pb-210 to the sediment, default \code{phi.shape=2}.
#' @param phi.mean Mean parameter of the prior gamma distribution used for the influx of Pb-210 to the sediment, default \code{phi.mean=50}.
#' @param s.shape Shape parameter of the prior gamma distribution used for the supported Pb-210 to the sediment, default \code{s.shape=5}.
#' @param s.mean Mean parameter of the prior gamma distribution used for the supported Pb-210 to the sediment, default \code{s.mean=10}.
#' @param Al Parameter used to limit the chronologies described in Aquino-Lopez et al. (2018) for the minimum distinguishable unsupported activity; default \code{Al=0.1}.
#' @param date.sample Date (in calendar years, e.g., AD 2023) at which the core was measured for Pb-120. This date will be used as a surface date and is assumed to have no uncertainty. If the date is not provided (in the .csv file or as \code{date.sample}), Plum will ask for it.
#' @param n.supp This value will delete n.supp data points from the deepest part of the core, and these points will then be used exclusively to estimate the supported activity. If this option is used, a constant supported Pb-210 will be assumed, \code{n.supp=-1}.
#' @param remove.tail Whether or not to remove the tail measurements when plotting. Sometimes automated removal might go wrong, or additional dates exist further down, so then this option can be used to avoid removing the tail 210Pb measurements. Is set to FALSE if there are non-210Pb data further down the core.
#' @param ra.case How to use radium-226 measurements if they are provided in the core's .csv file. 1 = assume constant radium, 2 = assume varying radium and use the radium measurements as individual estimates of supported Pb-210. If no radium measurements are present, use \code{ra.case=0}.
#' @param Bqkg This variable indicates whether total Pb-210 is expressed in Bq/kg (default; \code{Bqkg=TRUE}) or dpm/g if set to FALSE.
#' @param seed Seed used for C++ executions; if it is not assigned then the seed is set by system. Default \code{seed=NA}.
#' @param prob Confidence interval to report. This should lie between 0 and 1, default \code{prob=0.95} (95 \%).
#' @param d.min Minimum depth of age-depth model (use this to extrapolate to depths higher than the top dated depth).
#' @param d.max Maximum depth of age-depth model (use this to extrapolate to depths below the bottom dated depth).
#' @param d.by Depth intervals at which ages are calculated. Defaults to \code{d.by=1}.
#' @param depth.unit Units of the depths. Defaults to \code{depth.unit="cm"}.
#' @param age.unit Units of the ages. Defaults to \code{age.unit="yr"}.
#' @param unit Deprecated and replaced by \code{depth.unit}.
#' @param depths By default, Plum will calculate the ages for the depths \code{d.min} to \code{d.max} in steps of \code{d.by}.
#' Alternative depths can be provided as, e.g., \code{depths=seq(0, 100, length=500)} or as a file, e.g., \code{depths=read.table("CoreDepths.txt"}. See also \code{depths.file}
#' @param depths.file By default, Plum will calculate the ages for the depths \code{d.min} to \code{d.max} in steps of \code{d.by}.
#' If \code{depths.file=TRUE}, Plum will read a file containing the depths for which you require ages.
#' This file, containing the depths in a single column without a header, should be stored within \code{coredir},
#' and its name should start with the core's name and end with `_depths.txt'. Then specify \code{depths.file=TRUE} (default \code{FALSE}). See also \code{depths}.
#' @param acc.shape The prior for the accumulation rate consists of a gamma distribution with two parameters.
#' Its shape is set by acc.shape (default \code{acc.shape=1.5}; higher values result in more peaked shapes).
#' @param acc.mean The accumulation rate prior consists of a gamma distribution with two parameters. Its mean is set by acc.mean (default \code{acc.mean=10} yr/cm (or whatever age or depth units are chosen),
#' which can be changed to, e.g., 5, 10 or 50 for different kinds of deposits). Multiple values can be given in case of hiatuses or boundaries, e.g., Plum(hiatus.depths=23, acc.mean=c(5,20))
#' @param mem.strength The prior for the memory (dependence of accumulation rate between neighbouring depths) is a beta distribution, which looks much like the gamma distribution.
#'  but its values are always between 0 (no assumed memory) and 1 (100\% memory). Its default settings of \code{mem.strength=10}
#'  (higher values result in more peaked shapes) allow for a large range of posterior memory values.
#' @param mem.mean The prior for the memory is a beta distribution, which looks much like the gamma distribution but
#' its values are always between 0 (no assumed memory) and 1 (100\% memory). Its default settings of \code{mem.mean=0.5}
#' allow for a large range of posterior memory values.
#' @param boundary The assumed depths of any boundary, which divides sections of different accumulation rate regimes (e.g., as indicated by major change in the stratigraphy). No hiatus is assumed between these sections, and memory is reset crossing the boundary. Different accumulation priors can be set for the sections above and below the boundary, e.g., \code{acc.mean=c(5, 20)}. See also \code{hiatus.depths}, \code{mem.mean}, \code{acc.mean} and \code{acc.shape}. Setting many boundaries might not work, and having more than one boundary per model section (see \code{'thick'}) might not work either.
#' @param hiatus.depths The assumed depths for any hiatus should be provided as, e.g.,
#' \code{hiatus.depths=20} for one at 20cm depth, and \code{hiatus.depths=c(20,40)} for two hiatuses at 20 and 40 cm depth.
#' @param hiatus.max The prior for the maximum length of the hiatus. Hiatus length is a uniform distribution, with equal probabilities between 0 and \code{hiatus.max} yr (or whatever other \code{age.unit} is chosen).
#' @param add Add a value to the maximum hiatus length if a boundary is chosen. Defaults to 100 yr (or whatever other age unit is chosen). Can be adapted if Plum complains that the parameters are out of support.
#' @param after Sets a short section above and below hiatus.depths within which to calculate ages. For internal calculations - do not change.
#' @param cc Calibration curve for C-14 dates: \code{cc=1} for IntCal20 (northern hemisphere terrestrial), \code{cc=2} for Marine20 (marine),
#' \code{cc=3} for SHCal20 (southern hemisphere terrestrial). For dates that are already on the cal BP scale use \code{cc=0}.
#' @param cc1 For northern hemisphere terrestrial 14C dates (IntCal20).
#' @param cc2 For marine 14C dates (Marine20).
#' @param cc3 For southern hemisphere 14C dates (SHCal20).
#' @param cc4 Use an alternative curve (3 columns: cal BP, 14C age, error, separated by white spaces and saved as a plain-text file). See \code{cc.dir}.
#' @param cc.dir Directory where the calibration curves for C14 dates \code{cc} are located. By default \code{cc.dir=""} since they are loaded into R's memory.
#' For example, use \code{cc.dir="."} to choose current working directory, or \code{cc.dir="Curves/"} to choose sub-folder \code{Curves/}. Note that all calibration curves should reside in the same directory. If you want to add a custom-built curve, put it in the directory where the default calibration curves are (probably \code{list.files(paste0(.libPaths(), "/IntCal/extdata/"))}).
#' Alternatively produce a new folder, and add your curve as well as the default calibration curves there (cc1, cc2 and cc3; e.g., \code{write.table(copyCalibrationCurve(1), "./3Col_intcal20.14C", sep="\t")}.)
#' @param postbomb Use a postbomb curve for negative (i.e. postbomb) 14C ages. \code{0 = none, 1 = NH1, 2 = NH2, 3 = NH3, 4 = SH1-2, 5 = SH3}
#' @param F14C Radiocarbon ages can be provided as F14C values. If doing so, please indicate here which dates were entered as F14C (e.g., if the first 4 dates are in F14C, write \code{F14C=1:4}). The F14C values in your .csv file will then be replaced by their corresponding C14 ages.
#' @param pMC Radiocarbon ages can be provided as pMC values. If doing so, please indicate here which dates were entered as pMC (e.g., if the first 4 dates are in pMC, write \code{pMC=1:4}). The pMC values in your .csv file will then be replaced by their corresponding C14 ages.
#' @param delta.R Mean of core-wide age offsets (e.g., regional marine offsets).
#' @param delta.STD Error of core-wide age offsets (e.g., regional marine offsets).
#' @param t.a The dates are treated using the t distribution by default (\code{normal=FALSE}).
#' The t-distribution has two parameters, t.a and t.b, set at 3 and 4 by default (see Christen and Perez, 2009).
#' If you want to assign narrower error distributions (more closely resembling the normal distribution), set t.a and t.b at for example 33 and 34 respectively (e.g., for specific dates in your .csv file).
#' For symmetry reasons, t.a must always be equal to t.b-1.
#' @param t.b The dates are treated using the t distribution by default (\code{normal=FALSE}).
#' The t-distribution has two parameters, t.a and t.b, set at 3 and 4 by default (see Christen and Perez, 2010).
#' If you want to assign narrower error distributions (more closely resembling the normal distribution), set t.a and t.b at for example 33 and 34 respectively (e.g., for specific dates in your .csv file).
#' For symmetry reasons, t.a must always be equal to t.b-1.
#' @param normal By default, Plum uses the t-distribution to treat the dates. Use \code{normal=TRUE} to use the normal/Gaussian distribution. This will generally give higher weight to the dates.
#' @param suggest If initial analysis of the data indicates abnormally slow or fast accumulation rates, Plum will suggest to change the prior.
#'  Also, if the length of the core would cause too few or too many sections with the default settings, Plum will suggest an alternative section thickness \code{thick}, and it will suggest approaches to estimating supported Pb-120. 
#'  Accept these suggested alternative settings by typing "y" (or "yes please" if you prefer to be polite), or leave as is by typing "n" (or anything else, really). To get rid of these suggestions, use \code{suggest=FALSE}.
#' @param reswarn Plum will warn you if the number of sections lies outside the safe range (default between 10 and 200 sections;
#' \code{reswarn=c(10,200)}). Too few sections could lead to an `elbowy' model while with too many sections the modelling process can get lost,
#'  resulting in age-models far away from the dated depths.
#' @param remember Plum will try to remember which settings you have applied to your cores (default \code{remember=TRUE}). If you run into inconsistencies or other problems,
#' try running your core again with \code{remember=FALSE}, or, start cleanly by typing \code{Plum.cleanup()}.
#' @param ask By default Plum will ask you to confirm that you want to run the core with the provided settings. Disable this using \code{ask=FALSE} (e.g., for batch runs).
#' @param run In order to load an existing Plum run instead of producing a new one, you can use \code{run=FALSE}.
#' @param defaults Name of the file containing settings for the core. For internal use only - do not change.
#' @param sep Separator between the fields of the plain text file containing the dating information. Default \code{sep=","}.
#' @param dec Character for decimal points. Default to \code{dec="."}.
#' @param runname Text to add to the corename for specific runs, e.g., \code{runname="MyCore_Test1"}.
#' @param slump Upper and lower depths of any sections of assumed abrupt accumulation, that require excising before age-modelling (and adding after age-modelling). Requires pairs of depths, e.g., \code{slump=c(10,15,60,67)} for slumps at 67-60 and 15-10 cm core depth.
#' @param BCAD The calendar scale of graphs and age output-files is in cal BP (calendar or calibrated years before the present, where the present is AD 1950) by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param ssize The approximate amount of iterations to store at the end of the MCMC run. Default 2000; decrease for faster (but less reliable) runs or increase for cores where the MCMC mixing (panel at upper-left corner of age-model graph) appears problematic.
#' @param th0 Starting years for the MCMC iterations.
#' @param burnin Amount of initial, likely sub-optimal MCMC iterations that will be removed.
#' @param MinAge Deprecated - use youngest.age instead.
#' @param youngest.age Minimum age limit for Bacon runs, default at current year in cal BP. To set plot limits, use \code{age.min} instead.
#' @param MaxAge Deprecated - use oldest.age instead.
#' @param oldest.age Maximum age limit for Bacon runs, default at 1,000,000 cal BP. To set plot limits, use \code{age.max} instead.
#' @param cutoff Avoid plotting very low probabilities of date distributions (default \code{cutoff=0.001}).
#' @param rounded Rounding of calendar years. Defaults to 1 decimal. 
#' @param plot.pdf Produce a pdf file of the age-depth plot. Defaults to \code{plot.pdf=TRUE} after a Plum run.
#' @param dark Darkness of the greyscale age-depth model. The darkest grey value is \code{dark=1} by default.
#' Lower values will result in lighter grey but values >1 are not allowed.
#' @param date.res Date distributions are plotted using \code{date.res=100} segments by default.
#' @param age.res Resolution or amount of greyscale pixels to cover the age scale of the age-model plot. Default \code{yr.res=200}.
#' @param close.connections Internal option to close connections after a run. Default \code{close.connections=TRUE}.
#' @param save.info By default, a variable called `info' with relevant information about the run (e.g., core name, priors, settings, ages, output) is saved into the working directory. Note that this will overwrite any existing variable with the same name - as an alternative, one could run, e.g., \code{myvar <- Bacon()}, followed by supplying the variable \code{myvar} in any subsequent commands.
#' @param older.than an option to enable dates at the limit of C-14 dating. If there are older.than dates (works only for non-210Pb data), they tell us that the core should be older than a certain age at that depth. For example, if the 7th and 8th dates in the core's 'otherdates' .csv file are older-than dates, use as \code{older.than=c(7,8)}. The MCMC run could be problematic if the older-than ages do not fit with the other information.
#' @param younger.than an option to provide younger-than ages, for example a historical pollen marker. If there are younger-than dates (works only for non-210Pb data), they tell us that the core should be younger than a certain age at that depth. For example, if the 7th and 8th dates in the core's 'otherdates' .csv file are younger.than dates, use as \code{younger.than=c(7,8)}. The MCMC run could be problematic if the younger.than ages do not fit with the other information.
#' @param save.elbowages If you want to have a file with the MCMC-derived ages for all the age-depth model's elbows, set \code{save.elbowages=TRUE} and a file with the ages will be saved in the core's folder, ending in "_elbowages.txt".
#' @param verbose Provide feedback on what is happening (default \code{verbose=TRUE}).
#' @param ... options for the age-depth graph. See the agedepth and calib.plot functions. 
#' @author Maarten Blaauw, J. Andres Christen, Marco A. Aquino L.
#' @return An age-depth model graph, its age estimates, a summary, and the info variable which contains all relevant information.
#' @examples
#' \donttest{
#'   Plum(ask=FALSE, ssize=1000, coredir=tempfile(), date.sample=2018.5, ra.case=0, n.supp=3)
#' }
#' @references
#' Aquino-Lopez, M.A., Blaauw, M., Christen, J.A., Sanderson, N., 2018. Bayesian analysis of 210Pb dating. Journal of Agricultural, Biological, and Environmental Statistics 23, 317-333.
#'
#' Blaauw, M. and Christen, J.A., 2011. Flexible paleoclimate age-depth models using an autoregressive gamma process. Bayesian Analysis 6, 457-474.
#'
#' Christen, J.A., Perez E.S., 2010. A new robust statistical model for radiocarbon data. Radiocarbon 51, 1047-1059.
#'
#' Hogg et al., 2020. SHCal20 Southern Hemisphere calibration, 0-55,000 years cal BP. Radiocarbon 62, 759-778.
#'
#' Hua et al., 2022. Atmospheric radiocarbon for the period 1950-2019. Radiocarbon 64(4), 723-745, \doi{10.1017/RDC.2021.95}
#'
#' Hughen et al., 2020. Marine20-the marine radiocarbon age calibration curve (0-55,000 cal BP). Radiocarbon 62, 779-820.
#'
#' Jones, V.J., Stevenson, A.C., Battarbee, R.W., 1989. Acidification of lakes in Galloway, south west Scotland - a diatom and pollen study of the post-glacial history of the Round Loch of Glenhead.
#' Journal of Ecology 77, 1-23.
#'
#' Reimer et al., 2020. The IntCal20 Northern Hemisphere radiocarbon age calibration curve (0–55 cal kBP). Radiocarbon 62, 725-757.
#'
#' @export
Plum <- function(core="HP1C", thick=1, otherdates=NA, coredir="", phi.shape=2, phi.mean=50, s.shape=5, s.mean=10, Al=0.1, date.sample=c(), n.supp=c(), remove.tail=TRUE, ra.case=c(), Bqkg=TRUE, seed=NA, prob=0.95, d.min=0, d.max=NA, d.by=1, depths.file=FALSE, depths=c(), depth.unit="cm", age.unit="yr", unit=depth.unit, acc.shape=1.5, acc.mean=10, mem.strength=10, mem.mean=0.5, boundary=NA, hiatus.depths=NA, hiatus.max=10000, add=c(), after=.0001/thick, cc=1, cc1="IntCal20", cc2="Marine20", cc3="SHCal20", cc4="ConstCal", cc.dir="", postbomb=0, F14C=c(), pMC=c(), delta.R=0, delta.STD=0, t.a=3, t.b=4, normal=FALSE, suggest=TRUE, reswarn=c(10,200), remember=TRUE, ask=TRUE, run=TRUE, defaults="defaultPlum_settings.txt", sep=",", dec=".", runname="", slump=c(), BCAD=FALSE, ssize=4000, th0=c(), burnin=min(500, ssize), MinAge=c(), youngest.age=c(), MaxAge=c(), oldest.age=c(), cutoff=.001, rounded=1, plot.pdf=TRUE, dark=1, date.res=100, age.res=200, close.connections=TRUE, save.info=TRUE, older.than=c(), younger.than=c(), save.elbowages=FALSE, verbose=TRUE, ...) {
  # Check coredir and if required, copy example file in core directory
  coredir <- assign_coredir(coredir, core, ask, isPlum=TRUE)
  if(core == "HP1C" || core == "LL14") {
    dir.create(paste0(coredir, core, "/"), showWarnings = FALSE, recursive = TRUE)
    fileCopy <- system.file(paste0("extdata/Cores/", core), package="rplum")
    file.copy(fileCopy, coredir, recursive = TRUE, overwrite=FALSE)
  }

  # set the calibration curve
  if(cc.dir=="")
    cc.dir <- system.file("extdata", package="rintcal")
  cc.dir <- validateDirectoryName(cc.dir)

  # default_settings.txt is located within system.file
  defaults <- system.file("extdata", defaults, package=packageName())
  # read in the data, adapt settings from defaults if needed
  tmp <- read.dets.plum(core=core, coredir=coredir, n.supp=n.supp, date.sample=date.sample, sep=sep, dec=dec, cc=cc, Bqkg=Bqkg, ra.case=ra.case, suggest=suggest)

  dets <- tmp[[1]]
  supportedData <- tmp[[2]]
  ra.case <- tmp[[3]]
  date.sample <- tmp[[4]]
  detsOrig <- tmp[[5]]
  n.supp <- tmp[[6]]
  Bqkg <- tmp[[7]]

#  if(is.na(date.sample)) {
#    ans <- readline("Please provide a date (in AD) for the sample: ")

#    if( grepl("^[0-9.][0-9]*[.]?[0-9]*[0-9.]$",ans) == FALSE ) 
#      if( grepl("^[0-9]+$", ans) == FALSE )
#        stop("The date must be a numeric value", call.=FALSE)
#    date.sample = as.numeric(ans)
#  }
  theta0 <- 1950 - date.sample # date.sample is on AD scale, but internally we work with cal BP
  if(length(youngest.age) == 0)
    youngest.age <- min(theta0, 1950 - as.integer(format(Sys.time(), "%Y"))) else # was max
      if(youngest.age > theta0)
        message("setting theta0 to be the youngest of youngest.age and date.sample")
  theta0 <- min(youngest.age, theta0)

  if(length(oldest.age) == 0)
    oldest.age <- max(1e6, round(dets[,2] + (5*dets[,3])))

  # use the correct radiometric units # Mar 2022: check if correct. Should also allow for other conversions.
  if(Bqkg)
    dets[,6] <- dets[,6]*10.0 else {
      dets[,6] <- dets[,6]*500./3.
      Al <- Al*500./3.
    }
	
  detsBacon <- c()
  if(!is.na(otherdates)) { # core also has cal BP or C-14 dates
    csv.file <- paste0(coredir, core, "/", otherdates)
	detsBacon <- read.dets(core, coredir, otherdates, sep=sep, dec=dec, cc=cc)

    if(length(F14C) > 0) { # April 2025
  	  if(min(detsBacon[F14C,2]) < 0 || max(detsBacon[F14C,2]) > 3) 
        stop("The F14C values cannot be negative and are unlikely to be >3. Are you sure these values are in F14C?")		
      asC14 <- rice::F14CtoC14(detsBacon[F14C,2], detsBacon[F14C,3])
  	  detsBacon[F14C,2] <- round(asC14[,1],0)
  	  detsBacon[F14C,3] <- round(asC14[,2],0)
  	  rbacon:::fastwrite(as.data.frame(detsBacon), csv.file, sep=sep, dec=dec, row.names=FALSE, quote=FALSE) 
  	  message(paste("replaced F14C values with C14 ages in", csv.file))  
    }
    if(length(pMC) > 0) { # April 2025
  	  if(min(detsBacon[pMC,2]) < 0 || max(detsBacon[pMC,2]) > 300) 
        stop("The pMC values cannot be negative and are unlikely to be >300. Are you sure these values are in pMC?")		
      asC14 <- rice::pMCtoC14(detsBacon[pMC,2], detsBacon[pMC,3])
  	  detsBacon[pMC,2] <- round(asC14[,1])
  	  detsBacon[pMC,3] <- round(asC14[,2])
  	  rbacon:::fastwrite(as.data.frame(detsBacon), csv.file, sep=sep, dec=dec, row.names=FALSE, quote=FALSE) 
  	  message(paste("replaced pMC values with C14 ages in", csv.file))  
    }

    detsPlum <- dets
    # merge radiocarbon and 210Pb dates into the same variable dets
    dets <- merge_dets(dets, detsBacon, delta.R, delta.STD, t.a, t.b, cc)
  } else {
    detsPlum <- dets
    for(i in (ncol(dets)+1):9) {
      if(i == 5) {
        col <- array(delta.R, dim=c(nrow(dets),1))
      } else if(i == 6) {
          col <- array(delta.STD, dim=c(nrow(dets),1))
      } else if(i == 7) {
          col <- array(t.a, dim=c(nrow(dets),1))
      } else if(i == 8) {
          col <- array(t.b, dim=c(nrow(dets),1))
      } else if(i==9) {
          col <- array(5, dim=c(nrow(dets),1))
      }
      dets <- cbind(dets, col)
    }
    colnames(dets) <- c("labID", "X210Pb.Bq.kg.", "sd.210Pb.", "depth.cm.", "thickness.cm.", "density.g.cm.3.",  "t.a", "t.b", "cc")  
  }
  
  if(!is.na(otherdates)) { # if other, non-210Pb dates are to be included
    # give feedback about calibration curves used
    if(ncol(detsBacon) > 4 && length(cc) > 0) {
      cc.csv <- unique(detsBacon[,5])
      if(verbose) {
        if(length(cc.csv) == 1) {
          if(cc.csv != cc)
            message(" Using calibration curve specified within the .csv file,", cc[cc.csv], "\n")
        } else
          if(min(cc.csv) == 0)
            message(" Using a mix of cal BP and calibrated C-14 dates\n")
          else
            message(" Using several C-14 calibration curves\n")
	  }
    }

    if(suggest) { # adapt prior for mean accumulation rate?
      sugg <- sapply(c(1,2,5), function(x) x*10^(-1:2)) # some suggested "round" values
      if(nrow(detsBacon) < 2) # if there's only 1 entry, e.g. a Cs peak...
        sugg <- acc.mean else { # then don't use this but simply use the prior
          ballpacc <- lm(detsBacon[,2]*1.1 ~ detsBacon[,4])$coefficients[2] # very rough acc.rate estimates, uncalibrated dates
          ballpacc <- abs(sugg - ballpacc) # absolute differences between given acc.mean and suggested ones
          ballpacc <- ballpacc[ballpacc > 0] # do not suggest 0
          sugg <- sugg[order(ballpacc)[1]] # suggest rounded acc.rate with lowest absolute difference
        }
      if(!sugg %in% acc.mean) {
        ans <- readline(message(" Ballpark estimates suggest changing the prior for acc.mean to ", sugg, " ", age.unit, "/", depth.unit, ". OK? (y/N) "))
        if(tolower(substr(ans,1,1)) == "y")
          acc.mean <- sugg else
            message(" No problem, using the provided prior\n")
      }
    }
  }

  if(!is.na(boundary[1]))
    boundary <- sort(unique(boundary))
  if(!is.na(hiatus.depths[1])) {
    hiatus.depths <- sort(unique(hiatus.depths))
    if(length(acc.mean) == 1)
      acc.mean <- rep(acc.mean, length(hiatus.depths)+1)
  }

  #n.supp <<- n.supp # tmp?
  #set$n.supp <- n.supp

  if(nrow(detsPlum) <= 7) {
    message("Warning! Very few data points. Setting the bottom one to be background - scary.")

    # replace with request for chosen number of background
    bg <- min(1, n.supp)
  } else 
      bg <- check.equi(detsPlum, FALSE)

  if(suggest) {
    # check if the depths in the det file are bottom depths, and not, say, midpoints
    # it does this by calculating the top depths and ensuring they are not above d.min
    if(min(detsPlum[,4] - detsPlum[,5]) < d.min) # then we have a problem
      stop(paste0("The depths in ", core, ".csv should be the bottom depths of the measured slices. Not the midpoints! Or adapt d.min?\n"), call.=TRUE)

    # check if accrates might need adaptation
    # only check if n.supp not provided

    drange <- detsPlum[1:(nrow(detsPlum)-bg),4] # range of depths with unsupported Pb
    accrate <- (max(drange) - min(drange)) # assuming 100 years as fixed 210Pb limit, ugly
    agelim <- (1/0.03114) * log(phi.mean/Al) # Eq. 7 from Aquino et al. 2018
    message("\nPrior for acc.mean set at ", acc.mean, " ", age.unit, "/", depth.unit,
      ", ballpark estimate ", abs(round(agelim/accrate, 1))," ", age.unit, "/", depth.unit,
      if(abs(agelim/accrate) < 1)
        ", which seems quite fast. Adapt acc.mean?",
      if(abs(agelim/accrate) > 20)
        ", which seems quite slow. Adapt acc.mean?", "\n")
  }

  info <- .plum.settings(core=core, coredir=coredir, dets=dets, detsPlum=detsPlum, detsBacon=detsBacon, thick=thick, remember=remember, d.min=d.min, d.max=d.max, d.by=d.by, depths.file=depths.file, slump=slump, acc.mean=acc.mean, acc.shape=acc.shape, mem.mean=mem.mean, mem.strength=mem.strength, boundary=boundary, hiatus.depths=hiatus.depths, hiatus.max=hiatus.max, BCAD=BCAD, cc=cc, postbomb=postbomb, cc1=cc1, cc2=cc2, cc3=cc3, cc4=cc4, depth.unit=depth.unit, normal=normal, t.a=t.a, t.b=t.b, delta.R=delta.R, delta.STD=delta.STD, prob=prob, defaults=defaults, runname=runname, ssize=ssize, dark=dark, youngest.age=youngest.age, oldest.age=oldest.age, cutoff=cutoff, age.res=age.res, after=after, age.unit=age.unit, supportedData=supportedData, date.sample=date.sample, Al=Al, phi.shape=phi.shape, phi.mean=phi.mean, s.shape=s.shape, s.mean=s.mean, ra.case=ra.case, Bqkg=Bqkg, n.supp=n.supp)

  info$n.supp <- n.supp
  
  # optionally, make the info variable available in the working environment (default, but will overwrite any existing variable with the name 'info')
  info$save.info <- save.info
  if(save.info)
    assign_to_global("info", info)

  info$isplum <- TRUE # identify this as a plum run
  info$th0 <- th0
  info$cc.dir <- cc.dir
  info$seed <- seed

  if(!is.na(otherdates)) {
    info$hasBaconData <- TRUE
    info$detsBacon <- detsBacon
    info$detsPlum <- detsPlum
    info$detsOrig <- detsOrig
  } else {
    info$hasBaconData <- FALSE
    info$detsPlum <- detsPlum
    info$detsOrig <- detsOrig
  }

  ### check for initial mistakes
  if(length(MinAge) > 0)
    warning("please do not use MinAge, instead use the option youngest.age", call.=FALSE)
  if(length(MaxAge) > 0)
    warning("please do not use MaxAge, instead use the option oldest.age", call.=FALSE)
  if(any(info$acc.shape == info$acc.mean))
    stop("acc.shape cannot be equal to acc.mean", call.=FALSE)
  if(round(info$t.b - info$t.a) != 1)
    stop("t.b - t.a should always be 1, check the manual", call.=FALSE)
  if(min(acc.shape) < 1)
    message("\nWarning, using values <1 for acc.shape might cause unexpected results\n")

  ## calibrate dates
  if(info$hasBaconData && info$cc > 0) # confirm we are using radiocarbon dates
    if(info$postbomb == 0 && ((ncol(info$detsBacon) == 4 && min(info$detsBacon[,2]) < 0) ||
      ncol(info$detsBacon)>4 && max(info$detsBacon[,5]) > 0 && min(info$detsBacon[info$detsBacon[,5] > 0,2]) < 0))
        stop("you have negative C14 ages so should select a postbomb curve", call.=FALSE)
  if(info$hasBaconData)  # only calibrate radiocarbon dates
    info$calib <- bacon.calib(info$detsBacon, info, date.res, cc.dir=cc.dir)

  ### find some relevant values
  info$rng <- c()
  if(info$hasBaconData)
    for(i in 1:length(info$calib$probs)) {
      tmp <- info$calib$probs[[i]]
      info$rng <- range(info$rng, tmp[which(tmp[,2]>cutoff),1])
    }

  if(length(th0) == 0) # provide two ball-park/initial age estimates
    info$th0 <- c(min(youngest.age,theta0, na.rm=TRUE)+.1, min(youngest.age, theta0, na.rm=TRUE)+.2)
  info$th0[info$th0 < info$youngest.age] <- info$youngest.age # otherwise twalk will not start

  ### assign depths
  if(length(depths) == 0)
    depths <- seq(info$d.min, info$d.max, by=d.by) # was info$d.by
  if(depths.file) {
    dfile <- paste0(info$coredir, info$core, "/", info$core, "_depths.txt")
    if(!file.exists(dfile))
      stop("I cannot find the file ", paste0(info$coredir, info$core, "/", info$core, "_depths.txt"), call.=FALSE)
    depths <- read.table(dfile, header=FALSE)[,1]
    if(!is.numeric(depths[1]))
      stop("File should contain numbers only, no headers", call.=FALSE)
  }
  info$depths <- depths
  if(min(depths) < info$d.min)
    info$d.min <- min(depths)
  if(max(depths) > info$d.max)
    info$d.max <- max(depths)+thick # added +thick as per David Hatton's bug report

  info$elbows <- seq(info$d.min, info$d.max, by=thick)
  info$K <- length(info$elbows)
  info$cK <- info$d.min+(info$thick*info$K) # the maximum depth to be used by the bacon model

  # stop and warn if hiatus.depths conflict with other parameters
  if(!is.na(info$hiatus.depths[1]) || !is.na(info$boundary[1])) {
    ifelse(is.na(info$boundary[1]), hd <- info$hiatus.depths, hd <- info$boundary)
    if(min(hd) < info$d.min) # hiatus above core top
      stop("cannot have hiatus above the core's top depth. Adapt hiatus.depths or d.min.", call.=FALSE)
    if(max(hd)+info$thick > info$d.max)
      stop("the age-depth model should have at least one section below the one containing the deepest hiatus. Adapt thick or d.max?", call.=FALSE)
    if(length(hd) > 1) { # then check for how far separated hiatuses are
      above <- c()
      for(i in hd)
        above <- c(above, max(which(info$elbows <= i))) # find the section top of each hiatus
        if(any(diff(above) < 2)) # stop if fewer than 2 section elbows separating hiatuses
          stop("we need at least 2 section elbows between hiatuses. Choose fewer hiatuses, different depths, more sections (decrease thick) or a different d.min.\n ", call.=FALSE)
    }
  }

  ans <- "n"
  if(suggest)
    if(length(reswarn) == 2)
      if(info$K < min(reswarn)) {
        sugg <- pretty(thick*(info$K/min(reswarn)), 10)
        sugg <- min(sugg[sugg>0])
        ans <- readline(message(" Warning, the current value for thick, ", thick, ", will result in very few age-model sections (", info$K, ", not very flexible). Suggested maximum value for thick: ", sugg, " OK? (y/n) "))
      } else
        if(info$K > max(reswarn)) {
          sugg <- max(pretty(thick*(info$K/max(reswarn))))
          ans <- readline(message(" Warning, the current value for thick, ", thick, ", will result in very many age-model sections (", info$K, ", possibly hard to run). Suggested minimum value for thick: ", sugg, " OK? (y/n) "))
        }
    if(tolower(substr(ans, 1, 1)) == "y") {
      message(" OK, setting thick to ", sugg, "\n")
      thick <- sugg
      info$thick = thick
      info$elbows <- seq(floor(info$d.min), ceiling(info$d.max), by=thick)
      if(length(info$slump) > 0) # why here, and not a few lines later?
        info$elbows <- seq(floor(info$d.min), toslump(ceiling(info$d.max), info$slump), by=thick)
      info$K <- length(info$elbows)
      info$cK <- info$d.min+(info$thick*info$K) # the maximum depth to be used by the bacon model
    }

  ### prepare for any slumps
  if(length(slump) > 0) {
    if(length(slump) %% 2 == 1)
      stop("slumps need both upper and lower depths. Please check the manual", call.=FALSE)
    slump <- matrix(sort(slump), ncol=2, byrow=TRUE)
    info$slump <- slump

    slumpdmax <- toslump(ceiling(info$d.max), slump)
    info$elbows <- seq(floor(info$d.min), slumpdmax, by=thick)
    info$K <- length(info$elbows)
    info$cK <- info$d.min+(info$thick*info$K) # the maximum depth to be used by the bacon model

    info$slumpfree <- toslump(depths, slump)
    info$slumphiatus <- toslump(info$hiatus.depths, slump) # check
    if(!is.na(info$boundary[1])) {
      info$slumpboundary <- toslump(info$boundary, slump) # check
      info$slumphiatus <- info$slumpboundary
    }
    slumpdets <- info$dets
    slumpdets[,4] <- toslump(slumpdets[,4], slump, remove=FALSE) # dates within slumps are not removed
    info$slumpdets <- slumpdets[!is.na(slumpdets[,4]),]
  }

  ### produce files
  info$prefix <- paste0(coredir, core, "/", core, runname, "_", info$K)
  info$coredir <- coredir
  info$plum.name <- paste0(core, runname, "_", info$K,".plum")
  info$plum.file <- paste0(coredir, core, "/",core, runname, "_", info$K,".plum")
  info$bacon.file <- paste0(info$prefix, ".bacon")
  if(!file.exists(outfile <- paste0(info$prefix, ".out")))
    file.create(outfile)
  #create output file for data values
  if(!file.exists(plumOutFile <- paste0(info$prefix, "_plum.out") ))
    file.create(plumOutFile)
  if(!file.exists(baconTmpOutFile <- paste0(info$prefix, "_bacon.out") ))
    file.create(baconTmpOutFile)

  ### store values (again) for future manipulations
  if(BCAD)
    info$BCAD <- TRUE
  if(!is.na(boundary[1])) {
    if(length(slump) > 0)
      boundary <- info$slumpboundary
    info$hiatus.depths <- boundary
    if(length(add) == 0)
      add <- info$acc.mean # then add a short (max)hiatus, large enough not to crash Bacon but not affect the chronology much. Needs more work
    info$hiatus.max <- add
  }
  
  if(save.info)
    assign_to_global("info", info)

  prepare <- function() {
    oldpar <- par(mar=c(3,3,1,1), mgp=c(1.5,.7,.0), bty="l", xaxs="i")
    on.exit(par(oldpar))         	
    ### plot initial data and priors
    if( !info$hasBaconData){
      pn <- c(1:4, rep(5,4))
      if(!is.na(info$hiatus.depths[1])) # was ...hiatus.depths)[1])
        if(is.na(info$boundary[1]))
          pn <- c(1:5, rep(6,5))
      layout(matrix(pn, nrow=2, byrow=TRUE), heights=c(.3,.7))
    } else {
      pn <- c(1,2,3,4,5,5,6,6)
      if(!is.na(info$hiatus.depths[1])) # was ...hiatus.depths)[1])
        if(is.na(info$boundary[1]))
          pn <- c(1,1,2,2,3,3,4,4,5,5, rep(6,5), rep(7,5))
      layout(matrix(pn, nrow=2, byrow=TRUE), heights=c(.3,.7))
    }

    PlotAccPrior(info$acc.shape, info$acc.mean, info, depth.unit=depth.unit, age.unit=age.unit)
    PlotMemPrior(info$mem.strength, info$mem.mean, thick, info)

    if(!is.na(info$hiatus.depths)[1])
      if(is.na(info$boundary)[1])
        PlotHiatusPrior(info$hiatus.max, info$hiatus.depths, info)

    PlotPhiPrior(, , info)
    PlotSuppPrior(info)

    if(info$hasBaconData)
      rbacon::calib.plot(info, dets=info$detsBacon, BCAD=BCAD, new.plot=TRUE, plot.dists=TRUE, height=1)
    draw.pbmeasured(info)
    legend("top", core, bty="n", cex=1.5)
  }

  cook <- function() {
    plum.its(ssize, info) # new June 2021
    txt <- paste0(info$prefix, ".bacon")
    ssize <- as.integer(ssize)
    bacon(txt, outfile, ssize, cc.dir)
    rbacon::scissors(burnin, info, save.info=save.info)

    agedepth(info, BCAD=BCAD, depths.file=depths.file, depths=depths, verbose=TRUE, age.unit=age.unit, depth.unit=depth.unit, remove.tail=remove.tail, save.info=save.info, ...)
    #Plum.agedepth(info, BCAD=BCAD, depths.file=depths.file, depths=depths, verbose=TRUE, age.unit=age.unit, depth.unit=depth.unit, ...) # tmp May 21

    if(plot.pdf)
      if(interactive())
        if(length(dev.list()) > 0)
          dev.copy2pdf(file=paste0(info$prefix, ".pdf")) else {
            pdf(file=paste0(info$prefix, ".pdf"))
            rbacon::agedepth(info, BCAD=BCAD, depths.file=depths.file, depths=depths, verbose=FALSE, age.unit=age.unit, depth.unit=depth.unit, rounded=rounded, ...)
            # Plum.agedepth(info, BCAD=BCAD, depths.file=depths.file, depths=depths, verbose=FALSE, age.unit=age.unit, depth.unit=depth.unit, rounded=rounded, ...) # tmp May 21
            dev.off()
          }
  }

  # no need to cut off tails if we have dates below the above-background 210Pb data
  if(info$hasBaconData)
    if(max(info$detsBacon[,4]) > max(info$dets[which(info$dets$cc==5),2])) 
      remove.tail <- FALSE 

  ### run plum if initial graphs seem OK; run automatically, not at all, or only plot the age-depth model
  write.plum.file(info, younger.than=younger.than, older.than=older.than, save.info=save.info)
  if(!run)
    prepare() else
      if(!ask)
        cook() else {
          prepare()
          ans <- readline(message("Run ", core, " with ", info$K, " sections? (Y/n) "))
          ans <- tolower(substr(ans,1,1))[1]
          if(ans=="y" || ans=="")
            cook() else
              message("  OK. Please adapt settings.\n")
        }

  #if(close.connections)
  #  close(outfile)  
  
  if(save.elbowages) {
    saved <- sapply(info$elbows, Bacon.Age.d)
    write.table(saved, paste0(info$prefix, "_",info$d.min, "_", thick, "_elbowages.txt"), row.names=FALSE, col.names=FALSE)
  }

  invisible(info) # MB April 2024
}
