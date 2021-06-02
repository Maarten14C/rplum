# rplum 0.2.2
* removed closeAllConnections() option as requested by CRAN
* if there is only one date in otherdates (e.g., a historical Cs peak), acc.rate isn't set to overly low values any more
* added a function to estimate how many iterations will be run and returned

# rplum 0.2.1
* further separation of rbacon and rplum
* added LL14, a core with Pb-210, radon and C14 dates (core kindly provided by Dr Lysanna Anderson, USGS)
* many sundry bug fixes related mainly to plotting

# rplum 0.2.0
* now depends on the rbacon c++ code instead of carrying a duplicate of the code. Same for many of rbacon's internal functions for plotting etc.
* now depends on IntCal package for radiocarbon calibration curves
* now only calculates the age-depth model for Pb data that are above background level
* enhanced how to take into account different scenarios with radon.case and n.supp

# rplum 0.1.5.1
* removed vignettes as they threw errors which made the package unacceptable by CRAN

# rplum 0.1.5
* updated to the latest radiocarbon calibration curves IntCal20
* repaired a bug with plotting when using BCAD=TRUE
* reduced sizes of files of example cores
* removed inset option for now as it caused many post-plotting issues

# rplum 0.1.4
* clarified author names
* replaced cat() by message() or warning() where possible
* lots of minor bug fixes related to plotting and post-run analysis

# rplum 0.1.3
* added single quotes to package names as requested by CRAN maintainer
* further clarified copyright statements

# rplum 0.1.2
* clarified copyright statements

# rplum 0.1
* first version of the 'rplum' package, developed from 'pyplum' at https://github.com/maquinolopez
