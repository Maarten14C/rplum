## Test environments
* local Fedora install, R 4.0.2
* local Mac OS, R 4.0.2
* rhub win-builder (devel and release)

## R CMD check results

* This is a resubmission. I have removed the vignette for now, to avoid the notes.

0 errors | 0 warnings | 1 note

* 'et' and 'al' are not misspelled words.

* on some rhub installs, a note is produced regarding a sub-directory extdata of >1MB. This is required to: 
  1) enable examples to be performed quickly by loading existing MCMC runs, and 
  2) to include recently released calibration curves for radiocarbon dates, which are more detailed and thus larger than previous ones

* some examples take >5s to run as these are based on MCMC runs. Examples have been set such as to make these runs as short as possible while still resulting in reasonable outputs.
  


