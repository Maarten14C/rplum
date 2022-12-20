# to enable running rbacon functions from within rplum

# these functions run c++ code and are as yet unexported (not sure how to export them)
events <- utils::getFromNamespace("events", "rbacon")
# bacon <- utils::getFromNamespace("bacon", "rbacon") # since already called like this within the Plum function

# these functions are exported by rbacon
agedepth <- rbacon::agedepth
proxy.ghost <- rbacon::proxy.ghost

# these are R functions in rbacon which are not exported (because they are internal)
draw.pbmeasured <- utils::getFromNamespace("draw.pbmeasured", "rbacon")
bacon.calib <- utils::getFromNamespace("bacon.calib", "rbacon")
calib.plot <- utils::getFromNamespace("calib.plot", "rbacon")
read.dets <- utils::getFromNamespace("read.dets", "rbacon")
toslump <- utils::getFromNamespace("toslump", "rbacon")
PlotAccPrior <- utils::getFromNamespace("PlotAccPrior", "rbacon")
PlotMemPrior <- utils::getFromNamespace("PlotMemPrior", "rbacon")
PlotHiatusPrior <- utils::getFromNamespace("PlotHiatusPrior", "rbacon")
PlotPhiPrior <- utils::getFromNamespace("PlotPhiPrior", "rbacon")
PlotSuppPrior <- utils::getFromNamespace("PlotSuppPrior", "rbacon")
PlotAccPost <- utils::getFromNamespace("PlotAccPost", "rbacon")
PlotMemPost <- utils::getFromNamespace("PlotMemPost", "rbacon")
PlotHiatusPost <- utils::getFromNamespace("PlotHiatusPost", "rbacon")
PlotPhiPost <- utils::getFromNamespace("PlotPhiPost", "rbacon")
PlotSuppPost <- utils::getFromNamespace("PlotSuppPost", "rbacon")
PlotLogPost <- utils::getFromNamespace("PlotLogPost", "rbacon")
Plum.AnaOut <- utils::getFromNamespace("Plum.AnaOut", "rbacon")
agedepth.ghost <- utils::getFromNamespace("agedepth.ghost", "rbacon")
Bacon.rng <- utils::getFromNamespace("Bacon.rng", "rbacon")
hiatus.slopes <- utils::getFromNamespace("hiatus.slopes", "rbacon")
assign_to_global <- utils::getFromNamespace("assign_to_global", "rbacon")
assign_coredir <- utils::getFromNamespace("assign_coredir", "rbacon")
validateDirectoryName <- utils::getFromNamespace("validateDirectoryName", "rbacon")
