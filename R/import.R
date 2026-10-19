# to enable running rbacon functions from within rplum

require(rbacon)

# these functions run c++ code and are unexported
events <- utils::getFromNamespace("events", "rbacon")
#bacon <- utils::getFromNamespace("bacon", "rbacon") # since already called like this within the Plum function

# these functions are loaded through import rbacon
#agedepth <- rbacon::agedepth
#proxy.ghost <- rbacon::proxy.ghost
#importFrom(rbacon, agedepth)
#importFrom(rbacon, proxy.ghost)

rbacon_internal <- function(fun)
  getFromNamespace(fun, "rbacon")
draw.pbmeasured <- rbacon_internal("draw.pbmeasured")
bacon.calib <- rbacon_internal("bacon.calib")
calib.plot <- rbacon_internal("calib.plot")
fastwrite <- rbacon_internal("fastwrite")
read.dets <- rbacon_internal("read.dets")
toslump <- rbacon_internal("toslump")
PlotAccPrior <- rbacon_internal("PlotAccPrior")
PlotMemPrior <- rbacon_internal("PlotMemPrior")
PlotHiatusPrior <- rbacon_internal("PlotHiatusPrior")
PlotPhiPrior <- rbacon_internal("PlotPhiPrior")
PlotSuppPrior <- rbacon_internal("PlotSuppPrior")
PlotAccPost <- rbacon_internal("PlotAccPost")
PlotMemPost <- rbacon_internal("PlotMemPost")
PlotHiatusPost <- rbacon_internal("PlotHiatusPost")
PlotPhiPost <- rbacon_internal("PlotPhiPost")
PlotSuppPost <- rbacon_internal("PlotSuppPost")
PlotLogPost <- rbacon_internal("PlotLogPost")
Plum.AnaOut <- rbacon_internal("Plum.AnaOut")
agedepth.ghost <- rbacon_internal("agedepth.ghost")
Bacon.rng <- rbacon_internal("Bacon.rng")
hiatus.slopes <- rbacon_internal("hiatus.slopes")
assign_to_global <- rbacon_internal("assign_to_global")
assign_coredir <- rbacon_internal("assign_coredir")
validateDirectoryName <- rbacon_internal("validateDirectoryName")
