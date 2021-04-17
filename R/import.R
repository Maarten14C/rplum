# to enable running rbacon functions from within rplum
events <- utils::getFromNamespace("events", "rbacon")
bacon <- utils::getFromNamespace("bacon", "rbacon")
#agedepth <- rbacon::agedepth # currently using tailed-made version for rplum
#calib.plumbacon.plot <- utils::getFromNamespace("calib.plumbacon.plot", "rbacon") # currently using rplum version instead of the one from rbacon
proxy.ghost <- rbacon::proxy.ghost # is this a good idea?
calib.plot <- utils::getFromNamespace("calib.plot", "rbacon")
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
#Plum.AnaOut <- utils::getFromNamespace("Plum.AnaOut", "rbacon")
agedepth.ghost <- utils::getFromNamespace("agedepth.ghost", "rbacon")
Bacon.rng <- utils::getFromNamespace("Bacon.rng", "rbacon")
hiatus.slopes <- utils::getFromNamespace("hiatus.slopes", "rbacon")
assign_to_global <- utils::getFromNamespace("assign_to_global", "rbacon")
