library("iSEEde")
library("airway")
library("DESeq2")
library("iSEE")

airway <- readRDS("dds.rds")

iSEE(airway, initial = list(
  DETable(ContrastName="dex: trt vs untrt", HiddenColumns = c("baseMean", 
    "lfcSE", "stat"), PanelWidth = 4L),
  VolcanoPlot(ContrastName="dex: trt vs untrt", PanelWidth = 4L),
  MAPlot(ContrastName="dex: trt vs untrt", PanelWidth = 4L)
))

# shiny::runApp(app)
