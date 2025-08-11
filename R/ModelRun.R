# Run functions

source("R/Functions/ReadData.R")
source("R/Functions/ResidualFunction.R")
source("R/Functions/PredictXenm.R")
source("R/Functions/UptakeModelFunction.R")
source("R/Functions/PlotModel.R")

pcb.ind <- "PCB8"
data <- read_data(pcb.ind)
model_out <- run_uptake_model(data$enm, data$cair, pcb.ind)
plot_uptake_model(model_out, pcb.ind)

