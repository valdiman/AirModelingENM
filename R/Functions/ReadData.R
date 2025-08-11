# Read Data

read_data <- function(pcb.ind) {
  # Read the data
  enm.data <- read.csv("Data/ENM.csv", check.names = FALSE)
  cair.data <- read.csv("Data/Ca.csv")
  
  # Extract relevant columns from enm.data
  enm <- enm.data[, c("sample", "time", pcb.ind)]
  
  # Extract relevant columns from enm.data
  cair <- cair.data[, c("sample", "time", pcb.ind)]
  
  return(list(enm = enm, cair = cair))
}

