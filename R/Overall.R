#
# Packages and libraries --------------------------------------------------
# Install packages
install.packages("dplyr")
install.packages("reshape2")
install.packages("ggplot2")
install.packages("deSolve")
install.packages("tidyr")
install.packages("gridExtra")

# Load libraries
{
  library(dplyr) # organize data
  library(reshape2) # organize data
  library(ggplot2) # plotting
  library(deSolve) # solving differential equations
  library(tidyr)
  library(gridExtra)
  library(rlang)
}

# Read Raw Data -----------------------------------------------------------
# Read the data
enm.rawdata <- read.csv("Data/ENM.csv", check.names = FALSE)
cair.rawdata <- read.csv("Data/Ca.csv")

# Extract relevant columns from enm.data
pcb.ind <- "PCB8"
enm <- enm.rawdata[, c("sample", "time", pcb.ind)]

# Extract relevant columns from enm.data
cair <- cair.rawdata[, c("sample", "time", pcb.ind)]


# Plot Raw Data -----------------------------------------------------------
plot.enm <- ggplot(enm, aes(x = time, y = !!sym(pcb.ind))) +
  geom_point(shape = 21, size = 2.5, color = "blue") +
  theme_bw() +
  labs(x = expression(bold("Time (h)")),
       y = bquote(bold(.(pcb.ind) ~ "PAN (ng/g)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10, vjust = 0.1),
        axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10),
        aspect.ratio = 3/2)

plot.cair <- ggplot(cair, aes(x = time, y = !!sym(pcb.ind))) +
  geom_point(shape = 21, size = 2.5, color = "red") +
  theme_bw() +
  labs(x = expression(bold("Time (day)")),
       y = bquote(bold(.(pcb.ind) ~ "Cair (ng/m3)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10, vjust = 0.1),
        axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10),
        aspect.ratio = 3/2)

# Combine both plots
gridExtra::grid.arrange(plot.enm, plot.cair, ncol = 2)



