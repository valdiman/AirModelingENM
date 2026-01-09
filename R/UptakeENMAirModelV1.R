#
# Packages and libraries --------------------------------------------------
# Install packages
install.packages("dplyr")
install.packages("reshape2")
install.packages("ggplot2")
install.packages("deSolve")
install.packages("tidyr")
install.packages("gridExtra")
install.packages("minpack.cl")

# Load libraries
{
  library(dplyr) # organize data
  library(reshape2) # organize data
  library(ggplot2) # plotting
  library(deSolve) # solving differential equations
  library(tidyr)
  library(gridExtra)
  library(rlang)
  library(minpack.lm)
}

# Read Raw Data -----------------------------------------------------------
{
  rawdata <- read.csv("Data/Data.csv", check.names = FALSE)
  sr <- read.csv("Data/WBSamplingRateStat.csv", check.names = FALSE) # [m/d]
  logKoa <- read.csv("Data/logKoa.csv")
}

# Organize data -----------------------------------------------------------
wb.data <- rawdata[rawdata$sample == 'WB', ]
pan.data <- rawdata[rawdata$sample == 'PAN', ]
wb.sr <- data.frame(congener = sr$congener, ko = sr$ko)

# Calculate air concentration from WB -------------------------------------
# Calculate logKwb [m3/m3]
# Regression created with data from Tromp et al 2019 (Table 2, Wristband)
# & Frederiksen et al 2022 (Table 3)
logKwb <- data.frame(
  congener = logKoa$congener,
  logKwb = 0.6156 * logKoa$logKoa + 2.161) # R2 = 0.96

# remove metadata from wb.data
wb.data.2 <- wb.data[, 6:177]

# Use effective volume
# Add time
wb.data.2 <- cbind(wb.data$time, wb.data.2)
# Add WB volume
wb.data.2 <- cbind(wb.data$vol, wb.data.2)
# Add WB area
wb.data.2 <- cbind(wb.data$area, wb.data.2)
# Change column names
names(wb.data.2)[names(wb.data.2) == 'wb.data$time'] <- 'time.hr'
names(wb.data.2)[names(wb.data.2) == 'wb.data$vol'] <- 'vol.WB'
names(wb.data.2)[names(wb.data.2) == 'wb.data$area'] <- 'area.WB'

n_samples <- nrow(wb.data.2)
n_pcbs <- nrow(logKwb)

vol_matrix  <- matrix(rep(wb.data.2$vol.WB, times = n_pcbs),
                      nrow = n_samples, ncol = n_pcbs)
area_matrix <- matrix(rep(wb.data.2$area.WB, times = n_pcbs),
                      nrow = n_samples, ncol = n_pcbs)
time_matrix <- matrix(rep(wb.data.2$time.hr, times = n_pcbs),
                      nrow = n_samples, ncol = n_pcbs)
logK_matrix <- matrix(rep(logKwb$logKwb, each = n_samples),
                      nrow = n_samples, ncol = n_pcbs)
ko_matrix  <- matrix(rep(wb.sr$ko, each = n_samples),
                     nrow = n_samples, ncol = n_pcbs)

# Calculate veff for samples as a 5 x 172 matrix [m3]
veff <- 10^logK_matrix * vol_matrix *
  (1 - exp(-ko_matrix * area_matrix / vol_matrix / 10^logK_matrix * time_matrix / 24))

# Compute concentration [ng/m3]
conc.wb <- wb.data.2[, 4:175] / veff
conc.wb$time <- wb.data$time
tPCB.conc.wb <- as.data.frame(rowSums(conc.wb[, 1:172], na.rm = TRUE))

# Select PCBi -------------------------------------------------------------
pcb.ind <- "PCB16"

# Remove first observation and use constant Cair
conc.wb.i2 <- conc.wb[, c("time", pcb.ind)][-1, ]
Cair_const <- mean(conc.wb.i2[[pcb.ind]])

# PAN data
pan.i <- pan.data[, c("sample", "time", "vol", "area", pcb.ind)]
dpan <- 0.06 / pan.i$vol[1]  # PAN density
# Include time 0
pan_times <- c(0, pan.i$time)
pan_vals  <- c(0, pan.i[[pcb.ind]])

# Solve uptake equation ---------------------------------------------------
# Analytical solution function
pan_model <- function(time, ku, ke) {
  (ku * Cair_const) / (ke * dpan) * (1 - exp(-ke * time))
}

# Residuals for fitting

residuals_fun <- function(pars) {
  ku <- pars["ku"]
  ke <- pars["ke"]
  sim <- pan_model(pan_times, ku, ke)
  sim - pan_vals
}

# Starting parameters
start_pars <- c(ku = 80000, ke = 0.01)

# Fit using nonlinear least squares (Levenberg-Marquardt)
fit <- nls.lm(par = start_pars, fn = residuals_fun)

# Get optimized parameters
fitted_ku <- fit$par["ku"]
fitted_ke <- fit$par["ke"]

# Predicted PAN and goodness-of-fit
sim <- pan_model(pan_times, fitted_ku, fitted_ke)
res <- sim - pan_vals
RSS  <- sum(res^2)
TSS  <- sum((pan_vals - mean(pan_vals))^2)
R2   <- 1 - RSS / TSS
RMSE <- sqrt(mean(res^2))

cat("\nModel performance:\n")
cat("R2   =", R2, "\n")
cat("RMSE =", RMSE, "\n")

# Sampler parameters ------------------------------------------------------
# Kpan calculations
ku <- fit$par["ku"]
ke <- fit$par["ke"]
logKpan <- log10(ku / ke * 1000^2 / dpan) # [L/kg]
# Sampling rate
Rs <- ku * pan.i$vol[1] * 24 # [m3/d]
# Time to reach 90% equilibrium
t90.h <- log(10) / ke # [h]
t90.d <- log(10) / ke  / 24 # [d]

# Create a summary data frame
model_summary <- data.frame(
  pcb.ind = pcb.ind,
  ku = ku,
  ke = ke,
  logKpan = logKpan,
  Rs = Rs,
  t90 = t90.h,
  t90 = t90.d,
  RMSE = RMSE,
  R2 = R2,
  row.names = NULL)

# Add units to column names
colnames(model_summary) <- c(
  "congener",
  "ku (1/h)",
  "ke (1/h)",
  "logKpan (L/kg)",
  "Rs (m3/d)",
  "t90 (h)",
  "t90 (d)",
  "RMSE",
  "R2")

# Print the summary for verification
print(model_summary)

# Save the summary to a CSV file
summary_filename <- paste0("Output/Data/Model/Summary/Summary_", pcb.ind, ".csv")
write.csv(model_summary, file = summary_filename, row.names = FALSE)

# Plot observations & modeling predictions --------------------------------
# Finer time grid for plotting
time_grid <- seq(min(pan_times), max(pan_times), by = 0.1)

# ODE for plotting
ode_func_plot <- function(t, state, pars) {
  ku <- pars["ku"]
  ke <- pars["ke"]
  X <- state[1]
  dX <- ku / dpan * Cair_const - ke * X
  return(list(c(dX)))
}

state0 <- c(X = 0)
out_plot <- ode(y = state0, times = time_grid, func = ode_func_plot,
                parms = fit$par, method = "rk4")
smooth_df <- data.frame(time = out_plot[, "time"], Predicted = out_plot[, "X"])

# Observed data
obs_df <- data.frame(
  time = pan.i$time,
  Observed = pan.i[[pcb.ind]])

# Plot
plot.uptake <- ggplot() +
  geom_point(data = obs_df, aes(x = time, y = Observed), 
             shape = 21, color = "black", size = 2.5) +
  geom_line(data = smooth_df, aes(x = time, y = Predicted),
            color = "black", linewidth = 0.5) +
  theme_bw() +
  labs(x = expression(bold("Time (hours)")),
       y = bquote(bold("ng" ~ .(pcb.ind) ~ "/g PAN"))) +
  theme(aspect.ratio = 3/5,
        axis.text = element_text(face = "bold", size = 10),
        axis.title = element_text(face = "bold", size = 11),
        panel.grid = element_blank())

plot.uptake

# Save plot in folder
plot_filename <- paste0("Output/Plots/Model/Plot_", pcb.ind, ".png")
ggsave(plot_filename, plot = plot.uptake, width = 5,
       height = 3, dpi = 500)

