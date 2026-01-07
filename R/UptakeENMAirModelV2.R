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

# Organize data
wb.data <- rawdata[rawdata$sample == 'WB', ]
pan.data <- rawdata[rawdata$sample == 'PAN', ]
wb.sr <- data.frame(congener = sr$congener, ko = sr$ko)

# Calculate air concentration from WB
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

# Extract relevant columns from pan.data
pcb.ind <- "PCB52"
pan.i <- pan.data[, c("sample", "time", "vol", "area", pcb.ind)]
# Extract relevant columns from conc.wb
conc.wb.i <- conc.wb[, c("time", pcb.ind)]

# Pan density calculation
m <- 0.06 # [g]
pan.i$den <- m / pan.i$vol  # [g/m3]

# Measured data
pan_times <- pan.i$time # PAN measurement times
pan_vals  <- pan.i[[pcb.ind]] # PAN observed concentrations

cair_times <- conc.wb.i$time # measured Cair times
cair_vals  <- conc.wb.i[[pcb.ind]] # measured Cair values
dpan <- pan.i$den[1] # PAN density [g/m3]

# ODE derivative function
ode_func <- function(t, state, pars) {
  ku <- pars["ku"] # [1/h]
  ke <- pars["ke"] # [1/h])
  
  # Use exact measured Cair at nearest previous time
  idx <- max(which(cair_times <= t))
  Cair <- cair_vals[idx]
  
  X <- state[1]
  dX <- ku / dpan * Cair - ke * X
  return(list(c(dX)))
}

# Function to predict Xpan at PAN measurement times
predict_Xpan <- function(pars) {
  state0 <- c(X = 0)  # initial amount
  out <- ode(y = state0, times = pan_times, func = ode_func, parms = pars)
  Xpan <- out[, "X"]
  return(Xpan)
}

# Residual function for fitting
residuals_fun <- function(pars) {
  pars <- setNames(pars, names(start_pars))  # ensures names
  sim <- predict_Xpan(pars)
  return(sim - pan_vals)
}

# Starting parameters
start_pars <- c(ku = 1000000, ke = 0.1) # [1/h]

# Fit using Levenberg-Marquardt
fit <- nls.lm(par = start_pars, fn = residuals_fun)

# Fitted parameters
fit$par

# Predicted values using the fitted parameters
sim <- predict_Xpan(fit$par)
# Residuals
res <- sim - pan_vals
# Residual sum of squares
RSS <- sum(res^2)
# Total sum of squares
TSS <- sum((pan_vals - mean(pan_vals))^2)
# R-squared
R2 <- 1 - RSS / TSS
# RMSE (Root Mean Square Error)
RMSE <- sqrt(mean(res^2))
# MAE (Mean Absolute Error)
MAE <- mean(abs(res))

# Print metrics
cat("Fitted parameters:\n")
print(fit$par)
cat("\nGoodness-of-fit:\n")
cat("R2:", R2, "\nRMSE:", RMSE, "\nMAE:", MAE, "\n")

# Sampler parameters
# Kpan calculations
ku <- fit$par[1]
ke <- fit$par[2]
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
summary_filename <- paste0("Output/Data/Model/Summary2/Summary_", pcb.ind, ".csv")
write.csv(model_summary, file = summary_filename, row.names = FALSE)

# Finer time grid for plotting
time_grid <- seq(min(pan_times), max(pan_times), by = 1)  # hours

# Linear interpolation (with rule=2 to avoid NA)
cair_interp <- approxfun(cair_times, cair_vals, rule = 2)

ode_func_plot <- function(t, state, pars) {
  ku <- pars["ku"]
  ke <- pars["ke"]
  X <- state[1]
  Cair <- cair_interp(t)
  dX <- ku / dpan * Cair - ke * X
  return(list(c(dX)))
}

state0 <- c(X = 0)
out_plot <- ode(y = state0, times = time_grid, func = ode_func_plot, parms = fit$par)
smooth_df <- data.frame(time = out_plot[, "time"], Predicted = out_plot[, "X"])

obs_df <- data.frame(
  time = pan.i$time,
  Observed = pan.i[[pcb.ind]])

plot.uptake <- ggplot() +
  geom_point(data = obs_df, aes(x = time, y = Observed), 
             shape = 21, color = "black", size = 2.5) +
  geom_line(data = smooth_df, aes(x = time, y = Predicted), 
            color = "black", linewidth = 0.4) +
  theme_bw() +
  labs(x = expression(bold("Time (hours)")),
       y = bquote(bold("ng" ~ .(pcb.ind) ~ "/g PAN"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10, vjust = 0.1),
        axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10),
        aspect.ratio = 1,
        panel.grid = element_blank())

plot.uptake

# Save plot in folder
plot_filename <- paste0("Output/Plots/Model2/Plot_", pcb.ind, ".png")
ggsave(plot_filename, plot = plot.uptake, width = 5,
       height = 5, dpi = 500)



