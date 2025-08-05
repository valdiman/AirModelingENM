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
# Read the data
rawdata <- read.csv("Data/Data.csv", check.names = FALSE)
sr <- read.csv("Data/WBSamplingRateStat.csv", check.names = FALSE)
logKoa <- read.csv("Data/logKoa.csv")

# Organize data
wb.data <- rawdata[rawdata$sample == 'WB', ]
pan.data <- rawdata[rawdata$sample == 'PAN', ]
wb.sr <- data.frame(congener = sr$congener, ko = sr$ko)

# Calculate air concentration from wb
# Calculate logKws
# Regression created with data from Tromp et al 2019 (Table 2, Wristband)
# & Frederiksen et al 2022 (Table 3)
logKwb <- data.frame(
  congener = logKoa$congener,
  logKwb = 0.6156 * logKoa$logKoa + 2.161) # R2 = 0.96

# Use effective volume. Adult WBs
Kwb <- 10^logKwb$logKwb
ko  <- wb.sr$ko
Vwb <- wb.data$vol
Awb <- wb.data$area
t <- wb.data$time

veff_wb <- matrix(NA, nrow = 6, ncol = 172)

for (i in 1:length(Vwb)) {
  # Loop over 173 congeners (columns)
  for (j in 1:length(Kwb)) {
    veff_wb[i, j] <- Kwb[j] * Vwb[i] * (1 - exp(-ko[j] * Awb[i] / Vwb[i] / Kwb[j] * t[i] / 24))
  }
}

# Compute concentration
conc.wb <- wb.data[, 5:176] / veff_wb

# Add congener names to the columns
colnames(conc.wb) <- logKwb$congener
# Add sample and time
conc.wb$sample <- wb.data$sample
conc.wb$time <- wb.data$time

# Extract relevant columns from pan.data
pcb.ind <- "PCB56"
pan.i <- pan.data[, c("sample", "time", pcb.ind)]

# Extract relevant columns from conc.wb
cair.i <- conc.wb[, c("sample", "time", pcb.ind)]

# Model uptake of ENM -----------------------------------------------------
# Extract times and values
pan_times <- pan.data$time
pan_vals <- pan.data[[pcb.ind]]

# Constant average Cair from last 4 samples
cair.i_avg <- mean(tail(cair.i[[pcb.ind]], 5))

# Pan density calculation
A <- 38 / 100^2 # [m2]
th <- 0.015 / 100 # [m]
V <- A * th # [m3]
m <- 0.2 # [g]
denm <- m / V  # [g/m3]

# Define ODE with constant Cair
predict_Xpan_constCair <- function(pars, times, cair.i_avg, denm) {
  ku <- pars["ku"]
  ke <- pars["ke"]
  
  ode_func <- function(t, state, parameters) {
    Xpan <- state[1]
    dXpan <- ku / denm * cair.i_avg - ke * Xpan
    return(list(c(dXpan)))
  }
  
  state0 <- c(Xpan = 0)  # Xpan @ time 0 = 0
  
  ode_out <- ode(y = state0, times = times, func = ode_func, parms = NULL)
  Xpan_pred <- ode_out[, "Xpan"]
  return(Xpan_pred)
}

# Residuals function for fitting
fit_model_constCair <- function(pars, times, Xpan_obs, cair.i_avg, denm) {
  sim_Xpan <- predict_Xpan_constCair(pars, times, cair.i_avg, denm)
  return(sim_Xpan - Xpan_obs)
}

# Starting parameters
start_pars <- c(ku = 1000, ke = 0.1) # [1/h]

# Fit using constant Cair
fit <- nls.lm(
  par = start_pars,
  fn = fit_model_constCair,
  lower = c(0, 0),
  control = nls.lm.control(maxiter = 200),
  times = pan_times,
  Xpan_obs = pan_vals,
  cair.i_avg = cair.i_avg,
  denm = denm
)

summary(fit)

# Predicted values at original time points
Xpan_fitted <- predict_Xpan_constCair(fit$par, pan_times, cair.i_avg, denm)

# Model evaluation metrics (performance metrics)
# Root mean squared error (RMSE)
RMSE <- sqrt(mean((pan_vals - Xpan_fitted)^2))
# Residual sum of squares (RSS)
RSS <- sum((pan_vals - Xpan_fitted)^2) # [ng/g]
# Total sum of squares (TSS)
TSS <- sum((pan_vals - mean(pan_vals))^2)
# R-squared
R2 <- 1 - RSS / TSS

# Sampler parameters
# Kenm calculations
ku <- fit$par[1]
ke <- fit$par[2]
logKpan <- log10(ku / ke * 1000^2 / denm) # [L/kg]
# Sampling rate
Rs <- ku * V *24 # [m3/d]
# Time to reach 90% equilibrium
t90 <- log(10) / ke #[h]

# Create a summary data frame
model_summary <- data.frame(
  pcb.ind = pcb.ind,
  ku = ku,
  ke = ke,
  logKpan = logKpan,
  Rs = Rs,
  t90 = t90,
  RMSE = RMSE,
  R2 = R2,
  row.names = NULL
)

# Add units to column names
colnames(model_summary) <- c(
  "congener",
  "ku (1/h)",
  "ke (1/h)",
  "logKpan (L/kg)",
  "Rs (m3/d)",
  "t90 (h)",
  "RMSE",
  "R2"
)

# Print the summary for verification
print(model_summary)

# Save the summary to a CSV file
summary_filename <- paste0("Output/Data/Model/Summary_", pcb.ind, ".csv")
write.csv(model_summary, file = summary_filename, row.names = FALSE)

# Predict for plotting
pars_fit <- fit$par
time_smooth <- seq(min(pan_times), max(pan_times), by = 0.1)
Xpan_pred_const <- predict_Xpan_constCair(pars_fit, time_smooth, cair.i_avg, denm)

# Plot
obs_df <- data.frame(time = pan_times, Observed = pan_vals)
smooth_df <- data.frame(time = time_smooth, Predicted = Xpan_pred_const)

plot.uptake <- ggplot() +
  geom_point(data = obs_df, aes(x = time, y = Observed), shape  = 21,
             color = "black", size = 2.5) +
  geom_line(data = smooth_df, aes(x = time, y = Predicted),
            color = "black", linewidth = 0.4) +
  theme_bw() +
  labs(x = expression(bold("Time (hours)")),
       y = bquote(bold("ng" ~.(pcb.ind) ~ "/g PAN"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10, vjust = 0.1),
        axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10),
        aspect.ratio = 1,
        panel.grid = element_blank())

plot.uptake

# Save plot in folder
plot_filename <- paste0("Output/Plots/Model/Plot_", pcb.ind, ".png")
ggsave(plot_filename, plot = plot.uptake, width = 5,
       height = 5, dpi = 500)

