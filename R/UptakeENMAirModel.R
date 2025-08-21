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
rawdata <- read.csv("Data/Data.csv", check.names = FALSE)
sr <- read.csv("Data/WBSamplingRateStat.csv", check.names = FALSE)
logKoa <- read.csv("Data/logKoa.csv")

# Organize data
wb.data <- rawdata[rawdata$sample == 'WB', ]
# Remove 4 hrs wb - outlier
wb.data <- wb.data[-1,]
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

veff_wb <- matrix(NA, nrow = 5, ncol = 172)

for (i in 1:length(Vwb)) {
  # Loop over 173 congeners (columns)
  for (j in 1:length(Kwb)) {
    veff_wb[i, j] <- Kwb[j] * Vwb[i] * (1 - exp(-ko[j] * Awb[i] / Vwb[i] / Kwb[j] * t[i] / 24))
  }
}

# Compute concentration
conc.wb <- wb.data[, 6:177] / veff_wb
# Average 5 sampling concentrations
conc.wb.av <- as.data.frame(t(colMeans(conc.wb)))
# Add sample and time
conc.wb.av$sample <- c('conc.ave')

# Extract relevant columns from pan.data
pcb.ind <- "PCB11"
pan.i <- pan.data[, c("sample", "time", pcb.ind)]

# Extract relevant columns from conc.wb
cair.i <- conc.wb.av[, c(pcb.ind)]

# Pan density calculation
A <- pan.data$area # [m2]
V <- pan.data$vol[1] # [m3]
m <- 0.06 # [g]
denm <- m / V  # [g/m3]

# Model uptake of ENM -----------------------------------------------------
# Extract times and values
pan_times <- pan.data$time
pan_vals <- pan.data[[pcb.ind]]

# Define ODE with constant Cair
predict_Xpan <- function(pars, times, cair.i, denm) {
  ku <- pars["ku"]
  ke <- pars["ke"]
  
  ode_func <- function(t, state, parameters) {
    Xpan <- state[1]
    dXpan <- ku / denm * cair.i - ke * Xpan
    return(list(c(dXpan)))
  }
  
  state0 <- c(Xpan = 0)  # Xpan @ time 0 = 0
  
  ode_out <- ode(y = state0, times = times, func = ode_func, parms = NULL)
  Xpan_pred <- ode_out[, "Xpan"]
  return(Xpan_pred)
}

# Residuals function for fitting
fit_model_constCair <- function(pars, times, Xpan_obs, cair.i, denm) {
  sim_Xpan <- predict_Xpan(pars, times, cair.i, denm)
  return(sim_Xpan - Xpan_obs)
}

# Starting parameters
start_pars <- c(ku = 100000, ke = 0.1) # [1/h]

# Fit using constant Cair
fit <- nls.lm(
  par = start_pars,
  fn = fit_model_constCair,
  lower = c(0, 0),
  control = nls.lm.control(maxiter = 200),
  times = pan_times,
  Xpan_obs = pan_vals,
  cair.i = cair.i,
  denm = denm
)

summary(fit)

# Predicted values at original time points
Xpan_fitted <- predict_Xpan(fit$par, pan_times, cair.i, denm)

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
summary_filename <- paste0("Output/Data/Model/Summary/Summary_", pcb.ind, ".csv")
write.csv(model_summary, file = summary_filename, row.names = FALSE)

# Predict for plotting
pars_fit <- fit$par
time_smooth <- seq(0, max(pan_times), by = 0.1) # Start @ time 0
Xpan_pred <- predict_Xpan(pars_fit, time_smooth, cair.i, denm)

obs_df <- data.frame(time = pan_times, Observed = pan_vals)
smooth_df <- data.frame(time = time_smooth, Predicted = Xpan_pred)

# Save predicted data
pred_filename <- paste0("Output/Data/Model/Prediction/Pred_", pcb.ind, ".csv")
write.csv(smooth_df, file = pred_filename, row.names = FALSE)

# Plot
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

# Display plot
plot.uptake

# Save plot in folder
plot_filename <- paste0("Output/Plots/Model/Plot_", pcb.ind, ".png")
ggsave(plot_filename, plot = plot.uptake, width = 5,
       height = 5, dpi = 500)

