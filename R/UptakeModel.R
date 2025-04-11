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
enm.rawdata <- read.csv("Data/ENM.csv", check.names = FALSE)
cair.rawdata <- read.csv("Data/Ca.csv", check.names = FALSE)

# Extract relevant columns from enm.data
pcb.ind <- "PCB86+97+109+119"
enm <- enm.rawdata[, c("sample", "time", pcb.ind)]

# Extract relevant columns from enm.data
cair <- cair.rawdata[, c("sample", "time", pcb.ind)]

# Model uptake of ENM -----------------------------------------------------
# Remove outlier
enm <- enm[-5, ]

# Extract times and values
enm_times <- enm$time
enm_vals <- enm[[pcb.ind]]

# Constant average Cair from last 4 samples
cair_avg <- mean(tail(cair[[pcb.ind]], 4))

# Pan density calculation
A <- 9.5 * 9.5 / 100^2 # [m2]
th <- 0.015 / 100 # [m]
V <- A * th # [m3]
m <- 0.2 # [g]
denm <- m / V  # [g/m3]

# Define ODE with constant Cair
predict_Xenm_constCair <- function(pars, times, cair_avg, denm) {
  ku <- pars["ku"]
  ke <- pars["ke"]
  
  ode_func <- function(t, state, parameters) {
    Xenm <- state[1]
    dXenm <- ku / denm * cair_avg - ke * Xenm
    return(list(c(dXenm)))
  }
  
  state0 <- c(Xenm = 0)  # Xenm @ time 0 = 0
  
  ode_out <- ode(y = state0, times = times, func = ode_func, parms = NULL)
  Xenm_pred <- ode_out[, "Xenm"]
  return(Xenm_pred)
}

# Residuals function for fitting
fit_model_constCair <- function(pars, times, Xenm_obs, cair_avg, denm) {
  sim_Xenm <- predict_Xenm_constCair(pars, times, cair_avg, denm)
  return(sim_Xenm - Xenm_obs)
}

# Starting parameters
start_pars <- c(ku = 1000, ke = 0.1)

# Fit using constant Cair
fit <- nls.lm(
  par = start_pars,
  fn = fit_model_constCair,
  lower = c(0, 0),
  control = nls.lm.control(maxiter = 200),
  times = enm_times,
  Xenm_obs = enm_vals,
  cair_avg = cair_avg,
  denm = denm
)

summary(fit)

# Predicted values at original time points
Xenm_fitted <- predict_Xenm_constCair(fit$par, enm_times, cair_avg, denm)

# Root mean squared error (RMSE)
RMSE <- sqrt(mean((enm_vals - Xenm_fitted)^2))

# Residual sum of squares (RSS)
RSS <- sum((enm_vals - Xenm_fitted)^2) # [ng/g]

# Total sum of squares (TSS)
TSS <- sum((enm_vals - mean(enm_vals))^2)

# R-squared
R2 <- 1 - RSS / TSS

# Kenm calculations
ku <- fit$par[1]
ke <- fit$par[2]
logKenm <- log10(ku / ke * 1000^2 / denm) # [L/kg]
# Time to reach 90% equilibrium
t90 <- log(10) / ke #[h]

# Create a summary data frame
model_summary <- data.frame(
  ku = ku,
  ke = ke,
  logKenm = logKenm,
  t90 = t90,
  RMSE = RMSE,
  R2 = R2
)

# Save the summary to a CSV file
summary_filename <- paste0("Output/Data/Model/Summary_", pcb.ind, ".csv")
write.csv(model_summary, file = summary_filename, row.names = FALSE)

# Print the summary for verification
print(model_summary)

# Predict for plotting
pars_fit <- fit$par
time_smooth <- seq(min(enm_times), max(enm_times), by = 0.1)
Xenm_pred_const <- predict_Xenm_constCair(pars_fit, time_smooth, cair_avg, denm)

# Plot
obs_df <- data.frame(time = enm_times, Observed = enm_vals)
smooth_df <- data.frame(time = time_smooth, Predicted = Xenm_pred_const)

plot.uptake <- ggplot() +
  geom_point(data = obs_df, aes(x = time, y = Observed), shape  = 21,
             color = "black", size = 2.5) +
  geom_line(data = smooth_df, aes(x = time, y = Predicted),
            color = "black", size = 0.4) +
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

