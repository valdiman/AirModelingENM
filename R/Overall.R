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

# Model uptake of ENM -----------------------------------------------------

# Step 1: Prepare the data
cair_times <- cair$time
cair_vals <- cair$PCB8

enm <- enm[-5, ] # remove sample 5
enm_times <- enm$time
enm_vals <- enm$PCB8

# Known constant
denm <- 11080 # [g/m3]

# Step 2: Define the ODE system
model <- function(t, state, parameters, data) {
  with(as.list(c(state, parameters, data)), {
    Cair_t <- approx(cair_times, cair_vals, xout = t, rule = 2)$y  # interpolated Cair
    dXenm <- (ku / denm) * Cair_t - ke * Xenm
    list(c(dXenm))
  })
}

# Step 3: Residual function for model fitting
fit_model <- function(pars, times, Xenm_obs, cair_times, cair_vals, denm) {
  state <- c(Xenm = Xenm_obs[1])
  out <- ode(
    y = state,
    times = times,
    func = model,
    parms = pars,
    method = "lsoda",
    data = list(cair_times = cair_times, cair_vals = cair_vals, denm = denm)
  )
  sim_Xenm <- out[, "Xenm"]
  return(sim_Xenm - Xenm_obs)  # residuals to minimize
}

# Step 4: Fit the model using nls.lm (NOT nlsLM)
start_pars <- c(ku = 0.1, ke = 0.01)

fit <- nls.lm(
  par = start_pars,
  fn = fit_model,
  lower = c(0, 0),  # optional: bounds to keep parameters positive
  control = nls.lm.control(maxiter = 200),
  times = enm_times,
  Xenm_obs = enm_vals,
  cair_times = cair_times,
  cair_vals = cair_vals,
  denm = denm
)

# Step 5: Results
summary(fit)

predict_Xenm <- function(pars, times, cair_times, cair_vals, denm) {
  ku <- pars["ku"]
  ke <- pars["ke"]
  
  # Interpolate Cair over time
  Cair_interp <- approxfun(cair_times, cair_vals, rule = 2)
  
  # ODE system
  ode_func <- function(t, state, parameters) {
    Xenm <- state[1]
    Cair <- Cair_interp(t)
    dXenm <- ku / denm * Cair - ke * Xenm
    return(list(c(dXenm)))
  }
  
  # Initial condition (first value of Xenm)
  state0 <- c(Xenm = 0)  # or use enm_vals[1] if you prefer
  
  # Solve ODE
  ode_out <- ode(y = state0, times = times, func = ode_func, parms = NULL)
  
  # Extract predicted Xenm
  Xenm_pred <- ode_out[, "Xenm"]
  return(Xenm_pred)
}

pars_fit <- fit$par

time_smooth <- seq(min(enm_times), max(enm_times), by = 0.1)

Xenm_pred_smooth <- predict_Xenm(pars_fit, time_smooth, cair_times, cair_vals, denm)

smooth_df <- data.frame(
  time = time_smooth,
  Predicted = Xenm_pred_smooth
)

# Original observed data
obs_df <- data.frame(
  time = enm_times,
  Observed = enm_vals
)

ggplot() +
  geom_point(data = obs_df, aes(x = time, y = Observed),
             color = "blue", size = 3) +
  geom_line(data = smooth_df, aes(x = time, y = Predicted),
            color = "red", size = 1) +
  labs(x = "Time (hours)", y = "Xenm (ng/g)") +
  theme_minimal()


cair_avg <- mean(tail(cair$PCB8, 4))

predict_Xenm_constCair <- function(pars, times, cair_avg, denm) {
  ku <- pars["ku"]
  ke <- pars["ke"]
  
  # ODE system with constant Cair
  ode_func <- function(t, state, parameters) {
    Xenm <- state[1]
    dXenm <- ku / denm * cair_avg - ke * Xenm
    return(list(c(dXenm)))
  }
  
  # Initial condition
  state0 <- c(Xenm = 0)  # or enm_vals[1] if preferred
  
  # Solve ODE
  ode_out <- ode(y = state0, times = times, func = ode_func, parms = NULL)
  
  Xenm_pred <- ode_out[, "Xenm"]
  return(Xenm_pred)
}

# Predict with constant Cair
Xenm_pred_const <- predict_Xenm_constCair(pars_fit, time_smooth, cair_avg, denm)

# Plot
ggplot() +
  geom_point(data = obs_df, aes(x = time, y = Observed), color = "blue", size = 3) +
  geom_line(aes(x = time_smooth, y = Xenm_pred_const), color = "darkgreen", size = 1) +
  labs(
    title = "Observed vs Predicted Xenm (Constant Cair)",
    x = "Time (hours)",
    y = "Xenm (ng/g)"
  ) +
  theme_minimal()


