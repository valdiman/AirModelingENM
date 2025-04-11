# Predict Xenm

# Define ODE with constant Cair
predict_Xenm_constCair <- function(pars, times, cair_avg, denm) {
  ku <- pars["ku"]
  ke <- pars["ke"]
  
  ode_func <- function(t, state, parameters) {
    Xenm <- state[1]
    dXenm <- ku / denm * cair_avg - ke * Xenm
    return(list(c(dXenm)))
  }
  
  state0 <- c(Xenm = 0)
  
  ode_out <- ode(y = state0, times = times, func = ode_func, parms = NULL)
  Xenm_pred <- ode_out[, "Xenm"]
  return(Xenm_pred)
}

