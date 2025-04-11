# Residual function

fit_model_constCair <- function(pars, times, Xenm_obs, cair_avg, denm) {
  sim_Xenm <- predict_Xenm_constCair(pars, times, cair_avg, denm)
  return(sim_Xenm - Xenm_obs)
}

