# Fitting and summary function

run_uptake_model <- function(enm, cair, pcb.ind) {
  enm <- enm[-5, ]  # Remove outlier
  
  enm_times <- enm$time
  enm_vals <- enm[[pcb.ind]]
  cair_avg <- mean(tail(cair[[pcb.ind]], 4))
  
  A <- 9.5 * 9.5 / 100^2 # [m2]
  th <- 0.2 / 100        # [m]
  V <- A * th            # [m3]
  m <- 0.2               # [g]
  denm <- m / V          # [g/m3]
  
  start_pars <- c(ku = 1000, ke = 0.1)
  
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
  
  Xenm_fitted <- predict_Xenm_constCair(fit$par, enm_times, cair_avg, denm)
  
  RMSE <- sqrt(mean((enm_vals - Xenm_fitted)^2))
  RSS <- sum((enm_vals - Xenm_fitted)^2)
  TSS <- sum((enm_vals - mean(enm_vals))^2)
  R2 <- 1 - RSS / TSS
  
  ku <- fit$par["ku"]
  ke <- fit$par["ke"]
  logKenm <- log10(ku / ke * 1000^2 / denm)
  
  return(list(
    fit = fit,
    RMSE = RMSE,
    R2 = R2,
    ku = ku,
    ke = ke,
    logKenm = logKenm,
    fitted_vals = Xenm_fitted,
    cair_avg = cair_avg,
    denm = denm,
    enm_times = enm_times,
    enm_vals = enm_vals
  ))
}
