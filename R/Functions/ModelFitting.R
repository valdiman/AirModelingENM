# Model fitting

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

# Fit summary
summary(fit)
