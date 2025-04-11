# Plotting predictions and observations

plot_uptake_model <- function(model_out, pcb.ind) {
  # Extract the required components from model_out
  enm_times <- model_out$enm_times
  enm_vals <- model_out$enm_vals
  fitted_vals <- model_out$fitted_vals
  
  # Create a smooth sequence of time points
  time_smooth <- seq(min(enm_times), max(enm_times), by = 0.1)
  
  # Assuming a simple linear interpolation for predicted values over the smoothed time points
  smooth_fitted_vals <- approx(enm_times, fitted_vals, xout = time_smooth)$y
  
  # Create a data frame for observed and smoothed predicted values
  obs_df <- data.frame(time = enm_times, Observed = enm_vals)
  smooth_df <- data.frame(time = time_smooth, Predicted = smooth_fitted_vals)
  
  # Plot
  plot_uptake <- ggplot() +
    geom_point(data = obs_df, aes(x = time, y = Observed), shape = 21,
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
  
  print(plot_uptake)
  
  # Optionally save the plot if needed
  # ggsave("Output/Plots/Uptake/PCB.png", plot = plot_uptake, width = 5,
  #        height = 5, dpi = 500)
}

# Now you can call the function like this:
plot_uptake_model(model_out, pcb.ind)
