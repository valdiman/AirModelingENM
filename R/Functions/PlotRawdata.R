# Plot raw data

plot_obs <- function(format_obs_data, pcb_name, save_path = "Output/Plots/RawData") {
  plot.spme <- ggplot(format_obs_data, aes(x = time, y = mf_Obs)) +
    geom_point(shape = 21, size = 2.5, color = "red") +
    theme_bw() +
    labs(x = expression(bold("Time (day)")),
         y = bquote(bold(.(pcb_name) ~ "SPME (ng/cm)"))) +
    theme(axis.text.y = element_text(face = "bold", size = 10),
          axis.title.y = element_text(face = "bold", size = 10, vjust = 0.1),
          axis.text.x = element_text(face = "bold", size = 10),
          axis.title.x = element_text(face = "bold", size = 10),
          aspect.ratio = 3/2)
  
  plot.puf <- ggplot(format_obs_data, aes(x = time, y = mpuf_Obs)) +
    geom_point(shape = 21, size = 2.5, color = "red") +
    theme_bw() +
    labs(x = expression(bold("Time (day)")),
         y = bquote(bold(.(pcb_name) ~ "PUF (ng/PUF)"))) +
    theme(axis.text.y = element_text(face = "bold", size = 10),
          axis.title.y = element_text(face = "bold", size = 10, vjust = 0.1),
          axis.text.x = element_text(face = "bold", size = 10),
          axis.title.x = element_text(face = "bold", size = 10),
          aspect.ratio = 3/2)
  
  # Combine both plots
  combined_plot <- gridExtra::grid.arrange(plot.spme, plot.puf, ncol = 2)
  
  # Define the file path for saving
  file_name <- paste0(save_path, "/", pcb_name, "_plot.png")
  
  # Save the combined plot directly
  ggsave(file_name, plot = combined_plot, width = 10, height = 5, dpi = 300)
  
  return(combined_plot)
}