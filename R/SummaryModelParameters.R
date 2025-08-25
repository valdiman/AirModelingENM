# Compile model data generated

# Packages and libraries --------------------------------------------------
# Install packages
install.packages("dplyr")
install.packages("ggplot2")
install.packages("readr")
install.packages("purrr")

# Load libraries
{
  library(dplyr) # organize data
  library(ggplot2) # plotting
  library(readr)
  library(purrr)
}

# Read generated data -----------------------------------------------------
# List all files
files <- list.files(path = "Output/Data/Model/Summary", 
                    pattern = "^Summary_PCB.*\\.csv$", 
                    full.names = TRUE)

# Read them all into a list of data frames
df_list <- map(files, read_csv)

# Combine data frame
all_data <- bind_rows(df_list, .id = "source_file")

# Read logKoa
logKoa <- read.csv("Data/logKoa.csv")

# PCB diffusivity in air
T.air <- 22 # [C]
P <- 1013 #[mbar]
D.water.air <- (10^(-3)*1013.25*((273.15+T.air)^1.75*((1/28.97) +
                                                        (1/18.0152))^(0.5))/P/(20.1^(1/3)
                                                                               + 9.5^(1/3))^2) # [cm2/s]
D.PCB.air <- D.water.air*(logKoa$MW/18.0152)^(-0.5) # [cm2/s]

D.PCB.air <- data.frame(
  congener = logKoa$congener,
  d.PCB.air = D.PCB.air
)

# Remove PCBs with R2 <= 0.85
all_data2 <- all_data[!is.infinite(all_data$`logKpan (L/kg)`), ]
all_data2 <- all_data2[all_data2$R2 >= 0.85, ]

# Select data -------------------------------------------------------------
all_data_sel <- all_data2 %>%
  select(congener,
         ku = `ku (1/h)`,
         ke = `ke (1/h)`,
         logKpan = `logKpan (L/kg)`,
         Rs = `Rs (m3/d)`,
         t90 = `t90 (h)`)

# Join with logKoa by congener
combined_df <- all_data_sel %>%
  left_join(logKoa, by = "congener") %>%
  left_join(D.PCB.air, by = "congener")

# Plots -------------------------------------------------------------------
# See plot
# (1) log Koa vs log ku
ggplot(combined_df, aes(x = logKoa, y = log10(ku))) +
  geom_point(shape = 21, color = "black", size = 2.5) +
  geom_text(aes(label = congener), vjust = -0.5, size = 3) +
  theme_bw() +
  labs(x = expression(bold("log Koa")),
       y = bquote(bold("log ku (1/h)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10, vjust = 0.1),
        axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10),
        aspect.ratio = 1,
        panel.grid = element_blank())

# Fit linear model
lm_fit <- lm(log10(ku) ~ logKoa, data = combined_df)

# Extract coefficients
coef_fit <- coef(lm_fit)
intercept <- coef_fit[1]
slope <- coef_fit[2]

# Extract summary to get R² and p-value
lm_sum <- summary(lm_fit)
R2 <- lm_sum$r.squared
p_slope <- lm_sum$coefficients["logKoa", "Pr(>|t|)"]  # p-value of slope

eq_text <- sprintf("y = %.2f x + %.2f", slope, intercept)
R2_text <- sprintf("R² = %.2f", R2)
p_text <- sprintf("p = %.1e", p_slope)

plot.ku <- ggplot(combined_df, aes(x = logKoa, y = log10(ku))) +
  geom_point(shape = 21, color = "black", size = 2.5) +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  annotate("text", x = min(combined_df$logKoa), 
           y = max(log10(combined_df$ku)) + 0.4, 
           label = eq_text, hjust = 0, vjust = 1, size = 4) +
  annotate("text", x = min(combined_df$logKoa), 
           y = max(log10(combined_df$ku)) + 0.3, 
           label = R2_text, hjust = 0, vjust = 1, size = 4) +
  annotate("text", x = min(combined_df$logKoa), 
           y = max(log10(combined_df$ku)) + 0.2, 
           label = p_text, hjust = 0, vjust = 1, size = 4) +
  theme_bw() +
  labs(x = expression(bold("log Koa")),
       y = bquote(bold("log ku (1/h)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10, vjust = 0.1),
        axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10),
        aspect.ratio = 1,
        panel.grid = element_blank())

# See plot
plot.ku
# Save plot in folder
ggsave("Output/Plots/Regressions/logKoa_logku.png", plot = plot.ku,
       width = 5, height = 5, dpi = 500)

# (2) d.PCB.air vs log ku
ggplot(combined_df, aes(x = d.PCB.air, y = log10(ku))) +
  geom_point(shape = 21, color = "black", size = 2.5) +
  geom_text(aes(label = congener), vjust = -0.5, size = 3) +
  theme_bw() +
  labs(x = expression(bold("Air Diffusivity (cm2/s)")),
       y = bquote(bold("log ku (1/h)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10, vjust = 0.1),
        axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10),
        aspect.ratio = 1,
        panel.grid = element_blank())

# Fit linear model
lm_fit <- lm(log10(ku) ~ d.PCB.air, data = combined_df)

# Extract coefficients
coef_fit <- coef(lm_fit)
intercept <- coef_fit[1]
slope <- coef_fit[2]

# Extract summary to get R² and p-value
lm_sum <- summary(lm_fit)
R2 <- lm_sum$r.squared
p_slope <- lm_sum$coefficients["d.PCB.air", "Pr(>|t|)"]  # p-value of slope

eq_text <- sprintf("y = %.2f x + %.2f", slope, intercept)
R2_text <- sprintf("R² = %.2f", R2)
p_text <- ifelse(p_slope < 0.001, "p < 0.001",
                 sprintf("p = %.3f", p_slope))

plot.ku2 <- ggplot(combined_df, aes(x = d.PCB.air, y = log10(ku))) +
  geom_point(shape = 21, color = "black", size = 2.5) +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  annotate("text", x = min(combined_df$d.PCB.air), 
           y = max(log10(combined_df$ku)) + 0.4, 
           label = eq_text, hjust = 0, vjust = 1, size = 4) +
  annotate("text", x = min(combined_df$d.PCB.air), 
           y = max(log10(combined_df$ku)) + 0.3, 
           label = R2_text, hjust = 0, vjust = 1, size = 4) +
  annotate("text", x = min(combined_df$d.PCB.air), 
           y = max(log10(combined_df$ku)) + 0.2, 
           label = p_text, hjust = 0, vjust = 1, size = 4) +
  theme_bw() +
  labs(x = expression(bold("log Koa")),
       y = bquote(bold("log ku (1/h)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10, vjust = 0.1),
        axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10),
        aspect.ratio = 1,
        panel.grid = element_blank())

# See plot
plot.ku2
# Save plot in folder
ggsave("Output/Plots/Regressions/DPCBAir_logku.png", plot = plot.ku2,
       width = 5, height = 5, dpi = 500)

# (3) log Koa vs log ke
ggplot(combined_df, aes(x = logKoa, y = log10(ke))) +
  geom_point(shape = 21, color = "black", size = 2.5) +
  geom_text(aes(label = congener), vjust = -0.5, size = 3) +
  theme_bw() +
  labs(x = expression(bold("log Koa")),
       y = bquote(bold("log ke (1/h)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10, vjust = 0.1),
        axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10),
        aspect.ratio = 1,
        panel.grid = element_blank())

# Fit linear model
lm_fit <- lm(log10(ke) ~ logKoa, data = combined_df)

# Extract coefficients
coef_fit <- coef(lm_fit)
intercept <- coef_fit[1]
slope <- coef_fit[2]

# Extract summary to get R² and p-value
lm_sum <- summary(lm_fit)
R2 <- lm_sum$r.squared
p_slope <- lm_sum$coefficients["logKoa", "Pr(>|t|)"]  # p-value of slope

eq_text <- sprintf("y = %.2f x + %.2f", slope, intercept)
R2_text <- sprintf("R² = %.2f", R2)
p_text <- ifelse(p_slope < 0.001, "p < 0.001",
                 sprintf("p = %.3f", p_slope))

plot.ke <- ggplot(combined_df, aes(x = logKoa, y = log10(ke))) +
  geom_point(shape = 21, color = "black", size = 2.5) +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  annotate("text", x = min(combined_df$logKoa), 
           y = max(log10(combined_df$ke)) + 0.5, 
           label = eq_text, hjust = 0, vjust = 1, size = 4) +
  annotate("text", x = min(combined_df$logKoa), 
           y = max(log10(combined_df$ke)) + 0.4, 
           label = R2_text, hjust = 0, vjust = 1, size = 4) +
  annotate("text", x = min(combined_df$logKoa), 
           y = max(log10(combined_df$ke)) + 0.3, 
           label = p_text, hjust = 0, vjust = 1, size = 4) +
  theme_bw() +
  labs(x = expression(bold("log Koa")),
       y = bquote(bold("log ke (1/h)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10, vjust = 0.1),
        axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10),
        aspect.ratio = 1,
        panel.grid = element_blank())

# See plot
plot.ke
# Save plot in folder
ggsave("Output/Plots/Regressions/logKoa_logke.png", plot = plot.ke,
       width = 5, height = 5, dpi = 500)

# (4) log Koa vs log Kpan
ggplot(combined_df, aes(x = logKoa, y = logKpan)) +
  geom_point(shape = 21, color = "black", size = 2.5) +
  geom_text(aes(label = congener), vjust = -0.5, size = 3) +
  theme_bw() +
  labs(x = expression(bold("log Koa")),
       y = bquote(bold("log Kpan"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10, vjust = 0.1),
        axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10),
        aspect.ratio = 1,
        panel.grid = element_blank())
  
# Fit linear model
lm_fit <- lm(logKpan ~ logKoa, data = combined_df)

# Extract coefficients
coef_fit <- coef(lm_fit)
intercept <- coef_fit[1]
slope <- coef_fit[2]

# Extract summary to get R² and p-value
lm_sum <- summary(lm_fit)
R2 <- lm_sum$r.squared
p_slope <- lm_sum$coefficients["logKoa", "Pr(>|t|)"]  # p-value of slope

eq_text <- sprintf("y = %.2f x + %.2f", slope, intercept)
R2_text <- sprintf("R² = %.2f", R2)
p_text <- ifelse(p_slope < 0.001, "p < 0.001",
                 sprintf("p = %.3f", p_slope))

plot.kpan <- ggplot(combined_df, aes(x = logKoa, y = logKpan)) +
  geom_point(shape = 21, color = "black", size = 2.5) +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  annotate("text", x = min(combined_df$logKoa), 
           y = max(combined_df$logKpan), 
           label = eq_text, hjust = 0, vjust = 1, size = 4) +
  annotate("text", x = min(combined_df$logKoa), 
           y = max(combined_df$logKpan) - 0.05, 
           label = R2_text, hjust = 0, vjust = 1, size = 4) +
  annotate("text", x = min(combined_df$logKoa), 
           y = max(combined_df$logKpan) - 0.1, 
           label = p_text, hjust = 0, vjust = 1, size = 4) +
  theme_bw() +
  labs(x = expression(bold("log Koa")),
       y = bquote(bold("log Kpan"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10, vjust = 0.1),
        axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10),
        aspect.ratio = 1,
        panel.grid = element_blank())

# See plot
plot.kpan
# Save plot in folder
ggsave("Output/Plots/Regressions/logKao_logKpan.png", plot = plot.kpan,
       width = 5, height = 5, dpi = 500)

# (4) log Koa vs Rs
ggplot(combined_df, aes(x = logKoa, y = Rs)) +
  geom_point(shape = 21, color = "black", size = 2.5) +
  geom_text(aes(label = congener), vjust = -0.5, size = 3) +
  theme_bw() +
  labs(x = expression(bold("log Koa")),
       y = bquote(bold("Rs (m3/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10, vjust = 0.1),
        axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10),
        aspect.ratio = 1,
        panel.grid = element_blank())

# (4) log Koa vs t90
ggplot(combined_df, aes(x = logKoa, y = t90/24)) +
  geom_point(shape = 21, color = "black", size = 2.5) +
  geom_text(aes(label = congener), vjust = -0.5, size = 3) +
  theme_bw() +
  labs(x = expression(bold("log Koa")),
       y = bquote(bold("t90 (d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10, vjust = 0.1),
        axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10),
        aspect.ratio = 1,
        panel.grid = element_blank())

ggplot(combined_df, aes(x = logKoa, y = t90/24)) +
  geom_point(shape = 21, color = "black", size = 2.5) +
  theme_bw() +
  labs(x = expression(bold("log Koa")),
       y = bquote(bold("t90 (d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10, vjust = 0.1),
        axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10),
        aspect.ratio = 1,
        panel.grid = element_blank())

