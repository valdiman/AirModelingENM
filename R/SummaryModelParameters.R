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

T.air <- 22 # [C]
P <- 1013 #[mbar]
D.water.air <- (10^(-3)*1013.25*((273.15+T.air)^1.75*((1/28.97) +
                                                        (1/18.0152))^(0.5))/P/(20.1^(1/3)
                                                                               + 9.5^(1/3))^2) # [cm2/s]
# PCB diffusivity in air
D.PCB.air <- D.water.air*(logKoa$MW/18.0152)^(-0.5) # [cm2/s]

D.PCB.air <- data.frame(
  congener = logKoa$congener,
  d.PCB.air = D.PCB.air
)

# Select data -------------------------------------------------------------
all_data_sel <- all_data %>%
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
