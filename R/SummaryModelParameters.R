# Compile model data generated

library(dplyr)
library(readr)
library(purrr)


# Read generated data


# list all files that start with "Summary_PCB"
files <- list.files(path = "Output/Data/Model/Summary", 
                    pattern = "^Summary_PCB.*\\.csv$", 
                    full.names = TRUE)

# read them all into a list of data frames
df_list <- map(files, read_csv)

# Combine data frame
all_data <- bind_rows(df_list, .id = "source_file")

logKoa <- read.csv("Data/logKoa.csv")

# Select only the needed columns from all_data
all_data_sel <- all_data %>%
  select(congener,
         ku = `ku (1/h)`,
         ke = `ke (1/h)`,
         logKpan = `logKpan (L/kg)`,
         Rs = `Rs (m3/d)`,
         t90 = `t90 (h)`)

# Join with logKoa by congener
combined_df <- all_data_sel %>%
  left_join(logKoa, by = "congener")

plot(x = combined_df$logKoa, combined_df$t90)

