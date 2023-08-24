# Simulation of experiment

# Libraries
library(tidyverse) 

# The goal of this script is to take the values outputted by the Baranyi-model earlier and see if they fit with the change in abundance observed both in NGS and in competition experiments

# Numerical Baranyi-Roberts model parameters

# # INPUT
# - A data frame containing the names of the peptides (peptide), their respective µmax (mumax) and their starting OD (y0)
# - A vector of time steps at which to compute the ODs
# - The constant across all strains: h0, K

# # OUTPUT
# A data frame containing the OD of each strain at each time point, according to the Baranyi-Roberts model (see Baranyi, 1994)

# Exemple Data
K = 2
h0 = 5

time = seq(0,75000, by = 250)
input_df <- data.frame(peptide = c("EP1", "EP2", "DP1", "DP2", "DP3"),
                       mumax = c(26,25,22,17,16)*1e-5,
                       y0 = rep(1e-3,5))

competition_simulation <- function(input_df, 
                                   time = seq(0,75000, by = 250),
                                   K = 2,
                                   h0 = 5){
  z0 = 1 # My understanding of the math is that this value doesn't matter since Kz is computed from it
  Kz = (exp(h0) - 1) * z0 # Compute Kz
  
  # Create the names we will use for OD and z (metabolic states in the DF)
  rownames_OD = paste0(input_df$peptide,"_OD")
  rownames_z = paste0(input_df$peptide,"_z")
  
  # Extract the vector of mumaxes
  mumax_vector = input_df$mumax
  
  # Pre-allocate storage_df
  num_times = length(time)
  num_peptides = nrow(input_df)
  storage_df = data.frame(matrix(nrow = num_times, ncol = num_peptides * 2 + 1))
  colnames(storage_df) = c(rownames_OD, rownames_z, "time")
  storage_df$time = time
  storage_df[1, rownames_OD] = input_df$y0
  storage_df[1, rownames_z] = z0
  
  # Loop through time points
  for(i in 2:num_times) {
    dt = time[i] - time[i-1]
    ut_dt = 1 - sum(storage_df[i-1, rownames_OD]) / K
    dz_vector = mumax_vector * storage_df[i-1, rownames_z]
    zt_vector = storage_df[i-1, rownames_z] + dz_vector * dt
    dy_vector = ut_dt * storage_df[i-1, rownames_OD] * mumax_vector * zt_vector / (Kz + zt_vector)
    yt_vector = storage_df[i-1, rownames_OD] + dy_vector * dt
    
    storage_df[i, rownames_OD] = yt_vector
    storage_df[i, rownames_z] = zt_vector
  }
  
  return(storage_df)
}
# 
# ptm <- Sys.time()
# x <- competition_simulation(input_df, time = seq(from = 1, to = 5000, by = 1))
# print(Sys.time()-ptm)

# This function just takes the DF produced by growth_simulation() and edit it to make it suitable for ggplots

pivot_simulation_df <- function(df){
  # Grab the names of the columns with _OD
  colnames <- colnames(df)
  rownames_OD <- colnames[grepl("_OD",colnames)]
  
  
  df = df %>% 
    select(all_of(c(rownames_OD, "time"))) %>%
    rename_with(~gsub("_OD", "", .x))
  
  graph_df = df %>%
    pivot_longer(-time, names_to = "peptide", values_to = "OD")
  
  return(graph_df)
}



