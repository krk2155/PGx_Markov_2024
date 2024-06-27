# Name:
# Date            Accomplishment
# 6/27/2024       - Create basic Markov model   

#Mac
setwd("/Users/kun-wookim/Library/CloudStorage/OneDrive-VUMC/Research_discrete-event-simulation/r/PGx_Markov_2024")

# Load packages -----------------------------------------------------------
load.lib<-c("flexsurv", "msm", "dplyr",
            "ggplot2", "reshape2")
install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)

