# Date          Task
# 7/2/2024        - remove PGx & No-PGx; only 1 trace remains


#Mac
setwd("/Users/kun-wookim/Library/CloudStorage/OneDrive-VUMC/Research_discrete-event-simulation/r/PGx_Markov_2024")

# Load packages -----------------------------------------------------------
load.lib<-c("flexsurv", "msm", "dplyr", "markovchain", "heemod", "dampack",
            "devtools", "ggplot2", "reshape2")
install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)

### Define functions
source("check_transition_probabilities.R")

  
## General setup

cycle_length <- 1 # cycle length equal one year (use 1/12 for monthly)

n_pop <- 1000
n_age_init <- 40 # age at baseline 
n_age_max <- 80 # maximum age of follow up 
n_cycles <- (n_age_max - n_age_init)/cycle_length # time horizon, number of cycles 
v_names_states <- c("A","D") # the 2 health states of the model: # Alive (A), Dead (D)
n_states <- length(v_names_states) # number of health states 
d_e <- 0.03 # annual discount rate for QALYs of 3% 
d_c <- 0.03 # annual discount rate for costs of 3% 
v_names_str <- c("No PGx") # store the strategy names

t <- 0

# Transition probabilities (annual), and hazard ratios (HRs)  -------------
shape_gompertz <- 0.101 # shape parameter from Gompertz model of secular death for 40-year-old woman
rate_gompertz <- 0.001 # rate parameter from Gompertz model of secular death for 40-year-old woman
r_a_pct <- 0.1 # condition indication percentage
r_a_dur <- 10 # condition indication time period
r_HC <- -log(1-r_a_pct)/r_a_dur # rate of developing condition
p_HC <- 1- exp(-r_HC) # probability of developing condition

# State rewards -----------------------------------------------------------
## Costs -------------------------------------------------------------------

c_H <-0 # annual cost of being Healthy 
c_C <- 0 # annual cost of having a condition 
c_D <- 0 # annual cost of being dead


## Utilities ---------------------------------------------------------------
u_H <- 1 # annual utility of being Healthy 
u_C <- 1 # annual utility of having a condition
u_D <- 0 # annual utility of being dead

# Transition Probabilities2 -----------------------------------------------
## Mortality rates
## Probabilities of dying (Cycle-specific) ----------------------------------
p_background_death_gompertz <- function(t, rate_gompertz, shape_gompertz) {
  1-exp(-rate_gompertz*exp(shape_gompertz*t))
}


# Create an empty dataframe for the lifetable
lifetable <- data.frame(
  Age = n_age_init:(n_age_init + 41),  # Let's consider up to age 80
  Population = numeric(42),
  Prob_Death_Annual = numeric(42)
)
lifetable$Population[1] <- n_pop

# create a lifetable
for (i in 2:nrow(lifetable)) {
  age <- lifetable$Age[i]
  prev_age <- lifetable$Age[i - 1]
  lifetable$Prob_Death_Annual[i-1] <- p_background_death_gompertz(i-1, rate_gompertz, shape_gompertz)
  lifetable$Population[i] <- lifetable$Population[i - 1] * lifetable$Prob_Death_Annual[i-1]
}


# Extract age-specific mortality rate for all cycles
v_p_mort_by_age <- lifetable %>%
  dplyr::filter(Age >= n_age_init & Age <=n_age_max)%>%
  dplyr::select(Prob_Death_Annual)%>%
  as.matrix()

v_p_HDage <- rep(v_p_mort_by_age, each =1/cycle_length)

# Starting Population -----------------------------------------------------
v_m_init <- c(H = 1, D = 0) # initial state vector


## Initialize cohort trace for SoC -----------------------------------------
## m_M: trace matrix
## SoC: Standard of Care
m_M <- matrix(NA, 
              nrow = (n_cycles + 1), ncol = n_states, 
              dimnames = list(0:n_cycles, v_names_states)) 
# Store the initial state vector in the first row of the cohort trace 
m_M[1, ] <- v_m_init 

## Initialize transition probability matrix for strategy SoC ---------------
a_P_SoC <- array (0, dim = c(n_states, n_states, n_cycles+1),
               dimnames = list(v_names_states, v_names_states, 0:(n_cycles))) # row and column names


### Fill in matrix ----------------------------------------------------------
# From H
a_P_SoC["A", "A", ] <- (1 - v_p_HDage) 
a_P_SoC["A", "D", ] <- v_p_HDage
a_P_SoC["D", "A", ] <- 0
a_P_SoC["D", "D", ] <- 1

### Check if transition probability matrices are valid
a_P_SoC[, , 1]

## Check that all rows sum to 1 
rowSums(a_P_SoC[, , 1]) == 1 

# Iterative solution of time-independent cSTM
for(t in 1:n_cycles){
  # For SoC 
  m_M[t + 1, ] <- m_M[t, ] %*% a_P_SoC[, , t]
}



# Cost and effectiveness outcomes -----------------------------------------
## State rewards -----------------------------------------------------------

# Vector of state utilities under SOC
v_u_SoC <- c(H = u_H, D = u_D) * cycle_length

# Vector of state costs under SoC
v_c_SoC <- c(H = c_H, D = c_D) * cycle_length

# Vector of QALYs under SoC
v_qaly_SoC <- m_M %*% v_u_SoC 

# Vector of costs under SoC
v_cost_SoC <- m_M %*% v_c_SoC 


## Within-cycle correction -------------------------------------------------
# First, we define two functions to identify if a number is even or odd

is_even <- function(x) x %% 2 == 0 
is_odd <- function(x) x %% 2 != 0 
## Vector with cycles 
v_cycles <- seq(1, n_cycles + 1) 
length(v_cycles)
## Generate 2/3 and 4/3 multipliers for even and odd entries, respectively 
v_wcc <- is_even(v_cycles)*(2/3) + is_odd(v_cycles)*(4/3) 
length(v_wcc)
## Substitute 1/3 in first and last entries 
v_wcc[1] <- v_wcc[n_cycles + 1] <- 1/3


## Discounting future rewards ----------------------------------------------

# Discount weight for effects
v_dwe <- 1 / ((1 + (d_e * cycle_length)) ^ (0:n_cycles))

# Discount weight for costs
v_dwc <- 1 / ((1 + (d_c * cycle_length)) ^ (0:n_cycles))

## Expected discounted QALYs under SoC
n_tot_qaly_SoC <- t(v_qaly_SoC) %*% (v_dwe * v_wcc)

n_tot_qaly_SoC

## Expected discounted costs under SoC
n_tot_cost_SoC <- t(v_cost_SoC) %*% (v_dwc * v_wcc)

n_tot_cost_SoC

