#Mac
setwd("/Users/kun-wookim/Library/CloudStorage/OneDrive-VUMC/Research_discrete-event-simulation/r/PGx_Markov_2024")

# Load packages -----------------------------------------------------------
load.lib<-c("flexsurv", "msm", "dplyr", "markovchain", "heemod", "dampack",
            "ggplot2", "reshape2")
install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)

## General setup

cycle_length <- 1 # cycle length equal one year (use 1/12 for monthly)

n_age_init <- 40 # age at baseline 
n_age_max <- 120 # maximum age of follow up 
n_cycles <- (n_age_max - n_age_init)/cycle_length # time horizon, number of cycles 
v_names_states <- c("A","D") # the 2 health states of the model: # Alive (A), Dead (D)
n_states <- length(v_names_states) # number of health states 
d_e <- 0.03 # annual discount rate for QALYs of 3% 
d_c <- 0.03 # annual discount rate for costs of 3% 
v_names_str <- c("Standard of care") # store the strategy names



# Transition probabilities (annual), and hazard ratios (HRs)  -------------
r_HD <- 0.002 # constant annual rate of dying when Healthy (all-cause mortality rate) 
r_HS1 <- 0.15 # constant annual rate of becoming Sick when Healthy 
r_S1H <- 0.5 # constant annual rate of becoming Healthy when Sick 
r_S1S2 <- 0.105 # constant annual rate of becoming Sicker when Sick 
hr_S1 <- 3 # hazard ratio of death in Sick vs Healthy 
hr_S2 <- 10 # hazard ratio of death in Sicker vs Healthy

## Constant transition probability of becoming Sick when Healthy 
# transform rate to probability and scale by the cycle length 
p_HS1 <- 1 - exp(-r_HS1 * cycle_length) 
## Constant transition probability of becoming Healthy when Sick 
# transform rate to probability and scale by the cycle length 
p_S1H <- 1 - exp(-r_S1H * cycle_length) 
## Constant transition probability of becoming Sicker when Sick 
# transform rate to probability and scale by the cycle length 
p_S1S2 <- 1 - exp(-r_S1S2 * cycle_length)

# Effectiveness of treatment B 
hr_S1S2_trtB <- 0.6 # hazard ratio of becoming Sicker when Sick under treatment B


# State rewards -----------------------------------------------------------
## Costs -------------------------------------------------------------------

c_H <-2000 # annual cost of being Healthy 
c_S1 <- 4000 # annual cost of being Sick 
c_S2 <- 15000 # annual cost of being Sicker 
c_D <- 0 # annual cost of being dead

c_trtA <- 12000 # annual cost of receiving treatment A 
c_trtB <- 13000 # annual cost of receiving treatment B 


## Utilities ---------------------------------------------------------------
u_H <- 1 # annual utility of being Healthy 
u_S1 <- 0.75 # annual utility of being Sick 
u_S2 <- 0.5 # annual utility of being Sicker 
u_D <- 0 # annual utility of being dead 
u_trtA <- 0.95 # annual utility when receiving treatment A


# Transition Probabilities2 -----------------------------------------------
## Mortality rates

r_S1D <- r_HD * hr_S1 # annual rate of dying when Sick 
r_S2D <- r_HD * hr_S2 # annual rate of dying when Sicker 

## Probabilities of dying (Cycle-specific) ----------------------------------
cycle_length <- 1 
p_HD <- 1 - exp(-r_HD * cycle_length) # annual background mortality risk (i.e., probability) 
p_S1D <- 1 - exp(-r_S1D * cycle_length) # annual probability of dying when Sick 
p_S2D <- 1 - exp(-r_S2D * cycle_length) # annual probability of dying when Sicker



## Treatment Effect on Transition Sick->Sicker -----------------------------


## Transition probability of becoming Sicker when Sick for treatment B
# apply hazard ratio to rate to obtain transition rate of becoming Sicker when Sick # for treatment B 
r_S1S2_trtB <- r_S1S2 * hr_S1S2_trtB 
# transform rate to probability 
# probability to become Sicker when Sick under treatment B 
p_S1S2_trtB <- 1 - exp(-r_S1S2_trtB * cycle_length)


# Starting Population -----------------------------------------------------
v_m_init <- c(H = 1, S1 = 0, S2 = 0, D = 0) # initial state vector


## Initialize cohort trace for SoC -----------------------------------------
## Recall: "SoC" -> "Standard of Care"
m_M <- matrix(NA, 
              nrow = (n_cycles + 1), ncol = n_states, 
              dimnames = list(0:n_cycles, v_names_states)) 
# Store the initial state vector in the first row of the cohort trace 
m_M[1, ] <- v_m_init 

## Initialize cohort trace for strategies A, B, and AB  --------------------
#Structure and initial states are the same as for SoC 
m_M_strA <- m_M # Strategy A 
m_M_strB <- m_M # Strategy B 
m_M_strAB <- m_M # Strategy AB


## Initialize transition probability matrix for strategy SoC ---------------
m_P <- matrix (0, nrow = n_states, ncol = n_states,
               dimnames = list(v_names_states, v_names_states)) # row and column names

mode(v_names_states)
### Fill in matrix ----------------------------------------------------------
# From H

m_P["H", "H"] <- (1 - p_HD) * (1 - p_HS1) 
m_P["H", "S1"] <- (1 - p_HD) * p_HS1 
m_P["H", "D"] <- p_HD # From S1 
m_P["S1", "H"] <- (1 - p_S1D) * p_S1H

m_P["S1", "S1"] <- (1 - p_S1D) * (1 - (p_S1H + p_S1S2))

m_P["S1", "S2"] <- (1 - p_S1D) * p_S1S2 
m_P["S1", "D"] <- p_S1D # From S2 
m_P["S2", "S2"] <- 1 - p_S2D 
m_P["S2", "D"] <- p_S2D # From D 
m_P["D", "D"] <- 1


## Initialize transition probability matrix for strategy A as a cop --------

m_P_strA <- m_P


## Initialize transition probability matrix for strategy B -----------------
m_P_strB <- m_P ## Update only transition probabilities from S1 involving p_S1S2

m_P_strB["S1", "S1"] <- (1 - p_S1D) * (1 - (p_S1H + p_S1S2_trtB))

m_P_strB["S1", "S2"] <- (1 - p_S1D) * p_S1S2_trtB


##Initialize transition probability matrix for strategy AB as a co --------
m_P_strAB <- m_P_strB

### Check if transition probability matrices are valid
## Check that transition probabilities are [0, 1]

m_P >= 0 && m_P <= 1 
m_P_strA >= 0 && m_P_strA <= 1 
m_P_strB >= 0 && m_P_strB <= 1 
m_P_strAB >= 0 && m_P_strAB <= 1 

## Check that all rows sum to 1 
rowSums(m_P) == 1 
rowSums(m_P_strA) == 1 
rowSums(m_P_strB) == 1 
rowSums(m_P_strAB) == 1

# Iterative solution of time-independent cSTM

for(t in 1:n_cycles){
  
  # For SoC 
  m_M[t + 1, ] <- m_M[t, ] %*% m_P 
  # For strategy A 
  m_M_strA[t + 1, ] <- m_M_strA[t, ] %*% m_P_strA 
  # For strategy B 
  m_M_strB[t + 1, ] <- m_M_strB[t, ] %*% m_P_strB 
  # For strategy AB 
  m_M_strAB[t + 1, ] <- m_M_strAB[t, ] %*% m_P_strAB
  
}


# Cost and effectiveness outcomes -----------------------------------------
## State rewards -----------------------------------------------------------

# Vector of state utilities under SOC
v_u_SoC <- c(H = u_H, S1 = u_S1, S2 = u_S2, D = u_D) * cycle_length

# Vector of state costs under SoC
v_c_SoC <- c(H = c_H, S1 = c_S1, S2 = c_S2, D = c_D) * cycle_length

# Vector of state utilities for strategy A
v_u_strA <- c(H = u_H, S1 = u_trtA, S2 = u_S2, D = u_D) * cycle_length

# Vector of state utilities for strategy B
v_u_strB <- c(H = u_H, S1 = u_S1, S2 = u_S2, D = u_D) * cycle_length

# Vector of state utilities for strategy AB
v_u_strAB <- c(H = u_H, S1 = u_trtA, S2 = u_S2, D = u_D) * cycle_length

# Vector of state costs for strategy A

v_c_strA <- c(H = c_H,
              S1 = c_S1 +c_trtA,
              S2 = c_S2 + c_trtA, 
              D = c_D) * cycle_length 

# Vector of state costs for strategy B 
v_c_strB <- c(H = c_H, 
                S1 = c_S1 + c_trtB, 
                S2 = c_S2 + c_trtB, 
                D = c_D) * cycle_length 
# Vector of state costs for strategy AB 
v_c_strAB <- c(H = c_H, 
               S1 = c_S1 + (c_trtA + c_trtB), 
               S2 = c_S2 + (c_trtA + c_trtB), 
               D = c_D) * cycle_length


# Vector of QALYs under SoC
v_qaly_SoC <- m_M %*% v_u_SoC 

# Vector of costs under SoC
v_cost_SoC <- m_M %*% v_c_SoC 

# Vector of QALYs for strategy A
v_qaly_strA <- m_M_strA %*% v_u_strA 

# Vector of costs for strategy A
v_cost_strA <- m_M_strA %*% v_c_strA 

# Vector of QALYs for strategy B
v_qaly_strB <- m_M_strB %*% v_u_strB 

# Vector of costs for strategy B
v_cost_strB <- m_M_strB %*% v_c_strB 

# Vector of QALYs for strategy AB
v_qaly_strAB <- m_M_strAB %*% v_u_strAB 

# Vector of costs for strategy AB
v_cost_strAB <- m_M_strAB %*% v_c_strAB


## Within-cycle correction -------------------------------------------------
# First, we define two functions to identify if a number is even or odd

is_even <- function(x) x %% 2 == 0 
is_odd <- function(x) x %% 2 != 0 
## Vector with cycles 
v_cycles <- seq(1, n_cycles + 1) 

## Generate 2/3 and 4/3 multipliers for even and odd entries, respectively 
v_wcc <- is_even(v_cycles)*(2/3) + is_odd(v_cycles)*(4/3) 

## Substitute 1/3 in first and last entries 
v_wcc[1] <- v_wcc[n_cycles + 1] <- 1/3


## Discounting future rewards ----------------------------------------------

# Discount weight for effects
v_dwe <- 1 / ((1 + (d_e * cycle_length)) ^ (0:n_cycles))

# Discount weight for costs
v_dwc <- 1 / ((1 + (d_c * cycle_length)) ^ (0:n_cycles))

## Expected discounted QALYs under SoC
n_tot_qaly_SoC <- t(v_qaly_SoC) %*% (v_dwe * v_wcc)

## Expected discounted costs under SoC
n_tot_cost_SoC <- t(v_cost_SoC) %*% (v_dwc * v_wcc)

## Expected discounted QALYs for strategy A
n_tot_qaly_strA <- t(v_qaly_strA) %*% (v_dwe * v_wcc)

## Expected discounted costs for strategy A
n_tot_cost_strA <- t(v_cost_strA) %*% (v_dwc * v_wcc)

## Expected discounted QALYs for strategy B
n_tot_qaly_strB <- t(v_qaly_strB) %*% (v_dwe * v_wcc)

## Expected discounted costs for strategy B
n_tot_cost_strB <- t(v_cost_strB) %*% (v_dwc * v_wcc)

## Expected discounted QALYs for strategy AB
n_tot_qaly_strAB <- t(v_qaly_strAB) %*% (v_dwe * v_wcc)

## Expected discounted costs for strategy AB
n_tot_cost_strAB <- t(v_cost_strAB) %*% (v_dwc * v_wcc)


# Incremental cost-effectiveness ratios (ICERs) ---------------------------
### Vector of costs

v_cost_str <- c(n_tot_cost_SoC, 
                n_tot_cost_strA, 
                n_tot_cost_strB, 
                n_tot_cost_strAB) 

### Vector of effectiveness 
v_qaly_str <- c(n_tot_qaly_SoC, 
                n_tot_qaly_strA, 
                n_tot_qaly_strB, 
                n_tot_qaly_strAB) 

### Calculate incremental cost-effectiveness ratios (ICERs) 
df_cea <- dampack::calculate_icers(cost = v_cost_str,
                                   effect = v_qaly_str, 
                                   strategies = v_names_str)


# Probabilistic sensitivity analysis --------------------------------------



