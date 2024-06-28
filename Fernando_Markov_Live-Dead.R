#Mac
setwd("/Users/kun-wookim/Library/CloudStorage/OneDrive-VUMC/Research_discrete-event-simulation/r/PGx_Markov_2024")

# Load packages -----------------------------------------------------------
load.lib<-c("flexsurv", "msm", "dplyr", "markovchain", "heemod",
            "ggplot2", "reshape2")
install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)

## General setup

cycle_length <- 1 # cycle length equal one year (use 1/12 for monthly)

n_age_init <- 25 # age at baseline 
n_age_max <- 100 # maximum age of follow up 
n_cycles <- (n_age_max - n_age_init)/cycle_length # time horizon, number of cycles 
v_names_states <- c("H", "S1", "S2", "D") # the 4 health states of the model: # Healthy (H), Sick (S1), Sicker (S2), Dead (D) 
n_states <- length(v_names_states) # number of health states 
d_e <- 0.03 # annual discount rate for QALYs of 3% 
d_c <- 0.03 # annual discount rate for costs of 3% 
v_names_str <- c("Standard of care", # store the strategy names
                "Strategy A", 
                "Strategy B", 
                "Strategy AB")


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



