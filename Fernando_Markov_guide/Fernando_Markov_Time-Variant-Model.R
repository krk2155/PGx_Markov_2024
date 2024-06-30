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
survival_gompertz <- function(t, rate, shape) {
  exp(-rate / shape * (exp(shape * t) - 1))
}


# Create an empty dataframe for the lifetable
lifetable <- data.frame(
  Age = n_age_init:(n_age_init + 40),  # Let's consider up to age 80
  Population = numeric(41)
)
lifetable$Population[1] <- n_pop

# create a lifetable
for (i in 2:nrow(lifetable)) {
  age <- lifetable$Age[i]
  prev_age <- lifetable$Age[i - 1]
  S_prev <- survival_gompertz(prev_age, rate_gompertz, shape_gompertz)
  S_current <- survival_gompertz(age, rate_gompertz, shape_gompertz)
  lifetable$Population[i] <- lifetable$Population[i - 1] * (S_current / S_prev)
}

# Print the lifetable
print(lifetable)

# Add survival probability and hazard rate columns to the lifetable
lifetable$Survival_Probability <- numeric(41)
lifetable$Hazard_Rate <- numeric(41)

# Calculate survival probabilities and hazard rates
for (i in 1:nrow(lifetable)) {
  age <- lifetable$Age[i]
  if (i > 1) {
    lifetable$Survival_Probability[i] <- lifetable$Population[i] / lifetable$Population[i - 1]
  } else {
    lifetable$Survival_Probability[i] <- 1  # Initial survival probability is 1
  }
  lifetable$Hazard_Rate[i] <- rate_gompertz * exp(shape_gompertz * age)
}

v_r_HDage <- rep(v_r_mort_by_age, each = 1/cycle_length)

# NOTE (20240630): 
# 1. Double-check if the lifetable was created correctly
# 2. add the values from the Hazard Rate column to "v_r_HDage"
# 3. convert v_r_HDage into v_p_HDage
# 4. use v_p_HDage in the matrix for Alive -> Dead transition probability

# Print the extended lifetable
print(lifetable)



## Probabilities of dying (Cycle-specific) ----------------------------------
r_HD <- rate_gompertz*(exp(shape_gompertz*t)) # rate of background mortality
p_HD <- 1-exp(-r_HD) # probability of background mortality


# Starting Population -----------------------------------------------------
v_m_init <- c(H = 1, D = 0) # initial state vector


## Initialize cohort trace for SoC -----------------------------------------
## Recall: "SoC" -> "Standard of Care"
m_M <- matrix(NA, 
              nrow = (n_cycles + 1), ncol = n_states, 
              dimnames = list(0:n_cycles, v_names_states)) 
# Store the initial state vector in the first row of the cohort trace 
m_M[1, ] <- v_m_init 

## Initialize cohort trace for strategies A, B, and AB  --------------------
#Structure and initial states are the same as for SoC 
m_M_No_PGx <- m_M # Strategy Without PGx Testing
m_M_PGx <- m_M # Strategy With PGx Testing



## Initialize transition probability matrix for strategy SoC ---------------
a_P_SoC <- array (0, dim = c(n_states, n_states, n_cycles),
               dimnames = list(v_names_states, v_names_states, 0:(n_cycles-1))) # row and column names

### Fill in matrix ----------------------------------------------------------
# From H

m_P["A", "A", ] <- (1 - p_HD) 
m_P["A", "D", ] <- p_HD
m_P["D", "A", ] <- 0
m_P["D", "D", ] <- 1

## Initialize transition probability matrix for strategy A as a cop --------

m_M_No_PGx <- m_P


## Initialize transition probability matrix for strategy B -----------------
m_M_PGx <- m_P ## Update only transition probabilities from S1 involving p_S1S2

### Check if transition probability matrices are valid
## Check that transition probabilities are [0, 1]

## Check that all rows sum to 1 
rowSums(m_P) == 1 
rowSums(m_M_No_PGx) == 1 
rowSums(m_M_PGx) == 1 

# Iterative solution of time-independent cSTM

for(t in 1:n_cycles){
  
  # For SoC 
  m_M[t + 1, ] <- m_M[t, ] %*% m_P[, , t]
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


# Incremental cost-effectiveness ratios (ICERs) ---------------------------
### Vector of costs

v_cost_str <- c(n_tot_cost_SoC) 

### Vector of effectiveness 
v_qaly_str <- c(n_tot_qaly_SoC) 

calculate_icers <- function(cost, effect, strategies) {
  # checks on input
  n_cost <- length(cost)
  n_eff <- length(effect)
  n_strat <- length(strategies)
  if (n_cost != n_eff | n_eff != n_strat) {
    stop("cost, effect, and strategies must all be vectors of the same length", call. = FALSE)
  }
  
  # coerce to character, in case they are provided as numeric
  char_strat <- as.character(strategies)
  
  # create data frame to hold data
  df <- data.frame("Strategy" = char_strat,
                   "Cost" = cost,
                   "Effect" = effect,
                   stringsAsFactors = FALSE)
  nstrat <- nrow(df)
  
  # if only one strategy was provided, return df with NAs for incremental
  if (nstrat == 1) {
    df[, c("ICER", "Inc_Cost", "Inc_Effect")] <- NA
    return(df)
  }
  
  # three statuses: dominated, extended dominated, and non-dominated
  d <- NULL
  
  # detect dominated strategies
  # dominated strategies have a higher cost and lower effect
  df <- df %>%
    arrange(.data$Cost, desc(.data$Effect))
  
  # iterate over strategies and detect (strongly) dominated strategies
  # those with higher cost and equal or lower effect
  for (i in 1:(nstrat - 1)) {
    ith_effect <- df[i, "Effect"]
    for (j in (i + 1):nstrat) {
      jth_effect <- df[j, "Effect"]
      if (jth_effect <= ith_effect) {
        # append dominated strategies to vector
        d <- c(d, df[j, "Strategy"])
      }
    }
  }
  
  # detect weakly dominated strategies (extended dominance)
  # this needs to be repeated until there are no more ED strategies
  ed <- vector()
  continue <- TRUE  # ensure that the loop is run at least once
  while (continue) {
    # vector of all dominated strategies (strong or weak)
    dom <- union(d, ed)
    
    # strategies declared to be non-dominated at this point
    nd <- setdiff(strategies, dom)
    
    # compute icers for nd strategies
    nd_df <- df[df$Strategy %in% nd, ] %>%
      compute_icers()
    
    # number non-d
    n_non_d <- nrow(nd_df)
    
    # if only two strategies left, we're done
    if (n_non_d <= 2) {
      break
    }
    
    # strategy identifiers for non-d
    nd_strat <- nd_df$Strategy
    
    # now, go through non-d strategies and detect any
    # with higher ICER than following strategy
    ## keep track of whether any ED strategies are picked up
    # if not, we're done - exit the loop
    new_ed <- 0
    for (i in 2:(n_non_d - 1)) {
      if (nd_df[i, "ICER"] > nd_df[i + 1, "ICER"]) {
        ed <- c(ed, nd_strat[i])
        new_ed <- new_ed + 1
      }
    }
    if (new_ed == 0) {
      continue <- FALSE
    }
  }
  
  # recompute icers without weakly dominated strategies
  nd_df_icers <- nd_df[!(nd_df$Strategy %in% dom), ] %>%
    mutate(Status = "ND") %>%
    compute_icers()
  
  # dominated and weakly dominated
  d_df <- df[df$Strategy %in% d, ] %>%
    mutate(ICER = NA, Status = "D")
  
  ed_df <- df[df$Strategy %in% ed, ] %>%
    mutate(ICER = NA, Status = "ED")
  
  # when combining, sort so we have ref,ND,ED,D
  results <- bind_rows(d_df, ed_df, nd_df_icers) %>%
    arrange(desc(.data$Status), .data$Cost, desc(.data$Effect))
  
  # re-arrange columns
  results <- results %>%
    select(.data$Strategy, .data$Cost, .data$Effect,
           .data$Inc_Cost, .data$Inc_Effect, .data$ICER, .data$Status)
  
  # declare class of results
  class(results) <- c("icers", "data.frame")
  return(results)
}

### Calculate incremental cost-effectiveness ratios (ICERs) 
df_cea <- calculate_icers(cost = v_cost_str,
                                   effect = v_qaly_str, 
                                   strategies = v_names_str)
df_cea

# Probabilistic sensitivity analysis --------------------------------------



