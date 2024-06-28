# Name: Kun-Woo Kim
# Date                        #Task
#------------     ---------------------------
# 6/27/2024       - Create basic Markov model   

#Mac
setwd("/Users/kun-wookim/Library/CloudStorage/OneDrive-VUMC/Research_discrete-event-simulation/r/PGx_Markov_2024")

# Load packages -----------------------------------------------------------
load.lib<-c("flexsurv", "msm", "dplyr",
            "ggplot2", "reshape2")
install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)

# Define Parameters --------------------------------------------------------------
p_Die_background <- 1.0  # Expression: 1
shape_gompertz <- 0.101  # Expression: 0.101
rate_gompertz <- 0.001  # Expression: 0.001
r_a_pct <- 0.1  # Expression: 0.1
r_a_dur <- 10.0  # Expression: 10
r_condition_annual <- -log(1-r_a_pct)/r_a_dur  # Expression: 0.01053605
p_condition_annual <- 1-exp(-r_condition_annual)  # Expression: 0.01048074

# Discount Rate
discountRates <- c(0.03, 0.03)
startDiscountCycle <- 0
halfCycle <- TRUE


# Define Variables --------------------------------------------------------------
# time variable
t <- 0

# cohort size
cohortSize <- 1000

# Create Live-Dead status variable


# Initialize Prevalence --------------------------------------------------------------

# Define the transition matrix
Q <- matrix(c(0, 1, 0, 0), nrow = 2, byrow = TRUE)
rownames(Q) <- colnames(Q) <- c("live", "dead")


curPrev <- c()
curPrev[1] <- cohortSize * 1  # Live
curPrev[2] <- cohortSize * 0  # Condition
curPrev[3] <- cohortSize * 0  # Die
newPrev <- curPrev  # copy inital prev

# Define the survival function for the Gompertz model --------------------------------------------------------------
# Fit a Gompertz model using flexsurvreg
survival_gompertz <- function(t) {
  exp(-rate_gompertz / shape_gompertz * (exp(shape_gompertz * t) - 1))
}

survival_gompertz <- function(t) {
  1-shape_gompertz * (exp(rate_gompertz * t) - 1)
}


########## Markov Chain: Live-Die ##########
print("Running Markov Chain: Live-Die ")
numStates <- 2
colNames <- c("Cycle", "Live", "Die", "Cycle_Cost", "Cum_Cost", "Cycle_Dis_Cost", "Cum_Dis_Cost", "Cycle_QALE", "Cum_QALE", "Cycle_Dis_QALE", "Cum_Dis_QALE")
trace <- data.frame(matrix(nrow=0, ncol=11))
names(trace) <- colNames

# Initialize prevalence
curPrev <- c()
curPrev[1] <- cohortSize * 1  # Live
curPrev[2] <- cohortSize * 0  # Die
newPrev <- curPrev  # copy inital prev

# Initialize variables
t <- 0  # initialize cycle
r_Die_bg <- 1-exp(-shape_gompertz*(exp(rate_gompertz * t) - 1))
p_Die_bg <- 1-exp(-r_Die_bg)
             
# Run chain
# Initialize outcomes
Live_minusDie_Cost <- 0; Live_minusDie_Dis_Cost <- 0
Live_minusDie_QALE <- 0; Live_minusDie_Dis_QALE <- 0

terminate <- FALSE
while (terminate == FALSE) {
 # Update progress
 if (t %% 10 == 0){
   cat(t, sep="")
 } else {
   cat(".", sep="")
 }
 
 # Cycle outcomes
 cycleCost <- 0; cycleCost_dis <- 0
 cycleQALE <- 0; cycleQALE_dis <- 0
 
 # Update prevalence
 curPrev <- newPrev
 # Simulate state transitions
 
 # Live
 # Update rewards
 cycleCost <- cycleCost+curPrev[1] * 0
 cycleQALE <- cycleQALE+curPrev[1] * 1
 
 # Calculate child probs
 sumProb <- 0
 childProbs <- c()
 
 # P(New Die|Prev Survived)
 childProbs[2] <- p_Die_bg; sumProb <- sumProb + childProbs[2]  # Prob Background Death
 # P(New Survive|Prev Survived)
 childProbs[1] <- 1.0 - sumProb  # Complementary prob
 
 # N of New Surviving Pop (prev. Surviving Pop x % of Surviving)
 prev_Survive <- curPrev[1] * childProbs[1]
 # N of New Dead Pop (prev. Surviving Pop x % of Dead)
 prev_BackgroundDeath <- curPrev[1] * childProbs[2]
 # Survive
 # Transition to Live
 newPrev[1] <- newPrev[1] - prev_Survive
 newPrev[1] <- newPrev[1] + prev_Survive
 # Background Death
 # Transition to Die
 newPrev[1] <- newPrev[1] - prev_BackgroundDeath
 newPrev[2] <- newPrev[2] + prev_BackgroundDeath
 
 # Die
 # Update rewards
 cycleCost <- cycleCost+curPrev[2] * 0
 cycleQALE <- cycleQALE+curPrev[2] * 0
 
 # Calculate child probs
 sumProb <- 0
 childProbs <- c()
 childProbs[1] <- 1; sumProb <- sumProb + childProbs[1]  # Prob Dead
 prev_Dead <- curPrev[2] * childProbs[1]
 # Dead
 # Transition to Die
 newPrev[2] <- newPrev[2] - prev_Dead
 newPrev[2] <- newPrev[2] + prev_Dead
 
 # Update outcomes
 discount <- ifelse(t >= startDiscountCycle, 1 / ((1 + discountRates[1]) ^ (t - startDiscountCycle + 1)), 1)
 cycleCost_dis <- cycleCost * discount
 if (t == 0 & halfCycle) {  # half-cycle correction
   cycleCost <- 0.5 * cycleCost; cycleCost_dis <- 0.5 * cycleCost_dis
 }
 Live_minusDie_Cost <- Live_minusDie_Cost + cycleCost
 Live_minusDie_Dis_Cost <- Live_minusDie_Dis_Cost + cycleCost_dis
 discount <- ifelse(t >= startDiscountCycle, 1 / ((1 + discountRates[2]) ^ (t - startDiscountCycle + 1)), 1)
 cycleQALE_dis <- cycleQALE * discount
 if (t == 0 & halfCycle) {  # half-cycle correction
   cycleQALE <- 0.5 * cycleQALE; cycleQALE_dis <- 0.5 * cycleQALE_dis
 }
 Live_minusDie_QALE <- Live_minusDie_QALE + cycleQALE
 Live_minusDie_Dis_QALE <- Live_minusDie_Dis_QALE + cycleQALE_dis
 
 # Update trace
 row <- c(t, curPrev, cycleCost, Live_minusDie_Cost, cycleCost_dis, Live_minusDie_Dis_Cost, cycleQALE, Live_minusDie_QALE, cycleQALE_dis, Live_minusDie_Dis_QALE)
 row <- data.frame(matrix(row, nrow=1))
 names(row) <- colNames
 trace <- rbind(trace, row)
 
 # Check termination condition
 terminate <- (t==40)
 if (terminate & halfCycle){  # half cycle-correction, update last trace row
   lastRow <- nrow(trace)
   Live_minusDie_Cost <- Live_minusDie_Cost - cycleCost * 0.5
   trace[lastRow, 4] <- trace[lastRow, 4] * 0.5  # Cost
   trace[lastRow, 5] <- trace[lastRow - 1, 5] + trace[lastRow, 4]  # cum Cost
   Live_minusDie_Dis_Cost <- Live_minusDie_Dis_Cost - cycleCost_dis * 0.5
   trace[lastRow, 6] <- trace[lastRow, 6] * 0.5  # Cost discounted
   trace[lastRow, 7] <- trace[lastRow - 1, 7] + trace[lastRow, 6]  # cum Cost discounted
   Live_minusDie_QALE <- Live_minusDie_QALE - cycleQALE * 0.5
   trace[lastRow, 8] <- trace[lastRow, 8] * 0.5  # QALE
   trace[lastRow, 9] <- trace[lastRow - 1, 9] + trace[lastRow, 8]  # cum QALE
   Live_minusDie_Dis_QALE <- Live_minusDie_Dis_QALE - cycleQALE_dis * 0.5
   trace[lastRow, 10] <- trace[lastRow, 10] * 0.5  # QALE discounted
   trace[lastRow, 11] <- trace[lastRow - 1, 11] + trace[lastRow, 10]  # cum QALE discounted
 }
 
 t <- t + 1  # next cycle
 
}  # end cycle while loop
setwd("/Users/kun-wookim/Library/CloudStorage/OneDrive-VUMC/Research_discrete-event-simulation/r/PGx_Markov_2024")
write.csv(trace, "Live-Die_Trace_20240627.csv")  # write trace
cat("done!")

# Report totals
print(paste("Cost:", Live_minusDie_Cost))
print(paste("Cost (Discounted):", Live_minusDie_Dis_Cost))
print(paste("QALE:", Live_minusDie_QALE))
print(paste("QALE (Discounted):", Live_minusDie_Dis_QALE))
print("")

             
