###############################################################################
# START FOR INTERACTIVE SESSION 1                                              #
###############################################################################
###############################################################################
# Gillespie exact                                                             #
###############################################################################
for(jj in 1:100){
###############################################################################
# Parameters and initial conditions                                            #
###############################################################################

S <- 998           # number susceptible
I <- 1              # number infected
R <- 1              # number recovered
time <- 0

beta <- .5           # transmission rate
gamma <- 1/7        # recovery rate

###############################################################################
# Gillespie Step                                                              #
###############################################################################
gilstep <- function(SIR, beta, gamma){
  #S
  #I
  #R 
  # time
  # beta = transmission rate
  # gamma = recovery rate
  
  times <- rexp(2, c(beta * SIR[1] * SIR[2]/sum(SIR[1:3]),
                     SIR[2] * gamma))
  
  return(list(change_state = which.min(times),
              timestep = min(times)))
}

###############################################################################
# Simulate over time                                                          #
###############################################################################
counter <- 0
while(all(I>0)){        #continue until I is depleted
  counter <- counter + 1    #counter for number of transitions
# object that is SIR is throwing this
# 
  SIRtmp <- c(S[counter], I[counter], R[counter], time[counter])       #current SIR states
  step <- gilstep(SIRtmp, beta, gamma)
  if(step$change_state ==1){    # if transition is an infection, reduce S and increase I
    SIRtmp[1] <- SIRtmp[1] - 1
    SIRtmp[2] <- SIRtmp[2] + 1 
  }
  if(step$change_state ==2){    # if transition is an infection, reduce S and increase R
    SIRtmp[2] <- SIRtmp[2] - 1
    SIRtmp[3] <- SIRtmp[3] + 1 
  }
  SIRtmp[4] <- SIRtmp[4] + step$timestep      # increment time
    
#Append changes
  S <- c(S,SIRtmp[1])
  I <- c(I,SIRtmp[2])
  R <- c(R,SIRtmp[3])
  time <- c(time,SIRtmp[4])

  #cat(S[counter],"-",I[counter],"-",R[counter],"-",time[counter], ".\n")
}

###############################################################################
# Plotting                                                                    #
###############################################################################

if(jj==1){plot(time,I, xlab="time", ylab="prevalence of infection", type="l", lty=1,col=rgb(0,0,0,.2),xlim=c(0,100),ylim=c(0,450))}
if(jj>1){lines(time,I, xlab="time", ylab="prevalence of infection", lty=1,col=rgb(0,0,0,.1))}
}
###############################################################################
# STOP FOR INTERACTIVE SESSION 1                                              #
###############################################################################
###############################################################################
###############################################################################
# START FOR INTERACTIVE SESSION 2                                             #
###############################################################################
###############################################################################
# Tau Leaping                                                                 #
###############################################################################
###############################################################################
# Parameters and initial conditions                                           #
###############################################################################

S <- 998           # number susceptible
I <- 1              # number infected
R <- 1              # number recovered
time <- 0

beta <- .5           # transmission rate
gamma <- 1/7        # recovery rate

###############################################################################
# Single time step                                                            #
###############################################################################
sir_step <- function (sims, S, I, R, beta,gamma, delta.t, ...) {
  # adapted from Aaron King's code
  # sims = number of simulations
  # S = initial susceptible population
  # I = initial infected population
  # R = initial recovered population
  # beta = transmission rate
  # gamma = recovery rate
  
  N <- S+I+R    # total population size
  dSI <- rpois(n=sims,beta*S*(I/N)*delta.t)  # new incident infections
  dIR <- rpois(n=sims,gamma*I*delta.t)       # recoveries
  # dSI <- rbinom(n=sims,size=S,prob=1-exp(-beta*(I/N)*delta.t))  # new incident infections
  # dIR <- rbinom(n=sims,size=I,prob=1-exp(-gamma*delta.t))       # recoveries
  
  S <- S - dSI            # change in S
  I <- I + dSI - dIR      # change in I
  R <- R + dIR            # change in R
  cbind(S, I, R, dSI) # note that dSI are the new incident infections
}

###############################################################################
# set up and storage for states                                               #
###############################################################################
T <- 100        # time to simulate over, we only care about start
sims <- 1000    # number of simulations

Smat <- matrix(S,1,sims)      # storage item for S for all simulations
Imat <- matrix(I,1,sims)      # storage item for I for all simulations
Rmat <- matrix(R,1,sims)      # storage item for R for all simulations
new_cases <- matrix(0,1,sims) # storage item for N for all simulations
Nmat <- S+I+R

for(ts in 2:T){     #loop over time, ts is the index
  
  out <- sir_step(sims, Smat[ts-1,], Imat[ts-1,], Rmat[ts-1,], beta, gamma, delta.t=1) # call to SIR step function above

  Smat <- rbind(Smat,out[,1])  # update state
  Imat <- rbind(Imat,out[,2])  # update state
  Rmat <- rbind(Rmat,out[,3])  # update state
  Nmat <- rbind(Nmat,out[,1]+out[,2]+out[,3])  # update state -- note population size isn't changing, but this could be updated with births/deaths
  new_cases <- rbind(new_cases,out[,4])  # update state
}

###############################################################################
# plotting                                                                   #
###############################################################################
matplot(Imat, type="l", lty=1,col=rgb(0,0,1,.1),add=T)
###############################################################################
# STOP FOR INTERACTIVE SESSION 2                                              #
###############################################################################
###############################################################################
###############################################################################
# START FOR INTERACTIVE SESSION 3                                             #
###############################################################################
###############################################################################
# Chain Binomial                                                              #
###############################################################################
###############################################################################
# Parameters and initial conditions                                           #
###############################################################################
S <- 998           # number susceptible
I <- 1              # number infected
R <- 1              # number recovered
time <- 0

beta <- 3.5           # transmission rate

###############################################################################
# Single time step                                                            #
###############################################################################
cb_step <- function (sims, S, I, R, beta, ...) {
  # adapted from Aaron King's code
  # S = initial susceptible population
  # I = initial infected population
  # R = initial recovered population
  # beta is transmission rate
  
  N <- S+I+R    # total population size
  newI <- rbinom(n=sims, S, 1-exp(-beta * I/N))
  newS <- S - newI
  newR <- I

  cbind(newS, newI, newR) 
}

###############################################################################
# Run over many steps                                                         #
###############################################################################

T <- 20 # note here that the time step is one infectious generation time, so 7 days from gamma above
Smat <- matrix(S,1,sims)      # storage item for S for all simulations
Imat <- matrix(I,1,sims)      # storage item for I for all simulations
Rmat <- matrix(R,1,sims)      # storage item for R for all simulations
Nmat <- S+I+R

for(ts in 2:T){
  out <- cb_step(sims, Smat[ts-1,], Imat[ts-1,], Rmat[ts-1,], beta)
  Smat <- rbind(Smat,out[,1])  # update state
  Imat <- rbind(Imat,out[,2])  # update state
  Rmat <- rbind(Rmat,out[,3])  # update state
  Nmat <- rbind(Nmat,out[,1]+out[,2]+out[,3])  # update state -- note population size isn't changing, but this could be updated with births/deaths
  
}

###############################################################################
# plotting                                                                   #
###############################################################################
matplot(seq(1,7*T,by=7), Imat,type="l", lty=1, col=rgb(1,0,0,.1),add=T)
###############################################################################
# STOP FOR INTERACTIVE SESSION 3                                              #
###############################################################################



