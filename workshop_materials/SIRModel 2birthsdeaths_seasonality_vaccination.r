#SIR model for IIT Bombay measles workshop 2024
#Version 2: Model + births/deaths

# a function to calculate the seasonal variation in transmission rate
# timescale here is assumed to be daily timesteps, with 365 in a year
# mean.beta is as before
# seas is the amplitude of the seasonal forcing -- bigger values give bigger differences
# between the high and low season
beta.fx<-function(t,mean.beta,seas,timescale=365){mean.beta*(1-seas*cos(2*pi*(t)/timescale))}

model <- function(brate=0.0001,drate=0.0001,gamma=0.1,mean_beta=.5, seasonality = 0.2,vaccov=0.8,vactime=500,vacdegree=1,tmax=2000){
  # brate = birth rate
  # drate = death rate
  # gamma = recovery rate
  # mean_beta = mean of transmission rate
  # seasonality = amplitude of seasonality
  # vaccov = coverage with routine vaccination of children
  # vactime = the time step at which vaccination is introduced
  # vacdegree = the effectiveness of vacciantion (1 is perfect vaccination)
  # tmax = maximum time step; duration of simulation
  
  S=rep(0,tmax+1)                # storage for susceptibles
  I=rep(0,tmax+1)                # storage for infecteds
  R=rep(0,tmax+1)                # storage for recovereds
  V=rep(0,tmax+1)                # storage for vaccinated
  I[1]=0.01                      # initial proportion infected  
  S[1]=1-I[1]                    # initial proportion susceptible
  R[1]=0                         # initial proportion recovered
  V[1]=0                         # initial proportion vaccinated
  
  for(t in 1:vactime){
    beta.t = beta.fx(t,mean.beta = mean_beta, seas=seasonality)
    S[t+1]=S[t]-beta.t*S[t]*I[t]-drate*S[t]+brate*(S[t]+I[t]+R[t]+V[t])
    I[t+1]=I[t]+beta.t*S[t]*I[t]+beta.t*(1-vacdegree)*V[t]*I[t]-gamma*I[t]-drate*I[t]
    R[t+1]=R[t]+gamma*I[t]-drate*R[t]
    V[t+1]=0
  }
  for(t in (vactime+1):tmax){
    beta.t = beta.fx(t,mean.beta = mean_beta, seas=seasonality)
    S[t+1]=S[t]-beta.t*S[t]*I[t]-drate*S[t]+brate*(S[t]+I[t]+R[t]+V[t])*(1-vaccov)
    I[t+1]=I[t]+beta.t*S[t]*I[t]+beta.t*(1-vacdegree)*V[t]*I[t]-gamma*I[t]-drate*I[t]
    R[t+1]=R[t]+gamma*I[t]-drate*R[t]
    V[t+1]=V[t]+brate*(S[t]+I[t]+R[t]+V[t])*vaccov-beta.t*(1-vacdegree)*V[t]*I[t]-drate*V[t]
  }
  
  cbind(S, I, R, V)
  
}

#Plot graph
plot.model <- function(tmax=2000){

  plot(x=0:tmax,y=model.out[1:(tmax+1),"S"],type="l",ylim=c(0,1.1),xlab="time in years",ylab="proportion in class",axes=F)
  lines(x=0:tmax,y=model.out[1:(tmax+1),"I"],col="red")
  lines(x=0:tmax,y=model.out[1:(tmax+1),"R"],col="blue")
  legend("top",legend=c("S","I","R"), col=c("black","red","blue"),
         lty=c(1,1,1),bty="n",ncol=3)

}

#Example output
model.out <- model(tmax=100000, brate=0.001,drate=0.001, mean_beta=.5, seasonality=.5, vaccov=0.0) # run at baseline parameters
#plot.model(tmax=100)
plot.model(tmax=10000)    # plot the time series
axis(2)
axis(1,at = seq(365,10000, by=365), labels = 1:floor(10000/365))
abline(v=seq(1,10000,by=365),col=grey(.75))  # add vertical lines to indicate years to 

################################################################################
# explore the consequences of changing the birth rate and the strength of seasonality
# 1. change the vaccination coverage to 0.2
# 2. change the vaccination coverage to 0.4
################################################################################

q
