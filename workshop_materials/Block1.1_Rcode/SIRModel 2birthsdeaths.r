#SIR model for IIT Bombay measles workshop 2024
#Version 2: Model + births/deaths

#Definitions of variables and parameters
#S = proportion of people who are susceptible to the infection
#I = proportion of people who are infected (and infectious)
#R = proportion of people who are resistant to the infection
#tmax = maximum time to run simulation
#rec_rate = rate (in person-days) of recovery from the infection
#beta = probability of transmission between any two individuals
#brate = rate (in person-days) of births
#drate = rate (in person-days) of deaths

model <- function(brate=0.0001,drate=0.0001,rec_rate=0.1,beta=1,tmax=2000){

  S=rep(0,tmax+1)                                                                           #susceptibles
  I=rep(0,tmax+1)                                                                            #infecteds
  R=rep(0,tmax+1)                                                                           #recovereds
  I[1]=0.01
  S[1]=1-I[1]
  R[1]=0

  for(t in 1:tmax){
    S[t+1]=S[t]-beta*S[t]*I[t]-drate*S[t]+brate*(S[t]+I[t]+R[t])
    I[t+1]=I[t]+beta*S[t]*I[t]-rec_rate*I[t]-drate*I[t]
    R[t+1]=R[t]+rec_rate*I[t]-drate*R[t]
  }

  cbind(S, I, R)
}

#Plot graph
plot.model <- function(tmax=2000){

  plot(x=0:tmax,y=model.out[1:(tmax+1),"S"],type="l",ylim=c(0,1.1),xlab="time",ylab="")
  lines(x=0:tmax,y=model.out[1:(tmax+1),"I"],col="red")
  lines(x=0:tmax,y=model.out[1:(tmax+1),"R"],col="blue")
  legend("top",legend=c("S","I","R"), col=c("black","red","blue"),
         lty=c(1,1,1),bty="n",ncol=3)

}

#Example output
model.out <- model(tmax=10000, brate=0.0001,drate=0.0001)
plot.model(tmax=100)
plot.model(tmax=10000)