#SIR model for IIT Bombay measles workshop 2024
#Version 3: Model + births/deaths + vaccination

#Definitions of variables and parameters
#S = proportion of people who are susceptible to the infection
#I = proportion of people who are infected (and infectious)
#R = proportion of people who are resistant to the infection
#tmax = maximum time to run simulation
#rec_rate = rate (in person-days) of recovery from the infection
#beta = probability of transmission between any two individuals
#brate = rate (in person-days) of births
#drate = rate (in person-days) of deaths
#vaccov = vaccine coverage
#vactime = time (in days) when vaccination starts
#vacdegree = vaccine degree (reduction in risk of infection to a vaccinated individual)

model <- function(brate=0.00005,drate=0.00005,rec_rate=0.1,beta=1,vaccov=0.8,vactime=500,vacdegree=1,tmax=2000){

  S=rep(0,tmax+1)                                                                           #susceptibles
  I=rep(0,tmax+1)                                                                            #infecteds
  R=rep(0,tmax+1)                                                                           #recovereds
  V=rep(0,tmax+1)                                                                           #vaccinated
  I[1]=0.01
  S[1]=1-I[1]
  R[1]=0
  V[1]=0
 
  for(t in 1:vactime){
    S[t+1]=S[t]-beta*S[t]*I[t]-drate*S[t]+brate*(S[t]+I[t]+R[t]+V[t])
    I[t+1]=I[t]+beta*S[t]*I[t]+beta*(1-vacdegree)*V[t]*I[t]-rec_rate*I[t]-drate*I[t]
    R[t+1]=R[t]+rec_rate*I[t]-drate*R[t]
    V[t+1]=0
  }
  for(t in (vactime+1):tmax){
    S[t+1]=S[t]-beta*S[t]*I[t]-drate*S[t]+brate*(S[t]+I[t]+R[t]+V[t])*(1-vaccov)
    I[t+1]=I[t]+beta*S[t]*I[t]+beta*(1-vacdegree)*V[t]*I[t]-rec_rate*I[t]-drate*I[t]
    R[t+1]=R[t]+rec_rate*I[t]-drate*R[t]
    V[t+1]=V[t]+brate*(S[t]+I[t]+R[t]+V[t])*vaccov-beta*(1-vacdegree)*V[t]*I[t]-drate*V[t]
  }

  cbind(S, I, R, V)
}

#Plot graph
plot.model <- function(tmax=2000){

  plot(x=0:tmax,y=model.out[1:(tmax+1),"S"],type="l",ylim=c(0,1.1),xlab="time",ylab="")
  lines(x=0:tmax,y=model.out[1:(tmax+1),"I"],col="red")
  lines(x=0:tmax,y=model.out[1:(tmax+1),"R"],col="blue")
  lines(x=0:tmax,y=model.out[1:(tmax+1),"V"],col="blue", lty=2)
  legend("top",legend=c("S","I","R","V"), col=c("black","red","blue","blue"),
         lty=c(1,1,1,2),bty="n",ncol=4)

}

#Example output
model.out <- model(vaccov=0,tmax=10000, brate=0.0001,drate=0.0001)
plot.model(tmax=100)
plot.model(tmax=10000)