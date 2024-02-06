## Block 2.1 Estimating FOI by Fitting Catalytic Model to IgG Serological Data
## Fitting with OLS


# Read in data
df <- read.csv(file="rubella.sero.csv")

#Look over a range of lambda values and extract the sum of square residuals for each lambda
lambda.range <- seq(0.01, 0.25, 0.01)
ss.residuals <- rep(NA, length(lambda.range)) #storate for sum of square residuals
for (i in 1:length(lambda.range)){
  ages <- df$Age
  prop.seropos.obs <- df$Prop.Pos
  prop.seropos.est <- 1-exp(-(lambda.range[i]*ages))
  ss.residuals[i] <- sum((prop.seropos.obs-prop.seropos.est)^2)
}
plot(lambda.range, ss.residuals, ylab="sum of square residuals", xlab="lambda")


#Lets look at model fit of a few different lambda
points(lambda.range[3], ss.residuals[3], col="blue", pch=16, cex=1.5)
points(lambda.range[10], ss.residuals[10], col="aquamarine2", pch=16, cex=1.5)
points(lambda.range[which(ss.residuals==min(ss.residuals))], 
       ss.residuals[which(out==min(ss.residuals))], col="red", pch=16, cex=1.5)

plot(df$Age, df$Prop.Pos, cex=0.025*df$Tot, 
     pch=16, xlab="age (year)", xlim=c(1,15), ylim=c(0,1), ylab="proportion seropositive")
lines(df$Age, 1-exp(-(lambda.range[3]*df$Age)), col="blue", lwd=2)
lines(df$Age, 1-exp(-(lambda.range[10]*df$Age)), col="aquamarine2", lwd=2)
lines(df$Age, 1-exp(-(lambda.range[which(ss.residuals==min(ss.residuals))]*df$Age)), col="red", lwd=2)


# FYI - you can easily fit catalytic model using glm  
fit=glm(cbind(Tot-Pos,Pos)~-1+Age, family=binomial(link="log"), data=df)
summary(fit) #lambda is 0.064363
residuals(fit)
