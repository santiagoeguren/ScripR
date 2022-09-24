################################################################################
#Cargar librerias
###############################################################################

#CTS
library(SymTS)
#CTS
library(MixedTS)



################################################################################
#Generar Datos
################################################################################


dat=rCTS(1000, 1.5, c = 0.01, ell = 0.084, mu = 0)
hist(dat,probability = T)
curve(dCTS(x, 1.5, c = 0.01, ell = 0.084, mu = 0), -0.5, 0.5,add=TRUE,col="red")
curve(dnorm(x, 0,sd(dat)), add=TRUE,col="Blue")


mean(dat)

sd(dat)



################################################################################
#Ver consistencia es VAR[X]
################################################################################



dat=rCTS(2, 1.5, c = 0.01, ell = 0.084, mu = 0)

sdr_estimate=sd(dat)


i=2

while (i<=1000) {
  
  
  
  
  dat[i+1]=rCTS(1, 1.5, c = 0.01, ell = 0.084, mu = 0)
  
  
  sdr_estimate[i]=sd(dat)
  
  i=i+1
  
}



plot(sdr_estimate,type="l")





################################################################################
#Ver teorema del limite central
################################################################################



mu_hat=mean(rCTS(100, 1.5, c = 0.01, ell = 0.084, mu = 0))



i=2

while (i<=1000) {
  
  
  
  mu_hat[i]=mean(rCTS(100, 1.5, c = 0.01, ell = 0.084, mu = 0))
  
  
  
  
  i=i+1
  
}



hist(mu_hat,probability = T)
curve(dnorm(x, mean(mu_hat),sd(mu_hat)),add=TRUE,col="Blue")



shapiro.test(mu_hat)

qqnorm(mu_hat)




