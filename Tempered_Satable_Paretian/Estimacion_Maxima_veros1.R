################################################################################
#Cargar librerias
################################################################################

# para la función mle
library(stats4) 
# para la función mle2
library(bbmle) 
#CTS
library(SymTS)


################################################################################
#Fuentes
################################################################################

##https://rpubs.com/SRifo/415830


################################################################################
#Función mle del paquete stats4
################################################################################


dat=rCTS(100, alpha=1.72, c = 0.01, ell = 0.084, mu = 10)
hist(dat,probability = T)
curve(dCTS(x, alpha=1.72, c = 0.01, ell = 0.084, mu = 0),add=TRUE,col="red")
curve(dnorm(x, 0,sd(dat)), add=TRUE,col="Blue")


NegLogLik = function(alpha,c,ell,mu){-sum(log(dCTS(x, alpha,c,ell,mu)))}




EMV1 = mle(NegLogLik, start = list(alpha=1.72, c = 0.01, ell = 0.084, mu = 1))


EMV1

summary(EMV1)

confint(EMV1)



################################################################################
#Función mle2 del paquete bbmle
################################################################################

#1 Forma

EMV3 = mle2(NegLogLik,start = list(mu=90,sigma=5 ), data = list(x))


























NLL = function(pars, data) {
  # Extract parameters from the vector
  alpha = pars[1]
  c = pars[2]
  ell = pars[3]
  mu = pars[4]
  # Calculate Negative Log-LIkelihood
  -sum(log(dCTS(x=data, alpha, c, mu)))
}




mle = optim(par = c(alpha=1.72, c = 0.01, ell = 0.084, mu =1), fn = NLL, data = dat,
            control = list(parscale = c(alpha=1.72, c = 0.01, ell = 0.084, mu = 1)))

help(optim)

mle$par
mle$value
mle$counts
mle$convergence
mle$message











