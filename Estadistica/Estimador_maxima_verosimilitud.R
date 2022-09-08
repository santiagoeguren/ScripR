################################################################################
#Cargar librerias
################################################################################

# para la función mle
library(stats4) 
# para la función mle2
library(bbmle) 
#EstimationTools
library(EstimationTools)


################################################################################
#Fuentes
################################################################################

##https://rpubs.com/SRifo/415830
#https://www.r-bloggers.com/2019/08/maximum-likelihood-estimation-from-scratch/


################################################################################
#Función mle del paquete stats4
################################################################################

#-------------------------------------------------------------------------------
#Distribucion normal
#-------------------------------------------------------------------------------


x = rnorm(n = 100,mean = 98,sd = 10)


#"""
#Ahora construimos el negativo de la función de log-verosimilitud.
#"""


NegLogLik = function(mu,sigma){-sum(dnorm(x,mu,sigma,log = TRUE))}


#"""
#Luego usamos la función mle del paquete stats4 para realizar 
#la estimación por máxima verosimilitud:
#"""


EMV1 = mle(NegLogLik, start = list(mu=1, sigma=1))


EMV1

summary(EMV1)

confint(EMV1)



################################################################################
#Función mle2 del paquete bbmle
################################################################################

#1 Forma

EMV3 = mle2(NegLogLik,start = list(mu=97.41412,sigma=10.28713 ), data = list(x))

summary(EMV3)

#2 Forma
x1=rnorm(1000, mean=10, sd=5)

EMV4 = mle2(x1~dnorm(mu,sigma),start = list(mu=0,sigma=1), data=data.frame(x1))

summary(EMV4)




################################################################################
#EstimationTools
################################################################################



z=rnorm(n = 1000, mean = 10, sd = 2)
fit1 = maxlogL(x = z, dist = "dnorm", start=c(2, 3),
                    lower=c(-15, 0), upper=c(15, 10))
summary(fit1)


################################################################################
#Optimun
################################################################################



NLL = function(pars, data) {
  # Extract parameters from the vector
  mu = pars[1]
  sigma = pars[2]
  # Calculate Negative Log-LIkelihood
  -sum(dnorm(x = data, mean = mu, sd = sigma, log = TRUE))
}







#sample = rnorm(100,mean = 10, sd=2)
#prod(dnorm(sample))
## [1] 2.23626e-58


mle = optim(par = c(mu =1 , sigma = 1), fn = NLL, data = x,
            control = list(parscale = c(mu = 1, sigma = 1)))



mle$par
mle$value
mle$counts
mle$convergence
mle$message

##          mu       sigma 
## -0.07332745  0.90086176















