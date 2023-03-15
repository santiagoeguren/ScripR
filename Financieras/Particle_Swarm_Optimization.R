

################################################################################
#Cargar librerias
################################################################################


library(pso)
library(qrmtools)
library(rugarch)
library(nvmix)

################################################################################
#Bibliografia
################################################################################

#https://rpubs.com/argaadya/intro-pso



################################################################################
#Ecuación Cuadratica
################################################################################





fitness = function(x){
  
  x=x[1]
  
  fitness = x**2
  
    return(fitness)
  
}



pso_finance = psoptim(par = rep(NA,1), fn = function(x){fitness(x)}, 
                      lower = -10, upper = 10, 
                      control = list(maxit = 10000, s = 100))




pso_finance





################################################################################
#Regresión lineal
################################################################################


#-------------------------------------------------------------------------------
#Crear datos

x_m=c(0:100)

y=NULL

for (i in 0:101){
  
   y[i]=10+1*x_m[i]+rnorm(1,0,1)
  
}






#B0=10
#B1=1

#sum((y-B0-B1*x_m)**2)



fitness = function(x){
  
  B0=numeric()
  B1=numeric()
  
  B0=x[1]
  B1=x[2]
  
  fitness = sum((y-B0-B1*x_m)**2)
  
  return(fitness)
  
}



pso_finance = psoptim(par = rep(NA,2), fn = function(x){fitness(x)}, 
                      lower = c(-10,-10), upper = c(10,10), 
                      control = list(maxit = 10000, s = 100))


pso_finance


################################################################################
#Portfolio
################################################################################



data("SP500_const", package = "qrmdata")

time=c("2010-01-01", "2012-12-31")


#Arma un data frame con 15 acciones

#"""
#Sera precios ajustado
#"""

x=SP500_const[paste0(time, collapse = "/"),   SP500_const_info$Subsector == "REITs"]

#Selecciono 2

x=x[,c(1,2,3)]

X = -returns(x)


#Parece que definis el modela ARIMA - GARCH, para cada accion

help(ugarchspec)

#std: s-studen
#“norm”: normal

uspec=rep(list(ugarchspec(distribution.model = "std")), ncol(X))

uspec

#Correr el modelo

fit.ARMA.GARCH=fit_ARMA_GARCH(X, ugarchspec.list = uspec, verbose = FALSE)

#Modelo ARIMA-GARCH

fit.ARMA.GARCH

fits=fit.ARMA.GARCH$fit 

summary(fits)

#Extraer los residuos


#Dudo si lo divide por el sigma para cada periodo
resi=lapply(fits, residuals, standardize = TRUE)

resi[[1]]

sd(resi[[1]])


X=as.matrix(do.call(merge, resi))

#Cambia el nombre

colnames(X)=colnames(x)
n = nrow(X)


#Hacer una lista con los modelos mixed


qmix_ =list(constant = "constant",
            inverse.gamma = "inverse.gamma",
            #inverse.burr = function(u, nu) (u^(-1/nu[2]) - 1)^(-1/nu[1]),
            pareto = "pareto")

#Limites de parametros para estimar

m.p.b_ = list(constant = c(0, 1e8), 
              inverse.gamma = c(1, 8),
              #inverse.burr = matrix(c(0.1, 0.1, 8, 8), ncol = 2),
              pareto = c(1, 8))


#Estimacion de parametros

#Ejecutar este comando se demora

fit.results = lapply(1:3, function(i) fitnvmix(X, qmix = qmix_[[i]], mix.param.bounds = m.p.b_[[i]]))

#Da los parametros estimados para cada accion

fit.results


################################################################################
#Particle Swarm Optimization
################################################################################

#Los sacas de la funcion anterior:"Fijate bien donde"

mean_stock=c(-0.001154   , -0.000893, -0.000851 )
nyse_cov=matrix(c(0.5080,0.4093,0.2210,0.4093,0.5122,0.2300,0.2210,0.2300,0.4913), ncol = 3)
nyse_cov
rf=0.04/10000000000000000
#rf=0


fitness = function(x){

# Assign weight for each stocks
weight_stock=numeric()

#Dos cantidad de acciones

for (i in 1:3) {
  
  weight_stock[i] = x[i]

}



# Calculate the numerator

f1=numeric()


for (i in 1:3) {
  
  f1[i] = weight_stock[i]*mean_stock[i]

  }



# Calculate the denominator
f2=numeric()

for (i in 1:3) {
  
  f3=numeric()
  
  for (j in 1:3) {
    
    f3[j]= weight_stock[i]*weight_stock[j]*nyse_cov[i,j]
  }
  
  f2[i] = sum(f3)
}


# Calculate Fitness Value
fitness = (sum(f1)-rf)/sum(f2)

# Penalize Constraint Violation ??

#fitness = fitness - 1e9 * (round(sum(weight_stock),10)-1)^2


return(fitness)

}


#Penalize

weight_stock=c(0.5,0.5)

round(sum(weight_stock),10)-1


- 1e9 *(round(sum(weight_stock),10)-1)^2


#"""
#Por definicion pso es min, le poner negativo para que sea max
#"""


#lower=Como minimo puede ser 0


pso_finance = psoptim(par = rep(NA,3), fn = function(x){-fitness(x)}, 
                       lower = rep(0.01,3), upper = rep(0.5,3), 
                       control = list(maxit = 10000, s = 100))


help(psoptim)

pso_finance

#Par pesos

#Value: Sharpe

