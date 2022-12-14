
########################################################################
#Cargar librerias
#########################################################################

#Descargar datos
library(quantmod)  
#Graficos
library(ggplot2) 
#Test de Dickey-Fuller
library(tseries)
#Zivot & Andrews Unit Root Test
library(urca) 
#Estimar correlograma  
library(astsa)
#Extraer stock sp500
library(rvest)
#Ordenar datos funcios selec
library(dplyr)
#Funcion tiempo
library(lubridate)
#Ordenar datos para graficos
library(tidyr)
#Garch
require(fGarch)






#################################################################################
#Descargar datos
################################################################################





qmix=function(u, nu) (u^(-1/nu[2]) - 1)^(-1/nu[1])
set.seed(274)
nu=c(2, 5)
d=1
n=500
scale=cov2cor(tcrossprod(matrix(runif(d * d), ncol = d)))
scale
m.p.b=matrix(c(0.1, 0.1, 8, 8), ncol = 2)
x=rnvmix(n, qmix = qmix, nu = nu, scale = scale)
fit.burr=fitnvmix(x, qmix = qmix, mix.param.bounds = m.p.b)
fit.burr





#-------------------------------------------------------------------------------
#1 Accion

qmix=function(u, nu) (u^(-1/nu[2]) - 1)^(-1/nu[1])
set.seed(274)
nu=c(2, 5)
d=1
n=500
scale=cov2cor(tcrossprod(matrix(runif(d * d), ncol = d)))
scale
m.p.b=matrix(c(0.1, 0.1, 8, 8), ncol = 2)

x=rnvmix(n, qmix = qmix, nu = nu, scale = scale)

hist(x)

fit.burr=fitnvmix(x, qmix = qmix, mix.param.bounds = m.p.b)
summary(fit.burr)

