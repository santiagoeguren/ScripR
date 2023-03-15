
################################################################################
#Cargar librerias
################################################################################


library(nvmix)
library(qrmtools)
library(rugarch)


################################################################################
#Descargar datos
################################################################################


set.seed(123)
data("SP500_const", package = "qrmdata")

time=c("2010-01-01", "2012-12-31")


#Arma un data frame con 15 acciones

#"""
#Sera precios ajustado
#"""

x=SP500_const[paste0(time, collapse = "/"),   SP500_const_info$Subsector == "REITs"]

#Selecciono 2

x=x[,c(1,4)]

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
            inverse.burr = function(u, nu) (u^(-1/nu[2]) - 1)^(-1/nu[1]),
            pareto = "pareto")

#Limites de parametros para estimar

m.p.b_ = list(constant = c(0, 1e8), 
              inverse.gamma = c(1, 8),
              inverse.burr = matrix(c(0.1, 0.1, 8, 8), ncol = 2),
              pareto = c(1, 8))


#Estimacion de parametros

#Ejecutar este comando se demora

fit.results = lapply(1:4, function(i) fitnvmix(X, qmix = qmix_[[i]], mix.param.bounds = m.p.b_[[i]]))

#Da los parametros estimados para cada accion

fit.results



qq.results=lapply(1:4, function(i) qqplot_maha(fitnvmix_object = fit.results[[i]]))
  

#Genera alphas

alpha = 1 - 1/10^seq(0.5, 3.5, by = 0.05)
alpha

#Genera  Matriz vacia

VaRs=matrix(NA, ncol = 5, nrow = length(alpha))
VaRs

ESs=matrix(NA, ncol = 5, nrow = length(alpha))
ESs



a=matrix(c(1,1), ncol = 1)
a



for (i in 1:4) {
  
    VaRs[, i]=VaR_nvmix(alpha, qmix = qmix_[[i]],
                        loc = as.numeric(t(a)%*%fit.results[[i]]$loc),
                        scale = as.numeric(t(a)%*%fit.results[[i]]$scale%*%a),
                        nu = fit.results[[i]]$nu)
    
    
    ESs[, i] = ES_nvmix(alpha, qmix = qmix_[[i]],
                             loc = as.numeric(t(a)%*%fit.results[[i]]$loc),
                             scale = as.numeric(t(a)%*%fit.results[[i]]$scale%*%a),
                             nu = fit.results[[i]]$nu)
}



help(VaR_nvmix)


VaRs
ESs


#Estimacion no parametrica



#sum.obs = rowSums(X)
sum.obs=X%*%a

VaRs[, 5] = quantile(sum.obs, probs = alpha)
ESs[, 5] = sapply(seq_along(alpha),  function(i) mean(sum.obs[sum.obs > VaRs[i, 5]]))

ESs
                       
