################################################################################
#Cargar librerias
################################################################################


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

library(plotly)

library(nvmix)

################################################################################
#Crear Modelo
################################################################################

getSymbols("CTAS",
           from = '2012-10-01',
           to = '2022-11-22',
           periodicity = "daily")


#-------------------------------------------------------------------------------
#Transformar

y_t=as.numeric(CTAS$CTAS.Adjusted)


x_t=log(y_t)


dx_t=diff(x_t,lag=1)


plot(dx_t,type="l")


#-------------------------------------------------------------------------------
#Modelar ARIMA-GARCH

model=garchFit(  ~ arma(1,3) + garch(1, 1),data=dx_t, trace = F)
summary(model)

acf2(dx_t)

#-------------------------------------------------------------------------------
#Estimacion



#-------------------------------------------------------------------------------
#Crear la matriz A


#Crear serie v_t

v_t=model@residuals/model@sigma.t
sd(v_t)

v_t=matrix(v_t, ncol=1)


#Modelar "inverse.gamma"

fit_ig=fitnvmix(v_t, qmix =  "inverse.gamma", mix.param.bounds = c(0.5, 10))
summary(fit_ig)

                            
qq.par=qqplot_maha(v_t, qmix = "inverse.gamma",loc = -0.004208, scale = 0.4758, df = 3.782, plot = T)
plot(qq.par)
plot(qq.par, plot.pars = list(log = "xy"))                                     

qq.par


#Modelar "Parero"

fit_p=fitnvmix(v_t, qmix =  "pareto", mix.param.bounds = c(0.5, 10))
summary(fit_p)



qq.par=qqplot_maha(v_t, qmix = "pareto",loc = -0.007864, scale = 0.3196, alpha = 1.367, plot = FALSE)
plot(qq.par)
plot(qq.par, plot.pars = list(log = "xy"))                                     

qq.par





#Modelar "Burr"


qmix=function(u, nu) (u^(-1/nu[2]) - 1)^(-1/nu[1])
#set.seed(274)
#nu=c(2, 5)
d=1
#n=500
scale=cov2cor(tcrossprod(matrix(runif(d * d), ncol = d)))
scale
m.p.b=matrix(c(0.1, 0.1, 8, 8), ncol = 2)
#x=rnvmix(n, qmix = qmix, nu = nu, scale = scale)
fit.burr=fitnvmix(v_t, qmix = qmix, mix.param.bounds = m.p.b)
fit.burr









#-------------------------------------------------------------------------------
#Graficar 

#Inverse gamma

d_ig=NULL
d_p=NULL
d_burr=NULL

x=seq(from=-8,to=8,by=0.01)

for (i in 1:length(x)){
  
  d_ig[i]=dnvmix(x[i], qmix = "inverse.gamma", loc = -0.004208, scale = 0.4758, df = 3.782)
  d_p[i]=dnvmix(x[i], qmix = "pareto", loc = -0.007864, scale = 0.3196, alpha = 1.367)
  d_burr[i]=dnvmix(x[i], qmix = qmix, nu = c(1.709,2.634), scale =0.2868)
  
}



dat=data.frame(v_t=v_t)
dat_ig=data.frame(d_ig=d_ig,x=x)
dat_p=data.frame(d_p=d_p,x=x)
dat_burr= data.frame(d_burr=d_burr,x=x)
 

g=ggplot(dat,aes(x = v_t))
g=g+geom_histogram(aes( y = ..density..), colour = 1, fill = "white",bins = 30)
#g=g+geom_line(data = dat_ig, aes(x = x, y = d_ig), color = "red")
#g=g+geom_line(data = dat_p, aes(x = x, y = d_p), color = "green")
g=g+geom_line(data = dat_burr, aes(x = x, y = d_burr), color = "blue")
#g=g+  xlim(-8, -2)
g                  
                   
#ggplotly(g)                   


