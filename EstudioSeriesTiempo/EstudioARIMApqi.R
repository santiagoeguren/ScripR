################################################################################
#Cargar librerias
################################################################################

#Graficos
library(ggplot2) 
#ACP/PACF
library(tseries)
#forecast
library(forecast)
#Atsa
library(astsa)
#Transormacion Yeo-Jhonson
#library(bestNormalize)
#Libreria test B-P
library(lmtest)
#Test de Dickey-Fuller-Zivot & Andrews Unit Root Test
library(urca)



###############################################################################
#Armar serie x[t]
##############################################################################



alpa_ordenada=10

delta=0.01


psi=0
pi_1=0
pi_2=0





z_t=rnorm(1,0,1)
z_t[2]=rnorm(1,0,1)
x_t=alpa_ordenada
x_t[2]=alpa_ordenada



t=c(0:9999)

i=3

while (i<=10000) {
  
  z_t[i]=rnorm(1,0,1)
  
  x_t[i]=alpa_ordenada+delta*t[i]+pi_1*x_t[i-1]+pi_2*x_t[i-2]+psi*z_t[i-1]+z_t[i]
  
  i=i+1
  
}




###############################################################################
#Notas
#Si pi=0 y psi=0 dXt=delta+Zt-1Zt --->ARIMA(0,1,1) --->tita1 = -1
#Si pi=0 y psi=0.5 dXt=delta+Zt-1Zt --->ARIMA(0,1,2)--->tita1=(-0.5) y tita2=-0.5
#Si pi=0.5 y psi=0 dXt=delta + 0.5dXt-1 + Zt - Zt-1 --->ARIMA(1,1,1) ---> tita= -1 y fi=0.5
# si delta=0.01 ---> mu = 0.01/(1-0.5)=0.02
#Si pi=0.5 y psi=0.5 dXt=delta + 0.5dXt-1 + Zt - 0.5Zt-1 - 0.5Zt-2 --->ARIMA(1,1,2) --->
#tita1= -0.5,  tita2= -0.5 y fi=0.5. Si delta=0.01 ---> mu = 0.01/(1-0.5)=0.0
###############################################################################







#Graficar

plot(x_t,type="l",col="blue")


#Diferenciar

dx_t=diff(x_t)

plot(dx_t,type="l",col="blue")

acf2(dx_t)

mean(dx_t)
var(dx_t)

###############################################################################
#Estimación de los parametros
##############################################################################



#------------------------------------------------------------------------------
#Estimación Maxima verosimilitud


#Estimación x_t

fit= Arima(x_t, order=c(0,1,1), include.drift =F,include.mean = T,include.constant = T)
summary(fit)

help(Arima)


#Estimación dx_t

fit_d= Arima(dx_t, order=c(0,0,0), include.drift = F,include.mean = T,include.constant = T)
summary(fit_d)


###############################################################################
#Estudio Predicción
##############################################################################

#--------------------------------------------------------------------------------


plot(forecast(fit, h=500),ylab="x_t", xlab="t",col="blue", xlim = c(9500,10500))
ARIMA_prediccion=forecast(fit, h=40)
ARIMA_prediccion

acf2(ARIMA_prediccion$residuals)
#--------------------------------------------------------------------------------


plot(forecast(fit_d, h=1000),ylab="dx_t", xlab="t",col="blue", xlim = c(9990,10040))
ARIMA_prediccion_d=forecast(fit_d, h=40)
ARIMA_prediccion_d

acf2(ARIMA_prediccion_d$residuals)
mean(ARIMA_prediccion_d$residuals)
var(ARIMA_prediccion_d$residuals)




##############################################################################
#Test de Dickey Fuller
##############################################################################

df=ur.df(x_t,type="trend",lags=1)
summary(df)











