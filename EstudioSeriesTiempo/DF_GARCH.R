################################################################################
#Cargar librerias
################################################################################

#Descargar datos
library(quantmod)  
#Graficos
library(ggplot2) 
#Complemeto graficos
require(tidyr)
#Test de Dickey-Fuller
library(tseries)
#Zivot & Andrews Unit Root Test
library(urca) 
# Estimar correlograma  
library(astsa)
#Extraer stock sp500
library(rvest)
#Garch
require(fGarch)
#forecast
library(forecast)
#Complemento grafici
library(plotly)


################################################################################
#Crear modelo
################################################################################


#----------------------------------------------------------------------------------------
#ARIMA(p,d,q)
#----------------------------------------------------------------------------------------

#AR

phi_1=0.5
phi_2=0

#MA
tita_1=0.2
tita_2=0.2

#mu
mu=10

#pendiente
beta=0.5


#Garch(p,q)
#Omega

alpa_0=2

#Z

alpa_1=0.5

#Sigma^2

beta_1=0.5


#---------------------------------------------------------------------------------------
#Crear serie de tiempo


z=NULL
x=NULL
sd=NULL

z[1]=rnorm(1,0,1)+beta*1+mu
z[2]=rnorm(1,0,1)+beta*2+mu

x[1]=z[1]
x[2]=z[2]

sd[1]=alpa_0
sd[2]=alpa_0




t=3


while (t<=500) {
  
  
  
  #GARCH
  
  sd[t]=sqrt(alpa_0+alpa_1*z[t-1]^2+beta_1*sd[t-1]^2)
  
  #ARMA
  
  z[t]=rnorm(1,0,sd[t])
  
  
  
  #Xt
  
  x[t]=mu+beta*t+phi_1*x[t-1]+phi_2*x[t-2]+tita_1*z[t-1]+tita_2*z[t-2]+z[t]
  
  t=t+1
  
}




#----------------------------------------------------------------------------------------
#Graficar Normal

plot(x,type="l",col="blue")



#----------------------------------------------------------------------------------------
#Diferenciar

dx_t=diff(x,lag=1)

plot(dx_t,type="l",col="blue")



model=garchFit( ~ arma(0, 0) + garch(1, 1),data =dx_t, trace = F)
summary(model)



help("garchFit")


model@residuals


w_t=cumsum((model@residuals/model@sigma.t)+  1.04202   )



plot(w_t,type="l",col="blue")


df=ur.df(x,type="trend",lags=3)
summary(df)



df=ur.df(w_t,type="trend",lags=3)
summary(df)



























