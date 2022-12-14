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

f_init_test='1962-05-01'
f_final_test='2022-09-01'

#----------------------------------------------------------------------
#Cargar simbolo

eq_1="DOW"


#------------------------------------------------------------------
#Descargar datos de yahoo


data_eq_1=new.env()
suppressWarnings(try(getSymbols(eq_1, src = 'yahoo', from = f_init_test, 
                                to=f_final_test,env = data_eq_1, auto.assign = T,
                                periodicity = "m"),silent = T))
suppressWarnings(try(for(i in ls(data_eq_1)) data_eq_1[[i]] =
                       adjustOHLC(data_eq_1[[i]], use.Adjusted=T),silent = TRUE)) 




#-------------------------------------------------------------------
#Extraer Precio ajustado

y_t=as.numeric(data_eq_1[[eq_1]][,6])


y_t_12=y_t[1]


i=2


while (i <= floor(length(y_t)/12)) {
  
  y_t_12[i]=y_t[1+12*i]
  
  i=i+1
  
}




x_t=log(y_t_12)


dx_t=diff(x_t)

plot(dx_t, type="l")



shapiro.test(dx_t)


hist(dx_t)

################################################################################
#consistencia var directamente
################################################################################
sdt_dat=NULL

sdt_dat[1]=sd(dx_t[1:10])


i=1


while (i<=length(dx_t)) {
  
  sdt_dat[i]=sd(dx_t[1:(10+i)])
  
  i=i+1
  
  
}


plot(sdt_dat, type="l")


###############################################################################
#Estimación ARIMA(p,d,q)-GARCH(p,q)
##############################################################################

#------------------------------------------------------------------------------
#Estimación omega




omega_dat=NULL

model=garchFit( ~ arma(1,3) + garch(1, 1),data=dx_t[1:100], trace = F)
omega_dat[1]=model@fit$coef[6]


i=1


while (i<=length(dx_t)) {
  
  model=garchFit( ~ arma(1,3) + garch(1, 1),data=dx_t[1:(100+i)], trace = F)
  omega_dat[i]=model@fit$coef[6]
  
  i=i+1
  
  
}


plot(omega_dat, type="l")

omega_dat



#------------------------------------------------------------------------------
#Estimación var residuos




std_vector=NULL

model=garchFit( ~ arma(1,3) + garch(1, 1),data=dx_t[1:10], trace = F)
std_vector[1]=sd((model@residuals)/(model@sigma.t))




i=1


while (i<=length(dx_t)) {
  
  model=garchFit( ~ arma(1,3) + garch(1, 1),data=dx_t[1:(10+i)], trace = F)
  std_vector[i]=sd((model@residuals)/(model@sigma.t))
  
  i=i+1
  
  
}


summary(model)


plot(std_vector, type="l")



model@sigma.t




hist(model@residuals)


shapiro.test(model@residuals)




#https://www.idrisstsafack.com/post/garch-models-with-r-programming-a-practical-example-with-tesla-stock



