###############################################################################
#Librerias
##############################################################################

#no cargar ya que GEvStable no funciona:

#Descargar datos
library(quantmod)  
#Graficos
library(ggplot2) 
#forecast
library(forecast)
#ACP/PACF
library(tseries)
#Atsa
library(astsa)
#Zivot & Andrews Unit Root Test
library(urca) 
#Garch
require(fGarch)
#Stable Paretian
library(stabledist)
#Librerias para GEVStable
library(Rsolnp)
library(skewt)

#-----------------------------------------------------------------------------
#Establecer direccion
#-----------------------------------------------------------------------------

#setwd("/home/santiago/...")
#getwd()


###############################################################################
#Descargar series x[t]
##############################################################################


#------------------------------------------------------------------------------
#Cargar Fechas
#-------------------------------------------------------------------------------
#Primero año--->mes--->dia  

f_init='2009-01-01'
f_final='2022-01-01'

#-------------------------------------------------------------------------------
#Cargar simbolo
#-------------------------------------------------------------------------------

eq="META"	


#-------------------------------------------------------------------------------
#Descargar datos de yahoo
#-------------------------------------------------------------------------------
#Primero año--->mes--->dia

data_eq=new.env()

getSymbols(eq, src = 'yahoo', from = f_init, to=f_final,env = data_eq, 
           auto.assign = T, periodicity = "d")

suppressWarnings(try(for(i in ls(data_eq)) data_eq[[i]] = adjustOHLC(data_eq[[i]],
                                                                     use.Adjusted=T),silent = TRUE)) 




#-------------------------------------------------------------------------------
#Extraer Precio ajustado

y_t=as.numeric(data_eq[[eq]][,6])



#-------------------------------------------------------------------------------
#Controlar los NA

i=1

while(i<length(y_t)+1){
  if(is.na(y_t[i])==TRUE){
    y_t[i]=y_t[i-1]
  }
  i=i+1
}



#-------------------------------------------------------------------------------
#Graficar

plot(y_t,type="l",col="blue")


#-------------------------------------------------------------------------------
#Transformar en log

x_t=log(y_t)

plot(x_t,type="l",col="blue")

#------------------------------------------------------------------------------
#Generar serie dx_t
#-------------------------------------------------------------------------------

#Diferenciar
dx_t=diff(x_t)



################################################################################
#Dividir train vs test
################################################################################



test=1500
train=length(dx_t)-test

dx_t_test=dx_t[c((length(dx_t)-(test)+1):length(dx_t))]


prediction=NULL
dx_t_arbitrage=NULL

i=0


while (i<=(test-1)) {
  




dx_t_train=dx_t[c((1+i):(train+i))]









#-------------------------------------------------------------------------------
#Graficar
#plot(dx_t_train,type="l",col="blue")


#-------------------------------------------------------------------------------
#ACF y PACF
#acf2(dx_t)




###############################################################################
#Estimación ARIMA(p,d,q)-GARCH(p,q)
##############################################################################

#------------------------------------------------------------------------------
#Estimación dx_t

model=garchFit(~ arma(1,3) + garch(1, 1),data=dx_t_train, trace = F)
#summary(model)  



#-------------------------------------------------------------------------------
#Estudio Prediccion
#-------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#Prediccion para dx_t

prediccion=predict(model,n.ahead=2,plot=F,conf=.95)


prediction[1+i]=prediccion$meanForecast[1]






if(prediccion$meanForecast[1]>=0){
  
         dx_t_arbitrage[1+i]=dx_t_test[1+i]
  
  
}else{
          
        dx_t_arbitrage[1+i]=-1*dx_t_test[1+i]
     
        }


i=i+1

}



retorno_real_arbitrage=cumsum(dx_t_arbitrage)
retorno_real_buy=cumsum(dx_t_test)
t=c(1:length(retorno_real_arbitrage))

prediction


df=data.frame(retorno_real_arbitrage=retorno_real_arbitrage,retorno_real_buy=retorno_real_buy,
              t=t)

g=ggplot(df, aes(x=t))                  
g=g +  geom_line(aes(y=retorno_real_arbitrage), colour="red")
g=g +  geom_line(aes(y=retorno_real_buy), colour="green")
g



#test=1500
#f_init='2009-01-01'
#f_final='2022-01-01'

tick=c("UNH","V","CTAS","CRM","AAPL","MCD","XOM","BRK-A","MSFT",META	)
gano=c(1,1,1,0,1,1,1,1,0,1)






