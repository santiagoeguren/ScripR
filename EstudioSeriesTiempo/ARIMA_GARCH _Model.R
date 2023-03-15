#-----------------------------------------------------------------------------------------
#Cargar librerias
#----------------------------------------------------------------------------------------

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


#---------------------------------------------------------------------------------------
#Crear modelo
#-------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------
#ARIMA(2,1,2)
#----------------------------------------------------------------------------------------

#AR

phi_1=0
phi_2=0

#MA
tita_1=0
tita_2=0
  
#mu
mu=0

#pendiente
beta=0

  
#Garch(p,q)
#Omega

alpa_0=0.5

#Z

alpa_1=0.3

#Sigma^2

beta_1=0.3



#---------------------------------------------------------------------------------------
#training vs test
training=700
test=300

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


while (t<=training+test) {
  
  
  
  #GARCH
  
  sd[t]=sqrt(alpa_0+alpa_1*z[t-1]^2+beta_1*sd[t-1]^2)
  
  #ARMA
  
  z[t]=rnorm(1,0,sd[t])
  
 
  
  #Xt
  
  x[t]=mu+beta*t+phi_1*x[t-1]+phi_2*x[t-2]+tita_1*z[t-1]+tita_2*z[t-2]+z[t]
  
  t=t+1
  
}

#Separar training vs test

x_training=x[c(1:training)]
x_test=c(rep(NA,training),x[c((training+1):(training+test))])






#----------------------------------------------------------------------------------------
#Graficar Normal

plot(x_training,type="l",col="blue")


#--------------------------------------------------------------------------------
#Graficar ggplot


t_training=c(1:training)


g=ggplot(data.frame(x_training=x_training,t_training=t_training), aes(x = t_training, y =x_training))
g=g+geom_line(aes(y=x_training),color = "blue", size=0.8,alpha=0.3) 
g=g+ geom_point(alpha=0.5, size=0.5)
g=g+xlab("tiempo") + ylab("x[t]")
g=g+theme(plot.title = element_text(hjust = 0.5))
g=g+geom_hline(yintercept = 0,size=0.8,alpha=0.3)
g=g+theme_minimal()
g


#---------------------------------------------------------------------------------
#Correlacion
#ACF --->MA
#PACF--->AR


acf2(x_training)


#----------------------------------------------------------------------------------------
#Modelar ARIMA (p,d,q)
#---------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------
#Ajustar Auto

#Siempre con drift

fit=auto.arima(x_training,seasonal = F,ic = c("aicc"), stepwise=T, allowdrift = T,allowmean = T)
summary(fit)
mean(x_training)
#Residuos
res=residuals(fit)
tsdisplay(res)


#--------------------------------------------------------------------------------
#Grafico Normal
#Prediccion

plot(forecast(fit, h=test),ylab="x", xlab="t",type="o",lwd=1,pch=1,cex=1)






#---------------------------------------------------------------------------------
#Notas
#Por defecto siempre se permite una constante, para d=1 o 0. Si Alowwdrift es falsa
#la contante se permite si d=0


#------------------------------------------------------------------------------------
#Ajustar manual

#Si d=0 ---> include.drift = F

fit= Arima(x_training, order=c(2,1,2), include.drift = F,include.mean = T,include.constant = T)
summary(fit)
mean(x_training)
#Residuos
res=residuals(fit)
tsdisplay(res)


#---------------------------------------------------------------------------------
#Notas

#https://robjhyndman.com/hyndsight/arimaconstants/

#Si hay tendencia: ARIMA(0,0,0) y Drift=T . Saca la tendencia en intercepto y pendiente y modela el proceso
#Si hay tendencia: ARIMA(0,1,0) y Drift=T o F . Modela el proceso como ARIMA

#Si y solo d=0 include.mean=F Force mu=0. Si d>0 no tiene efecto
#Estimar refrecion lineal siemple
#summary(lm(x~c(1:200)))


#-----------------------------------------------------------------------------
#Grafico Normal
#Prediccion

plot(forecast(fit, h=test),ylab="x", xlab="t",type="o",lwd=1,pch=1,cex=1)




#---------------------------------------------------------------------------
#Graficar ggplot
#https://stackoverflow.com/questions/28775036/ggplot-line-graph-with-na-values


prediction=forecast(fit, h=test)

#Armar datos
x_prediction=c(rep(NA,training),prediction$mean)
ci_up=c(rep(NA,training),prediction$upper[,2])
ci_low=c(rep(NA,training),prediction$lower[,2])
x_training_g=c(x_training,rep(NA,test))
t_training_g=c(1:(training+test))



#Base
g=ggplot(data.frame(x_training_g=x_training_g,t_training_g=t_training_g,x_prediction=x_prediction,ci_up=ci_up,ci_low=ci_low,x_test=x_test), aes(x = t_training_g, y =x_training_g))
g=g+xlab("tiempo") + ylab("x[t]")
g=g+theme(plot.title = element_text(hjust = 0.5))
g=g+geom_hline(yintercept = 0,size=0.8,alpha=0.3)
g=g+theme_minimal()


#Serie de tiempo observada
g=g+geom_line(aes(y=x_training_g),color = "blue", size=0.8,alpha=0.3,na.rm=TRUE) 
g=g+ geom_point(aes(y=x_training_g),alpha=0.5, size=0.5,na.rm=TRUE)

#Serie de tiempo prediccion
g=g+geom_line(aes(y=x_prediction),color = "red", size=0.8,alpha=0.3,na.rm=TRUE) 
g=g+ geom_point(aes(y=x_prediction),alpha=0.5, size=0.5,na.rm=TRUE)

#intervalos
g=g+geom_ribbon(aes(ymin=ci_low, ymax=ci_up), fill="red", alpha=0.2)


#Valores test

g=g+geom_line(aes(y=x_test),color = "green", size=0.8,alpha=0.3,na.rm=TRUE) 
g=g+ geom_point(aes(y=x_test),alpha=0.5, size=0.5,na.rm=TRUE)
g

#Graficar plotly

p=ggplotly(g)
p

#----------------------------------------------------------------------------------------
#Modelar ARIMA-GARCH
#---------------------------------------------------------------------------------------


#Si es integrado

#x=diff(x)

#Generar modelo

model=garchFit( x_training~ arma(0, 0) + garch(1, 1),data = x_training, trace = F)
summary(model)             

#help(garchFit)


#-----------------------------------------------------------------------------
#Grafico normal
#Prediccion

predict(model,n.ahead=test,plot=TRUE,conf=.95)
      

sqrt(0.64868/(1-0.17477-0.31171))

#plot(model@sigma.t, type="l")

#---------------------------------------------------------------------------
#Graficar ggplot
#https://stackoverflow.com/questions/28775036/ggplot-line-graph-with-na-values


prediction=predict(model,n.ahead=test,plot=T,conf=.95)





#Armar datos
x_prediction=c(rep(NA,training),prediction$meanForecast)
ci_up=c(rep(NA,training),prediction$upperInterval)
ci_low=c(rep(NA,training),prediction$lowerInterval)
x_training_g=c(x_training,rep(NA,test))
t_training_g=c(1:(training+test))



#Base
g=ggplot(data.frame(x_training_g=x_training_g,t_training_g=t_training_g,x_prediction=x_prediction,ci_up=ci_up,ci_low=ci_low,x_test=x_test), aes(x = t_training_g, y =x_training_g))
g=g+xlab("tiempo") + ylab("x[t]")
g=g+theme(plot.title = element_text(hjust = 0.5))
g=g+geom_hline(yintercept = 0,size=0.8,alpha=0.3)
g=g+theme_minimal()


#Serie de tiempo observada
g=g+geom_line(aes(y=x_training_g),color = "blue", size=0.8,alpha=0.3,na.rm=TRUE) 
g=g+ geom_point(aes(y=x_training_g),alpha=0.5, size=0.5,na.rm=TRUE)

#Serie de tiempo prediccion
g=g+geom_line(aes(y=x_prediction),color = "red", size=0.8,alpha=0.3,na.rm=TRUE) 
g=g+ geom_point(aes(y=x_prediction),alpha=0.5, size=0.5,na.rm=TRUE)

#intervalos
g=g+geom_ribbon(aes(ymin=ci_low, ymax=ci_up), fill="black", alpha=0.2)


#Valores test

g=g+geom_line(aes(y=x_test),color = "green", size=0.8,alpha=0.3,na.rm=TRUE) 
g=g+ geom_point(aes(y=x_test),alpha=0.5, size=0.5,na.rm=TRUE)
g

#Graficar plotly

p=ggplotly(g)
#p




sd(x_training)
sd(model@residuals)

sd(model@residuals/model@sigma.t)
































diffinv(x, xi = 1)







# Create specification for GARCH(1, 1)
spec <- garchSpec(model = list(omega = 0.05, alpha =  0.1, beta = 0.75), cond.dist = "norm")

# Simulate the model with n = 1000
sim <- garchSim(spec, n = 1000)

# Fit a GARCH (1, 1)
fit <- garchFit(formula = ~ garch(1, 1), data = sim, include.mean = F)

# Predict 40 steps ahead
pred <- predict(fit, n.ahead = 40)

# Concatenate the fitted model with the prediction, transform to time series
dat <- as.ts(c(sqrt(fit@h.t), pred = pred$standardDeviation))

# Create the plot
plot(window(dat, start = start(dat), end = 1000), col = "blue",
     xlim = range(time(dat)), ylim = range(dat),
     ylab = "Conditional SD", main = "Prediction based on GARCH model")

par(new=TRUE)

plot(window(dat, start = 1000), col = "red", axes = F, xlab = "", ylab = "", xlim = range(time(dat)), ylim = range(dat))
