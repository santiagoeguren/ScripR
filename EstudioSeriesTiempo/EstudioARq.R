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

fi_1=0.9
fi_2=0
fi_3=0
fi_4=0

sgma=1

c=2




x_t=NULL

mu=c/(1-fi_1-fi_2-fi_3-fi_4)


z_t[1]=rnorm(1,0,sgma)
z_t[2]=rnorm(1,0,sgma)
z_t[3]=rnorm(1,0,sgma)
z_t[4]=rnorm(1,0,sgma)


x_t[1]=mu+rnorm(1,0,sgma)
x_t[2]=mu+rnorm(1,0,sgma)
x_t[3]=mu+rnorm(1,0,sgma)
x_t[4]=mu+rnorm(1,0,sgma)




i=5

t=c(1:10000)

while (i<=10000) {
  
  z_t[i]=rnorm(1,0,sgma)
  
  x_t[i]=c+z_t[i]+fi_1*x_t[i-1]+fi_2*x_t[i-2]+fi_3*x_t[i-3]+fi_4*x_t[i-4]
  
  i=i+1
}



#sgma*sgma*(1+tita_1*tita_1+tita_2*tita_2+tita_3*tita_3+tita_4*tita_4)

#Graficar

plot(t,x_t,xlab="t",ylab="x_t",type="l",col="blue")



###############################################################################
#Estimación de los parametros
##############################################################################


#------------------------------------------------------------------------------
#Estimación Maxima verosimilitud




fit= Arima(x_t, order=c(1,0,0), include.drift = F,include.mean = T,include.constant = T)
summary(fit)




#-----------------------------------------------------------------------------
#Estimación regresión lineal para tita_1

x_t_1=x_t[-1]
z_t_1=z_t[-length(z_t)]
summary(fit)
summary(lm(x_t_1~z_t_1))

plot(z_t_1,x_t_1,xlab="z_t_1",ylab="x_t_1",col="blue")
cor_x_z_1=cor(x_t_1,z_t_1)
cor_x_z_1
mean(x_t_1)
var(x_t_1)


#-----------------------------------------------------------------------------
#Estimación regresión lineal tita_2

x_t_2=x_t[c(-1,-2)]
z_t_2=z_t[-length(z_t)]
z_t_2=z_t_2[-length(z_t_2)]
summary(fit)
summary(lm(x_t_2~z_t_2))

plot(z_t_2,x_t_2,xlab="z_t_2",ylab="x_t_2",col="blue")
cor_x_z_2=cor(x_t_2,z_t_2)
cor_x_z_2

#-----------------------------------------------------------------------------
#Estimación regresión lineal para tita_1 y tita_2

z_t_1=z_t_1[-1]

summary(fit)
summary(lm(x_t_2~z_t_1+z_t_2))


###############################################################################
#ACF
##############################################################################


#0
x_t_0=x_t[-c(1:10)]
length(x_t_0)

#-1
x_t_1=x_t[-length(x_t)]
x_t_1=x_t_1[-c(1:9)]
length(x_t_1)

#-2
x_t_2=x_t[-(length(x_t):(length(x_t)-1))]
x_t_2=x_t_2[-c(1:8)]
length(x_t_2)

#-3
x_t_3=x_t[-(length(x_t):(length(x_t)-2))]
x_t_3=x_t_3[-c(1:7)]
length(x_t_3)

#-4
x_t_4=x_t[-(length(x_t):(length(x_t)-3))]
x_t_4=x_t_4[-c(1:6)]
length(x_t_4)


#-5
x_t_5=x_t[-(length(x_t):(length(x_t)-4))]
x_t_5=x_t_5[-c(1:5)]
length(x_t_5)


#-6
x_t_6=x_t[-(length(x_t):(length(x_t)-5))]
x_t_6=x_t_6[-c(1:4)]
length(x_t_6)




#-7
x_t_7=x_t[-(length(x_t):(length(x_t)-6))]
x_t_7=x_t_7[-c(1:3)]
length(x_t_7)



#-8
x_t_8=x_t[-(length(x_t):(length(x_t)-7))]
x_t_8=x_t_8[-c(1:2)]
length(x_t_8)




#-9
x_t_9=x_t[-(length(x_t):(length(x_t)-8))]
x_t_9=x_t_9[-c(1:1)]
length(x_t_9)


#-9
x_t_10=x_t[-(length(x_t):(length(x_t)-9))]
length(x_t_10)

#-------------------------------------------------------------
#Xt-1

acf2(x_t)
plot(x_t_1,x_t_0)
summary(lm(x_t_0~x_t_1))

cor(x_t_0,x_t_1)

0.889563*(sd(x_t_0)/sd(x_t_1))



#-------------------------------------------------------------
#Xt-2

acf2(x_t)
plot(x_t_2,x_t_0)
summary(lm(x_t_0~x_t_2))

cor(x_t_0,x_t_2)

-0.004785*(sd(x_t_0)/sd(x_t_2))


###############################################################################
#Estudio Predicción
##############################################################################

#--------------------------------------------------------------------------
#Prediccion en base a la libreria

acf2(x_t)

plot(forecast(fit, h=12),ylab="x_t", xlab="t",col="blue", xlim = c(9950,10012))
ARIMA_prediccion=forecast(fit, h=12)
ARIMA_prediccion






