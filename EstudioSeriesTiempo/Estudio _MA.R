
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

tita_1=0.9
tita_2=0
tita_3=0
tita_4=0

sgma=1

c=2

z_t=NULL
x_t=NULL

z_t[1]=rnorm(1,0,sgma)
z_t[2]=rnorm(1,0,sgma)
z_t[3]=rnorm(1,0,sgma)
z_t[4]=rnorm(1,0,sgma)

x_t[1]=z_t[1]
x_t[2]=z_t[2]
x_t[3]=z_t[3]
x_t[4]=z_t[4]


i=5

t=c(1:10000)

while (i<=10000) {
  
  z_t[i]=rnorm(1,0,sgma)
  
  x_t[i]=c+z_t[i]+tita_1*z_t[i-1]+tita_2*z_t[i-2]+tita_3*z_t[i-3]+tita_4*z_t[i-4]
  
  i=i+1
}



sgma*sgma*(1+tita_1*tita_1+tita_2*tita_2+tita_3*tita_3+tita_4*tita_4)

#Graficar

plot(t,x_t,xlab="t",ylab="x_t",type="l",col="blue")


###############################################################################
#Estimación de los parametros
##############################################################################


#------------------------------------------------------------------------------
#Estimación Maxima verosimilitud




fit= Arima(x_t, order=c(0,0,1), include.drift = F,include.mean = T,include.constant = T)
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

-0.512181*(sd(x_t_0)/sd(x_t_1))



#-------------------------------------------------------------
#Xt-2

acf2(x_t)
plot(x_t_2,x_t_0)
summary(lm(x_t_0~x_t_2))

cor(x_t_0,x_t_2)

0.002482  *(sd(x_t_0)/sd(x_t_2))


###############################################################################
#Estudio Predicción
##############################################################################

#--------------------------------------------------------------------------
#Prediccion en base a la libreria

acf2(x_t)

plot(forecast(fit, h=12),ylab="x_t", xlab="t",col="blue", xlim = c(9950,10012))
ARIMA_prediccion=forecast(fit, h=12)
ARIMA_prediccion


#----------------------------------------------------------------------------------------
#Predicion real en base a z_t

#T+1
summary(fit)

-0.8966*z_t[length(z_t)] + 0*z_t[length(z_t)-1]+0*z_t[length(z_t)-2]+ 0 *z_t[length(z_t)-3]+ 1.9997




#T+2

-0.8966 *0 + 0*z_t[length(z_t)]+ 0*z_t[length(z_t)-1]+ 0 *z_t[length(z_t)-2] +  1.9997


#T+3

0.50637*0 + 0.4974*0+0.5034*z_t[length(z_t)]+ 0.5004 *z_t[length(z_t)-1] + 2.0079

#T+4

0.5063*0 + 0.4974*0+0.5034*0+0.5004 *z_t[length(z_t)] + 2.0079


#prediccion mediante estimación parametrica

qnorm(0.025, mean(x_t),sd(x_t))
qnorm(0.975, mean(x_t),sd(x_t))



#----------------------------------------------------------------------------------------
#Predicion real en base a X_t,.....





#---------------------------------------------------------------------------------
#Estimar a1,...,a10 en base al modelo lineal

ARIMA_prediccion


summary(lm(x_t_0~x_t_1+x_t_2+x_t_3+x_t_4+x_t_5+x_t_6+x_t_7+x_t_8+x_t_9+x_t_10))



#---------------------------------------------------------------------------------
#Estimar en forma teorica para xT+1

summary(fit)

tita_1_estimada= -0.8966 

sum_a=0

i=1

while(i<=200){
  
  sum_a=sum_a+((-1)**(i-1))*tita_1_estimada**i
  
  i=i+1
  
}



summary(fit)

mean_estimada=1.9997

a_0= mean_estimada*(1-sum_a)
a_0

valor_XT_1=0

i=1

while(i<=200){
  
  valor_XT_1=valor_XT_1+((-1)**(i-1))*(tita_1_estimada**i)*x_t[(length(x_t)-i+1)]
  
  i=i+1
  
}



valor_XT_1=valor_XT_1+a_0

valor_XT_1

ARIMA_prediccion





#---------------------------------------------------------------------------------
#Estimar en forma teorica para xT+2



valor_XT_2=0



x_t_XT_2=c(x_t,valor_XT_1)

i=1

while(i<=200){
  
  valor_XT_2=valor_XT_2+((-1)**(i-1))*(tita_1_estimada**i)*x_t_XT_2[(length(x_t_XT_2)-i+1)]
  
  i=i+1
  
}



valor_XT_2=valor_XT_2+a_0

valor_XT_2

ARIMA_prediccion






#-----------------------------------------------------------------------------
#Causalidad




x_t[length(z_t)]

z_t_cuasality=x_t[length(z_t)]

i=1

while (i<=200) {
  
  z_t_cuasality=z_t_cuasality+((-1)^i)*(tita_1_estimada^i)*x_t[(length(z_t)-i)]
  
  i=i+1
  
}


z_t_cuasality-a_0
z_t[length(z_t)]




###############################################################################
#Estudio PACF
##############################################################################


#---------------------------------------------------------------------
#Estudio de multicolinealidad

x1=rnorm(10000,0,1)
x2=rnorm(10000,0,1)+1*x1


y=NULL

i=1

while (i<=10000) {
  
  y[i]=10+0.5*x1[i]+0.5*x2[i]+rnorm(1,0,1)
  
  i=i+1
  
}

#---------------------------------------------
#beta=cov(y,x1)/var(x1)

cov(y,x1)/var(x1)


#------------------------------------------------------------------------------
#Para x1: el supuesto es que x1 en la independiente
#Mide el aporte real de x1
summary(lm(y~x1))


#Para x2


fit_x2_x1=lm(x2~x1)
summary(fit_x2_x1)

summary(lm(y~fit_x2_x1$residuals))






#Para x1 y x2
summary(lm(y~x1))
summary(lm(y~x1+x2))





#------------------------------------------------------------------------------
#Para x2: el supuesto es que x2 en la independiente
#Mide el aporte real de x2
summary(lm(y~x2))


#Para x2


fit_x2_x1=lm(x1~x2)
summary(fit_x2_x1)

summary(lm(y~fit_x2_x1$residuals))

#Para x2 y x1
summary(lm(y~x2))
summary(lm(y~x2+x1))










#-------------------------------------------------------------------------
#Estimación


acf2(x_t)

#para X -1
summary(lm(x_t_0~x_t_1))

#para X -2

#Una forma
summary(lm(x_t_0~x_t_1+x_t_2))

fit_x_t_2_x_t_1=lm(x_t_2~x_t_1)

summary(lm(x_t_0~fit_x_t_2_x_t_1$residuals))

#para X -3
summary(lm(x_t_0~x_t_1+x_t_2+x_t_3))

fit_x_t_2_x_t_1_x_t_3=lm(x_t_3~x_t_1+x_t_2)

summary(lm(x_t_0~fit_x_t_2_x_t_1_x_t_3$residuals))




