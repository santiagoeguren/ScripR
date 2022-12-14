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

fi_1=0.01
fi_2=0
fi_3=0
fi_4=0

sgma=0.001

c=0




x_t=NULL

mu=c/(1-fi_1-fi_2-fi_3-fi_4)


z_t=NULL
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
#Estimación para fi

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


#-10
x_t_10=x_t[-(length(x_t):(length(x_t)-9))]
length(x_t_10)







#-------------------------------------------------------------
#Xt-1

acf2(x_t)
plot(x_t_1,x_t_0)
summary(lm(x_t_0~x_t_1))

cor(x_t_0,x_t_1)

0.9787643*(sd(x_t_0)/sd(x_t_1))



#-------------------------------------------------------------
#Xt-2

acf2(x_t)
plot(x_t_2,x_t_0)
summary(lm(x_t_0~x_t_2))

cor(x_t_0,x_t_2)

9.582e-01*(sd(x_t_0)/sd(x_t_2))




summary(lm(x_t_0~x_t_1+x_t_2))


###############################################################################
#ACF
##############################################################################

#----------------------------------------------------------------------------
#Causality
#-------------------------------------------------------------------------------

acf2(x_t)

#-----------------------------------------------------------------------------
#Estimación regresión lineal para z_t_1

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
#Estimación regresión lineal z_t_2

x_t_2=x_t[c(-1,-2)]
z_t_2=z_t[-length(z_t)]
z_t_2=z_t_2[-length(z_t_2)]
summary(fit)
summary(lm(x_t_2~z_t_2))

plot(z_t_2,x_t_2,xlab="z_t_2",ylab="x_t_2",col="blue")
cor_x_z_2=cor(x_t_2,z_t_2)
cor_x_z_2

#-----------------------------------------------------------------------------
#Estimación regresión lineal para z_t_1 y z_t_2

z_t_1=z_t_1[-1]

summary(fit)
summary(lm(x_t_2~z_t_1+z_t_2))


#---------------------------------------------------------------
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


#-10
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




summary(lm(x_t_0~x_t_1+x_t_2))



###############################################################################
#Estudio Predicción
##############################################################################

#--------------------------------------------------------------------------
#Prediccion en base a la libreria

acf2(x_t)

plot(forecast(fit, h=40),ylab="x_t", xlab="t",col="blue", xlim = c(9950,10040))
ARIMA_prediccion=forecast(fit, h=40)
ARIMA_prediccion




#----------------------------------------------------------------------------------------
#Predicion en base al concepto de causalidad




#T+1
summary(fit)

i=1

x_T1=0

fi_estimado= 0.9787
constante=-6e-04*(1-fi_estimado)

i=0

while (i<=100) {
  

x_T1=x_T1+constante*(fi_estimado**i)+(fi_estimado**(i+1))*z_t[(length(z_t)-i)] 

i=i+1

}

x_T1




#T+2

z_t0=c(z_t,0)

x_T2=0

i=0

while (i<=100) {
  
  
  x_T2=x_T2+constante*(fi_estimado**i)+(fi_estimado**(i+1))*z_t0[(length(z_t0)-i)] 
  
  i=i+1
  
}

x_T2




#prediccion mediante estimación parametrica

qnorm(0.025, mean(x_t),sd(x_t))
qnorm(0.975, mean(x_t),sd(x_t))


#----------------------------------------------------------------------------------------
#Predicion real en base al concepto de invertibilidad





#---------------------------------------------------------------------------------
#Estimar a1,...,a10 en base al modelo lineal

ARIMA_prediccion


summary(lm(x_t_0~x_t_1+x_t_2+x_t_3+x_t_4+x_t_5+x_t_6+x_t_7+x_t_8+x_t_9+x_t_10))

summary(fit)

a_0= -1.281e-05*(1-  9.730e-01)


#Para x_T_1

x_T1= 9.730e-01 *x_t[length(x_t)]+ a_0
x_T1


#Para x_T_1

x_T2=0.8941024*x_T1+ a_0
x_T2


###############################################################################
#Estudio PACF
##############################################################################


#---------------------------------------------------------------------
#Estudio de multicolinealidad


x2=rnorm(10000,0,1)
e1=rnorm(10000,0,1)
x1=e1+1*x2

y=NULL

i=1

while (i<=10000) {
  
  y[i]=0+0.999*x1[i]+0.001*x2[i]+rnorm(1,0,1)
  
  i=i+1
  
}


#------------------------------------------------------------------------------
#Nota

# y = 0.5 *x1 + 0.5 *x1

#Si x1=e1+x2, reemplazamos: y=0.5*(e1+x2)+0.5*x2 = 0.5*e1+0.5*x2+0.5*x2
#y=0.5*e1+1x2 [1]

#Si x2=0.5*x1+e2, reemplazamos: y = 0.5*x1+0.5*(0.5*x1+e2)=0.5*x1+0.25*x1+0.5*e2
#y=0.75*x1+0.5*e2      [2]




#-------------------------------------------------------------------------------
#Para x1
#------------------------------------------------------------------------------

#Regresion equación [2]


#ACF 1
summary(lm(y~x1))



#-------------------------------------------------------------------------------
#Eliminar el efecto de x2 en x1





fit_x1_x2=lm(x1~x2)
summary(fit_x1_x2)


#Regresión equación [1]

#PACF
summary(lm(y~fit_x1_x2$residuals))

plot(fit_x1_x2$residuals,y)


plot(e1,y)
summary(lm(y~e1))



#-------------------------------------------------------------------------------
#Para x2
#-------------------------------------------------------------------------------

#Regresión equación [1]

#ACF 2
summary(lm(y~x2))



#-------------------------------------------------------------------------------
#Eliminar el efecto de x1 en x2


fit_x2_x1=lm(x2~x1)
summary(fit_x2_x1)

#Regresión equación [2]

#PACF
summary(lm(y~fit_x2_x1$residuals))

plot(x1,x2)


#-------------------------------------------------------------------------------
#Para x1 y x2

#PACF
summary(lm(y~x1+x2))

#-------------------------------------------------------------------------
#Estimación AFC


acf2(x_t)

#para X -1
summary(lm(x_t_0~x_t_1))

#para X -2

#Una forma
summary(lm(x_t_0~x_t_1+x_t_2))


fit_x1_x2=lm(x_t_1~x_t_2)
summary(lm(x_t_0~fit_x1_x2$residuals))


#para X -3
summary(lm(x_t_0~x_t_1+x_t_2+x_t_3))


