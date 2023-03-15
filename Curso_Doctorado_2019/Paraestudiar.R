


#-------------------------------------------------------------------------------------
#Crear modelo b?sico sin correlaci?n. normal e igual varianzas
#-------------------------------------------------------------------------------------

#Crear mu


mugeneral=500

mu1=-100
mu2=-100
mu3=200

#crear tau

tau1=mugeneral+mu1
tau2=mugeneral+mu2
tau3=mugeneral+mu3
  

#Crar des. est?ndar

sd1=20
sd2=20
sd3=20

#crear y1,y2,y3

y1=NULL
y2=NULL
y3=NULL



i=1

while(i<201){
  
  y1[i]=tau1+rnorm(1, mean=0, sd=sd1)
  y2[i]=tau2+rnorm(1, mean=0, sd=sd2)
  y3[i]=tau3+rnorm(1, mean=0, sd=sd3)
  
i=i+1  
}



#Crear ytotal y factor


y=c(y1,y2,y3)
x=c(rep(1,length(y1)),rep(2,length(y2)),rep(3,length(y3)))


#si no pones factor te va a tomar como si fuera un modelo lineal general

x=factor(x)



#Gr?ficar box-plot

boxplot(y~x,col="lightblue")

#si x es vector normal o factor gr?ficamente es igual


#-------------------------------------------------------------------------------------
#Modelo  sin media General
#-------------------------------------------------------------------------------------

modelo=lm(y~x-1)
summary(modelo)


#Anova

#H0: tau1=tau2=tau3    vs     H1: al menos uno es diferente

anova(modelo)

#-------------------------------------------------------------------------------------
#Modelo  con media General
#-------------------------------------------------------------------------------------


options(contrasts = c("contr.sum","contr.poly"))
#options(contrasts = c("contr.treatment","contr.poly"))  #por defecto 

modelo.media.general=lm(y~x)
summary(modelo.media.general)

#al mu general lo toma el del primero

anova(modelo.media.general)



#-------------------------------------------------------------------------------------
#tama?o de la muestra
#-------------------------------------------------------------------------------------

# se usa modelo con media general
#H0: ta1=tau2=tau3    vs   H1 la diferencia es 50





#ingresas datos
k=3                   #numero de factores
sigma2=20^2       #estimado del modelo con media general
Delta=50              #diferencia minima
piDelta=0.95          #potencia de la  prueba



#r0=round(100/k)
r0=1
alfa=0.05
a=(Delta^2)/(2*sigma2)
pot=0
r=r0


while(pot<piDelta){
  
  r=r+0.001
  f0=qf(1-alfa,k-1,k*(r-1))    # le resta los grados de liberta que son k
  delta2=a*r
  pot=(1-pf(f0,k-1,k*(r-1),delta2))
  
}

print(pot)

print(r)

ntotat=r*k
ntotat


#Gr?ficar bajo H0 verdadera

curve(df(x,k-1,k*(r-1)), col="blue", lwd=2,  yaxt="n",ylim=c(0,0.2), xlim=c(0,20))

#los grados de libertad se obtien del modelo con media general

#Gr?ficar bajo H1 verdadera

dd=((r*Delta^2)/(2*sigma2))^0.5

#el r es el mumero de poblaci?n de cada factor

curve(df(x,k-1,k*(r-1),delta2), col="red", lwd=2,add = T ,yaxt="n",ylim=c(0,0.2), xlim=c(0,20))

legend("topleft",legend=c("H0 Verdadera","H1 Verdadera"),col=c("blue","red"),pch=c(1,1))
       
abline(v=qf(0.95,k-1,k*(r-1)))
abline(v=qf(0.05,k-1,k*(r-1),delta2,lower.tail = T))














#-------------------------------------------------------------------------------------
#Contraste
#------------------------------------------------------------------------------------


boxplot(y~x,col="lightblue")

#uno solo (la sumatoria de los c es igual a cero)



contrastes=matrix(c(-1/2,1/4,1/4))

contrasts(x)=contrastes

ajuste=aov(y~x)

summary.aov(ajuste,split = list(x=list("Contraste 1"=1)))

# para dos


#el nrow es el numero de c, el ncol numero de contraste

contrastes=matrix(c(1/2,-1/4,-1/4,1/2,-1/2,0),nrow = 3,ncol = 2,byrow=F)

contrasts(x)=contrastes

ajuste=aov(y~x)

summary.aov(ajuste,split = list(x=list("Contraste 1 "=1,"Contraste 2"=2)))




#-------------------------------------------------------------------------------------
#Tukey
#------------------------------------------------------------------------------------

#modelo sin media general
library(multcomp)

modelo=lm(y~x-1)         #si pones el menos uno o no el resultado es igual
confint(modelo)

ajuste=aov(y~x-1)   #es igual que anova


compara=glht(ajuste,mcp(x="Tukey"))  #recorda el x es el nombre del factor

summary(compara)

plot(compara)



#-------------------------------------------------------------------------------------
#Dunnet
#------------------------------------------------------------------------------------

ajuste=aov(y~x-1)
confint(modelo)
ajuste=aov(y~x-1)

compara=glht(ajuste,mcp(x="Dunnett"))

summary(compara)

plot(compara)




#-------------------------------------------------------------------------------------
#Normalidad
#------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
#Reformular modelo con residuos no normales (F de fisher)
#------------------------------------------------------------------------------------



#Crear mu


mugeneral=100


mu1=0
mu2=0
mu3=0

#crear tau

tau1=mugeneral+mu1
tau2=mugeneral+mu2
tau3=mugeneral+mu3



#crear y1,y2,y3

y1=NULL
y2=NULL
y3=NULL



i=1

while(i<201){
  
  y1[i]=tau1+rf(1, 2, 30)
  y2[i]=tau2+rf(1, 2, 30)
  y3[i]=tau3+rf(1, 2, 30)
  
  i=i+1  
}



#Crear ytotal y factor


y=c(y1,y2,y3)
x=c(rep(1,length(y1)),rep(2,length(y2)),rep(3,length(y3)))

#si no pones factor te va a tomar como si fuera un modelo lineal general

x=factor(x)


#Gr?ficar box-plot

boxplot(y~x,col="lightblue")


#-------------------------------------------------------------------------------------
#Modelo  sin media General
#-------------------------------------------------------------------------------------

modelo=lm(y~x-1)
summary(modelo)




#-------------------------------------------------------------------------------------
#Modelo  con media General
#-------------------------------------------------------------------------------------


options(contrasts = c("contr.sum","contr.poly"))

modelo.media.general=lm(y~x)
summary(modelo.media.general)

#al mu general lo toma el del primero

anova(modelo.media.general)



#---------------------------------------------------------------------------------------
#Test de shapiro-Wilks
#-----------------------------------------------------------------------------------------

#H0: Es normal    vs    H1: no es normal




shapiro.test(rstandard(modelo))
shapiro.test(rstandard(modelo.media.general))






library(moments)


hist(modelo$residuals,freq = F)

kurtosis(modelo$residuals)
#leptoc?rtica, beta>3
#mesocurtica beta=3
#platicurtica beta<3

skewness(modelo$residuals) 

#0 sin sesgo
# + sesgo derecha
#- sesdp izquierda

#---------------------------------------------------------------------------------------
#Transformaci?n Box & Cox
#-----------------------------------------------------------------------------------------

library(MASS)
  
b=boxcox(modelo)
lambda.l <- b$x #lamda valores
lik <- b$y    #
bc <- cbind(lambda.l, lik)
sorted_bc <- bc[order(-lik),]
head(sorted_bc, n = 10)


lambda=sorted_bc[1,1]   #el lambda que maximiza la funci?n de verosimilitud
lambda



y.transformado=((y^lambda)-1)/lambda #si lambda distinto de cero
#y.transformado=log(y) #si lambda igual a cero


#Gr?ficar box-plot

boxplot(y.transformado~x,col="lightblue")



modelo.transformado=lm(y.transformado~x-1)
summary(modelo.transformado)



options(contrasts = c("contr.sum","contr.poly"))

modelo.transformado.media.general=lm(y.transformado~x)
summary(modelo.transformado.media.general)



#H0: Es normal    vs    H1: no es normal



shapiro.test(rstandard(modelo.transformado))
shapiro.test(rstandard(modelo.transformado.media.general))


#-------------------------------------------------------------------------------------
#Heroscedasticidad
#------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
#Reformular con una relaci?n lineal de las varianzas
#------------------------------------------------------------------------------------




#Crear mu


mugeneral=500

mu1=0
mu2=0
mu3=0

#crear tau

tau1=mugeneral+mu1
tau2=mugeneral+mu2
tau3=mugeneral+mu3

#Crar des. est?ndar

areal=0.2

sd1=10
sd2=10
sd3=10

#crear y1,y2,y3

y1=NULL
y2=NULL
y3=NULL



i=1

while(i<201){
  
  y1[i]=tau1+rnorm(1, mean=0, sd=sd1)
  y2[i]=tau2+rnorm(1, mean=0, sd=sd2)
  y3[i]=tau3+rnorm(1, mean=0, sd=sd3)
  
  i=i+1  
}



#Crear ytotal y factor


y=c(y1,y2,y3)
x=c(rep(1,length(y1)),rep(2,length(y2)),rep(3,length(y3)))

#si no pones factor te va a tomar como si fuera un modelo lineal general

x=factor(x)


#Gr?ficar box-plot

boxplot(y~x,col="lightblue")



#-------------------------------------------------------------------------------------
#Modelo  sin media General
#-------------------------------------------------------------------------------------

modelo=lm(y~x-1)
summary(modelo)

#-------------------------------------------------------------------------------------
#Modelo  con media General
#-------------------------------------------------------------------------------------


options(contrasts = c("contr.sum","contr.poly"))

modelo.media.general=lm(y~x)
summary(modelo.media.general)

#al mu general lo toma el del primero

anova(modelo.media.general)








#-----------------------------------------------------------------------------------
#test de bartlett
#----------------------------------------------------------------------------------
#H0: todas las varianzas son iguales vs H1: al menos una diferente


#1 forma- Residuos stunderizados. Mejor forma

bartlett.test(rstandard(modelo)~x-1)
bartlett.test(rstandard(modelo.media.general)~x)


#2 forma. solo si ni es igual para todos
bartlett.test(residuals(modelo)~x)

# 3 forma. solo para comparar los datos originales
bartlett.test(y~x)


#test Breusch-Pagan test
bptest(modelo,studentize = TRUE)


#-----------------------------------------------------------------------------------
#Minimo cuadrado ponderados
#----------------------------------------------------------------------------------

varianzas=tapply(y,x,var) 



pesos=c(rep(1/varianzas,each=length(y1)))



#sin media genera
modelo.p=lm(y~x-1,w=pesos)
summary(modelo.p)


#con media general


options(contrasts = c("contr.sum","contr.poly"))
modelo.p.media.general=lm(y~x,w=pesos)
summary(modelo.p.media.general)



# Test de Barlett


bartlett.test(rstandard(modelo.p)~x-1)
bartlett.test(rstandard(modelo.p.media.general)~x)



#-----------------------------------------------------------------------------------
#tranformaci?n que estabiliza la varianza
#----------------------------------------------------------------------------------


#Gr?ficar

vari=tapply(y,x,var)
medias=tapply(y,x,mean)


plot(medias,vari)


#calcular el a Real con modelo lineal

mod=lm(log(vari)~log(medias))
summary(mod)

q=coef(mod)[2]
q

a=exp(coef(mod)[1])
a

y.tr=y^(1-(q/2))


#modelo sin media general

modelo.tr=lm(y.tr~x-1)
summary(modelo.tr)


#modelo con media general



options(contrasts = c("contr.sum","contr.poly"))
modelo.tr.media.general=lm(y.tr~x)
summary(modelo.tr.media.general)




bartlett.test(rstandard(modelo.tr)~x-1)

bartlett.test(rstandard(modelo.tr.media.general)~x)



#-------------------------------------------------------------------------------------
#Autocorrelaci?n
#------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
#Reformular con autoccorelaci?n 
#------------------------------------------------------------------------------------


#Crear mu


mugeneral=500

mu1=0
mu2=0
mu3=0

#crear tau

tau1=mugeneral+mu1
tau2=mugeneral+mu2
tau3=mugeneral+mu3
#Crar des. est?ndar

sd1=20
sd2=20
sd3=20

#crear y1,y2,y3

y1=rnorm(1, mean=0, sd=sd1)+tau1
y2=rnorm(1, mean=0, sd=sd1)+tau2
y3=rnorm(1, mean=0, sd=sd1)+tau3



i=2

while(i<200){
  
  y1[i]=tau1+rnorm(1, mean=0, sd=sd1)+0.5*(y1[i-1]-tau1)
  y2[i]=tau2+rnorm(1, mean=0, sd=sd2)+0.5*(y2[i-1]-tau2)
  y3[i]=tau3+rnorm(1, mean=0, sd=sd3)+0.5*(y3[i-1]-tau3)
  
  i=i+1  
}



#Crear ytotal y factor


y=c(y1,y2,y3)
x=c(rep(1,length(y1)),rep(2,length(y2)),rep(3,length(y3)))

#si no pones factor te va a tomar como si fuera un modelo lineal general

x=factor(x)


#Gr?ficar box-plot

boxplot(y~x,col="lightblue")

#si x es vector normal o factor gr?ficamente es igual


#-------------------------------------------------------------------------------------
#Modelo  sin media General
#-------------------------------------------------------------------------------------

modelo=lm(y~x-1)
summary(modelo)




#-------------------------------------------------------------------------------------
#Modelo  con media General
#-------------------------------------------------------------------------------------



options(contrasts = c("contr.sum","contr.poly"))
modelo.media.general=lm(y~x)
summary(modelo.media.general)




#--------------------------------------------------------------------------------------
#Durbin-Wattson
#--------------------------------------------------------------------------------------
#H0: errores son independientes  vs H1: no son independientes


library(lmtest)

dwtest(modelo)



dwtest(modelo.media.general)

