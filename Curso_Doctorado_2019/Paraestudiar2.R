

#-------------------------------------------------------------------------------------
#Crear modelo b?sico sin correlaci?n. normal e igual varianzas
#-------------------------------------------------------------------------------------



#Modelo Poblacional

#Yij=mu+alpai+Betaj+(alpa*beta)ij+Eij       i=1,...,a   j=1,...,b


#Modelo muestral

#Yijt=mu+alpai+Betaj+(alpa*beta)ij+Eijt       i=1,...,a   j=1,...,b t=1,...,a*b



#Crear mu


mugeneral=500

#Crear alfa 

alfa1=200
alfa2=-100
alfa3=-100

#Crear beta 2

beta1=300
beta2=-300


#itentacciones

alfa1beta1=50
alfa2beta1=-100
alfa3beta1=50



alfa1beta2=-50
alfa2beta2=100
alfa3beta2=-50




#Crar des. est?ndar

sd1=20
sd2=30




y11=NULL
y12=NULL


i=1

while(i<201){
  
  y11[i]=mugeneral+alfa1+beta1+alfa1beta1+rnorm(1, mean=0, sd=sd1)
  y12[i]=mugeneral+alfa1+beta2+alfa1beta2+rnorm(1, mean=0, sd=sd2)
 
  
  i=i+1  
}






y21=NULL
y22=NULL


i=1

while(i<201){
  
  y21[i]=mugeneral+alfa2+beta1+alfa2beta1+rnorm(1, mean=0, sd=sd1)
  y22[i]=mugeneral+alfa2+beta2+alfa2beta2+rnorm(1, mean=0, sd=sd2)

  
  i=i+1  
}






y31=NULL
y32=NULL


i=1

while(i<201){
  
  y31[i]=mugeneral+alfa3+beta1+alfa3beta1+rnorm(1, mean=0, sd=sd1)
  y32[i]=mugeneral+alfa3+beta2+alfa3beta2+rnorm(1, mean=0, sd=sd2)
  
  
  i=i+1  
}


y=c(y11,y12,y21,y22,y31,y32)


#representa alfa son 3
alfa=c(rep("1",length(c(y11,y12))),rep("2",length(c(y21,y22))),rep("3",length(c(y31,y32))))

#representa beta son 2

beta=c(rep("1",length(y11)),rep("2",length(y12)),rep("1",length(y21)),rep("2",length(y22)),rep("1",length(y31)),rep("2",length(y32)))


alfa=factor(alfa)
levels(alfa)
beta=factor(beta)
levels(beta)

#-------------------------------------------------------------------------------------
#Gr?ficar boxplot
#-------------------------------------------------------------------------------------



#--------------------------ggplot--------------------------------------------------



library(ggplot2)
library(plotly)

g = ggplot(data.frame(y=y,alfa=alfa,beta=beta), aes(y=y,x=alfa)) 
g= g + geom_boxplot(aes(fill=alfa)) 
g= g+facet_grid(.~beta)
g=g+geom_jitter(alpha = 0.5, color = "tomato")
g=g+ylab("y")
g

p <- ggplotly(g)
p


#-------------------------------------------------------------------------------------
#Gr?ficar interaction plot
#-------------------------------------------------------------------------------------



#--------------------------B?sico R--------------------------------------------------


interaction.plot(alfa,beta,y,col = topo.colors(3),lty = c(2,3),lwd=2.5,xlab = "alfa",ylab = "y",
                 las=1,legend = F,ylim = c(min(y),max(y)))
abline(h = mean(y))
legend("topright",col = topo.colors(3),lty=c(2,3),legend = c("beta1","beta2"),lwd = 2.5)




#-------------------------------------------------------------------------------------
#Modelo  con media General
#-------------------------------------------------------------------------------------


options(contrasts = c("contr.sum","contr.poly"))

modelo.media.general=lm(y~alfa+beta+alfa*beta)
summary(modelo.media.general)





#------------------------------An?lisis Anova para cada uno de los factores---------------------------------------


anova(modelo.media.general)

#para dos factores
#MSA/MSE=F(a-1,n-ab) para factor 1 con a niveles
#MSB/MSE=F(b-1,n-ab) para factor 2 con b niveles
#MSAB/MSE=F((a-1)*(b-1),n-ab) intereccion
#
#H0: los alpa son iguales Vs   H1:los alpha no son iguales
#H0: los Beta son iguales Vs   H1:los Beta no son iguales
#H0: (alfa*beta)ij-(alfa*beta)iq-(alfa*beta)sj+(alfa*beta)sq=0 i disnto de s   j distinto de q     vs  (alfa*beta)ij-(alfa*beta)iq-(alfa*beta)sj+(alfa*beta)sq ditnto de cero


#--------------------------------------Gr?ficar residuos ggplot---------------------------------------


g = ggplot(data.frame(residuos=rstandard(modelo.media.general),index=c(1:length(rstandard(modelo.media.general))),alfa=alfa,beta=beta), aes(y=residuos,x=index)) 
g=g+geom_point(aes(colour = alfa),alpha = 0.5)
g=g+geom_hline(yintercept = 0)
g= g+facet_grid(.~beta)
g
p <- ggplotly(g)
p
                 
            
                 
#---------------------------------------------------------------------------------------
#An?lisis de los supuestos
#-----------------------------------------------------------------------------------------




#---------------------------------------------------------------------------------------
#Test de shapiro-Wilks
#-----------------------------------------------------------------------------------------

#H0: Es normal    vs    H1: no es normal



shapiro.test(rstandard(modelo.media.general))



#-----------------------------------------------------------------------------------
#test Breusch-Pagan test
#----------------------------------------------------------------------------------

#H0: todas las varianzas son iguales vs H1: al menos una diferente


#1 forma- Residuos stunderizados. Mejor forma


library(lmtest)

bptest(modelo.media.general)

bptest(modelo.media.general,studentize = TRUE)




#-------------------------------------------------------------------------------------
#Minimos cuadrados generalizados
#------------------------------------------------------------------------------------


summary(modelo.media.general)

varianzas=tapply(modelo.media.general$residuals,beta,var)
varianzas


pesos=NULL

i=1

while(i<length(y)+1){
  
  if(beta[i]=="1"){
    
    pesos[i]=1/varianzas[1]
  }else{
    
    pesos[i]=1/varianzas[2]
    
  }
  
  i=i+1
}





modelo.p=lm(y~alfa+beta+alfa*beta, w=pesos)
summary(modelo.p)



#---------------------------------------------------------------------------------------
#Test de shapiro-Wilks
#-----------------------------------------------------------------------------------------


#H0: Es normal    vs    H1: no es normal



shapiro.test(rstandard(modelo.p))



#-----------------------------------------------------------------------------------
#test Breusch-Pagan test
#----------------------------------------------------------------------------------
#H0: todas las varianzas son iguales vs H1: al menos una diferente


#1 forma- Residuos stunderizados. Mejor forma


library(lmtest)

bptest(modelo.p,studentize = TRUE)




g = ggplot(data.frame(residuos=rstandard(modelo.p),index=c(1:length(rstandard(modelo.p))),alfa=alfa,beta=beta), aes(y=residuos,x=index)) 
g=g+geom_point(aes(colour = alfa),alpha = 0.5)
g=g+geom_hline(yintercept = 0)
g= g+facet_grid(.~beta)
g
p <- ggplotly(g)
p








#--------------------------------------------------------------------------------------
#Durbin-Wattson
#--------------------------------------------------------------------------------------
#H0: errores son independientes  vs H1: no son independientes


library(lmtest)

dwtest(modelo.media.general)






#--------------------------------------------------------------------------------------
#Predicci?n
#--------------------------------------------------------------------------------------


# Para los valores ya datos de 

fitted.values(modelo.media.general)

prediciones=data.frame(cbind(alfa,beta,fitted.values(modelo.media.general)))


# paraa un valor que yo quiero

nuevos=data.frame(alfa=factor("3"),beta=factor("2"))
predict(modelo.media.general, newdata =nuevos,interval = c("confidence"))















###########################################################################################3
#modelo lineal un solo factor
############################################################################################



#Modelo poblacional

#y/x=B0+B1*x1+E

#Modelo muestral

#y/xi=B0+B1*x1i+Ei   i=1,...,n

x=c(1:200)
y=NULL


sd1=20
sd2=20
beta0=50

i=1

while(i<length(x)/2){

y[i]=x[i]+rnorm(1, mean=0, sd=sd1)+beta0

i=i+1

}

while(i<100+(length(x)/2)+1){
  
  y[i]=x[i]+rnorm(1, mean=0, sd=sd2)+beta0
  
  i=i+1
  
}



#----------------------------------Gr?ficar lineal-----------------------------------------

library(ggplot2)


g=ggplot(data.frame(y=y,x=x),aes(x=x,y=y))
g=g+geom_smooth(method=lm)
g=g+geom_point(shape=1)
g


#-------------------------------------------------------------------------------------------------
#Modelo 
#------------------------------------------------------------------------------------------------



modelo=lm(y~x)
summary(modelo)



#-----------------------------------------criterios AIC-------------------------


step(modelo)




#---------------------------------------------------------------------------------------
#An?lisis de los supuestos
#-----------------------------------------------------------------------------------------



#---------------------------------------------------------------------------------------
#Test de shapiro-Wilks
#-----------------------------------------------------------------------------------------

#H0: Es normal    vs    H1: no es normal

shapiro.test(rstandard(modelo))
shapiro.test(modelo$residuals) #????




#---------------------------------------------------------------------------------------
#Test de Breusch-Pagan
#-----------------------------------------------------------------------------------------


library(lmtest)

bptest(modelo)



y.transformado=log(y) 




modelo.transformado=lm(y.transformado~x)
summary(modelo.transformado)






#---------------------------------------------------------------------------------------
#Test de shapiro-Wilks
#-----------------------------------------------------------------------------------------

#H0: Es normal    vs    H1: no es normal

shapiro.test(rstandard(modelo.transformado))
shapiro.test(modelo.transformado$residuals) #????




#---------------------------------------------------------------------------------------
#Test de Breusch-Pagan
#-----------------------------------------------------------------------------------------


library(lmtest)

bptest(modelo.transformado)




#-------------------------------Gr?ficar residuos-------------------------------------------


g = ggplot(data.frame(residuos=rstandard(modelo.transformado),index=c(1:length(rstandard(modelo.transformado)))),aes(y=residuos,x=index))
g=g+geom_hline(yintercept = 0)
g=g+geom_point(alpha = 0.5)
g
p <- ggplotly(g)
p




#--------------------------------------------------------------------------------------
#Durbin-Wattson
#--------------------------------------------------------------------------------------
#H0: errores son independientes  vs H1: no son independientes


library(lmtest)

dwtest(modelo)





  





#--------------------------------------------------------------------------------------
#Predicci?n
#--------------------------------------------------------------------------------------


# Para los valores ya datos de 

fitted.values(modelo)

prediciones=data.frame(cbind(x,fitted.values(modelo)))


# paraa un valor que yo quiero

nuevos=data.frame(x=c(3,4.5))
predict(modelo, newdata =nuevos,interval = c("confidence"))





























###########################################################################################3
#modelo covariable
############################################################################################

#Modelo poblacional

#Yi=mu+mui+B1*x1*i+B2*x2i+...+Ei i=1,..,k

#Modelo Muestral


#Yit=mu+mui+B1*x1*it+B2*x2it+...+Eit i=1,..,k  t=1,...,ni


x=c(1:200)
y=NULL


sd1=10
beta0=30


fac1=20
fac2=-20

#Crear el factor
fact=factor(sample(1:2,200,replace=T))




i=1

while(i<length(x)+1){
  
  y[i]=x[i]+rnorm(1, mean=0, sd=sd1)+beta0
  
  i=i+1
  
}

i=1

while(i<length(x)+1){
  
     if(fact[i]=="1"){
       
       y[i]=y[i]+fac1
     }else{
       
       y[i]=y[i]+fac2
     }
    i=i+1

  }





#-----------------------------------------------------------------------------------------
#Graficar 
#-----------------------------------------------------------------------------------------

library(ggplot2)

g = ggplot(data.frame(y=y,x=x,fact=fact), aes(y=y,x=x,colour=fact)) 
g=g+geom_smooth(method=lm)
g=g+geom_point(shape=1)
g=g+labs(color="SExo")
g=g+xlab("Notas")+ ylab("Rendimiento")
g
 










#-----------------------------------------------------------------------------------------
#Modelo 
#-----------------------------------------------------------------------------------------



options(contrasts = c("contr.sum","contr.poly"))
#options(contrasts = c("contr.treatment","contr.poly"))  #por defecto

modelo=lm(y~x+fact)
summary(modelo)




#Gr?fico multiple

par(mfrow=c(2,2))
plot(x,residuals(modelo),col="red",pch=19,
     main="con la covariable",xlab="x",ylab="residuos")
#Residuos con orden
plot(x,residuals(modelo),col="blue",pch=19,
     main="con el orden",xlab="orden",ylab="residuos")
#Residuos con valores predichos
plot(fitted.values(modelo),residuals(modelo),pch=19,
     main="con valores ajustados",xlab="valores ajustados",ylab="residuos")
#Grafico cuantil-cuantil
qqnorm(residuals(modelo),col="red",pch=19,main="Grafico cuantil-cuantil",
       xlab="cuantiles normales",ylab="cuantiles muestrales")
qqline(residuals(modelo),col="blue",lwd=2)






#para saber si los cofactores y los factores sin significativos en el modelo


anova(modelo)


#---------------------------------------------------------------------------------------
#Test de shapiro-Wilks
#-----------------------------------------------------------------------------------------


#H0: Es normal    vs    H1: no es normal



shapiro.test(rstandard(modelo))

#---------------------------------------------------------------------------------------
#Test de Breusch-Pagan
#-----------------------------------------------------------------------------------------


#H0:sigma^2_{1}=sigma^2_{2}=...=sigma^2_{1,5} vs sigma^2_{i} ditinto sigma^2_{j} i distinto a j 

library(lmtest)

bptest(modelo)



#--------------------------------------------------------------------------------------
#Durbin-Wattson
#--------------------------------------------------------------------------------------
#H0: errores son independientes  vs H1: no son independientes


library(lmtest)


dwtest(modelo)
















###########################################################################################3
#modelo no normales
############################################################################################

#-------------------------------------------------------------------------------------
# Para comparar medias 
#-------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------
#Poblaciones  dependientes
#-------------------------------------------------------------------------------------

#Modelo poblaci?n
#Y=(Y1,Y2)^t    E(Y)=(tau1,tau2)^t

#      [ var(Y1)           Cov(Y1,Y2)  ]
#   E= [                               ]
#      [ Cov(Y1,Y2)        Var(Y2)     ]

#------------------------------------------------------------------------------------
#Test de Wilcoxon 
#------------------------------------------------------------------------------------


#H0: y1-y2=d=0                    vs            #H1:y1-y2=d  >  0

guarde=c(82,69,73,43,58,56,76,65)
casa=c(63,42,74,37,51,43,80,62)

d=guarde-casa


wilcox.test(d,alternative = "greater",conf.int = T)



#-------------------------------------------------------------------------------------
#Poblaciones  independientes - no normales
#-------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------
#Comparar la dispercion - Test de Fligner y Killeen
#--------------------------------------------------------------------------------------


#H0: Var1 = Var2      vs       Var1 distinto de Var2


A=c(37,49,51,62,74,89,44,53,17)
B=c(40,39,62,62,73,40,38,33)


library(npsm)

fk.test(A,B)

#lo mismo podes mirar el 1



#--------------------------------------------------------------------------------------
#Comparar las mediana - Disperci?n igual -  Wilcoxon y Mann-Whitney
#--------------------------------------------------------------------------------------

#H0: diferencias de medianas =delta=0   vs H1: diferencias de medianas =delta distinto de 0

A=c(37,49,51,62,74,89,44,53,17)
B=c(40,39,62,62,73,40,38,33)


wilcox.test(A,B,exact = F,conf.int = T,alternative = c("two.sided"),conf.level = 0.95)







#--------------------------------------------------------------------------------------
#Comparar las mediana - Dispercion desigual - Fligner-Policello
#--------------------------------------------------------------------------------------

#H0: mediana1 = mediana2         vs          H1: mediana1 distinto de mediana2


y1=c(0.3,1,0.5,0.9,1.8,3.4,0.5,2.1,4.8,0.8)
y2=c(4.3,14.4,18.2,4.9,7.4,9.4,15.1,6.4,11.8,1)


library(npsm)


fp.test(y1,y2)  #las medianas son diferentes





#--------------------------------------------------------------------------------------
#Modelo de regresion lineal
#--------------------------------------------------------------------------------------

library(Rfit)




#Modelo poblacional

#y/x=B0+B1*x1+E

#Modelo muestral

#y/xi=B0+B1*x1i+Ei   i=1,...,n

#--------------------------------Generar Modelo con F-fisher------------------------------



x=c(1:200)
y=NULL



beta0=50

i=1

while(i<length(x)/2){
  
  y[i]=2*x[i]+rf(1, 2, 20)*10+beta0
  
  i=i+1
  
}

while(i<100+(length(x)/2)+1){
  
  y[i]=2*x[i]+rf(1, 2, 20)*10+beta0
  
  i=i+1
  
}



#----------------------------------Gr?ficar lineal-----------------------------------------

library(ggplot2)


g=ggplot(data.frame(y=y,x=x),aes(x=x,y=y))
g=g+geom_smooth(method=lm)
g=g+geom_point(shape=1)
g


#-------------------------------------------Ajustar------------------------------


ajuste=rfit(y~x)
summary(ajuste)





#--------------------------------------------------------------------------------------
#Prueba similar a la bondad de ajute
#--------------------------------------------------------------------------------------

#H0: B0=B1=...=Bk=0     vs  al menos un beta es distinto de cero



#d

#H0: El modelo no difiere de una constante vs H1: El modelo difiere de una constante




drop.test(ajuste)







#--------------------------------------------------------------------------------------
#Un factor 
#--------------------------------------------------------------------------------------


#modelo pobacional

#Yi=taui+Ei, i=1,..,k

#modelo poblacional

#Yij=tauij+Eij, i=1,..,k   j=1,...,nk


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
  
  y1[i]=tau1+rf(1, 2, 20)*10
  y2[i]=tau2+rf(1, 2, 20)*10
  y3[i]=tau3+rf(1, 2, 20)*10
  
  i=i+1  
}



#Crear ytotal y factor


y=c(y1,y2,y3)
x=c(rep(1,length(y1)),rep(2,length(y2)),rep(3,length(y3)))


#si no pones factor te va a tomar como si fuera un modelo lineal general

x=factor(x)



#Gr?ficar box-plot

boxplot(y~x,col="lightblue")






#--------------------------------------------------------------------------------------
#Armar el modelo
#--------------------------------------------------------------------------------------

options(contrasts = c("contr.sum","contr.poly"))

ajuste=rfit(y~x)
summary(ajuste)





#--------------------------------------------------------------------------------------
# prueba bondad de ajuste
#--------------------------------------------------------------------------------------

#H0= tau1=tau2=...=tauk      vs   H1: al menos un tau es diferente

#
modelo=oneway.rfit(y,x)
modelo


#--------------------------------------------------------------------------------------
# Ver las diferencias de los tau
#--------------------------------------------------------------------------------------

summary.oneway.rfit(modelo,method = "tukey") #compora 1 con 2, 1 con 3 y 3 con 2. Te muestra la diferencia
#no me dice cual es el valor estimado




#--------------------------------------------------------------------------------------
#Mas de un factor
#--------------------------------------------------------------------------------------


#Modelo Poblacional

#Yij=mu+alpai+Betaj+(alpa*beta)ij+Eij       i=1,...,a   j=1,...,b


#Modelo muestral

#Yijt=mu+alpai+Betaj+(alpa*beta)ij+Eijt       i=1,...,a   j=1,...,b t=1,...,nab




#-------------------------------------------------------------------------------------
#Crear modelo b?sico sin correlaci?n
#-------------------------------------------------------------------------------------

#Crear mu


mugeneral=500

#Crear alfa 

alfa1=200
alfa2=-100
alfa3=-100

#Crear beta 2

beta1=300
beta2=-300


#itentacciones

alfa1beta1=50
alfa2beta1=-100
alfa3beta1=50



alfa1beta2=-50
alfa2beta2=100
alfa3beta2=-50




#Crar des. est?ndar

sd1=20
sd2=20




y11=NULL
y12=NULL


i=1

while(i<201){
  
  y11[i]=mugeneral+alfa1+beta1+alfa1beta1+rf(1, 2, 20)*10
  y12[i]=mugeneral+alfa1+beta2+alfa1beta2+rf(1, 2, 20)*10
  
  
  i=i+1  
}






y21=NULL
y22=NULL


i=1

while(i<201){
  
  y21[i]=mugeneral+alfa2+beta1+alfa2beta1+rf(1, 2, 20)*10
  y22[i]=mugeneral+alfa2+beta2+alfa2beta2+rf(1, 2, 20)*10
  
  
  i=i+1  
}






y31=NULL
y32=NULL



i=1

while(i<201){
  
  y31[i]=mugeneral+alfa3+beta1+alfa3beta1+rf(1, 2, 20)*10
  y32[i]=mugeneral+alfa3+beta2+alfa3beta2+rf(1, 2, 20)*10
  
  
  i=i+1  
}


y=c(y11,y12,y21,y22,y31,y32)


#representa alfa son 3
alfa=c(rep("1",length(c(y11,y12))),rep("2",length(c(y21,y22))),rep("3",length(c(y31,y32))))

#representa beta son 2

beta=c(rep("1",length(y11)),rep("2",length(y12)),rep("1",length(y21)),rep("2",length(y22)),rep("1",length(y31)),rep("2",length(y32)))


alfa=factor(alfa)
levels(alfa)
beta=factor(beta)
levels(beta)



#--------------------------------Gr?ficar------------------------------------------------

library(ggplot2)
library(plotly)

g = ggplot(data.frame(y=y,alfa=alfa,beta=beta), aes(y=y,x=alfa)) 
g= g + geom_boxplot(aes(fill=alfa)) 
g= g+facet_grid(.~beta)
g=g+geom_jitter(alpha = 0.5, color = "tomato")
g=g+ylab("y")
g

p <- ggplotly(g)
p




#-------------------------------------------------------------------------------------------
#ajustar modelo
#-----------------------------------------------------------------------------------------




modelo.rangos=rfit(y~alfa+beta+alfa*beta)
summary(modelo.rangos)



#-------------------------------------------------------------------------------------------
#Anova
#-----------------------------------------------------------------------------------------

#H0: los alpa son iguales Vs   H1:los alpha no son iguales
#H0: los Beta son iguales Vs   H1:los Beta no son iguales
#H0: (alfa*beta)ij-(alfa*beta)iq-(alfa*beta)sj+(alfa*beta)sq=0    vs  (alfa*beta)ij-(alfa*beta)iq-(alfa*beta)sj+(alfa*beta)sq ditnto de cero

raov(y~alfa+beta+alfa*beta)  #no me dice cual es el valor estimado





















#-------------------------------------------------------------------------------------------
#Analysis de la covarianza
#-----------------------------------------------------------------------------------------


#Modelo poblacion

#Yi = mu+mui+B1*x1i+B2*x2i+...+ei   i=1,..,k


#Modelo muestral

#Yij = mu+mui+B1*x1ij+B2*x2ij+...+eij   i=1,..,k  j=1,...,nk

#-------------------------------------------------------------------------------------------
#para un solo covarianza
#-----------------------------------------------------------------------------------------





x=c(1:200)
y=NULL



beta0=0


fac1=10
fac2=-10

fact=sample(0:1,200,replace=T)






i=1

while(i<length(x)+1){
  
  y[i]=1*x[i]+rnorm(1, mean=0, sd=1)+beta0
  
  i=i+1
  
}




i=1

while(i<length(x)+1){
  
  if(fact[i]==1){
    
    y[i]=y[i]+fac1
  }else{
    
    y[i]=y[i]+fac2
  }
  i=i+1
  
}




#-----------------------------------------------------------------------------------------
#Graficar 
#-----------------------------------------------------------------------------------------

library(ggplot2)

g = ggplot(data.frame(y=y,x=x,fact=fact), aes(y=y,x=x,colour=factor(fact))) 
g=g+geom_smooth(method=lm)
g=g+geom_point(shape=1)
g



#-------------------------------------------------------------------------------------------
#ajustar modelo
#-----------------------------------------------------------------------------------------


#Para el Rfit no corre la funcion contrast: Toma en la ordenada el primer factor

modelo.rangos=rfit(y~fact+x)
summary(modelo.rangos)


#-------------------------------------------------------------------------------------------
#analisis anova
#-----------------------------------------------------------------------------------------



library(npsm)


y.group = cbind(y,fact)  #es la y con los dos factores
xcov = x               # es el cofactor



#para uno
#mide si las pendientes para cada grupo son iguales o no (Homog Slopes)

onecovaheter(2,y.group,xcov)



#Groups mide la covarianza si p <0.05 hay diferencia en los tau

#xint1 es el fact1+beta0
#xint2 es el fact2+beta0
  
#int3 beta 1


#-------------------------------------------------------------------------------------------
#para k covarianc
#-----------------------------------------------------------------------------------------





x=c(1:200)
y=NULL



beta0=0

#El supuesto es que tenes de un factor en este caso 3, cada factor tiene dos opciones o nive: es como alfa y beta

fac1=50
fac2=-50

fact=sample(0:1,200,replace=T)

fact1=sample(0:1,200,replace=T)

fact2=sample(0:1,200,replace=T)

#Supuesto el factor 1 y 2, con dos niveles cada uno no impactan en y



i=1

while(i<length(x)+1){
  
  y[i]=2*x[i]+rf(1, 10, 50)*10+beta0
  
  i=i+1
  
}




i=1

while(i<length(x)+1){
  
  if(fact[i]==1){
    
    y[i]=y[i]+fac1
  }else{
    
    y[i]=y[i]+fac2
  }
  i=i+1
  
}






#-----------------------------------------------------------------------------------------
#Graficar 
#-----------------------------------------------------------------------------------------

library(ggplot2)

g = ggplot(data.frame(y=y,x=x,fact=fact), aes(y=y,x=x,colour=factor(fact))) 
g=g+geom_smooth(method=lm)
g=g+geom_point(shape=1)
g


library(npsm)


#-------------------------------------------------------------------------------------------
#ajustar modelo
#-----------------------------------------------------------------------------------------




modelo.rangos=rfit(y~fact+fact1+fact2+x)
summary(modelo.rangos)

#-------------------------------------------------------------------------------------------
#anova
#-----------------------------------------------------------------------------------------



levels = c(2,2,2)      #numeros de niveles o de cada factor, supongamos si el factor2 tienen 3 niveles entonces levels = c(2,2,3) 
y.group = cbind(y,fact,fac1,fact2)  #es la y con el vector completo
xcov = x               # es el cofactor


#para k 

ajuste=kancova(levels,y.group,xcov,print.table=T)
ajuste




#xfull29 es la pendinte


#para uno








        