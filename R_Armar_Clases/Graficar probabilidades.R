
################################################################################
#library
################################################################################


library(bayesAB)
library(ggplot2)
library(Pareto)

#ver
#https://cran.r-project.org/web/packages/Pareto/vignettes/Pareto.html


################################################################################
#Normal
################################################################################

#-------------------------------------------------------------------------------
#Cola izquierda/Derecha


#media y sd
mean=0
sd=1

#intervalo donde va a ser pintada el area
lb=-10
ub=-1.565

#Crear valores funcion
x=seq(-4,4,length=100)*sd + mean
hx=dnorm(x,mean,sd)


#Graficar
plot(x, hx, type="n", xlab="z", ylab="",
     main="", axes=FALSE)

i <- x >= lb & x <= ub
lines(x, hx)
polygon(c(lb,x[i],ub), c(0,hx[i],0), col="blue")

area=pnorm(ub, mean, sd) - pnorm(lb, mean, sd)

#colocar Prob

#result=paste("P(X<490)=0,16")
#mtext(result,3)

#Configurar ejes

#Para que salga la linea
axis(1, at=c(-1000000,1000000), pos=0) 
axis(1, at=0, pos=0) 
axis(1, at=-1.565, pos=0) 








################################################################################
#Chi Cuadrado
################################################################################



#-------------------------------------------------------------------------------
#Cola izquierda/derecha


#media y sd
grado_libertad=199

#intervalo donde va a ser pintada el area
lb=0
ub=179.597

#Crear valores funcion
x=seq(120,280,length=100)
hx=dchisq(x,grado_libertad)


#Graficar
plot(x, hx, type="n", xlab="u", ylab="",
     main="", axes=FALSE)

i <- x >= lb & x <= ub
lines(x, hx)
polygon(c(lb,x[i],ub), c(0,hx[i],0), col="blue")


#colocar Prob

#result=paste("P(X<490)=0,16")
#mtext(result,3)

#Configurar ejes

#Para que salga la linea
axis(1, at=c(0,1000000), pos=0) 
axis(1, at=ub, pos=0) 








################################################################################
#t stundet
################################################################################



#-------------------------------------------------------------------------------
#Cola  izquierda/derecha


#media y sd
grado_libertad=199

#intervalo donde va a ser pintada el area
lb=-4
ub=-1.428

#Crear valores funcion
x=seq(-4,4,length=100)
hx=dt(x,grado_libertad)


#Graficar
plot(x, hx, type="n", xlab="t", ylab="",
     main="", axes=FALSE)

i <- x >= lb & x <= ub
lines(x, hx)
polygon(c(lb,x[i],ub), c(0,hx[i],0), col="blue")


#colocar Prob

#result=paste("P(X<490)=0,16")
#mtext(result,3)

#Configurar ejes

#Para que salga la linea
axis(1, at=c(-100,100), pos=0) 
axis(1, at=0, pos=0) 
axis(1, at=ub, pos=0) 





################################################################################
#F Fisher
################################################################################



#-------------------------------------------------------------------------------
#Cola izquierda/derecha


#media y sd
n_1=19
m_1=19

#intervalo donde va a ser pintada el area

lb=1
ub=10

#Crear valores funcion
x=seq(0,3,length=100)
hx=df(x,n_1,m_1)


#Graficar
plot(x, hx, type="n", xlab="f", ylab="",
     main="", axes=FALSE)

i <- x >= lb & x <= ub
lines(x, hx)
polygon(c(lb,x[i],ub), c(0,hx[i],0), col="blue")


#colocar Prob

#result=paste("P(X<490)=0,16")
#mtext(result,3)

#Configurar ejes

#Para que salga la linea
axis(1, at=c(-100,100), pos=0) 
axis(1, at=1, pos=0) 
























#-------------------------------------------------------------------------
#1 forma


mean=6200; sd=529.15
lb=7000; ub=100000

x <- seq(-4,4,length=100)*sd + mean
hx <- dnorm(x,mean,sd)

plot(x, hx, type="n", xlab="w", ylab="",
     main="Normal Distribution", axes=FALSE)

i <- x >= lb & x <= ub
lines(x, hx)
polygon(c(lb,x[i],ub), c(0,hx[i],0), col="red")

area <- pnorm(ub, mean, sd) - pnorm(lb, mean, sd)

axis(1, at=6200, pos=0) 
axis(1, at=c(lb), pos=0) 





#-------------------------------------------------------------------------
#2 forma


mean=3000; sd=600
lb=3000; ub=4000

x <- seq(-4,4,length=100)*sd + mean
hx <- dnorm(x,mean,sd)

plot(x, hx, type="n", xlab="u", ylab="",
     main="Normal Distribution", axes=FALSE)

i <- x >= lb & x <= ub
lines(x, hx)
polygon(c(lb,x[i],ub), c(0,hx[i],0), col="red")



axis(1, at=6000, pos=0) 
axis(1, at=c(0,lb,ub,9000), pos=0) 






#-------------------------------------------------------------------------
#3 Cola izquiera








 #-------------------------------------------------------------------------
#3 Expontential negativz
  


#-------------------------------------------------------------------------
#3 Cola izquiera
  
  
  
mean=500; sd=10
lb=-10; ub=490
  
x <- seq(-4,4,length=100)*sd + mean
hx <- dnorm(x,mean,sd)
  
plot(x, hx, type="n", xlab="x", ylab="",
       main="Distribución Normal", axes=FALSE)
  
i <- x >= lb & x <= ub
lines(x, hx)
polygon(c(lb,x[i],ub), c(0,hx[i],0), col="blue")
  
area <- pnorm(ub, mean, sd) - pnorm(lb, mean, sd)
result <- paste("P(X<490)=0,16")
mtext(result,3)
  
  
axis(1, at=c(-10,0,490,1000), pos=0) 
axis(1, at=500, pos=0) 

  


  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  




#----------------------------------------------------------------------------------
#Para X
#----------------------------------------------------------------------------------

g=plotPareto(1, 3)+ggtitle('Función densidad para X')
g=g+theme_bw()
g=g+xlab("x")+ylab("fX(x)")
g=g+xlim(0, 3)
g


#----------------------------------------------------------------------------------
#Para Y
#----------------------------------------------------------------------------------

g=plotPareto(1, 4)+ggtitle('Función densidad para Y')
g=g+theme_bw()
g=g+xlab("y")+ylab("fY(y)")
g=g+xlim(0, 3)
g
  

#----------------------------------------------------------------------------------
#Para X
#----------------------------------------------------------------------------------


x=as.numeric(seq(0,3,by=0.001))
plot(x,dPareto(x, 1, 3),xlab = "x", ylab = "fX(x)"
     ,main="Función densidad para X"
     ,col="darkgreen"
     ,type="l"
     ,lwd=1)
points(1,0, col = "darkgreen",pch=19) 
points(1,3, col = "darkgreen") 
      
        
    

#----------------------------------------------------------------------------------
#Para Y
#----------------------------------------------------------------------------------


y=as.numeric(seq(0,3,by=0.001))
plot(y,dPareto(x, 1, 4),xlab = "y", ylab = "fY(y)"
     ,main="Función densidad para Y"
     ,col="blue"
     ,type="l"
     ,lwd=1)
points(1,0, col = "blue",pch=19) 
points(1,4, col = "blue") 





#----------------------------------------------------------------------------------
#Para X entre 1.1 y 2
#----------------------------------------------------------------------------------

lb=1.1; ub=2

x=as.numeric(seq(0,3,by=0.001))
hx=dPareto(x, 1, 3)


plot(x,hx,xlab = "x", ylab = "fX(x)"
     ,main="Función densidad para X"
     ,col="darkgreen"
     ,type="l"
     ,lwd=1
     ,axes=FALSE)
points(1,0, col = "darkgreen",pch=19) 
points(1,3, col = "darkgreen") 



i <- x >= lb & x <= ub
lines(x, hx,col = "darkgreen")
polygon(c(lb,x[i],ub), c(0,hx[i],0), col="red")

area=0.6260

result=paste("P(",lb,"< X <",ub,") =",
                signif(area, digits=3))
mtext(result,3)


axis(1, at=c(0,lb,ub,3), pos=0) 
axis(2, at=seq(0, 3, 1), pos=0) 



#----------------------------------------------------------------------------------
#Para Y entre 1.1 y 2
#----------------------------------------------------------------------------------

lb=1.1; ub=2

y=as.numeric(seq(0,3,by=0.001))
hx=dPareto(x, 1, 4)


plot(y,hx,xlab = "y", ylab = "fY(y)"
     ,main="Función densidad para Y"
     ,col="blue"
     ,type="l"
     ,lwd=1
     ,axes=FALSE)
points(1,0, col = "blue",pch=19) 
points(1,4, col = "blue") 



i <- x >= lb & x <= ub
lines(x, hx,col = "blue")
polygon(c(lb,x[i],ub), c(0,hx[i],0), col="red")

area=0.6205

result=paste("P(",lb,"< Y <",ub,") =",
             signif(area, digits=4))
mtext(result,3)


axis(1, at=c(0,lb,ub,3), pos=0) 
axis(2, at=seq(0, 4, 1), pos=0) 































