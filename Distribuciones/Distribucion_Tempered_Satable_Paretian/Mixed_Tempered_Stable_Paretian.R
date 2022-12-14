################################################################################
#Cargar librerias
################################################################################


library(MixedTS)


################################################################################
#Generar Datos
################################################################################




################################################################################
#Referencias
################################################################################


#https://www.rdocumentation.org/packages/MixedTS/versions/1.0.4/topics/dMixedTS-methods
#https://rdrr.io/rforge/MixedTS/man/mle.MixedTS.html


#-------------------------------------------------------------------------------
# Density of MixedTS with Gamma
#-------------------------------------------------------------------------------

ParamEx1=setMixedTS.param(mu0=0, mu=0, sigma=0.047, a=4,
                           alpha=1.4, lambda_p=1.4665, lambda_m=1.9103, 
                           Mixing="Gamma")



#help(setMixedTS.param)



#-------------------------------------------------------------------------------
#Generar valores aleatorios

dat=rMixedTS(10000,object=ParamEx1)

dat@Data


#-------------------------------------------------------------------------------
#Determinar funcion densidad

x=seq(-0.2,0.2,length=10000)

dens=dMixedTS(x=x,object=ParamEx1)

#help(dMixedTS)

#-------------------------------------------------------------------------------
#Graficar

plot(dens)
hist(dat@Data, breaks = 60,probability = TRUE,border = "white",
     col = "lightgray",add=T)




################################################################################
#Estimacion
################################################################################



ParamEx2=setMixedTS.param(mu0=0, mu=0, sigma=1, a=2,
                          alpha=1.5, lambda_p=1, lambda_m=1, 
                          Mixing="Gamma")



Rand1=rMixedTS(x=5000,object=ParamEx2, setSup=10,setInf=-10,N=2^9)


Rand1@Data=dat@Data


est1=mle.MixedTS(object=Rand1 , setSup=10,setInf=-10,N=2^9)


summary(est1)




ParamEstimate=setMixedTS.param(mu0=0.002536208, mu= -0.426604638, sigma= 0.035314697, a=6.153769541,
                          alpha= 1.230753892, lambda_p= 3.073330898, lambda_m= 3.177889888, 
                          Mixing="Gamma")




#-------------------------------------------------------------------------------
#Determinar funcion densidad

x=seq(-0.2,0.2,length=10000)

dens=dMixedTS(x=x,object=ParamEstimate)

#help(dMixedTS)

#-------------------------------------------------------------------------------
#Graficar

plot(dens)
hist(dat@Data, breaks = 60,probability = TRUE,border = "white",
     col = "lightgray",add=T)









################################################################################
#Estimacion acciones
################################################################################

#################################################################################
#Descargar datos
################################################################################

f_init_test='2012-05-01'
f_final_test='2022-09-01'

#----------------------------------------------------------------------
#Cargar simbolo

eq_1="V"


#------------------------------------------------------------------
#Descargar datos de yahoo


data_eq_1=new.env()
suppressWarnings(try(getSymbols(eq_1, src = 'yahoo', from = f_init_test, 
                                to=f_final_test,env = data_eq_1, auto.assign = T,
                                periodicity = "d"),silent = T))
suppressWarnings(try(for(i in ls(data_eq_1)) data_eq_1[[i]] =
                       adjustOHLC(data_eq_1[[i]], use.Adjusted=T),silent = TRUE)) 




#-------------------------------------------------------------------
#Extraer Precio ajustado

y_t=as.numeric(data_eq_1[[eq_1]][,6])

x_t=log(y_t)


dx_t=diff(x_t)

plot(dx_t, type="l")



model=garchFit( ~ arma(1,3) + garch(1, 1),data=dx_t, trace = F)

plot(model@residuals, type="l")





ParamEx2=setMixedTS.param(mu0=0, mu=0, sigma=1, a=2,
                          alpha=1.5, lambda_p=1, lambda_m=1, 
                          Mixing="Gamma")



Rand1=rMixedTS(x=5000,object=ParamEx2, setSup=10,setInf=-10,N=2^9)


Rand1@Data=model@residuals*100


est1=mle.MixedTS(object=Rand1 , setSup=10,setInf=-10,N=2^9)


summary(est1)




ParamEstimate=setMixedTS.param(mu0= 0.42732551 , mu=  -0.19005606, sigma=   0.71361240  , a=  4.89189348  ,
                               alpha=    1.52875518, lambda_p=  0.15281382 , lambda_m=0.09713054, 
                               Mixing="Gamma")







#-------------------------------------------------------------------------------
#Determinar funcion densidad

x=seq(-10,10,length=10000)

dens=dMixedTS(x=x,object=ParamEstimate)

#help(dMixedTS)

#-------------------------------------------------------------------------------
#Graficar

plot(dens)
hist(model@residuals*100, breaks = 200,probability = TRUE,border = "white",
     col = "lightgray",add=T)




################################################################################
#Teorema del limite central
################################################################################




#-------------------------------------------------------------------------------
# Density of MixedTS with Gamma
#-------------------------------------------------------------------------------

ParamEx1=setMixedTS.param(mu0=0, mu=0, sigma=0.047, a=4,
                          alpha=1.4, lambda_p=1.4665, lambda_m=1.9103, 
                          Mixing="Gamma")



#help(setMixedTS.param)



#-------------------------------------------------------------------------------
#Generar valores aleatorios


xBarra=NULL


i=1

while (i<=5000) {
  
  
  dat=rMixedTS(1000,object=ParamEx1)
  
  xBarra[i]=mean(dat@Data)
  
  i=i+1
  
}



hist(xBarra)

shapiro.test(xBarra)



